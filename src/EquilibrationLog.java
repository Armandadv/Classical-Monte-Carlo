import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

/*
 * Pour AVRO (ecrivains des spins, fusion de Output.java)
 */
import org.apache.avro.Schema;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;

/**
 * Journal CSV du diagnostic d'equilibration, une ligne par temperature et un fichier par
 * replique (couple fichier d'entree / champ).
 *
 * <p>Les colonnes de diagnostic de {@code MonteCarlo.observablesData}
 * ({@code MonteCarlo.COL_THERM} .. {@code MonteCarlo.COL_G_PROD}) ne sont pas ecrites dans un
 * fichier Avro d'observables : depuis la fusion du 12/09/2026, les colonnes d'observables
 * vivent fusionnees dans CE csv, au layout historique fige, compatible avec les depouillements
 * existants. Ce CSV est le seul endroit ou l'on retrouve, en clair, la longueur de
 * thermalisation reellement consommee, t0, tau, g, N_eff (prevu <i>et</i> obtenu), les z de
 * Geweke et le temps passe a chaque temperature. C'est ce qui permet de verifier <i>a
 * posteriori</i> que le controleur adaptatif a bien concentre le calcul autour de la
 * transition.</p>
 *
 * <p>Cette classe heberge aussi les ecrivains Avro des spins (fusion de Output.java, supprime) :
 * {@link #spinOutputAvro} ecrit l'avro final {@code spins_finaux_h<champ>.avro} et
 * {@link #spinsSnapshotAvro} les instantanes par temperature {@code spins_T<tag>_h<champ>.avro}
 * (option --spins-per-T). Les noms sont deterministes : chaque replique ecrit dans SON
 * repertoire de sortie (date-heure, fichier d'entree et champ, cf. AppParallel), le fichier
 * n'a donc plus besoin de les porter.</p>
 *
 * <p>Ecriture immediate (pas de tampon differe entre deux temperatures) : une simulation
 * interrompue laisse un fichier exploitable jusqu'a la derniere temperature terminee.</p>
 *
 * <p>Toutes les valeurs numeriques sont formatees en {@link Locale#US} (point decimal), de
 * facon a etre relisibles par pandas / numpy quelle que soit la locale de la machine.</p>
 */
public final class EquilibrationLog implements AutoCloseable {

    /**
     * En-tete du fichier. L'ordre est celui des composantes de {@link TemperaturePoint} ecrites
     * par {@link #write(TemperaturePoint)} — les deux se lisent cote a cote et se corrigent
     * ensemble.
     */
    /** Une colonne : son nom d'en-tete et la facon de l'extraire d'un {@link TemperaturePoint}. */
    public record Column(String name, java.util.function.Function<TemperaturePoint, String> format) { }

    /**
     * Definition <b>unique</b> des colonnes : l'en-tete et chaque ligne en sont derives, ils ne
     * peuvent donc pas diverger (voir {@code EquilibrationLogTest}).
     */
    public static final java.util.List<Column> COLUMNS = java.util.List.of(
            new Column("mcIdx",         p -> Integer.toString(p.mcIdx())),
            new Column("T",             p -> num(p.T())),
            new Column("thermSweeps",   p -> Long.toString(p.thermSweeps())),
            new Column("t0",            p -> Long.toString(p.t0Sweep())),
            new Column("tau",           p -> num(p.tauSweeps())),
            new Column("g",             p -> num(p.g())),
            new Column("neff",          p -> num(p.neff())),
            new Column("zE",            p -> num(p.zE())),
            new Column("zM",            p -> num(p.zM())),
            new Column("prodSweeps",    p -> Integer.toString(p.prodSweeps())),
            new Column("converged",     p -> p.converged() ? "1" : "0"),
            new Column("acceptance",    p -> num(p.acceptance())),
            new Column("sigma",         p -> num(p.sigma())),
            new Column("sigmaTherm",    p -> num(p.sigmaEndOfThermalization())),
            new Column("ePerSite",      p -> num(p.ePerSite())),
            new Column("eErr",          p -> num(p.eErr())),
            new Column("mNorm",         p -> num(p.mNorm())),
            new Column("seconds",       p -> num(p.seconds())),
            new Column("measureEvery",  p -> Integer.toString(p.measureEvery())),
            new Column("nMeasurements", p -> Integer.toString(p.nMeasurements())),
            new Column("neffAchieved",  p -> num(p.neffAchieved())),
            new Column("gProduction",   p -> num(p.gProduction())),
            // Observables fusionnees (12/09/2026) : l'ancien fichier Observables_*.csv est
            // supprime, ses colonnes uniques vivent ici — le champ et la serie de production
            // de laquelle VarE (matiere premiere de C_v) et les composantes de M sont tirees.
            new Column("Hx",     p -> num(p.hx())),
            new Column("Hy",     p -> num(p.hy())),
            new Column("Hz",     p -> num(p.hz())),
            new Column("Std_E",  p -> num(p.stdE())),
            new Column("VarE",   p -> num(p.varE())),
            new Column("Mx",     p -> num(p.mx())),
            new Column("Std_Mx", p -> num(p.stdMx())),
            new Column("My",     p -> num(p.my())),
            new Column("Std_My", p -> num(p.stdMy())),
            new Column("Mz",     p -> num(p.mz())),
            new Column("Std_Mz", p -> num(p.stdMz())));

    /** En-tete CSV effectif : les 33 noms de {@link #COLUMNS} joints par des virgules. */
    public static final String HEADER = COLUMNS.stream().map(Column::name)
            .collect(java.util.stream.Collectors.joining(","));

    private final Path file;
    private final BufferedWriter out;
    private boolean closed;

    /**
     * Ouvre (ou ecrase) le fichier et y ecrit l'en-tete.
     *
     * @throws UncheckedIOException si le fichier ne peut pas etre cree (les repertoires
     *         parents sont crees au besoin). Si l'ouverture reussit mais l'ecriture de
     *         l'en-tete echoue (disque plein, systeme de fichiers en lecture seule remonte
     *         tardivement), le writer est ferme avant que l'exception ne soit propagee : le
     *         constructeur echoue, personne ne recevra l'objet, donc personne n'appellera
     *         {@link #close()} — le descripteur serait perdu jusqu'a la fin du processus, et
     *         une production lance des milliers de repliques.
     */
    public EquilibrationLog(Path file) {
        if (file == null) throw new IllegalArgumentException("file == null");
        this.file = file.toAbsolutePath();
        BufferedWriter w;
        try {
            Path parent = this.file.getParent();
            if (parent != null) Files.createDirectories(parent);
            w = Files.newBufferedWriter(this.file, StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new UncheckedIOException("impossible d'ouvrir " + this.file, e);
        }
        try {
            w.write(HEADER);
            w.newLine();
            w.flush();
        } catch (IOException e) {
            try {
                w.close();
            } catch (IOException suppressed) {
                e.addSuppressed(suppressed);
            }
            throw new UncheckedIOException("impossible d'ecrire l'en-tete de " + this.file, e);
        }
        this.out = w;
    }

    /** Chemin effectif du fichier. */
    public Path file() { return file; }

    /**
     * Chemin par defaut d'un journal : {@code <outputDir>/equilibration.csv}.
     *
     * <p>Le nom est deterministe et le meme pour toutes les repliques : chaque replique ecrit
     * dans SON repertoire de sortie — la date-heure du lancement, le fichier d'entree et le
     * champ font partie du NOM DU REPERTOIRE (cf. AppParallel), plus besoin de les porter dans
     * le nom du CSV.</p>
     */
    public static Path defaultPath(String outputDir) {
        return Paths.get(outputDir == null ? "." : outputDir, "equilibration.csv");
    }

    /**
     * Ecrit une ligne de diagnostic.
     *
     * <p>Un seul parametre, et non les dix-huit valeurs d'avant : une signature faite de dix-huit
     * {@code double} consecutifs n'offre aucune protection contre l'interversion de deux
     * arguments, et le compilateur ne peut rien signaler. {@link TemperaturePoint} nomme chaque
     * champ, et l'ordre des colonnes est defini une fois pour toutes par {@link #HEADER}.</p>
     */
    public void write(TemperaturePoint p) {
        if (p == null) throw new IllegalArgumentException("p == null");
        if (closed) throw new IllegalStateException("journal deja ferme : " + file);
        StringBuilder sb = new StringBuilder(256);
        for (int i = 0; i < COLUMNS.size(); i++) {
            if (i > 0) sb.append(',');
            sb.append(COLUMNS.get(i).format().apply(p));
        }
        try {
            out.write(sb.toString());
            out.newLine();
            out.flush();
        } catch (IOException e) {
            throw new UncheckedIOException("ecriture impossible dans " + file, e);
        }
    }

    /** {@code %.10g} en {@link Locale#US}; NaN et les infinis sont ecrits tels quels. */
    private static String num(double v) {
        if (Double.isNaN(v)) return "NaN";
        if (Double.isInfinite(v)) return v > 0 ? "Inf" : "-Inf";
        return String.format(Locale.US, "%.10g", v);
    }

    /** Force l'ecriture immediate des lignes deja produites ; sans effet apres {@link #close()}. */
    public void flush() {
        if (closed) return;
        try {
            out.flush();
        } catch (IOException e) {
            throw new UncheckedIOException("flush impossible sur " + file, e);
        }
    }

    /** Ferme le fichier ; idempotent (sans effet si deja ferme). */
    @Override
    public void close() {
        if (closed) return;
        closed = true;
        try {
            out.close();
        } catch (IOException e) {
            throw new UncheckedIOException("fermeture impossible de " + file, e);
        }
    }

    // -------------------------------------------------- ecrivains Avro des spins
    // Fusion de Output.java (supprime) : memes schemas, meme snappy ; nommage deterministe
    // (le repertoire de sortie porte l'identite de la replique) ; IOException avalees
    // volontairement (printStackTrace), le controle de fichiers produits d'AppParallel
    // s'appuie dessus.

    /**
     * Avro final {@code spins_finaux_h<champ>.avro} de la replique (S_x, S_y, S_z par site,
     * snappy) : etat final des spins apres le recuit. Le champ est formate en {@link Locale#US}
     * avec six decimales, comme le nom du repertoire de sortie d'AppParallel. Les IOException
     * sont avalees volontairement ({@code printStackTrace}) : le controle de fichiers produits
     * d'AppParallel s'appuie dessus.
     */
    public static void spinOutputAvro(Lattice lattice, String outputDirPath, double field) {
        try {

            Schema schema = new Schema.Parser().parse("{" +
                    "  \"type\": \"record\"," +
                    "  \"name\": \"Spin\"," +
                    "  \"fields\": [" +
                    "    {\"name\": \"S_x\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_y\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_z\", \"type\": \"float\"}" +
                    "  ]" +
                    "}");

            File outputDir = new File(outputDirPath);
            if (!outputDir.exists()) {
                outputDir.mkdirs();  // Créer le répertoire s'il n'existe pas
            }
            File file = new File(outputDir,
                    "spins_finaux_h" + String.format(Locale.US, "%.6f", field) + ".avro");
            DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(schema);
            DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter);
            dataFileWriter.setCodec(CodecFactory.snappyCodec());
            dataFileWriter.create(schema, file);

            for (int i = 0; i < lattice.size; i++) {
                double[] S_xyz = lattice.getSpin1D(i);
                GenericRecord record = new GenericData.Record(schema);
                record.put("S_x", (double) S_xyz[0]);
                record.put("S_y", (double) S_xyz[1]);
                record.put("S_z", (double) S_xyz[2]);
                dataFileWriter.append(record);
            }
            dataFileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Instantane des spins a une temperature donnee (option --spins-per-T) : meme contenu que
     * l'avro final, nomme {@code spins_T<tag>_h<champ>.avro} — {@code tag} (typiquement la
     * temperature formatee en %.5f par l'appelant) pour etre triable par temperature, le champ
     * formate en {@link Locale#US} avec six decimales comme le repertoire de sortie. IOException
     * avalees volontairement ({@code printStackTrace}), comme {@link #spinOutputAvro}.
     */
    public static void spinsSnapshotAvro(Lattice lattice, String outputDirPath, String tag, double field) {
        try {
            Schema schema = new Schema.Parser().parse("{" +
                    "  \"type\": \"record\"," +
                    "  \"name\": \"Spin\"," +
                    "  \"fields\": [" +
                    "    {\"name\": \"S_x\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_y\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_z\", \"type\": \"float\"}" +
                    "  ]" +
                    "}");

            File outputDir = new File(outputDirPath);
            if (!outputDir.exists()) {
                outputDir.mkdirs();  // Créer le répertoire s'il n'existe pas
            }
            File file = new File(outputDir,
                    "spins_T" + tag + "_h" + String.format(Locale.US, "%.6f", field) + ".avro");
            DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(schema);
            DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter);
            dataFileWriter.setCodec(CodecFactory.snappyCodec());
            dataFileWriter.create(schema, file);

            for (int i = 0; i < lattice.size; i++) {
                double[] S_xyz = lattice.getSpin1D(i);
                GenericRecord record = new GenericData.Record(schema);
                record.put("S_x", (double) S_xyz[0]);
                record.put("S_y", (double) S_xyz[1]);
                record.put("S_z", (double) S_xyz[2]);
                dataFileWriter.append(record);
            }
            dataFileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
