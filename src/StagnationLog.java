import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

/**
 * Journal CSV du diagnostic de decorrelation et de stagnation, une ligne par temperature et un
 * fichier par replique (couple fichier d'entree / champ), calque sur {@link EquilibrationLog}.
 *
 * <p>Ce fichier repond a un besoin different de {@code equilibration.csv} : celui-ci dit si la
 * thermalisation a converge et combien de temps elle a pris ; celui-la dit si, en plus, la
 * configuration de spins retenue est <i>utilisable</i> pour du criblage rapide — decorrelee de la
 * temperature precedente (memoire des domaines, recouvrement q, cf. {@link StagnationPoint}) et
 * decorrelee d'elle-meme d'une mesure de production a l'autre — ou si l'algorithme a du
 * constater qu'il ne pouvait plus ameliorer l'etat courant et s'arreter par stagnation plutot
 * que par convergence.</p>
 *
 * <p>Ecriture immediate (pas de tampon differe entre deux temperatures) : une simulation
 * interrompue laisse un fichier exploitable jusqu'a la derniere temperature terminee.</p>
 *
 * <p>Toutes les valeurs numeriques sont formatees en {@link Locale#US} (point decimal), de
 * facon a etre relisibles par pandas / numpy quelle que soit la locale de la machine, avec les
 * memes conventions que {@link EquilibrationLog} (NaN et infinis ecrits tels quels, booleens en
 * 0/1) pour que les deux CSV soient depouilles par le meme script.</p>
 *
 * <p>Exceptions non verifiees, comme {@link EquilibrationLog} : les deux journaux sont ecrits
 * cote a cote dans la meme boucle (typiquement {@code options.log.write(point)} sans
 * {@code try/catch}, cf. {@code MonteCarlo.runParallel}) ; si l'un lancait une exception verifiee
 * et l'autre non, l'appelant devrait envelopper l'un des deux pour rien, et le code de la boucle
 * deviendrait asymetrique sans aucun gain — une panne d'ecriture (disque plein, systeme de
 * fichiers en lecture seule) est de toute facon fatale a la replique, verifiee ou non.</p>
 */
public final class StagnationLog implements AutoCloseable {

    /** Une colonne : son nom d'en-tete et la facon de l'extraire d'un {@link StagnationPoint}. */
    public record Column(String name, java.util.function.Function<StagnationPoint, String> format) { }

    /**
     * Definition <b>unique</b> des colonnes : l'en-tete et chaque ligne en sont derives, ils ne
     * peuvent donc pas diverger (voir {@code StagnationLogTest}). L'ordre est celui des
     * composantes de {@link StagnationPoint}.
     */
    public static final java.util.List<Column> COLUMNS = java.util.List.of(
            new Column("mcIdx",                    p -> Integer.toString(p.mcIdx())),
            new Column("T",                        p -> num(p.T())),
            new Column("thermSweeps",              p -> Long.toString(p.thermSweeps())),
            new Column("prodSweeps",               p -> Long.toString(p.prodSweeps())),
            new Column("qInf",                     p -> num(p.qInf())),
            new Column("tauOverlapSweeps",         p -> num(p.tauOverlapSweeps())),
            new Column("memoryRatio",              p -> num(p.memoryRatio())),
            new Column("temperaturesDecorrelated", p -> p.temperaturesDecorrelated() ? "1" : "0"),
            new Column("driftZ",                   p -> num(p.driftZ())),
            new Column("driftPerSweep",            p -> num(p.driftPerSweep())),
            new Column("stopReason",               p -> sanitize(p.stopReason())),
            new Column("converged",                p -> p.converged() ? "1" : "0"),
            new Column("measureEvery",             p -> Integer.toString(p.measureEvery())),
            new Column("tauProductionSweeps",      p -> num(p.tauProductionSweeps())),
            new Column("measurementsPerTau",       p -> num(p.measurementsPerTau())),
            new Column("measurementsDecorrelated", p -> p.measurementsDecorrelated() ? "1" : "0"),
            new Column("neffAchieved",             p -> num(p.neffAchieved())),
            new Column("neffSaturated",            p -> p.neffSaturated() ? "1" : "0"),
            new Column("seconds",                  p -> num(p.seconds())),
            new Column("Hx",                       p -> num(p.hx())),
            new Column("Hy",                       p -> num(p.hy())),
            new Column("Hz",                       p -> num(p.hz())));

    public static final String HEADER = COLUMNS.stream().map(Column::name)
            .collect(java.util.stream.Collectors.joining(","));

    private final Path file;
    private final BufferedWriter out;
    private boolean closed;

    /**
     * Ouvre (ou ecrase) le fichier et y ecrit l'en-tete.
     *
     * @throws UncheckedIOException si le fichier ne peut pas etre cree (les repertoires parents
     *         sont crees au besoin). Si l'ouverture reussit mais l'ecriture de l'en-tete echoue
     *         (disque plein, systeme de fichiers en lecture seule remonte tardivement), le writer
     *         est ferme avant que l'exception ne soit propagee : le constructeur echoue, personne
     *         ne recevra l'objet, donc personne n'appellera {@link #close()} — le descripteur
     *         serait perdu jusqu'a la fin du processus, et une production lance des milliers de
     *         repliques.
     */
    public StagnationLog(Path file) {
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
     * Chemin par defaut d'un journal : {@code <outputDir>/decorrelation.csv}.
     *
     * <p>Meme convention que {@link EquilibrationLog#defaultPath(String)} : le nom est
     * deterministe et identique pour toutes les repliques, car chaque replique ecrit dans SON
     * repertoire de sortie (date-heure du lancement, fichier d'entree et champ font partie du
     * nom du repertoire, cf. {@code AppParallel}). Seul le prefixe distingue, au tri
     * alphabetique, les deux journaux d'une meme replique, ecrits cote a cote.</p>
     */
    public static Path defaultPath(String outputDir) {
        return Paths.get(outputDir == null ? "." : outputDir, "decorrelation.csv");
    }

    /**
     * Ecrit une ligne de diagnostic.
     *
     * <p>Un seul parametre, et non les dix-huit valeurs d'avant : voir {@link StagnationPoint}
     * pour la justification, identique a celle de {@link TemperaturePoint}.</p>
     */
    public void write(StagnationPoint p) {
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

    /**
     * Assainit {@code stopReason} pour qu'il ne puisse jamais casser le CSV : une virgule
     * decalerait toutes les colonnes suivantes, un saut de ligne creerait une ligne fantome. On
     * choisit de remplacer plutot que de refuser (une raison d'arret inattendue ne doit pas faire
     * echouer toute une replique de calcul) : tout caractere hors {@code [A-Za-z0-9_-]} devient
     * {@code '_'}. Le choix est documente ici pour qu'il soit visible du depouillement : une
     * valeur avec des underscores a la place de ponctuation est un signal, pas un bug.
     */
    private static String sanitize(String stopReason) {
        return stopReason.replaceAll("[^A-Za-z0-9_-]", "_");
    }

    /** {@code %.10g} en {@link Locale#US}; NaN et les infinis sont ecrits tels quels. */
    private static String num(double v) {
        if (Double.isNaN(v)) return "NaN";
        if (Double.isInfinite(v)) return v > 0 ? "Inf" : "-Inf";
        return String.format(Locale.US, "%.10g", v);
    }

    public void flush() {
        if (closed) return;
        try {
            out.flush();
        } catch (IOException e) {
            throw new UncheckedIOException("flush impossible sur " + file, e);
        }
    }

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
}
