import java.io.BufferedWriter;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Locale;

/**
 * Journal CSV de la susceptibilite AC, une ligne par couple (temperature, periode) — copie
 * structurelle stricte de {@link EquilibrationLog} (memes conventions de format, de nettoyage de
 * nom de fichier et d'ecriture immediate) appliquee aux points {@link ACPoint} produits par
 * {@link ACSusceptibility#aggregate}.
 *
 * <p>Ecriture immediate (flush apres chaque ligne) : une exploration interrompue laisse un
 * fichier exploitable jusqu'au dernier point termine. Toutes les valeurs numeriques sont
 * formatees en {@link Locale#US} (point decimal), lisibles par pandas / numpy quelle que soit la
 * locale de la machine.</p>
 */
public final class ACSusceptibilityLog implements AutoCloseable {

    /** Une colonne : son nom d'en-tete et la facon de l'extraire d'un {@link ACPoint}. */
    public record Column(String name, java.util.function.Function<ACPoint, String> format) { }

    /**
     * Definition <b>unique</b> des colonnes : l'en-tete et chaque ligne en sont derives, ils ne
     * peuvent donc pas diverger (voir {@code ACSusceptibilityLogTest}).
     */
    public static final java.util.List<Column> COLUMNS = java.util.List.of(
            new Column("T",               p -> num(p.T())),
            new Column("omega",           p -> num(p.omega())),
            new Column("periodSweeps",    p -> Integer.toString(p.periodSweeps())),
            new Column("h0",              p -> num(p.h0())),
            new Column("hStatic",         p -> num(p.hStatic())),
            new Column("nReplicas",       p -> Integer.toString(p.nReplicas())),
            new Column("chiPrime",        p -> num(p.chiPrime())),
            new Column("chiPrimeErr",     p -> num(p.chiPrimeErr())),
            new Column("chiSecond",       p -> num(p.chiSecond())),
            new Column("chiSecondErr",    p -> num(p.chiSecondErr())),
            new Column("thermSweepsMean", p -> num(p.meanThermSweeps())),
            new Column("acceptanceMean",  p -> num(p.meanAcceptance())),
            new Column("secondsMean",     p -> num(p.meanSeconds())));

    public static final String HEADER = COLUMNS.stream().map(Column::name)
            .collect(java.util.stream.Collectors.joining(","));

    private final Path file;
    private final BufferedWriter out;
    private boolean closed;

    /**
     * Ouvre (ou ecrase) le fichier et y ecrit l'en-tete.
     *
     * @throws UncheckedIOException si le fichier ne peut pas etre cree (les repertoires parents
     *         sont crees au besoin). Si l'ouverture reussit mais l'ecriture de l'en-tete echoue,
     *         le writer est ferme avant que l'exception ne soit propagee : le constructeur
     *         echoue, personne ne recevra l'objet, donc personne n'appellera {@link #close()}.
     */
    public ACSusceptibilityLog(Path file) {
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
     * Chemin par defaut d'un journal : {@code <outputDir>/ac_susceptibility_<fileName>_h<champ>.csv}.
     *
     * <p>Memes regles que {@link EquilibrationLog#defaultPath} : le nom du fichier d'entree est
     * nettoye des caracteres qui posent probleme dans un nom de fichier (tout sauf
     * {@code [A-Za-z0-9._-]}, remplace par {@code _}), un eventuel suffixe {@code .txt} est
     * retire, et le champ est formate en {@link Locale#US} avec six decimales.</p>
     */
    public static Path defaultPath(String outputDir, String fileName, double field) {
        String safe = (fileName == null ? "replica" : fileName).replaceAll("[^A-Za-z0-9._-]", "_");
        if (safe.endsWith(".txt")) safe = safe.substring(0, safe.length() - 4);
        String h = String.format(Locale.US, "%.6f", field);
        return Paths.get(outputDir == null ? "." : outputDir,
                "ac_susceptibility_" + safe + "_h" + h + ".csv");
    }

    /**
     * Ecrit une ligne. Un seul parametre nomme plutot que treize {@code double} consecutifs, pour
     * la meme raison que {@link EquilibrationLog#write(TemperaturePoint)} : rien dans le
     * compilateur ne protege contre l'interversion silencieuse de deux arguments positionnels.
     */
    public void write(ACPoint p) {
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
