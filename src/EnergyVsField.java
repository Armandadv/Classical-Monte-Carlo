import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

import org.apache.avro.Schema;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;

import org.apache.commons.math3.util.FastMath;

/**
 * Outil autonome modernise d'apres {@code E_v_field} (supprime) : energie d'etats imposes
 * (tables de signes Ising et flop) en fonction du champ magnetique, pour l'identification des
 * phases. Une CLI en regle, une grille de champ en entier et une sortie Avro deterministe.
 *
 * <h2>Usage</h2>
 * <pre>
 * java -cp "bin:lib/*" EnergyVsField &lt;inputDir&gt; &lt;outputDir&gt; &lt;hInit&gt; &lt;hMax&gt; &lt;dH&gt;
 *                                  [--field-dir x,y,z] [--size La,Lb]
 * </pre>
 * <p>Cinq arguments positionnels exactement (sinon usage sur {@code stderr} et code 2) ; options
 * n'importe ou sur la ligne, le dernier lu gagne : {@code --field-dir x,y,z} (direction du champ
 * dans le tenseur g, defaut -1,-1,2) et {@code --size La,Lb} (reseau La x Lb, entiers &gt;= 1,
 * defaut 12,12). Pour chaque fichier {@code bcaoExplor*} du repertoire d'entree, le reseau est
 * construit par {@link ReplicaFactory}, puis le champ balaie {@code h_k = hInit + k*dH} ; a
 * chaque pas un record Avro est ecrit avec l'energie par site des etats Ising (E0..E3) et
 * l'energie minimale des etats flop sur (slot, phi) — E_flop, Phi, period. Sorties :
 * {@code outputDir/config_energy_H_<input>.avro} et son miroir CSV lisible directement
 * ({@code .csv}, memes valeurs apres arrondi float).</p>
 *
 * <h2>Etats imposes : tables de signes</h2>
 * <p>Les etats ne sont plus construits par cosinus mais par <b>tables de signes</b> — le double
 * zigzag valide par Monte-Carlo par l'utilisateur. Chaque site (i, j, t) porte le spin
 * {@code signe × B_KITAEV}, avec {@code signe = sigA[j]} sur le sous-reseau A (t = 0) et
 * {@code sigA[j + decalage]} sur B (decalage 1 : la rangee B suit la rangee A). E0 = FERRO
 * (tous les spins +, moment net +N) ; E1 = zigzag 1/2 {+,-}, moment net 0 ; E2 = uud 1/3
 * {+,+,-}, moment net +N/3 — l'origine du plateau d'aimantation a 1/3 — ; E3 = double zigzag
 * 1/4 {+,+,-,-} (l'ex-{@code PHASE_SIGN} de l'ancien code), moment net 0. Les etats flop
 * gardent la mecanique Sp/Sm de l'original (angle phi), mais avec le signe lu dans les tables
 * uud et double zigzag : deux slots seulement, de periodes 1/3 et 1/4.</p>
 *
 * <h2>Grille de champ</h2>
 * <p>Boucle en ENTIER : {@code h_k = hInit + k*dH} pour {@code k = 0..nSteps-1}, avec
 * {@code nSteps = floor((hMax - hInit)/dH + 1e-9) + 1} — hMax INCLUSIF, meme convention
 * qu'{@code AppParallel}. {@code h_{nSteps-1}} peut donc depasser hMax d'un ulp : grille voulue.
 * Rupture assumee avec la boucle historique {@code H_field += h_step}, dont l'erreur d'arrondi
 * cumulee perdait le dernier champ selon la valeur binaire de h_step. Au-dela de
 * {@link #MAX_STEPS} pas : rejet (EXIT_USAGE), aucun clamp.</p>
 *
 * <h2>Codes de retour</h2>
 * <p>{@link #run(String[])} renvoie 0 en cas de succes, 1 en cas d'echec d'execution, 2 pour une
 * erreur d'usage. Elle ne leve rien et n'appelle jamais {@code System.exit} : c'est
 * {@link #main(String[])} qui traduit le code en sortie de processus, ce qui permet aux tests
 * d'appeler {@code run} directement (comme {@code AppParallel}).</p>
 */
public final class EnergyVsField {

    private EnergyVsField() {}

    /** Code de retour : tout s'est bien passe. */
    public static final int EXIT_OK = 0;
    /** Code de retour : echec d'execution (le detail est imprime). */
    public static final int EXIT_FAILURE = 1;
    /** Code de retour : erreur d'usage (arguments invalides). */
    public static final int EXIT_USAGE = 2;

    /** Taille du reseau PAR DEFAUT (l'option --size la remplace). */
    private static final int LA = 12;
    /** Taille du reseau PAR DEFAUT (l'option --size la remplace). */
    private static final int LB = 12;

    /** Direction du champ dans le tenseur g (historique E_v_field). */
    private static final double[] FIELD_DIRECTION = { -1.d, -1.d, 2.d };
    /** Conversion en Tesla : 1 unite de h = 6.91 T (les spins sont normes a 1). */
    private static final double TESLA_PER_UNIT = 6.91;
    /** Garde-fou : nombre maximal de pas de champ (rejet au-dela, aucun clamp). */
    private static final int MAX_STEPS = 10_000;

    // ---- Geometrie des etats imposes. ----
    /** Axe des ondes flop. */
    private static final double[] AXIS_A = { -0.707107, 0.707107, 0. };
    /** AXIS des etats imposes (spins = signe × B_KITAEV, norme 1). */
    private static final double[] B_KITAEV = { -0.4082482904638631, -0.4082482904638631,
            0.816496580927726 };

    // ---- Tables de signes (double zigzag MC-valide de l'utilisateur). ----
    /** Zigzag simple 1/2 : {+,-} sur j, decalage 1 ; moment net nul. */
    private static final int[] SIG_ZIGZAG_12 = { 1, -1 };
    /** uud 1/3 : {+,+,-} sur j, decalage 1 — moment net +1/3 par site (plateau d'aimantation). */
    private static final int[] SIG_UUD_13 = { 1, 1, -1 };
    /**
     * Double zigzag 1/4 : {+,+,-,-} sur j, decalage 1 (ex-{@code PHASE_SIGN} de l'ancien code,
     * sigA du dessin utilisateur) ; moment net nul.
     */
    private static final int[] SIG_DZ_14 = { 1, 1, -1, -1 };
    /** Slots flop scannes : tables uud et double zigzag UNIQUEMENT (periodes 1/3 et 1/4). */
    private static final int[][] FLOP_SIGS = { SIG_UUD_13, SIG_DZ_14 };
    /** Periode associee a chaque slot flop (colonne period de la sortie). */
    private static final double[] FLOP_PERIOD = { 1. / 3., 0.25 };

    // ---- Grille d'angle phi du scan flop : PH_MIN .. PH_MAX par PH_STEP. ----
    private static final double PH_MIN = 0.1;
    private static final double PH_STEP = 0.01;
    private static final double PH_MAX = Math.PI / 4.;

    private static final String USAGE = String.join("\n",
        "Usage :",
        "  java -cp \"bin:lib/*\" EnergyVsField <inputDir> <outputDir> <hInit> <hMax> <dH>",
        "                                      [--field-dir x,y,z] [--size La,Lb]",
        "",
        "Energie d'etats imposes (Ising + flop) en fonction du champ, pour chaque fichier",
        "bcaoExplor* de <inputDir> ; ecrit config_energy_H_<input>.avro et son miroir CSV",
        "(memes valeurs) dans <outputDir>.",
        "",
        "Positionnels :",
        "  inputDir    repertoire contenant les fichiers bcaoExplor_*",
        "  outputDir   repertoire de sortie (cree au besoin)",
        "  hInit       premiere valeur du champ, >= 0 (unites du g-tenseur)",
        "  hMax        borne visee, > hInit (INCLUSE dans la grille, cf. ci-dessous)",
        "  dH          pas du champ, > 0",
        "",
        "Options :",
        "  --field-dir x,y,z  direction du champ dans le tenseur g (defaut -1,-1,2) ;",
        "                     vecteur non nul, normalise par l'outil",
        "  --size La,Lb       taille du reseau La x Lb, entiers >= 1 (defaut 12,12)",
        "",
        "Grille de champ en entier : h_k = hInit + k*dH pour k = 0..nSteps-1, avec",
        "nSteps = floor((hMax - hInit)/dH + 1e-9) + 1 (hMax inclus, meme convention",
        "qu'AppParallel) ; h_nSteps-1 peut depasser hMax d'un ulp (grille voulue).",
        "Plus de " + MAX_STEPS + " pas de champ : rejet (EXIT_USAGE).",
        "",
        "Codes de retour : 0 succes, 1 echec d'execution, 2 erreur d'usage.");

    /**
     * Point d'entree du processus : delegue a {@link #run(String[])} et traduit le code en
     * sortie de processus (y compris 0 : un outil CLI doit toujours fixer son statut).
     */
    public static void main(String[] args) {
        System.exit(run(args));
    }

    /**
     * Execute le balayage sur tous les fichiers bcaoExplor* du repertoire d'entree. Ne leve
     * rien : les erreurs d'usage renvoient {@link #EXIT_USAGE}, les echecs d'execution
     * {@link #EXIT_FAILURE}.
     *
     * @return {@link #EXIT_OK}, {@link #EXIT_FAILURE} ou {@link #EXIT_USAGE}.
     */
    public static int run(String[] args) {
        if (args == null || args.length == 0) {
            usage();
            return EXIT_USAGE;
        }
        // Options n'importe ou sur la ligne, positionnels ensuite (meme regle qu'AppParallel :
        // le dernier lu gagne). Validation immediate des options : avant tout traitement.
        double[] fieldDir = FIELD_DIRECTION.clone();
        int la = LA;
        int lb = LB;
        final List<String> positional = new ArrayList<>();
        for (int i = 0; i < args.length; i++) {
            final String s = args[i];
            if (s.equals("--field-dir")) {
                if (i + 1 >= args.length) {
                    System.err.println("--field-dir attend x,y,z");
                    return EXIT_USAGE;
                }
                fieldDir = parseFieldDir(args[++i]);
                if (fieldDir == null) return EXIT_USAGE; // message deja imprime
            } else if (s.equals("--size")) {
                if (i + 1 >= args.length) {
                    System.err.println("--size attend La,Lb");
                    return EXIT_USAGE;
                }
                // Valide ICI (et non par l'IllegalArgumentException de ReplicaFactory, qui
                // sortirait en EXIT_FAILURE) : taille rejettee avant tout traitement.
                final int[] size = parseSize(args[++i]);
                if (size == null) return EXIT_USAGE; // message deja imprime
                la = size[0];
                lb = size[1];
            } else if (s.startsWith("--")) {
                System.err.println("Option inconnue : " + s);
                usage();
                return EXIT_USAGE;
            } else {
                positional.add(s);
            }
        }
        if (positional.size() != 5) {
            usage();
            return EXIT_USAGE;
        }
        final File inputDir = new File(positional.get(0));
        final File outputDir = new File(positional.get(1));
        final double hInit, hMax, dH;
        try {
            hInit = Double.parseDouble(positional.get(2));
            hMax = Double.parseDouble(positional.get(3));
            dH = Double.parseDouble(positional.get(4));
        } catch (NumberFormatException e) {
            System.err.println("hInit, hMax et dH doivent etre des reels (recus : \""
                    + positional.get(2) + "\", \"" + positional.get(3) + "\", \"" + positional.get(4)
                    + "\")");
            return EXIT_USAGE;
        }
        // (NaN echoue aussi : toute comparaison avec NaN est fausse.)
        if (!(hInit >= 0.d)) {
            System.err.println("hInit doit etre >= 0 (recu " + positional.get(2) + ")");
            return EXIT_USAGE;
        }
        if (!(hMax > hInit)) {
            System.err.println("hMax doit etre > hInit (" + hInit + ") : recu " + positional.get(3));
            return EXIT_USAGE;
        }
        // Rejette aussi +Infini (sinon un seul champ serait scanne silencieusement).
        if (!(dH > 0.d) || Double.isInfinite(dH)) {
            System.err.println("dH doit etre > 0 et fini (recu " + positional.get(4) + ")");
            return EXIT_USAGE;
        }
        // Grille en entier : h_k = hInit + k*dH, hMax inclus (meme convention qu'AppParallel) ;
        // h_nSteps-1 peut depasser hMax d'un ulp. Comparaison en double AVANT le cast : le cast
        // direct saturerait a Integer.MAX_VALUE puis +1 deborderait vers MIN_VALUE, franchissant
        // faussement le garde ci-dessous (rejet par MAX_STEPS, jamais de clamp silencieux).
        final double ratio = Math.floor((hMax - hInit) / dH + 1e-9);
        if (ratio > MAX_STEPS - 1) {
            System.err.println("(hMax - hInit) / dH = " + ratio + " pas de champ, maximum "
                    + MAX_STEPS + " ; verifier dH (rejet, aucun clamp)");
            return EXIT_USAGE;
        }
        final int nSteps = (int) ratio + 1;

        if (!inputDir.isDirectory()) {
            System.err.println("Le dossier specifie n'existe pas ou n'est pas un dossier : "
                    + positional.get(0));
            return EXIT_USAGE;
        }
        final File[] files = inputDir.listFiles();
        if (files == null) {
            System.err.println("Impossible de lister " + positional.get(0));
            return EXIT_USAGE;
        }
        Arrays.sort(files); // ordre de traitement reproductible
        final List<File> inputs = new ArrayList<>();
        for (File f : files) {
            // isFile() : un sous-repertoire nomme bcaoExplor_* casserait Input (cf. AppParallel).
            if (f.isFile() && f.getName().startsWith("bcaoExplor")) inputs.add(f);
        }
        if (inputs.isEmpty()) {
            System.err.println("Aucun fichier bcaoExplor* dans " + positional.get(0) + " : rien a faire.");
            return EXIT_USAGE;
        }

        if (!outputDir.isDirectory() && !outputDir.mkdirs()) {
            System.err.println("Impossible de creer le repertoire de sortie : " + positional.get(1));
            return EXIT_FAILURE;
        }

        try {
            for (File file : inputs) {
                processOne(file, outputDir, hInit, nSteps, dH, fieldDir, la, lb);
            }
        } catch (Exception e) {
            // Echec d'execution (pas d'usage) : message + pile, code 1.
            System.err.println("Echec du balayage : " + e);
            e.printStackTrace(System.err);
            return EXIT_FAILURE;
        }
        return EXIT_OK;
    }

    /**
     * Analyse {@code "x,y,z"} en direction de champ : 3 reels finis, vecteur non nul (sinon
     * message sur {@code stderr} et {@code null} — l'appelant renvoie {@link #EXIT_USAGE}).
     */
    private static double[] parseFieldDir(String spec) {
        final String[] xyz = spec.split(",");
        if (xyz.length != 3) {
            System.err.println("--field-dir attend x,y,z (recu \"" + spec + "\")");
            return null;
        }
        final double[] dir = new double[3];
        try {
            for (int i = 0; i < 3; i++) dir[i] = Double.parseDouble(xyz[i]);
        } catch (NumberFormatException e) {
            System.err.println("--field-dir attend 3 reels (recu \"" + spec + "\")");
            return null;
        }
        if (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] <= 0.d
                || !Double.isFinite(dir[0]) || !Double.isFinite(dir[1]) || !Double.isFinite(dir[2])) {
            System.err.println("--field-dir doit etre un vecteur non nul a composantes finies : \""
                    + spec + "\"");
            return null;
        }
        return dir;
    }

    /**
     * Analyse {@code "La,Lb"} en taille de reseau : 2 entiers &gt;= 1 separes par une virgule
     * (sinon message sur {@code stderr} et {@code null} — l'appelant renvoie
     * {@link #EXIT_USAGE}).
     */
    private static int[] parseSize(String spec) {
        final String[] parts = spec.split(",");
        if (parts.length != 2) {
            System.err.println("--size attend La,Lb : 2 entiers >= 1 (recu \"" + spec + "\")");
            return null;
        }
        final int[] size = new int[2];
        try {
            size[0] = Integer.parseInt(parts[0]);
            size[1] = Integer.parseInt(parts[1]);
        } catch (NumberFormatException e) {
            System.err.println("--size attend La,Lb : 2 entiers >= 1 (recu \"" + spec + "\")");
            return null;
        }
        if (size[0] < 1 || size[1] < 1) {
            System.err.println("--size : La et Lb doivent etre >= 1 (recu \"" + spec + "\")");
            return null;
        }
        return size;
    }

    /** Un fichier d'entree : reseau, puis un Avro et son miroir CSV, un record par pas de champ. */
    private static void processOne(File file, File outputDir, double hInit, int nSteps, double dH,
            double[] fieldDir, int la, int lb) throws IOException, FileNotFoundException {
        final ReplicaFactory.Built built =
                ReplicaFactory.build(new Input(file.getAbsolutePath()), la, lb);
        final Scan scan = new Scan(built.lattice(), new Spin(built.lattice()), fieldDir);

        final Schema schema = createSchema();
        // Noms deterministes (remplacent horodatage + compteur statique de l'original) : relancer
        // l'outil avec les memes arguments ecrase proprement les fichiers precedents.
        final File out = new File(outputDir, "config_energy_H_" + file.getName() + ".avro");
        final File csvOut = new File(outputDir, "config_energy_H_" + file.getName() + ".csv");
        final DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(schema);
        try (DataFileWriter<GenericRecord> writer = new DataFileWriter<>(datumWriter);
                PrintWriter csv = new PrintWriter(csvOut, "UTF-8")) {
            writer.setCodec(CodecFactory.snappyCodec()); // codec historique
            writer.create(schema, out);
            csv.println(CSV_HEADER);
            scan.scan(hInit, nSteps, dH, schema, writer, csv);
        }
        System.out.println("Ecrit : " + out + " et " + csvOut.getName());
    }

    /** En-tete du miroir CSV (meme colonnes que le schema Avro, H en Tesla). */
    private static final String CSV_HEADER = "H_T,E0,E1,E2,E3,E_flop,Phi,period";

    /** Usage imprime sur stderr (l'erreur d'usage ne pollue pas stdout). */
    private static void usage() {
        System.err.println(USAGE);
    }

    /**
     * Schema Avro de l'original (record "Energy", champs float), TOUJOURS dans sa variante
     * flop : l'outil moderne ecrit systematiquement les 8 champs, tous finis (4 etats Ising
     * de tables + le trio flop).
     */
    private static Schema createSchema() {
        final String schemaString = "{" +
                "  \"type\": \"record\"," +
                "  \"name\": \"Energy\"," +
                "  \"fields\": [" +
                "    {\"name\": \"H\", \"type\": \"float\"}," +
                "    {\"name\": \"E0\", \"type\": \"float\"}," +
                "    {\"name\": \"E1\", \"type\": \"float\"}," +
                "    {\"name\": \"E2\", \"type\": \"float\"}," +
                "    {\"name\": \"E3\", \"type\": \"float\"}," +
                "    {\"name\": \"E_flop\", \"type\": \"float\"}," +
                "    {\"name\": \"Phi\", \"type\": \"float\"}," +
                "    {\"name\": \"period\", \"type\": \"float\"}" +
                "  ]" +
                "}";
        return new Schema.Parser().parse(schemaString);
    }

    /**
     * Balayage en champ des energies d'etats imposes sur un reseau deja construit : c'est le
     * corps de l'original {@code E_v_field.compute_Energy}, dont la construction cosinus des
     * etats est remplacee par les tables de signes (double zigzag MC-valide).
     */
    private static final class Scan {
        private final Lattice lattice;
        private final Spin spin;
        private final int Nsr;
        private final int Na;
        private final int Nb;
        private final int size;
        private final double[] fieldDir;

        Scan(Lattice lattice, Spin spin, double[] fieldDir) {
            this.lattice = lattice;
            this.spin = spin;
            this.Nsr = lattice.lengthSr;
            this.Na = lattice.lengthA;
            this.Nb = lattice.lengthB;
            this.size = lattice.size;
            this.fieldDir = fieldDir;
        }

        /**
         * Balaye {@code h_k = hInit + k*dH} (k = 0..nSteps-1) et append un record Avro plus sa
         * ligne CSV miroir (exactement les memes valeurs apres arrondi float) par pas.
         *
         * @throws IOException si l'ecriture Avro echoue (propagee vers {@link #EXIT_FAILURE}).
         */
        void scan(double hInit, int nSteps, double dH, Schema schema,
                DataFileWriter<GenericRecord> writer, PrintWriter csv) throws IOException {
            // Etats Ising independants du champ : construits UNE fois (comme l'original).
            final double[][] isingStates = { ferroState(), tableState(SIG_ZIGZAG_12, 1),
                    tableState(SIG_UUD_13, 1), tableState(SIG_DZ_14, 1) };

            for (int k = 0; k < nSteps; k++) {
                final double h = hInit + k * dH; // index entier : pas de derive d'arrondi
                // Champ dans la base de Kitaev via le tenseur g (cf. App / AppParallel).
                final double[] B = MathOps.scale(h,
                        MathOps.matrixXvector(lattice.getGtensor(), MathOps.normalize(fieldDir)));
                lattice.setInteractionField(B);

                final GenericRecord record = new GenericData.Record(schema);

                // ---- Etats Ising : une energie par table de signes (energie par site).
                final float[] eIsing = new float[4];
                for (int i = 0; i < isingStates.length; i++) {
                    lattice.setSpin1Dlattice(isingStates[i]);
                    eIsing[i] = (float) (spin.getEnergy() / lattice.size);
                    record.put("E" + i, eIsing[i]);
                }

                // ---- Etats flop : RESCAN complet a chaque champ, sur les DEUX slots seulement
                // (uud 1/3, double zigzag 1/4). L'original reinitialisait aussi ph a ph_min pour
                // chaque champ (E_v_field:280-281). Accumulation de ph conservée : PAS d'index
                // entier ici (bit-exact requis).
                double energyMin = Double.MAX_VALUE;
                double bestPh = PH_MIN;
                double bestPeriod = 0.d;
                for (double ph = PH_MIN; ph <= PH_MAX + PH_STEP; ph += PH_STEP) {
                    // Une paire (Sp, Sm) par valeur de ph (elle ne depend que de ph) : refactor
                    // pur, l'original reconstruisait par (ph, periode).
                    final double[][] pair = flopPair(ph);
                    for (int s = 0; s < FLOP_SIGS.length; s++) {
                        lattice.setSpin1Dlattice(flopState(FLOP_SIGS[s], pair[0], pair[1]));
                        final double e = spin.getEnergy(); // NON divisee ici : division apres le min
                        if (e < energyMin) {
                            energyMin = e;
                            bestPh = ph;
                            bestPeriod = FLOP_PERIOD[s];
                        }
                    }
                }

                // Remise en Tesla : 1 unite de h = 6.91 T (les spins sont normes a 1).
                final float hTesla = (float) (h * TESLA_PER_UNIT);
                final float eFlop = (float) (energyMin / lattice.size);
                final float phi = (float) bestPh;
                final float period = (float) bestPeriod;
                record.put("H", hTesla);
                record.put("E_flop", eFlop);
                record.put("Phi", phi);
                record.put("period", period);
                writer.append(record);

                // Miroir CSV : exactement les memes valeurs (apres arrondi float) que l'Avro.
                csv.println(String.format(Locale.US, "%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f,%.8f",
                        (double) hTesla, (double) eIsing[0], (double) eIsing[1], (double) eIsing[2],
                        (double) eIsing[3], (double) eFlop, (double) phi, (double) period));
            }
        }

        /** Index lineaire d'un site : {@code (floorMod(j,Lb)·La + floorMod(i,La))·Nsr + t}. */
        private int siteIdx(int i, int j, int t) {
            return (Math.floorMod(j, Nb) * Na + Math.floorMod(i, Na)) * Nsr + t;
        }

        /**
         * Signe du site (j, t) : {@code sigA[j]} sur le sous-reseau A (t == 0),
         * {@code sigA[j + decalage]} sur B — le decalage 1 du zigzag MC-valide (la rangee B
         * suit la rangee A).
         */
        private static int signAt(final int[] sigA, final int decalage, final int j, final int t) {
            return (t == 0) ? sigA[Math.floorMod(j, sigA.length)]
                    : sigA[Math.floorMod(j + decalage, sigA.length)];
        }

        /**
         * Etat impose par table de signes ({@code sigA}, {@code decalage}) : chaque spin vaut
         * {@code signe × B_KITAEV}. C'est l'etat double zigzag MC-valide de l'utilisateur.
         */
        private double[] tableState(final int[] sigA, final int decalage) {
            final double[] config = new double[size * 3];
            for (int i = 0; i < Na; i++) {
                for (int j = 0; j < Nb; j++) {
                    for (int t = 0; t < Nsr; t++) {
                        final int signe = signAt(sigA, decalage, j, t);
                        final int idx = siteIdx(i, j, t);
                        config[3 * idx] = signe * B_KITAEV[0];
                        config[3 * idx + 1] = signe * B_KITAEV[1];
                        config[3 * idx + 2] = signe * B_KITAEV[2];
                    }
                }
            }
            return config;
        }

        /** Etat FERRO : tous les spins le long de +B_KITAEV (cas degenerate, pas de table). */
        private double[] ferroState() {
            final double[] config = new double[size * 3];
            for (int s = 0; s < size; s++) {
                config[3 * s] = B_KITAEV[0];
                config[3 * s + 1] = B_KITAEV[1];
                config[3 * s + 2] = B_KITAEV[2];
            }
            return config;
        }

        /**
         * Etat flop d'un slot : meme mecanique Sp/Sm que l'original (paire a l'angle ph,
         * {@link #flopPair}), mais le signe +/- vient de la table ({@code sigA}, decalage 1)
         * au lieu d'un cosinus : spin = Sp si signe &gt; 0, Sm sinon.
         */
        private double[] flopState(final int[] sigA, final double[] Sp, final double[] Sm) {
            final double[] config = new double[size * 3];
            for (int i = 0; i < Na; i++) {
                for (int j = 0; j < Nb; j++) {
                    for (int t = 0; t < Nsr; t++) {
                        final double[] v = signAt(sigA, 1, j, t) > 0 ? Sp : Sm;
                        final int idx = siteIdx(i, j, t);
                        config[3 * idx] = v[0];
                        config[3 * idx + 1] = v[1];
                        config[3 * idx + 2] = v[2];
                    }
                }
            }
            return config;
        }

        /**
         * Paire d'ondes flop (Sp, Sm) a l'angle {@code ph} (ex
         * {@code createRepeatedPatternArray_flop}) : Sp croissant-cosinus, Sm oppose.
         */
        private static double[][] flopPair(final double ph) {
            final double th = Math.PI / 3.;

            final double[] Sp = MathOps.add(MathOps.scale(FastMath.sin(th - ph), AXIS_A),
                    MathOps.scale(FastMath.cos(th - ph), B_KITAEV));

            final double[] Sm = MathOps.subtract(MathOps.scale(-FastMath.sin(th + ph), AXIS_A),
                    MathOps.scale(FastMath.cos(th + ph), B_KITAEV));

            return new double[][] { Sp, Sm };
        }
    }
}
