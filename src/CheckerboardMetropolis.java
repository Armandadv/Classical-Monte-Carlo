import java.util.SplittableRandom;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.LockSupport;
import java.util.random.RandomGenerator;
import java.util.random.RandomGeneratorFactory;

/**
 * Balayage Metropolis parallele « damier » : les sites d'une meme classe de couleur (cf.
 * {@link LatticeColoring}) etant deux a deux non adjacents, leurs champs locaux ne dependent
 * que de spins figes pendant la phase — les mettre a jour simultanement est rigoureusement
 * equivalent a un balayage sequentiel et laisse invariante la mesure de Boltzmann (Ren &amp;
 * Orkoulas 2007). L'ordre des classes est retire au hasard a chaque balayage.
 *
 * <p>References :
 * T. Preis, P. Virnau, W. Paul, J. J. Schneider, <i>J. Comput. Phys.</i> <b>228</b>, 4468 (2009) ;
 * M. Weigel, <i>J. Comput. Phys.</i> <b>231</b>, 3064 (2012) ;
 * J. A. Anderson, E. Jankowski, T. L. Grubb, M. Engel, S. C. Glotzer,
 * <i>J. Comput. Phys.</i> <b>254</b>, 27 (2013) ;
 * R. Ren, G. Orkoulas, <i>J. Chem. Phys.</i> <b>126</b>, 211102 (2007).</p>
 *
 * <p><b>Confinee a un thread</b> : toutes les methodes publiques, {@code close()} compris,
 * doivent etre appelees depuis le thread proprietaire. Protocole : champs de phase ecrits puis
 * {@code generation} (volatile) publie ; chaque worker termine par {@code done} (volatile) —
 * une seule classe traitee a la fois, aucune lecture de voisin ne concurrence une ecriture.
 * Un generateur par (couleur, tranche), decoupage independant du nombre de threads :
 * trajectoire identique bit a bit pour tout {@code nThreads}.</p>
 *
 * <p>Les gabarits d'echange (2 = 2-PARAM a.I + b.A, 3 = diagonale, 5 = Kitaev, 9 = pleine) sont
 * aplatis puis <b>groupees par forme</b> a la construction : les liaisons de chaque sous-reseau
 * sont triees par tri stable selon l'ordre canonique FULL -&gt; DEUX_PARAM -&gt; DIAG -&gt; KIT01 -&gt;
 * KIT02 -&gt; KIT12 (table de rangs, jamais l'octet brut), puis compactees groupe dos a dos dans un
 * {@code double[]} par sous-reseau (9 doubles par liaison pleine, 2 pour 2-PARAM {a, b}, 3 pour
 * la diagonale, 4 pour Kitaev {d0, d1, d2, v}) ; les voisins de chaque site ({@code reorderedNb})
 * et le gamma du terme quartique — positionnel sur les liaisons brutes dans {@code Lattice} —
 * sont reordonnes a l'identique (table gamma nulle sans terme quartique). Les noyaux branchent
 * une fois par groupe ({@code switch} previsible), deroulent une boucle serree par forme et
 * lisent le triangle superieur (convention symetrique de {@code Spin.exchangeEnergySym}) — une
 * partie antisymetrique (Dzyaloshinskii-Moriya) serait silencieusement ignoree : la construction
 * la refuse (controle de symetrie). Sur une entree deja triee par forme (la fixture reelle) le
 * tri est l'identite : trajectoire bit-a-bit identique a celle du noyau historique ; sinon
 * l'ordre d'accumulation change (derive ~1e-16 acceptee).</p>
 *
 * <p>Paralleliser a l'interieur d'un balayage n'amortit la barriere (~1-10 us) qu'a partir de
 * ~10^3 mises a jour par thread et par phase ({@code N / (numColors * nThreads) >= 1000}) ; en
 * dessous, preferer plus de repliques independantes avec moins de threads chacune. Mesure sur
 * M-series 8 coeurs (4P + 4E), honeycomb J1J2J3, 6 couleurs : 96x96 donne 1,86x a 2 threads,
 * 3,34x a 4 et seulement 1,82x a 8 — ne pas depasser les coeurs rapides, la barriere attend le
 * thread le plus lent.</p>
 */
public final class CheckerboardMetropolis implements AutoCloseable {

    /** Statistiques d'un balayage. */
    public record SweepStats(long attempted, long accepted) {
        public double acceptanceRate() { return attempted == 0L ? 0.d : (double) accepted / (double) attempted; }
        public double rejectionRate() { return attempted == 0L ? 1.d : 1.d - acceptanceRate(); }
    }

    /** Codes d'operation d'une phase, publies aux workers via {@code phaseOp}. */
    private static final int OP_METROPOLIS = 0;
    private static final int OP_ENERGY = 1;
    private static final int OP_MAGNETIZATION = 2;
    private static final int OP_OVERRELAX = 3;
    private static final int OP_SHUTDOWN = 4;

    /** Tolerance relative du controle de symetrie des matrices d'echange. */
    private static final double SYMMETRY_TOL = 1e-12;

    /** Formes de gabarits reconnues par {@code bondForms} et branchees dans les noyaux. */
    private static final byte F_DEUX_PARAM = 0, F_DIAG = 1, F_KIT01 = 2, F_KIT02 = 3, F_KIT12 = 4,
            F_FULL = 5;

    /**
     * Backoff d'attente des workers : ~SPIN_LIMIT tours de {@link Thread#onSpinWait()} (~30 us),
     * puis {@link LockSupport#parkNanos} a doublement exponentiel de PARK_NANOS (2 us) jusqu'a
     * MAX_PARK_NANOS (1 ms), remis a zero a chaque phase. Calibre par mesure sur Apple silicon :
     * la valeur classique (4000 iterations, scrutation 20 us) effondrait le debit (0,70x le
     * mono-thread en 30x30, 4 threads) ; 30 000 iterations / 2 us passent le meme cas a 2,9x.
     * Ne pas sur-souscrire les coeurs : un worker actif brule du CPU pendant son spin.
     */
    private static final int SPIN_LIMIT = 300_000;
    private static final long PARK_NANOS = 2_000L;
    private static final long MAX_PARK_NANOS = 1_000_000L;

    /** Padding anti-faux-partage : 16 long = 128 octets (taille de ligne de cache Apple/M-series). */
    private static final int LSTRIDE = 16;
    /** Padding anti-faux-partage : 16 double = 128 octets. */
    private static final int DSTRIDE = 16;

    private final Lattice lattice;
    private final LatticeColoring coloring;
    private final int nThreads;
    private final int nChunksPerColor;
    private final int nColors;
    private final int nSites;
    private final int nSublattices;

    /** [sous-reseau] liaisons triees par forme (ordre canonique, cf. javadoc de classe) puis
     *  compactees groupe dos a dos : 9 doubles par liaison pleine, 2 {a, b} pour 2-PARAM,
     *  3 {d0, d1, d2} pour la diagonale, 4 {d0, d1, d2, v} pour Kitaev. */
    private final double[][] compactBySublattice;
    /** [sous-reseau] groupes contigus de liaisons, 4 ints par groupe : {forme, premiere liaison,
     *  nb de liaisons, debut dans le compact} (cf. {@link #compactBySublattice}). */
    private final int[][] groupsBySublattice;
    /** [site][liaison triee] voisins reordonnes comme les liaisons triees du sous-reseau du site
     *  (remplace {@code lattice.getNeighborSites(site)} dans les noyaux chauds). */
    private final int[][] reorderedNb;
    /** [sous-reseau][liaison triee] gamma quartique reordonne (le gamma de {@code Lattice}
     *  est positionnel sur les liaisons brutes, cf. {@code getQuarticGamma}) ; {@code null}
     *  sans terme quartique. */
    private final int[][] gammaReordered;

    private final int[][] chunkStart; // [couleur][tranche] index de debut dans classSites(couleur)
    private final int[][] chunkEnd;
    private final int[] globalStart;  // [tranche] decoupage de 0..N-1 (energie, aimantation)
    private final int[] globalEnd;

    private final RandomGenerator[][] rng; // [couleur][tranche]
    /** Generateur maitre : tirage de l'ordre des couleurs ({@link #shuffleColors}). */
    private final RandomGenerator master;
    private final int[] colorOrder;

    private final long[] counters;  // 2 compteurs (tentatives, acceptations) par tranche, paddes
    private final double[] partials; // 3 accumulateurs par tranche, paddes

    /** Workers persistants (le thread appelant est le worker 0) ; daemons arretes par {@link #close()}. */
    private final Thread[] workers;
    /** Acquittements des workers de la phase courante (cible : {@code workers.length}). */
    private final AtomicInteger done = new AtomicInteger();
    /** Numero de phase ; l'ecriture volatile publie les champs de phase (cf. javadoc de classe). */
    private volatile int generation;
    /** Verrou dedie a la seule ecriture de {@link #workerFailure} (premier arrive, premier garde). */
    private final Object failureLock = new Object();
    private Throwable workerFailure; // garde par failureLock
    /** Demande d'arret hors protocole (close pendant une phase est interdit) : sortie des workers. */
    private volatile boolean shutdown;
    private boolean closed;

    /** Cache de {@code transpose(inverse(Kitaev_Basis))}, calcule paresseusement (cf. magnetization()). */
    private double[][] magnetizationBasis;

    // Champs de phase : ecrits avant l'ecriture volatile de generation, lus apres sa lecture.
    private int phaseOp;
    private int phaseColor;
    private double phaseTemperature;
    private double phaseSigma;
    private double[] spins;
    private double fieldX, fieldY, fieldZ;

    /** Raccourci : {@code nChunksPerColor = nThreads}. */
    public CheckerboardMetropolis(Lattice lattice, LatticeColoring coloring, int nThreads, long seed) {
        this(lattice, coloring, nThreads, nThreads, seed);
    }

    /**
     * Construction complete : valide les arguments, fige les tables chaudes (liaisons triees
     * par forme, voisins et gamma reordonnes) et demarre les {@code nThreads - 1} workers
     * persistants.
     *
     * @param lattice reseau porte par le sweeper.
     * @param coloring coloration propre <b>du meme reseau</b> (verifiee : meme instance de
     *        {@link Lattice}, classes partitionnant exactement {@code lattice.size}).
     * @param nThreads nombre de threads de mise a jour, appelant compris.
     * @param nChunksPerColor tranches par classe de couleur (aussi utilisee pour le decoupage
     *        global energie/aimantation).
     * @param seed graine du tirage des generateurs, un par (couleur, tranche).
     * @throws IllegalArgumentException arguments invalides, coloration incoherente, ou matrice
     *         d'echange non symetrique — les noyaux lisent le triangle superieur (convention de
     *         {@code Spin.exchangeEnergySym}) : une partie antisymetrique (Dzyaloshinskii-Moriya)
     *         serait silencieusement ignoree. L'etat du reseau (gabarits, terme quartique) doit
     *         etre fige AVANT la construction : les tables internes (liaisons triees, gamma
     *         quartique reordonne) sont construites ici et ne suivent pas les mutations ulterieures.
     */
    public CheckerboardMetropolis(Lattice lattice, LatticeColoring coloring,
                                  int nThreads, int nChunksPerColor, long seed) {
        if (lattice == null) throw new IllegalArgumentException("lattice == null");
        if (coloring == null) throw new IllegalArgumentException("coloring == null");
        if (nThreads < 1) throw new IllegalArgumentException("nThreads doit etre >= 1 : " + nThreads);
        if (nChunksPerColor < 1) {
            throw new IllegalArgumentException("nChunksPerColor doit etre >= 1 : " + nChunksPerColor);
        }
        if (coloring.lattice() != lattice) {
            throw new IllegalArgumentException(
                    "la coloration ne porte pas sur ce reseau (LatticeColoring.lattice() != lattice)");
        }
        long colored = 0L;
        for (int[] cl : coloring.classes()) colored += cl.length;
        if (colored != lattice.size) {
            throw new IllegalArgumentException("la coloration couvre " + colored
                    + " sites, le reseau en compte " + lattice.size);
        }
        this.lattice = lattice;
        this.coloring = coloring;
        this.nThreads = nThreads;
        this.nChunksPerColor = nChunksPerColor;
        this.nColors = coloring.numColors();
        this.nSites = lattice.size;
        this.nSublattices = lattice.lengthSr;
        // Aplatit et verifie la symetrie (refus a la construction, cf. checkSymmetry) ; le
        // resultat est un temporaire : il alimente le controle de reciprocite puis sert de
        // source au compact groupe — il n'est PAS retenu en champ (le stockage chaud est
        // compactBySublattice, groupe par forme, cf. buildGroupedBondTables).
        final double[][] flat = flattenTemplates(lattice);
        final byte[][] forms = bondForms(lattice);
        rejectSelfLoops(lattice);
        warnIfNonReciprocal(lattice, flat);
        final GroupedBonds grouped = buildGroupedBondTables(lattice, flat, forms);
        this.compactBySublattice = grouped.compact();
        this.groupsBySublattice = grouped.groups();
        this.reorderedNb = grouped.nb();
        this.gammaReordered = grouped.gamma();

        this.chunkStart = new int[nColors][nChunksPerColor];
        this.chunkEnd = new int[nColors][nChunksPerColor];
        for (int c = 0; c < nColors; c++) {
            split(coloring.classSites(c).length, nChunksPerColor, chunkStart[c], chunkEnd[c]);
        }
        this.globalStart = new int[nChunksPerColor];
        this.globalEnd = new int[nChunksPerColor];
        split(nSites, nChunksPerColor, globalStart, globalEnd);

        // Un generateur par (couleur, tranche), graine par SplittableRandom : etat independant
        // du nombre de threads.
        RandomGeneratorFactory<RandomGenerator> factory =
                RandomGeneratorFactory.of("Xoshiro256PlusPlus");
        SplittableRandom seeder = new SplittableRandom(seed);
        this.rng = new RandomGenerator[nColors][nChunksPerColor];
        for (int c = 0; c < nColors; c++) {
            for (int k = 0; k < nChunksPerColor; k++) rng[c][k] = factory.create(seeder.nextLong());
        }
        this.master = factory.create(seeder.nextLong());
        this.colorOrder = new int[nColors];
        for (int c = 0; c < nColors; c++) colorOrder[c] = c;

        this.counters = new long[nChunksPerColor * LSTRIDE];
        this.partials = new double[nChunksPerColor * DSTRIDE];

        this.workers = new Thread[nThreads - 1];
        for (int w = 0; w < workers.length; w++) {
            final int id = w + 1; // le thread appelant est le worker 0
            Thread t = new Thread(() -> workerLoop(id), "checkerboard-" + id);
            t.setDaemon(true);
            workers[w] = t;
            t.start();
        }
    }

    /**
     * Aplatit les gabarits d'echange en un {@code double[9 * z]} contigu par sous-reseau :
     * 2-PARAM {@code {a, b}} -> a.I + b.A, diagonale -> matrice pleine, Kitaev -> diag + sa
     * paire k, pleine -> copie. Deduplique au passage (un tableau par sous-reseau au lieu d'une
     * indirection par voisin). Chaque matrice est verifiee symetrique (cf. {@link #checkSymmetry}).
     * Le resultat n'est pas retenu en champ : temporaire de construction, il alimente le controle
     * de reciprocite puis sert de source au compact groupe ({@link #buildGroupedBondTables}).
     */
    static double[][] flattenTemplates(Lattice lattice) {
        final int lsr = lattice.lengthSr;
        final double[][] out = new double[lsr][];
        for (int s = 0; s < lsr; s++) {
            // Le site d'indice s appartient au sous-reseau s : idx = (j*La + i)*Lsr + s avec i=j=0.
            final double[][] mats = lattice.getInteractionMatrices(s);
            final int z = mats.length;
            final double[] flat = new double[9 * z];
            for (int n = 0; n < z; n++) {
                final double[] m = mats[n];
                final int b = 9 * n;
                if (m.length == 2) {
                    // 2-PARAM a.I + b.A — a sur la diagonale, b hors-diagonale
                    flat[b] = m[0];
                    flat[b + 4] = m[0];
                    flat[b + 8] = m[0];
                    flat[b + 1] = m[1];
                    flat[b + 3] = m[1];
                    flat[b + 5] = m[1];
                    flat[b + 7] = m[1];
                } else if (m.length == 3) {
                    flat[b] = m[0];
                    flat[b + 4] = m[1];
                    flat[b + 8] = m[2];
                } else if (m.length == 5) {
                    // Kitaev : diag(d0, d1, d2) + paire k de valeur v
                    flat[b] = m[0];
                    flat[b + 4] = m[1];
                    flat[b + 8] = m[2];
                    final int k = (int) m[3];
                    final double v = m[4];
                    if (k == 0) {
                        flat[b + 1] = v;
                        flat[b + 3] = v;
                    } else if (k == 1) {
                        flat[b + 2] = v;
                        flat[b + 6] = v;
                    } else {
                        flat[b + 5] = v;
                        flat[b + 7] = v;
                    }
                } else if (m.length == 9) {
                    System.arraycopy(m, 0, flat, b, 9);
                } else {
                    throw new IllegalArgumentException("sous-reseau " + s + ", voisin " + n
                            + " : gabarit de longueur " + m.length + " (2, 3, 5 ou 9 attendus)");
                }
                checkSymmetry(flat, b, s, n);
            }
            out[s] = flat;
        }
        return out;
    }

    /**
     * Forme de chaque gabarit, derivee de sa longueur : 2 = 2-PARAM, 3 = diagonale, 5 = Kitaev,
     * 9 = pleine. Classification de construction uniquement : elle alimente le tri par forme des
     * liaisons ({@link #buildGroupedBondTables}) ; les noyaux chauds ne branchent plus forme par
     * forme mais une fois par groupe pre-trie (un {@code switch} par groupe).
     */
    private static byte[][] bondForms(Lattice lattice) {
        final int lsr = lattice.lengthSr;
        final byte[][] out = new byte[lsr][];
        for (int s = 0; s < lsr; s++) {
            final double[][] mats = lattice.getInteractionMatrices(s);
            out[s] = new byte[mats.length];
            for (int n = 0; n < mats.length; n++) {
                final double[] m = mats[n];
                if (m.length == 2) out[s][n] = F_DEUX_PARAM;
                else if (m.length == 3) out[s][n] = F_DIAG;
                else if (m.length == 5) {
                    final int k = (int) m[3];
                    if (k < 0 || k > 2) throw new IllegalArgumentException("sous-reseau " + s
                            + ", voisin " + n + " : indice Kitaev k = " + k + " hors {0,1,2}");
                    out[s][n] = (byte) (F_KIT01 + k);
                }
                else if (m.length == 9) out[s][n] = F_FULL;
                else throw new IllegalArgumentException("sous-reseau " + s + ", voisin " + n
                        + " : gabarit de longueur " + m.length + " inconnue");
            }
        }
        return out;
    }

    /** Rang canonique d'une forme pour le tri des liaisons : FULL -&gt; DEUX_PARAM -&gt; DIAG -&gt;
     *  KIT01 -&gt; KIT02 -&gt; KIT12. Volontairement PAS l'ordre des octets (dont la numeration
     *  ne suit pas cette convention de layout) : le tri passe toujours par cette table de rangs.
     */
    private static int rankOf(byte form) {
        return switch (form) {
            case F_FULL -> 0;
            case F_DEUX_PARAM -> 1;
            case F_DIAG -> 2;
            case F_KIT01 -> 3;
            case F_KIT02 -> 4;
            default -> 5; // F_KIT12
        };
    }

    /** Nb de doubles de compact par liaison de la forme donnee : 9 (pleine), 2 (2-PARAM),
     *  3 (diagonale), 4 (Kitaev). */
    private static int compactStride(byte form) {
        return switch (form) {
            case F_FULL -> 9;
            case F_DEUX_PARAM -> 2;
            case F_DIAG -> 3;
            default -> 4; // F_KIT01, F_KIT02, F_KIT12
        };
    }

    /** Tables groupees produites une fois a la construction (cf. {@link #buildGroupedBondTables}) ;
     *  {@code gamma} est {@code null} sans terme quartique. */
    private record GroupedBonds(double[][] compact, int[][] groups, int[][] nb, int[][] gamma) {}

    /**
     * Construit les tables chaudes groupees par forme de liaison (cf. javadoc de classe) :
     * <ul>
     *   <li>{@code compactBySublattice} / {@code groupsBySublattice} : liaisons triees par ordre
     *       canonique de forme (table de rangs via {@link #rankOf}, tri stable — les liaisons
     *       d'une meme forme gardent leur ordre brut relatif), puis compactees groupe dos a dos
     *       (FULL 9 valeurs copiees ; 2-PARAM {flat[9n], flat[9n+1]} = {a, b} ; diagonale
     *       {flat[9n], flat[9n+4], flat[9n+8]} ; Kitaev {d0, d1, d2, v} avec v = flat[9n+1]
     *       (KIT01), flat[9n+2] (KIT02) ou flat[9n+5] (KIT12)) ;</li>
     *   <li>{@code gammaReordered} : gamma quartique reordonne dans l'ordre trie — le gamma de
     *       {@link Lattice#getQuarticGamma(int, int)} est positionnel sur les liaisons BRUTES du
     *       sous-reseau, il faut donc le reordonner en meme temps que les liaisons, sinon un
     *       fichier d'entree non deja trie biaiserait l'echantillonnage quartique ;</li>
     *   <li>{@code reorderedNb} : voisins de chaque site reordonnes comme les liaisons triees de
     *       son sous-reseau (l'indice m des groupes, du compact et du gamma designe le meme
     *       voisin dans tous les noyaux).</li>
     * </ul>
     * Sur une entree deja triee par forme (la fixture reelle : 3 liaisons J1 pleines en tete,
     * puis les 15 liaisons 2-PARAM), l'ordre trie est l'identite et les sommes s'accumulent dans
     * le meme ordre qu'avant : trajectoire bit-a-bit identique au noyau historique.
     */
    private static GroupedBonds buildGroupedBondTables(Lattice lattice, double[][] flat, byte[][] forms) {
        final int lsr = lattice.lengthSr;
        final boolean quartic = lattice.hasQuartic();
        final int[][] orders = new int[lsr][];
        final double[][] compactBySublattice = new double[lsr][];
        final int[][] groupsBySublattice = new int[lsr][];
        final int[][] gammaBySublattice = quartic ? new int[lsr][] : null;

        for (int sr = 0; sr < lsr; sr++) {
            final byte[] f = forms[sr];
            final int z = f.length;
            final double[] fs = flat[sr];

            // Rang canonique de chaque liaison brute.
            final int[] ranks = new int[z];
            for (int n = 0; n < z; n++) ranks[n] = rankOf(f[n]);

            // Tri stable par rang (insertion) : les liaisons d'une meme forme gardent leur
            // ordre brut relatif.
            final int[] order = new int[z];
            for (int n = 0; n < z; n++) {
                final int r = ranks[n];
                int m = n;
                while (m > 0 && ranks[order[m - 1]] > r) { order[m] = order[m - 1]; m--; }
                order[m] = n;
            }
            orders[sr] = order;

            // Groupes contigus apres tri + taille du compact.
            int nGroups = z == 0 ? 0 : 1;
            int compactSize = 0;
            for (int m = 0; m < z; m++) {
                compactSize += compactStride(f[order[m]]);
                if (m > 0 && ranks[order[m]] != ranks[order[m - 1]]) nGroups++;
            }
            final double[] compact = new double[compactSize];
            final int[] groups = new int[4 * nGroups];
            final int[] gamma = quartic ? new int[z] : null;

            int gi = 0, cb = 0, m = 0;
            while (m < z) {
                final int start = m;
                final int rank = ranks[order[m]];
                final byte form = f[order[m]];
                int count = 1;
                while (m + count < z && ranks[order[m + count]] == rank) count++;
                groups[gi] = form;
                groups[gi + 1] = start;
                groups[gi + 2] = count;
                groups[gi + 3] = cb; // debut du groupe dans le compact
                gi += 4;
                for (int k = 0; k < count; k++, m++) {
                    final int n = order[m];
                    final int b9 = 9 * n;
                    switch (form) {
                        case F_FULL -> {
                            System.arraycopy(fs, b9, compact, cb, 9);
                            cb += 9;
                        }
                        case F_DEUX_PARAM -> {
                            compact[cb] = fs[b9];         // a
                            compact[cb + 1] = fs[b9 + 1]; // b
                            cb += 2;
                        }
                        case F_DIAG -> {
                            compact[cb] = fs[b9];
                            compact[cb + 1] = fs[b9 + 4];
                            compact[cb + 2] = fs[b9 + 8];
                            cb += 3;
                        }
                        case F_KIT01 -> {
                            compact[cb] = fs[b9];
                            compact[cb + 1] = fs[b9 + 4];
                            compact[cb + 2] = fs[b9 + 8];
                            compact[cb + 3] = fs[b9 + 1]; // v
                            cb += 4;
                        }
                        case F_KIT02 -> {
                            compact[cb] = fs[b9];
                            compact[cb + 1] = fs[b9 + 4];
                            compact[cb + 2] = fs[b9 + 8];
                            compact[cb + 3] = fs[b9 + 2]; // v
                            cb += 4;
                        }
                        default -> { // F_KIT12
                            compact[cb] = fs[b9];
                            compact[cb + 1] = fs[b9 + 4];
                            compact[cb + 2] = fs[b9 + 8];
                            compact[cb + 3] = fs[b9 + 5]; // v
                            cb += 4;
                        }
                    }
                    if (gamma != null) gamma[m] = lattice.getQuarticGamma(sr, n);
                }
            }
            compactBySublattice[sr] = compact;
            groupsBySublattice[sr] = groups;
            if (gamma != null) gammaBySublattice[sr] = gamma;
        }

        // Voisins reordonnes pour TOUT site : reorderedNb[site][m] = voisin brut d'indice
        // order[site % lsr][m]. Remplace lattice.getNeighborSites(site) dans les trois chunks.
        final int[][] reorderedNb = new int[lattice.size][];
        for (int site = 0; site < lattice.size; site++) {
            final int[] nb = lattice.getNeighborSites(site);
            final int[] order = orders[site % lsr];
            final int z = order.length;
            final int[] rn = new int[z];
            for (int m = 0; m < z; m++) rn[m] = nb[order[m]];
            reorderedNb[site] = rn;
        }
        return new GroupedBonds(compactBySublattice, groupsBySublattice, reorderedNb, gammaBySublattice);
    }

    /**
     * Refuse un site voisin de lui-meme (reseau de taille 1 avec un couplage intra sous-reseau) :
     * le champ local contiendrait le spin du site lui-meme, non representable par une coloration.
     */
    private static void rejectSelfLoops(Lattice lattice) {
        for (int site = 0; site < lattice.size; site++) {
            final int[] nb = lattice.getNeighborSites(site);
            if (nb == null) continue;
            for (int j : nb) {
                if (j == site) {
                    throw new IllegalArgumentException("le site " + site + " est son propre voisin "
                            + "(reseau trop petit pour la portee des interactions)");
                }
            }
        }
    }

    /**
     * Verifie la reciprocite M_ij == M_ji^T (sinon l'energie n'est pas un hamiltonien bien defini
     * et la sur-relaxation ne la conserve plus exactement). Cas herite des g-tenseurs differents
     * par sous-reseau : annonce avec l'ecart maximal, sans refus.
     */
    private static void warnIfNonReciprocal(Lattice lattice, double[][] flat) {
        final int lsr = lattice.lengthSr;
        double worst = 0.d;
        int worstSite = -1, worstNb = -1;
        for (int site = 0; site < lattice.size; site++) {
            final int[] nb = lattice.getNeighborSites(site);
            if (nb == null) continue;
            final double[] fi = flat[site % lsr];
            for (int n = 0; n < nb.length; n++) {
                final int j = nb[n];
                if (j < site) continue;               // chaque paire une fois
                final int[] nbj = lattice.getNeighborSites(j);
                if (nbj == null) continue;
                int back = -1;
                for (int k = 0; k < nbj.length; k++) if (nbj[k] == site) { back = k; break; }
                if (back < 0) continue;               // lien non symetrise : deja hors hypothese
                final double[] fj = flat[j % lsr];
                double scale = 1.d, diff = 0.d;
                for (int r = 0; r < 3; r++) {
                    for (int c = 0; c < 3; c++) {
                        final double mij = fi[9 * n + 3 * r + c];
                        final double mji = fj[9 * back + 3 * c + r];   // transposee
                        scale = Math.max(scale, Math.abs(mij));
                        diff = Math.max(diff, Math.abs(mij - mji));
                    }
                }
                final double rel = diff / scale;
                if (rel > worst) { worst = rel; worstSite = site; worstNb = j; }
            }
        }
        if (worst > SYMMETRY_TOL) {
            System.err.printf(java.util.Locale.US,
                    "AVERTISSEMENT CheckerboardMetropolis : liens non reciproques (M_ij != M_ji^T), "
                    + "ecart relatif max %.3e sur le lien %d-%d ; l'energie n'est pas un hamiltonien "
                    + "bien defini et la sur-relaxation ne la conserve pas exactement.%n",
                    worst, worstSite, worstNb);
        }
    }

    /** Controle {@code |M - M^T| <= 1e-12 * (1 + max|M|)} sur la matrice pleine ecrite en {@code b}. */
    private static void checkSymmetry(double[] flat, int b, int sublattice, int neighbor) {
        double maxAbs = 0.d;
        for (int k = 0; k < 9; k++) maxAbs = Math.max(maxAbs, Math.abs(flat[b + k]));
        final double tol = SYMMETRY_TOL * (1.d + maxAbs);
        final double d10 = Math.abs(flat[b + 3] - flat[b + 1]);
        final double d20 = Math.abs(flat[b + 6] - flat[b + 2]);
        final double d21 = Math.abs(flat[b + 7] - flat[b + 5]);
        if (d10 > tol || d20 > tol || d21 > tol) {
            throw new IllegalArgumentException(
                    "matrice d'echange non symetrique (sous-reseau " + sublattice + ", voisin "
                    + neighbor + ") : |m10-m01|=" + d10 + ", |m20-m02|=" + d20 + ", |m21-m12|=" + d21
                    + " > tolerance " + tol + ". La convention symetrique de CheckerboardMetropolis "
                    + "(comme Spin.exchangeEnergySym) ignorerait silencieusement la partie "
                    + "antisymetrique (terme de Dzyaloshinskii-Moriya).");
        }
    }

    /** Decoupe [0, n) en {@code parts} tranches contigues de tailles quasi egales. */
    private static void split(int n, int parts, int[] start, int[] end) {
        int base = n / parts, rem = n % parts, cursor = 0;
        for (int k = 0; k < parts; k++) {
            int len = base + (k < rem ? 1 : 0);
            start[k] = cursor;
            cursor += len;
            end[k] = cursor;
        }
    }

    // --------------------------------------------------------------- API publique

    /**
     * Un balayage Metropolis complet : chaque classe de couleur est mise a jour une fois, dans un
     * ordre retire au hasard (Fisher-Yates avec le generateur maitre).
     */
    public SweepStats sweep(double temperature, double sigma) {
        checkOpen();
        beginPhaseState();
        this.phaseTemperature = temperature;
        this.phaseSigma = sigma;
        java.util.Arrays.fill(counters, 0L);
        shuffleColors();
        for (int ci = 0; ci < nColors; ci++) {
            runPhase(OP_METROPOLIS, colorOrder[ci]);
        }
        long attempted = 0L, accepted = 0L;
        for (int k = 0; k < nChunksPerColor; k++) {
            attempted += counters[k * LSTRIDE];
            accepted += counters[k * LSTRIDE + 1];
        }
        return new SweepStats(attempted, accepted);
    }

    /**
     * Sur-relaxation microcanonique (Creutz 1987 ; Brown &amp; Woch 1987) : {@code S <- 2 (S.h^) h^ - S}.
     * L'energie locale
     * {@code -S.h} est exactement conservee (les voisins, d'une autre couleur, sont figes) — ne
     * thermalise pas : a entrelacer avec des balayages Metropolis. Requiert des matrices
     * symetriques ; refuse si le terme quartique est actif (la reflexion ne conserve que le terme
     * quadratique).
     */
    public void overRelaxationSweep() {
        if (lattice.hasQuartic()) {
            throw new IllegalStateException("sur-relaxation incompatible avec le terme quartique : "
                    + "la reflexion S <- 2(S.h^)h^ - S ne conserve que l'energie quadratique, "
                    + "l'echantillonnage serait biaise. Utiliser --no-overrelax (ou --overrelax 0) "
                    + "avec une ligne QUARTIC dans le fichier d'entree.");
        }
        checkOpen();
        beginPhaseState();
        for (int c = 0; c < nColors; c++) runPhase(OP_OVERRELAX, c);
    }

    /**
     * Energie totale, meme convention que {@code Spin.getEnergy()} :
     * {@code E = sum_i [ -S_i . (sum_{j in nb(i), j < i} M^sym S_j) - S_i . H ]}
     * (chaque paire comptee une seule fois).
     */
    public double energy() {
        checkOpen();
        beginPhaseState();
        runPhase(OP_ENERGY, 0);
        double e = 0.d;
        for (int k = 0; k < nChunksPerColor; k++) e += partials[k * DSTRIDE];
        if (lattice.hasQuartic()) {
            // Terme quartique NN, meme convention que Spin.getQuarticEnergy() : somme
            // sequentielle hors chunks. b = 0 -> sautee.
            final double[] S = spins;
            final double b = lattice.getQuarticB();
            for (int site = 0; site < lattice.size; site++) {
                final int[] nb = lattice.getNeighborSites(site);
                final int off = site * 3;
                for (int i = 0; i < nb.length; i++) {
                    if (site > nb[i]) {
                        final int g = lattice.getQuarticGamma(site, i);
                        final double p = S[off + g] * S[nb[i] * 3 + g];
                        e -= b * p * p;
                    }
                }
            }
        }
        return e;
    }

    /**
     * Aimantation totale : somme des spins en base de Kitaev, convertie en base a*bc par
     * {@code transpose(inverse(lattice.getKitaevBasis()))}. Ce changement de base est calcule
     * une seule fois, paresseusement au premier appel : la base de Kitaev est supposee fixee
     * apres la construction du sweeper — la modifier ensuite ({@code lattice.setKitaevBasis(...)})
     * n'aurait aucun effet sur cette methode.
     */
    public double[] magnetization() {
        checkOpen();
        beginPhaseState();
        runPhase(OP_MAGNETIZATION, 0);
        double mx = 0.d, my = 0.d, mz = 0.d;
        for (int k = 0; k < nChunksPerColor; k++) {
            mx += partials[k * DSTRIDE];
            my += partials[k * DSTRIDE + 1];
            mz += partials[k * DSTRIDE + 2];
        }
        if (magnetizationBasis == null) {
            magnetizationBasis = MathOps.transpose(MathOps.inverse(lattice.getKitaevBasis()));
        }
        return MathOps.matrixXvector(magnetizationBasis, new double[] { mx, my, mz });
    }

    /** Norme de {@link #magnetization()}. */
    public double magnetizationNorm() {
        return MathOps.norm(magnetization());
    }

    /**
     * Adaptation de sigma, regle historique du sequentiel : {@code sigma *= 0.5 / rejection},
     * borne a [1e-4, 60]. Un rejet nul ou NaN donne 60 (legacy : 0.5/0 = +Inf, borne haute).
     */
    public static double nextSigma(double sigma, double rejectionRate) {
        if (rejectionRate <= 0.d || Double.isNaN(rejectionRate)) return 60.d; // legacy : 0.5/0 = +Inf
        double sig = sigma * (0.5d / rejectionRate);
        if (sig >= 60.d) sig = 60.d;
        if (sig <= 0.0001d) sig = 0.0001d;
        return sig;
    }

    /** Nombre de threads de mise a jour, appelant compris. */
    public int nThreads() { return nThreads; }

    /** Nombre de tranches par classe de couleur. */
    public int nChunksPerColor() { return nChunksPerColor; }

    /** Coloration portee par ce sweeper, figee a la construction. */
    public LatticeColoring coloring() { return coloring; }

    /** Reseau porte par ce sweeper (gabarits aplaties a la construction : non transferable). */
    public Lattice lattice() { return lattice; }

    /**
     * Arrete les threads de service. Idempotent ; a appeler depuis le thread proprietaire, jamais
     * pendant une phase. Ne leve jamais : un worker qui ne finit pas en 2 s est interrompu et
     * signale sur System.err (threads daemons).
     */
    @Override
    public synchronized void close() {
        if (closed) return;
        closed = true;
        if (workers.length == 0) return;
        phaseOp = OP_SHUTDOWN;
        shutdown = true;      // ecriture volatile : reveille aussi les workers hors protocole
        done.set(0);
        generation = generation + 1;
        for (Thread t : workers) {
            try {
                t.join(2000L);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
            }
            if (t.isAlive()) {
                t.interrupt();
                System.err.println("CheckerboardMetropolis.close() : le thread " + t.getName()
                        + " ne s'est pas termine en 2 s, interruption demandee (thread daemon, "
                        + "la JVM peut s'arreter malgre tout).");
            }
        }
    }

    /** Garde d'utilisation apres {@link #close()}. */
    private void checkOpen() {
        if (closed) throw new IllegalStateException("CheckerboardMetropolis deja ferme");
    }

    /** Relit les references susceptibles d'avoir ete remplacees entre deux appels. */
    private void beginPhaseState() {
        this.spins = lattice.getSpin1Dlattice();
        double[] h = lattice.getInteractionField();
        this.fieldX = h[0];
        this.fieldY = h[1];
        this.fieldZ = h[2];
    }

    /** Fisher-Yates de {@code colorOrder} avec le generateur maitre : ordre des couleurs au hasard. */
    private void shuffleColors() {
        for (int i = nColors - 1; i > 0; i--) {
            int j = master.nextInt(i + 1);
            int tmp = colorOrder[i];
            colorOrder[i] = colorOrder[j];
            colorOrder[j] = tmp;
        }
    }

    // --------------------------------------------------------------- orchestration

    /**
     * Publie une phase, execute la part du thread appelant, attend les workers. L'attente a lieu
     * meme en cas d'echec local (sinon {@code done} serait fausse pour la phase suivante) ;
     * l'exception d'origine est relancee apres la barriere, les echecs de workers ajoutes en
     * {@code suppressed}.
     */
    private void runPhase(int op, int color) {
        if (workers.length == 0) {
            phaseOp = op;
            phaseColor = color;
            runChunks(op, color, 0);
            return;
        }
        phaseOp = op;
        phaseColor = color;
        done.set(0);
        generation = generation + 1; // ecriture volatile : publie tous les champs de phase
        Throwable mainFailure = null;
        try {
            runChunks(op, color, 0);
        } catch (Throwable t) {
            mainFailure = t;
        }
        awaitWorkers(mainFailure); // barriere systematique : protocole toujours laisse coherent
    }

    /** Le thread {@code workerId} traite les tranches workerId, workerId + nThreads, ... */
    private void runChunks(int op, int color, int workerId) {
        for (int k = workerId; k < nChunksPerColor; k += nThreads) {
            switch (op) {
                case OP_METROPOLIS -> metropolisChunk(color, k);
                case OP_OVERRELAX -> overRelaxChunk(color, k);
                case OP_ENERGY -> energyChunk(k);
                case OP_MAGNETIZATION -> magnetizationChunk(k);
                default -> { }
            }
        }
    }

    /**
     * Barriere de fin de phase. {@code shutdown} n'est teste que dans la branche lente : un
     * {@code close()} concurrent (usage interdit) leve {@code IllegalStateException} au lieu de
     * bloquer, sans ajouter une lecture volatile a la boucle la plus chaude du protocole.
     */
    private void awaitWorkers(Throwable mainFailure) {
        final int target = workers.length;
        int spins = 0;
        long park = PARK_NANOS;
        while (done.get() != target) {
            if (++spins < SPIN_LIMIT) {
                Thread.onSpinWait();
            } else {
                if (shutdown) throw new IllegalStateException("ferme pendant une phase");
                LockSupport.parkNanos(park);
                park = Math.min(park << 1, MAX_PARK_NANOS); // backoff exponentiel
            }
        }
        Throwable workerFail;
        synchronized (failureLock) {
            workerFail = workerFailure;
            workerFailure = null;
        }
        if (mainFailure != null) {
            if (workerFail != null && workerFail != mainFailure) mainFailure.addSuppressed(workerFail);
            throw rethrow(mainFailure); // on relance l'exception d'origine, pas une seconde
        }
        if (workerFail != null) {
            throw new IllegalStateException("echec d'un thread checkerboard", workerFail);
        }
    }

    /** Enregistre l'echec d'un worker ; le premier arrive est conserve, les suivants en {@code suppressed}. */
    private void recordWorkerFailure(Throwable t) {
        synchronized (failureLock) {
            if (workerFailure == null) {
                workerFailure = t;
            } else if (workerFailure != t) {
                workerFailure.addSuppressed(t);
            }
        }
    }

    /** Relance telle quelle une RuntimeException ou une Error, emballe le reste. */
    private static RuntimeException rethrow(Throwable t) {
        if (t instanceof RuntimeException re) return re;
        if (t instanceof Error e) throw e;
        return new IllegalStateException(t);
    }

    /**
     * Boucle d'un worker persistant : attend la phase suivante en scrutant {@code generation}
     * (spin puis park exponentiel), execute ses tranches, acquitte via {@code done}. Sort sur
     * {@code shutdown} ou {@code phaseOp == OP_SHUTDOWN}.
     */
    private void workerLoop(int id) {
        int localGen = 0;
        for (;;) {
            localGen++;
            int spins = 0;
            long park = PARK_NANOS; // remis a zero a chaque nouvelle phase
            while (generation - localGen < 0) { // comparaison robuste au debordement d'int
                if (++spins < SPIN_LIMIT) {
                    Thread.onSpinWait();
                } else {
                    LockSupport.parkNanos(park);
                    park = Math.min(park << 1, MAX_PARK_NANOS); // backoff exponentiel
                    if (shutdown) return; // teste a chaque reveil, independamment de phaseOp
                }
            }
            if (shutdown) return; // teste a chaque reveil, independamment de phaseOp
            int op = phaseOp; // lu apres la lecture volatile de generation
            if (op == OP_SHUTDOWN) return;
            try {
                runChunks(op, phaseColor, id);
            } catch (Throwable t) {
                recordWorkerFailure(t);
            } finally {
                done.incrementAndGet(); // ecriture volatile : publie les ecritures de spins
            }
        }
    }

    // --------------------------------------------------------------- noyaux

    /**
     * Champ local {@code h = H + sum_n M_n S_nb[n]} en convention symetrique (triangle superieur,
     * comme {@code Spin.exchangeEnergySym}). Les liaisons sont pre-triees par forme a la
     * construction (ordre canonique FULL -&gt; DEUX_PARAM -&gt; DIAG -&gt; KIT01 -&gt; KIT02 -&gt; KIT12,
     * tri stable, cf. javadoc de classe) et compactees par groupe dos a dos dans {@code compact} ;
     * {@code groups} porte 4 ints par groupe {forme, premiere liaison, nb de liaisons, debut
     * compact} : le noyau branche une seule fois par groupe puis deroule une boucle serree sur
     * les liaisons du groupe, chaque cas omettant les termes nuls (derive ~1e-16 assumee pour
     * 2-PARAM par reassociation). {@code nb} est la table reordonnee du site ({@code reorderedNb})
     * : le m du groupe designe le meme voisin que le m du compact et du gamma. Sur une entree non
     * deja triee, l'ordre d'accumulation differe du noyau historique plat : derive ~1e-16
     * acceptee. Statique et sans allocation : inlinee par le JIT dans les trois noyaux,
     * {@code out} scalarise par l'analyse d'echappement.
     */
    private static void localField(double[] S, int[] nb, double[] compact, int[] groups,
                                   double hx0, double hy0, double hz0, double[] out) {
        double hx = hx0, hy = hy0, hz = hz0;
        for (int g = 0; g < groups.length; g += 4) {
            final int mStart = groups[g + 1];
            final int mEnd = mStart + groups[g + 2];
            switch (groups[g]) {
                case F_DEUX_PARAM -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 2) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 1] * (s1 + s2);
                        hy += compact[b + 1] * s0 + compact[b] * s1 + compact[b + 1] * s2;
                        hz += compact[b] * s2 + compact[b + 1] * (s0 + s1);
                    }
                }
                case F_DIAG -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 3) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0;
                        hy += compact[b + 1] * s1;
                        hz += compact[b + 2] * s2;
                    }
                }
                case F_KIT01 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 3] * s1;
                        hy += compact[b + 3] * s0 + compact[b + 1] * s1;
                        hz += compact[b + 2] * s2;
                    }
                }
                case F_KIT02 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 3] * s2;
                        hy += compact[b + 1] * s1;
                        hz += compact[b + 3] * s0 + compact[b + 2] * s2;
                    }
                }
                case F_KIT12 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0;
                        hy += compact[b + 1] * s1 + compact[b + 3] * s2;
                        hz += compact[b + 3] * s1 + compact[b + 2] * s2;
                    }
                }
                default -> { // F_FULL
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 9) {
                        final int o = nb[m] * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 1] * s1 + compact[b + 2] * s2;
                        hy += compact[b + 1] * s0 + compact[b + 4] * s1 + compact[b + 5] * s2;
                        hz += compact[b + 2] * s0 + compact[b + 5] * s1 + compact[b + 8] * s2;
                    }
                }
            }
        }
        out[0] = hx;
        out[1] = hy;
        out[2] = hz;
    }

    /**
     * Variante de {@link #localField} restreinte aux voisins {@code j < site} (chaque paire
     * comptee une seule fois), sans champ exterieur : la somme de l'energie totale. Meme
     * structure par groupes pre-tries que {@link #localField} (cf. javadoc de classe et de
     * {@code buildGroupedBondTables}), le filtre {@code site <= j} etant applique par liaison
     * dans chaque boucle de groupe ; sur une entree non deja triee, l'ordre d'accumulation
     * differe du noyau historique plat : derive ~1e-16 acceptee. Boucle distincte pour garder le
     * filtre hors du noyau general.
     */
    private static void localFieldLowerPairs(double[] S, int[] nb, double[] compact, int[] groups,
                                             int site, double[] out) {
        double hx = 0.d, hy = 0.d, hz = 0.d;
        for (int g = 0; g < groups.length; g += 4) {
            final int mStart = groups[g + 1];
            final int mEnd = mStart + groups[g + 2];
            switch (groups[g]) {
                case F_DEUX_PARAM -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 2) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 1] * (s1 + s2);
                        hy += compact[b + 1] * s0 + compact[b] * s1 + compact[b + 1] * s2;
                        hz += compact[b] * s2 + compact[b + 1] * (s0 + s1);
                    }
                }
                case F_DIAG -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 3) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0;
                        hy += compact[b + 1] * s1;
                        hz += compact[b + 2] * s2;
                    }
                }
                case F_KIT01 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 3] * s1;
                        hy += compact[b + 3] * s0 + compact[b + 1] * s1;
                        hz += compact[b + 2] * s2;
                    }
                }
                case F_KIT02 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 3] * s2;
                        hy += compact[b + 1] * s1;
                        hz += compact[b + 3] * s0 + compact[b + 2] * s2;
                    }
                }
                case F_KIT12 -> {
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 4) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0;
                        hy += compact[b + 1] * s1 + compact[b + 3] * s2;
                        hz += compact[b + 3] * s1 + compact[b + 2] * s2;
                    }
                }
                default -> { // F_FULL
                    for (int m = mStart, b = groups[g + 3]; m < mEnd; m++, b += 9) {
                        final int j = nb[m];
                        if (site <= j) continue; // paire comptee une seule fois
                        final int o = j * 3;
                        final double s0 = S[o], s1 = S[o + 1], s2 = S[o + 2];
                        hx += compact[b] * s0 + compact[b + 1] * s1 + compact[b + 2] * s2;
                        hy += compact[b + 1] * s0 + compact[b + 4] * s1 + compact[b + 5] * s2;
                        hz += compact[b + 2] * s0 + compact[b + 5] * s1 + compact[b + 8] * s2;
                    }
                }
            }
        }
        out[0] = hx;
        out[1] = hy;
        out[2] = hz;
    }

    /**
     * Noyau Metropolis d'une tranche : champ local via les liaisons triees, proposition
     * gaussienne renormalisee, dE quartique inclus, acceptation Metropolis ; compteurs
     * de tranche (tentatives, acceptations).
     */
    private void metropolisChunk(int color, int chunk) {
        final int start = chunkStart[color][chunk];
        final int end = chunkEnd[color][chunk];
        if (start >= end) return;
        final int[] sites = coloring.classSites(color);
        final RandomGenerator r = rng[color][chunk];
        final double[] S = spins;
        final int lsr = nSublattices;
        final double temperature = phaseTemperature;
        final double sigma = phaseSigma;
        final double hx0 = fieldX, hy0 = fieldY, hz0 = fieldZ;
        final double[] h = new double[3]; // une seule fois par tranche, hors de la boucle chaude
        final boolean quartic = lattice.hasQuartic();
        final double quarticB = lattice.getQuarticB();
        long attempted = 0L, accepted = 0L;

        for (int t = start; t < end; t++) {
            final int site = sites[t];
            final int sr = site % lsr;
            final int[] nb = reorderedNb[site];
            localField(S, nb, compactBySublattice[sr], groupsBySublattice[sr], hx0, hy0, hz0, h);
            final double hx = h[0], hy = h[1], hz = h[2];

            final int off = site * 3;
            final double ox = S[off], oy = S[off + 1], oz = S[off + 2];
            double nx = ox + sigma * r.nextGaussian();
            double ny = oy + sigma * r.nextGaussian();
            double nz = oz + sigma * r.nextGaussian();
            final double inv = 1.d / Math.sqrt(nx * nx + ny * ny + nz * nz);
            nx *= inv;
            ny *= inv;
            nz *= inv;

            double dE = -((nx - ox) * hx + (ny - oy) * hy + (nz - oz) * hz);
            if (quartic) {
                // dE4 = -b * somme_j (S_j^g)^2 * [(n^g)^2 - (o^g)^2], g = gamma de la liaison
                // (site, j) ; le site porte toutes ses liaisons (pas de garde i<j ici, elle ne
                // sert qu'au comptage de l'energie totale, cf. Spin.getQuarticEnergy).
                // Gamma reordonne comme les liaisons triees (positionnel brut dans Lattice,
                // cf. buildGroupedBondTables) : le g designe la meme liaison que nb[i].
                final int[] gamma = gammaReordered[sr];
                double qx = 0.d, qy = 0.d, qz = 0.d;
                for (int i = 0; i < nb.length; i++) {
                    final int g = gamma[i];
                    final double sjg = S[nb[i] * 3 + g];
                    if (g == 0) qx += sjg * sjg;
                    else if (g == 1) qy += sjg * sjg;
                    else qz += sjg * sjg;
                }
                dE -= quarticB * (qx * (nx * nx - ox * ox)
                        + qy * (ny * ny - oy * oy)
                        + qz * (nz * nz - oz * oz));
            }
            attempted++;
            if (dE <= 0.d || r.nextDouble() < Math.exp(-dE / temperature)) {
                S[off] = nx;
                S[off + 1] = ny;
                S[off + 2] = nz;
                accepted++;
            }
        }
        counters[chunk * LSTRIDE] += attempted;
        counters[chunk * LSTRIDE + 1] += accepted;
    }

    /** Reflexion microcanonique {@code S <- 2(S.h^)h^ - S} sur une tranche ; champ nul : site saute. */
    private void overRelaxChunk(int color, int chunk) {
        final int start = chunkStart[color][chunk];
        final int end = chunkEnd[color][chunk];
        if (start >= end) return;
        final int[] sites = coloring.classSites(color);
        final double[] S = spins;
        final int lsr = nSublattices;
        final double hx0 = fieldX, hy0 = fieldY, hz0 = fieldZ;
        final double[] h = new double[3];

        for (int t = start; t < end; t++) {
            final int site = sites[t];
            final int sr = site % lsr;
            final int[] nb = reorderedNb[site];
            localField(S, nb, compactBySublattice[sr], groupsBySublattice[sr], hx0, hy0, hz0, h);
            final double hx = h[0], hy = h[1], hz = h[2];

            final double h2 = hx * hx + hy * hy + hz * hz;
            if (h2 <= 1e-24) continue; // |h| <= 1e-12
            final int off = site * 3;
            final double ox = S[off], oy = S[off + 1], oz = S[off + 2];
            final double f = 2.d * (ox * hx + oy * hy + oz * hz) / h2;
            S[off] = f * hx - ox;
            S[off + 1] = f * hy - oy;
            S[off + 2] = f * hz - oz;
        }
    }

    /** Energie d'une tranche globale : paires ({@code j < site}) et champ exterieur, dans {@code partials}. */
    private void energyChunk(int chunk) {
        final int start = globalStart[chunk];
        final int end = globalEnd[chunk];
        final double[] S = spins;
        final int lsr = nSublattices;
        final double hx0 = fieldX, hy0 = fieldY, hz0 = fieldZ;
        final double[] h = new double[3];
        double e = 0.d;

        for (int site = start; site < end; site++) {
            final int sr = site % lsr;
            final int[] nb = reorderedNb[site];
            localFieldLowerPairs(S, nb, compactBySublattice[sr], groupsBySublattice[sr], site, h);
            final int off = site * 3;
            final double sx = S[off], sy = S[off + 1], sz = S[off + 2];
            e -= sx * h[0] + sy * h[1] + sz * h[2];
            e -= sx * hx0 + sy * hy0 + sz * hz0;
        }
        partials[chunk * DSTRIDE] = e;
    }

    /** Somme des spins d'une tranche globale (base de Kitaev) dans {@code partials}. */
    private void magnetizationChunk(int chunk) {
        final double[] S = spins;
        double mx = 0.d, my = 0.d, mz = 0.d;
        for (int site = globalStart[chunk]; site < globalEnd[chunk]; site++) {
            final int off = site * 3;
            mx += S[off];
            my += S[off + 1];
            mz += S[off + 2];
        }
        partials[chunk * DSTRIDE] = mx;
        partials[chunk * DSTRIDE + 1] = my;
        partials[chunk * DSTRIDE + 2] = mz;
    }
}
