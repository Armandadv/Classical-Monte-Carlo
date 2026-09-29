import java.util.ArrayList;
import java.util.List;

/**
 * Geometrie du reseau nid d'abeille (honeycomb) a deux sites de base, en <b>un seul</b> endroit.
 *
 * <p>La meme construction — vecteurs primitifs a1 = (1, 0) et a2 = (-1/2, sqrt(3)/2), sites de
 * base A = (0, 0, 0) et B = (2/3, 1/3, 0) en coordonnees fractionnaires (3 elements, z nul :
 * reseau planaire de test), couches de voisins trouvees
 * par distance — etait recopiee dans {@code AppParallel.honeycomb} (pour {@code --bench}),
 * {@code TestInputGenerator.neighbours} et {@code TestLatticeFactory.honeycomb}. Trois copies
 * d'une meme boucle a triple imbrication sur {@code (da, db, s2)} avec trois seuils de distance
 * en dur : une divergence d'un seul {@code TOL} ou d'une seule borne de balayage aurait donne
 * deux hamiltoniens differents dans deux tests censes se comparer. Les distances et le balayage
 * vivent donc ici.</p>
 *
 * <h2>Couches de voisins</h2>
 * <table>
 *   <caption>Distances et coordinations</caption>
 *   <tr><th>couche</th><th>distance</th><th>voisins</th><th>sous-reseau cible</th></tr>
 *   <tr><td>1</td><td>{@code D1 = 1/sqrt(3)}</td><td>3</td><td>l'autre</td></tr>
 *   <tr><td>2</td><td>{@code D2 = 1}</td><td>6</td><td>le meme</td></tr>
 *   <tr><td>3</td><td>{@code D3 = 2/sqrt(3)}</td><td>3</td><td>l'autre</td></tr>
 * </table>
 *
 * <h2>Ordre de parcours</h2>
 * <p>{@link #neighbours(int, int)} renvoie les voisins tries par {@code da} croissant, puis
 * {@code db}, puis sous-reseau cible : exactement l'ordre des boucles historiques
 * {@code for (da = -3..3) for (db = -3..3) for (s2 = 0..1)}. L'energie ne depend evidemment pas
 * de cet ordre, mais le format de fichier lu par {@link Input} associe positionnellement la
 * i-eme matrice d'echange au i-eme voisin : conserver l'ordre garde les fichiers generes
 * comparables a ceux d'avant.</p>
 *
 * <p>TODO : {@code test/TestLatticeFactory.java} construit encore ses voisins avec sa propre
 * copie de cette boucle. Il appartient a un autre lot et n'est pas modifie ici ; il pourrait
 * etre migre sur {@link #neighbours(int, int)} sans changer son comportement (memes distances,
 * meme {@code TOL}, meme ordre de parcours), ce qui supprimerait la derniere copie.</p>
 */
public final class HoneycombGeometry {

    private HoneycombGeometry() {}

    // Les vecteurs primitifs et les positions de base sont exposes par des methodes qui
    // renvoient des copies : un tableau `public static final` reste modifiable par n'importe
    // quel appelant, et UnitCell conserve la reference du tableau qu'on lui passe.
    private static final double[] A1 = { 1.d, 0.d };
    private static final double[] A2 = { -1.d / 2.d, Math.sqrt(3.d) / 2.d };
    private static final double[][] BASIS_FRACTIONAL = { { 0.d, 0.d, 0.d }, { 2.d / 3.d, 1.d / 3.d, 0.d } };

    /** Premier vecteur primitif a1 = (1, 0) (copie). */
    public static double[] a1() { return A1.clone(); }

    /** Second vecteur primitif a2 = (-1/2, sqrt(3)/2) (copie). */
    public static double[] a2() { return A2.clone(); }

    /**
     * Position fractionnaire du site de base {@code sublattice} : A = (0, 0, 0), B = (2/3, 1/3, 0)
     * (copie) — 3 éléments {x, y, 0} comme le fichier d'entrée : le z fractionnaire hors-plan
     * vaut 0 pour le réseau honeycomb de test (le z réel vient de l'input, cf. {@code Input} ;
     * il alimente le delta 3D de {@code DynamicStructureFactor}).
     */
    public static double[] basisFractional(int sublattice) {
        checkSublattice(sublattice);
        return BASIS_FRACTIONAL[sublattice].clone();
    }

    /** Distance des premiers voisins (3 par site, sur l'autre sous-reseau). */
    public static final double D1 = 1.d / Math.sqrt(3.d);
    /** Distance des deuxiemes voisins (6 par site, sur le meme sous-reseau). */
    public static final double D2 = 1.d;
    /** Distance des troisiemes voisins (3 par site, sur l'autre sous-reseau). */
    public static final double D3 = 2.d / Math.sqrt(3.d);

    /** Tolerance de comparaison des distances (les positions sont exactes a ~1e-16 pres). */
    public static final double TOL = 1e-6;

    /**
     * Amplitude du balayage sur les offsets de maille. La couche la plus lointaine traitee est
     * D3 = 2/sqrt(3) ~ 1.155 : un balayage de +/- 3 mailles (longueur 1) la couvre trois fois.
     */
    private static final int OFFSET_RANGE = 3;

    /** Nombre de sous-reseaux. */
    public static final int SUBLATTICES = 2;

    /**
     * Position cartesienne du site {@code sublattice} de la maille {@code (da, db)} :
     * {@code frac_x a1 + frac_y a2 + da a1 + db a2}.
     *
     * @param sublattice 0 (site A) ou 1 (site B).
     * @param da         offset de maille suivant a1.
     * @param db         offset de maille suivant a2.
     * @return {@code {x, y}}.
     */
    public static double[] position(int sublattice, int da, int db) {
        checkSublattice(sublattice);
        final double[] f = BASIS_FRACTIONAL[sublattice];
        return new double[] {
            f[0] * A1[0] + f[1] * A2[0] + da * A1[0] + db * A2[0],
            f[0] * A1[1] + f[1] * A2[1] + da * A1[1] + db * A2[1]
        };
    }

    /** Position cartesienne du site de base {@code sublattice} de la maille d'origine. */
    public static double[] position(int sublattice) {
        return position(sublattice, 0, 0);
    }

    /** Distance de la couche {@code shell} (1, 2 ou 3). */
    public static double shellDistance(int shell) {
        return switch (shell) {
            case 1 -> D1;
            case 2 -> D2;
            case 3 -> D3;
            default -> throw new IllegalArgumentException("couche 1, 2 ou 3 attendue : " + shell);
        };
    }

    /** Coordination de la couche {@code shell} : 3, 6, 3. */
    public static int coordination(int shell) {
        return switch (shell) {
            case 1 -> 3;
            case 2 -> 6;
            case 3 -> 3;
            default -> throw new IllegalArgumentException("couche 1, 2 ou 3 attendue : " + shell);
        };
    }

    /**
     * Couche a laquelle appartient une distance, ou 0 si elle ne correspond a aucune des trois
     * (a {@link #TOL} pres).
     */
    public static int shellOf(double distance) {
        if (Math.abs(distance - D1) < TOL) return 1;
        if (Math.abs(distance - D2) < TOL) return 2;
        if (Math.abs(distance - D3) < TOL) return 3;
        return 0;
    }

    /**
     * Voisins de couche {@code shell} du site de base {@code sublattice}.
     *
     * @param sublattice 0 (site A) ou 1 (site B).
     * @param shell      1, 2 ou 3.
     * @return une liste de {@code {da, db, s2}} : offset de maille et sous-reseau cible, dans
     *         l'ordre {@code da} puis {@code db} puis {@code s2} croissants. Sa taille vaut
     *         toujours {@link #coordination(int)}, ce qui est verifie avant de renvoyer.
     * @throws IllegalArgumentException si {@code sublattice} ou {@code shell} est hors domaine.
     * @throws IllegalStateException    si la coordination trouvee n'est pas celle attendue (le
     *                                  balayage ou une constante de distance seraient faux).
     */
    public static List<int[]> neighbours(int sublattice, int shell) {
        checkSublattice(sublattice);
        final double target = shellDistance(shell);
        final double[] p = position(sublattice, 0, 0);

        final List<int[]> out = new ArrayList<>(coordination(shell));
        for (int da = -OFFSET_RANGE; da <= OFFSET_RANGE; da++) {
            for (int db = -OFFSET_RANGE; db <= OFFSET_RANGE; db++) {
                for (int s2 = 0; s2 < SUBLATTICES; s2++) {
                    if (da == 0 && db == 0 && sublattice == s2) continue; // le site lui-meme
                    final double[] q = position(s2, da, db);
                    if (Math.abs(Math.hypot(q[0] - p[0], q[1] - p[1]) - target) < TOL) {
                        out.add(new int[] { da, db, s2 });
                    }
                }
            }
        }
        final int expected = coordination(shell);
        if (out.size() != expected) {
            throw new IllegalStateException("couche " + shell + " du sous-reseau " + sublattice
                    + " : " + out.size() + " voisins trouves, " + expected + " attendus");
        }
        return out;
    }

    /**
     * Nombre de voisins par site pour un jeu de couplages : une couche dont le J est nul est
     * omise a la construction, exactement comme dans {@code TestLatticeFactory}.
     */
    public static int expectedNeighbors(double J1, double J2, double J3) {
        return (J1 != 0.d ? coordination(1) : 0)
             + (J2 != 0.d ? coordination(2) : 0)
             + (J3 != 0.d ? coordination(3) : 0);
    }

    /** Matrice 3x3 identite (g-tenseur par defaut des reseaux de test). */
    public static double[][] identity() {
        return new double[][] { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
    }

    /** Matrice d'echange diagonale {@code diag(J, J, J)}. */
    public static double[][] diag(double j) {
        return new double[][] { { j, 0, 0 }, { 0, j, 0 }, { 0, 0, j } };
    }

    private static void checkSublattice(int sublattice) {
        if (sublattice < 0 || sublattice >= SUBLATTICES) {
            throw new IllegalArgumentException("sous-reseau 0 ou 1 attendu : " + sublattice);
        }
    }
}
