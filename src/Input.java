import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * Lit le fichier d'entree du modele. Format attendu :
 * <ol>
 *   <li>lignes 1-3 : matrice de changement de base (3 lignes de 3 valeurs), stockee transposee ;</li>
 *   <li>par sous-reseau (0 puis 1) : position (3 valeurs : x, y fractionnaires (a1, a2) et z
 *       fractionnaire hors-plan), g-tenseur diagonal
 *       (3 valeurs), base de Kitaev (9 valeurs, stockee transposee), sites en interaction
 *       (dx dy dz sous-reseau par voisin, dz ignore — le nombre de voisins est deduit de la
 *       ligne et doit etre identique pour les deux sous-reseaux), matrices d'interaction
 *       (9 valeurs par voisin) ;</li>
 *   <li>ligne finale optionnelle : QUARTIC b g0 g1 g2.</li>
 * </ol>
 */

public class Input {

    /** Nombre de valeurs par voisin sur une ligne de sites d'interaction : dx dy dz sous-reseau. */
    private static final int VALUES_PER_NEIGHBOR = 4;

    /** Nombre de valeurs par voisin sur une ligne de matrices d'interaction : matrice 3x3 a plat. */
    private static final int VALUES_PER_MATRIX = 9;

    /** nombre de sous-réseau : 2 pour le nid d'abeille */
    private static final int SUB_LATTICE = 2;

    private double[][] changeOfBasis; // change of basis : axes (a, b, c) dans la base Kitaev, transposee
    private double[][] atomicPositions; // positions (x, y, z) fractionnaires, [sous-reseau][coord]
    private double[][][] gTensors; // g-tenseur diagonal par sous-reseau, base a*bc
    private double[][][] kitaevBasis; // base de Kitaev par sous-reseau, stockee transposee
    private byte[][][] interactionSites; // voisins en interaction, [sous-reseau][voisin][dx,dy,sub]
    private double[][][][] interactionMatrices; // matrices d'interaction, [sous-reseau][voisin][3][3]
    private int neighborCount; // nombre de voisins en interaction par site

    // Terme quartique optionnel (premiers voisins) : ligne "QUARTIC <b> <g0> <g1> <g2>".
    // g* = composante gamma (x|y|z ou 0|1|2) de la liaison (sous-reseau 0, voisin i), dans
    // l'ordre de la ligne des voisins du sous-reseau 0. Absente -> b = 0 (historique).
    private double quarticB = 0.d;
    private int[] quarticGammaSub0;
    private int[] quarticGammaSub1;

    public Input(String inputFileName) {
        readFile(inputFileName);
    }

    private void readFile(String inputFileName) {
        try (BufferedReader br = new BufferedReader(new FileReader(inputFileName))) {
            changeOfBasis = readChangeOfBasis(br);

            atomicPositions = new double[SUB_LATTICE][];
            gTensors = new double[SUB_LATTICE][][];
            kitaevBasis = new double[SUB_LATTICE][][];
            interactionSites = new byte[SUB_LATTICE][][];
            interactionMatrices = new double[SUB_LATTICE][][][];

            // sous-reseau 0
            atomicPositions[0] = readPosition(br);
            gTensors[0] = readGTensor(br);
            kitaevBasis[0] = readKitaevBasis(br);
            interactionSites[0] = readInteractionSites(br);
            neighborCount = interactionSites[0].length;
            interactionMatrices[0] = readInteractionMatrices(br, neighborCount);

            // sous-reseau 1
            atomicPositions[1] = readPosition(br);
            gTensors[1] = readGTensor(br);
            kitaevBasis[1] = readKitaevBasis(br);
            interactionSites[1] = readInteractionSites(br);
            if (interactionSites[1].length != neighborCount) {
                throw new IllegalArgumentException("sites d'interaction : les deux sous-reseaux "
                        + "doivent declarer le meme nombre de voisins (" + neighborCount + " puis "
                        + interactionSites[1].length + ")");
            }
            interactionMatrices[1] = readInteractionMatrices(br, neighborCount);

            readQuartic(br);
        } catch (IOException e) {
            throw new IllegalArgumentException("Impossible de lire le fichier d'entree : "
                    + inputFileName, e);
        }
    }

    /**
     * Sites en interaction d'un sous-reseau : (dx, dy, dz, sous-reseau d'arrivee) par voisin.
     * dz est ignore (reseau planaire) ; tableau retourne = [dx, dy, sous-reseau d'arrivee].
     */
    private static byte[][] readInteractionSites(BufferedReader br) throws IOException {
        final String[] tokens = br.readLine().trim().split("\\s+");
        if (tokens.length % VALUES_PER_NEIGHBOR != 0) {
            throw new IllegalArgumentException("sites d'interaction : " + tokens.length
                    + " valeurs, multiple de " + VALUES_PER_NEIGHBOR
                    + " attendu (dx dy dz sous-reseau par voisin)");
        }
        final byte[][] sites = new byte[tokens.length / VALUES_PER_NEIGHBOR][3];
        for (int i = 0; i < sites.length; i++) {
            // La valeur 2 (dz) est sautee, cf. javadoc.
            for (int j = 0, k = 0; j < VALUES_PER_NEIGHBOR; j++) {
                if (j != 2) {
                    sites[i][k++] = Byte.parseByte(tokens[i * VALUES_PER_NEIGHBOR + j]);
                }
            }
        }
        return sites;
    }

    /** Lit une ligne et decoupe ses valeurs */
    private static String[] readTokens(BufferedReader br, int expected, String what) throws IOException {
        final String[] tokens = br.readLine().trim().split("\\s+");
        if (tokens.length != expected) {
            throw new IllegalArgumentException(what + " : " + tokens.length
                    + " valeurs, " + expected + " attendues");
        }
        return tokens;
    }

    /** Matrice de changement de base : 3 lignes successives de 3 valeurs, stockee transposee. */
    private static double[][] readChangeOfBasis(BufferedReader br) throws IOException {
        final double[][] cob = new double[3][3];
        for (int row = 0; row < 3; row++) {
            final String[] tokens = readTokens(br, 3, "change de base, ligne " + (row + 1));
            for (int col = 0; col < 3; col++) {
                cob[row][col] = Double.parseDouble(tokens[col]);
            }
        }
        return cob;
    }

    /**
     * Position d'un sous-reseau : 1 ligne de 3 valeurs KEPT — x, y fractionnaires (a1, a2) et
     * z fractionnaire hors-plan. Les trois composantes alimentent {@code atomicPositions} ;
     * la composante z alimente le décalage de phase inter-sous-réseaux du S(Q,ω)
     * (deltaZ = z_B − z_A, cf. {@code DynamicStructureFactor}).
     */
    private static double[] readPosition(BufferedReader br) throws IOException {
        final String[] tokens = readTokens(br, 3, "position");
        return new double[] { Double.parseDouble(tokens[0]), Double.parseDouble(tokens[1]),
                Double.parseDouble(tokens[2]) };
    }

    /** g-tenseur diagonal d'un sous-reseau (base a*bc) : 1 ligne de 3 valeurs. */
    private static double[][] readGTensor(BufferedReader br) throws IOException {
        final String[] tokens = readTokens(br, 3, "g-tenseur");
        final double[][] g = new double[3][3];
        for (int axis = 0; axis < 3; axis++) {
            g[axis][axis] = Double.parseDouble(tokens[axis]);
        }
        return g;
    }

    /** Base de Kitaev d'un sous-reseau, stockee transposee : 1 ligne de 9 valeurs. */
    private static double[][] readKitaevBasis(BufferedReader br) throws IOException {
        final String[] tokens = readTokens(br, VALUES_PER_MATRIX, "base de Kitaev");
        final double[][] basis = new double[3][3];
        for (int row = 0; row < 3; row++) {
            for (int col = 0; col < 3; col++) {
                basis[row][col] = Double.parseDouble(tokens[row * 3 + col]);
            }
        }
        return basis;
    }

    /** Matrices d'interaction 3x3 d'un sous-reseau : 1 ligne de 9 valeurs par voisin. */
    private static double[][][] readInteractionMatrices(BufferedReader br, int count) throws IOException {
        final String[] tokens = readTokens(br, count * VALUES_PER_MATRIX, "matrices d'interaction");
        final double[][][] matrices = new double[count][3][3];
        for (int neighbor = 0; neighbor < count; neighbor++) {
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    matrices[neighbor][row][col] = Double.parseDouble(
                            tokens[neighbor * VALUES_PER_MATRIX + row * 3 + col]);
                }
            }
        }
        return matrices;
    }

    /** Terme quartique optionnel : ligne finale "QUARTIC <b> <g0> <g1> <g2>", absente -> b = 0. */
    private void readQuartic(BufferedReader br) throws IOException {
        String extra;
        while ((extra = br.readLine()) != null) {
            final String trimmed = extra.trim();
            if (trimmed.isEmpty()) continue;
            if (!trimmed.startsWith("QUARTIC")) {
                throw new IllegalArgumentException(
                        "ligne finale non reconnue (QUARTIC attendu, ou rien) : " + trimmed);
            }
            if (quarticB != 0.d || quarticGammaSub0 != null) {
                throw new IllegalArgumentException("ligne QUARTIC dupliquee");
            }
            final String[] q = trimmed.split("\\s+");
            if (q.length != 5) {
                throw new IllegalArgumentException(
                        "QUARTIC attendu : QUARTIC <b> <g0> <g1> <g2> (g = x|y|z)");
            }
            if (neighborCount != 3) {
                throw new IllegalArgumentException("QUARTIC : mode restreint aux premiers "
                        + "voisins (3) ; le fichier declare " + neighborCount + " voisins");
            }
            quarticB = Double.parseDouble(q[1]);
            quarticGammaSub0 = new int[neighborCount];
            for (int i = 0; i < 3; i++) {
                quarticGammaSub0[i] = switch (q[2 + i]) {
                    case "x", "X", "0" -> 0;
                    case "y", "Y", "1" -> 1;
                    case "z", "Z", "2" -> 2;
                    default -> throw new IllegalArgumentException(
                            "gamma x|y|z attendu : " + q[2 + i]);
                };
            }
            // Gamma du sous-reseau 1 : chaque liaison NN vue de sub1 est l'anti-deplacement
            // exact d'une liaison de sub0 ; le gamma est une propriete de la liaison.
            quarticGammaSub1 = new int[neighborCount];
            final boolean[] matched0 = new boolean[neighborCount];
            for (int k = 0; k < neighborCount; k++) {
                final int dx1 = interactionSites[1][k][0];
                final int dy1 = interactionSites[1][k][1];
                boolean matched = false;
                for (int j = 0; j < neighborCount && !matched; j++) {
                    if (matched0[j]) continue;
                    if (interactionSites[0][j][0] == -dx1
                            && interactionSites[0][j][1] == -dy1
                            && interactionSites[0][j][2] == 1) {
                        quarticGammaSub1[k] = quarticGammaSub0[j];
                        matched0[j] = true;
                        matched = true;
                    }
                }
                if (!matched) {
                    throw new IllegalArgumentException("QUARTIC : la liaison du sous-reseau 1 ("
                            + dx1 + "," + dy1 + ") n'a pas d'anti-deplacement dans le sous-reseau 0");
                }
            }
        }
    }

    /** Matrice de changement de base : axes (a, b, c) dans la base de Kitaev, stockee transposee. */
    double[][] getChangeOfBasis() {
        return changeOfBasis;
    }

    /** Positions (x, y, z) des sites en coordonnees fractionnaires, [sous-reseau][coordonnee]. */
    double[][] getAtomicPositions() {
        return atomicPositions;
    }

    /** Base de Kitaev par sous-reseau, stockee transposee. */
    double[][][] getKitaevBasis() {
        return kitaevBasis;
    }

    /** g-tenseur diagonal par sous-reseau, exprime dans la base a*bc. */
    double[][][] getGTensors() {
        return gTensors;
    }

    /** Voisins en interaction, [sous-reseau][voisin] = [dx, dy, sous-reseau d'arrivee]. */
    byte[][][] getInteractionSites() {
        return interactionSites;
    }

    /** Matrices d'interaction 3x3 par voisin, [sous-reseau][voisin][ligne][colonne]. */
    double[][][][] getInteractionMatrices() {
        return interactionMatrices;
    }

    /** Nombre de voisins en interaction par site (fixe par la ligne du sous-reseau 0). */
    int getNeighborCount() {
        return neighborCount;
    }

    /** Couplage quartique premiers voisins b de -b(S_i^g S_j^g)^2 ; 0 si la ligne est absente. */
    double getQuarticB() {
        return quarticB;
    }

    /** Gamma (0=x, 1=y, 2=z) des liaisons NN du sous-reseau 0 ; null si b = 0. */
    int[] getQuarticGammaSub0() {
        return quarticGammaSub0;
    }

    /** Gamma des liaisons NN du sous-reseau 1 (derive par anti-deplacement) ; null si b = 0. */
    int[] getQuarticGammaSub1() {
        return quarticGammaSub1;
    }

}
