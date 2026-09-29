import java.util.ArrayList;

/**
 * Reseau periodique de spins construit depuis la maille {@link UnitCell} (honeycomb, 2
 * sous-reseaux) : tables de voisins, gabarits compacts des matrices d'echange J_ij par
 * voisin, champ d'interaction (Zeeman) et terme quartique premiers voisins optionnel.
 * Spins unitaires stockes a plat {@code [site*3 + (Sx, Sy, Sz)]}.
 */
public class Lattice {

    int size;
    int lengthA; // nombre de cellules selon a
    int lengthB; // nombre de cellules selon b
    int lengthSr; // nombre de sous-réseaux
    int[][] interactionSites; // [site][voisin] : index du voisin
    private double[] spins1D; // configuration des spins [site * 3 + (Sx, Sy, Sz)]
    private double[][] Gtensor;
    double[] interactionField;
    private double[][] Kitaev_basis; // base de Kitaev
    private double[][] COB; // changement de base abc -> Kitaev

    private double[][][] templateMatrices; // [sous-reseau][liaison] gabarit : {a,a,c} (XY-Z) | {d0,d1,d2} | {a,b} (2-PARAM) | {d0,d1,d2,k,v} | 9 valeurs
    private byte[] siteSublattice; // [N_total] - 1 byte par site

    // Terme quartique NN : -b (S_i^g S_j^g)^2. quarticGammaTemplates[sub][i] = gamma
    // (0=x,1=y,2=z) de la liaison template i du sous-reseau sub. b = 0 (defaut) :
    // aucune physique ajoutee, comportement historique bit-identique.
    private double quarticB = 0.d;
    private int[][] quarticGammaTemplates;

    /** Alloue le reseau N = {cellules selon a, cellules selon b, sous-reseaux} ; spins a zero. */
    public Lattice(int[] N) {
        lengthA = N[0];
        lengthB = N[1];
        lengthSr = N[2];
        size = lengthA * lengthB * lengthSr;
        spins1D = new double[size * 3];
        interactionSites = new int[size][];
        interactionField = new double[3];
        Kitaev_basis = new double[3][3];
        Gtensor = new double[3][3];
        COB = new double[3][3];
    }

    /** Index global d'un site : cellule (i, j) repliée périodiquement, sous-réseau sr. */
    int siteIndex(int i, int j, int sr) {
        return (Math.floorMod(j, lengthB) * lengthA + Math.floorMod(i, lengthA)) * lengthSr + sr;
    }

    /** Installe la base de Kitaev (3 axes). */
    void setKitaevBasis(double basis[][]) {
        this.Kitaev_basis = basis;
    }

    /** Installe la matrice de changement de base abc -> Kitaev. */
    void setCOB(double basis[][]) {
        this.COB = basis;
    }

    /** Installe le tenseur g (3x3) qui construit le champ exterieur. */
    void setGtensor(double g[][]) {
        this.Gtensor = g;
    }

    /** Tenseur g (3x3). */
    double[][] getGtensor() {
        return this.Gtensor;
    }

    /** Ecrit les 3 composantes du spin du site. */
    void setSpin1D(int site, double[] xyz_vec) {
        int index = 3 * site;
        spins1D[index++] = xyz_vec[0];
        spins1D[index++] = xyz_vec[1];
        spins1D[index++] = xyz_vec[2];
    }

    /** Remplace le tableau plat des spins (reference conservee, pas de copie). */
    public void setSpin1Dlattice(double[] spins1Dlattice) {
        this.spins1D = spins1Dlattice;
    }

    /** Copie des 3 composantes du spin du site. */
    public double[] getSpin1D(int site) {
        int index = 3 * site;
        return new double[] { spins1D[index], spins1D[index + 1], spins1D[index + 2] };
    }

    /** Tableau plat complet des spins (reference, pas copie) : [site*3 + (Sx, Sy, Sz)]. */
    double[] getSpin1Dlattice() {
        return this.spins1D;
    }

    /** Base de Kitaev (3 axes). */
    double[][] getKitaevBasis() {
        return Kitaev_basis;
    }

    /** Matrice de changement de base abc -> Kitaev. */
    double[][] getChangeOfBasis() {
        return COB;
    }

    /** Indices globaux des voisins du site, dans l'ordre de declaration de la maille. */
    int[] getNeighborSites(int site) {
        return interactionSites[site];
    }

    /**
     * Gabarits compacts J_ij du sous-reseau du site : [voisin][gabarit]. Tableau partage par
     * tous les sites d'un meme sous-reseau (base du groupement par identite dans Spin.getdE2).
     */
    double[][] getInteractionMatrices(int site) {
        return templateMatrices[siteSublattice[site]];
    }

    /** Champ d'interaction B (Zeeman), reference directe. */
    double[] getInteractionField() {
        return interactionField;
    }

    /** Installe le champ d'interaction (reference conservee). */
    void setInteractionField(double[] B) {
        interactionField = B;
    }

    // ------------------------------------------------------------------ terme quartique NN

    /**
     * Installe le terme quartique {@code -b (S_i^g S_j^g)^2} sur les premiers voisins.
     *
     * @param b           couplage quartique ; 0 revient a desactiver le terme.
     * @param gammaSub0   gamma (0=x, 1=y, 2=z) de chaque liaison template du sous-reseau 0.
     * @param gammaSub1   idem sous-reseau 1 (derive par anti-deplacement dans {@code Input}).
     */
    void setQuartic(double b, int[] gammaSub0, int[] gammaSub1) {
        this.quarticB = b;
        this.quarticGammaTemplates = new int[][] { gammaSub0.clone(), gammaSub1.clone() };
    }

    /** Couplage quartique NN ; 0 = aucun terme (comportement historique). */
    double getQuarticB() {
        return quarticB;
    }

    /** {@code true} si le terme quartique est actif (b != 0). */
    boolean hasQuartic() {
        return quarticB != 0.d;
    }

    /**
     * Gamma (0=x, 1=y, 2=z) de la liaison (site, voisin i) pour le terme quartique.
     * Ne doit etre appele que si {@link #hasQuartic()} (la ligne QUARTIC fournit le gamma
     * de chaque voisin premier) — sinon le tableau n'existe pas.
     */
    int getQuarticGamma(int site, int neighborIdx) {
        return quarticGammaTemplates[siteSublattice[site]][neighborIdx];
    }

    /** Deroule la maille pour construire voisins et sous-reseau de chaque site. */
    void setInteractionSites(UnitCell uc) {
        int totalSites = lengthA * lengthB * lengthSr;
        siteSublattice = new byte[totalSites]; // Allouer une seule fois

        int idx = 0;
        for (int j = 0; j < lengthB; j++) {
            for (int i = 0; i < lengthA; i++) {
                for (int t = 0; t < lengthSr; t++) {
                    int[] neighborSites = new int[uc.interactions.get(t).size()];

                    ArrayList<double[]> interaction = uc.interactions.get(t);
                    for (int v = 0; v < interaction.size(); v++) {
                        double[] x = interaction.get(v);
                        int[] offset = { (int) x[1], (int) x[2], (int) x[3] };
                        neighborSites[v] = siteIndex(offset[0] + i, offset[1] + j, offset[2]);
                    }

                    interactionSites[idx] = neighborSites;
                    siteSublattice[idx] = (byte) t; // Stocker le sous-réseau (0 ou 1)
                    idx++;
                }
            }
        }
    }


    /** Extrait de la maille le gabarit compact de chaque liaison, par sous-reseau. */
    void createTemplates(UnitCell uc) {
        templateMatrices = new double[lengthSr][][];

        for (int sr = 0; sr < lengthSr; sr++) {
            ArrayList<double[]> interaction = uc.interactions.get(sr);
            templateMatrices[sr] = new double[interaction.size()][];

            for (int v = 0; v < interaction.size(); v++) {
                double[] x = interaction.get(v);

                // Stocker LE gabarit compact de la forme taggée
                final int tag = (int) x[4];
                if (tag == 0) {
                    // XY-Z : diag(a, a, c) -> gabarit diagonale {a, a, c}
                    templateMatrices[sr][v] = new double[] { x[5], x[5], x[6] };
                } else if (tag == 1) {
                    // Diagonale générale -> {d0, d1, d2}
                    templateMatrices[sr][v] = new double[] { x[5], x[6], x[7] };
                } else if (tag == 2) {
                    // 2-PARAM a.I + b.A -> gabarit compact {a, b}
                    templateMatrices[sr][v] = new double[] { x[5], x[6] };
                } else if (tag == 3) {
                    // Kitaev -> {d0, d1, d2, k, v}
                    templateMatrices[sr][v] = new double[] { x[5], x[6], x[7], x[8], x[9] };
                } else if (tag == 4) {
                    // Pleine -> 9 valeurs
                    templateMatrices[sr][v] = new double[] {
                            x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13]
                    };
                } else {
                    throw new IllegalArgumentException("sous-reseau " + sr + ", liaison " + v
                            + " : tag de forme inconnu " + tag);
                }
            }
        }
    }
}
