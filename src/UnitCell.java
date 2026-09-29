import java.util.ArrayList;

/**
 * Maille 2D du reseau : deux vecteurs primitifs ({@code primitive}, longueur 2), sites de
 * base en positions fractionnaires ({@code basis}, indice = sous-reseau), g-tenseurs et
 * interactions voisin par voisin. Chaque matrice d'echange est rotatee dans la base de
 * Kitaev puis rangee en gabarit compact taggue ({@link #addInteraction}) ; {@link Lattice}
 * n'en retient que ces gabarits et les offsets de voisins.
 */
class UnitCell {
    /** Tolérance de reconnaissance de la forme 2-PARAM (décision utilisateur explicite). */
    private static final double EPS_DEUX_PARAM = 1e-12;

    /** Vecteurs primitifs (a1, a2) : tableau de longueur 2. */
    double[][] primitive;
    /** Positions fractionnaires des sites de base (indice = sous-reseau). */
    ArrayList<double[]> basis;
    /** g-tenseurs 3x3 par sous-reseau ; seul le premier declare est lu. */
    ArrayList<ArrayList<double[][]>> interactionsGtensor;
    ArrayList<ArrayList<double[]>> interactions; // interaction sur la maille unité

    /** Maille unité 2D définie par ses deux vecteurs primitifs (a1, a2). */
    UnitCell(double[] a1, double[] a2) {
        this.primitive = new double[][] { a1, a2 };
        this.basis = new ArrayList<>();
        this.interactionsGtensor = new ArrayList<>();
        this.interactions = new ArrayList<>();
    }

    /**
     * Ajoute l'interaction de b1 vers b2 (offset = déplacement (dx, dy, dz) vers le voisin),
     * avec la matrice M exprimée dans la base a*bc. Elle est stockée transformée dans la
     * base de Kitaev (gtMg = gᵀ·M·g), sous la forme compacte taguée : un <b>octet de forme
     * en index 4</b>, après le préfixe {@code [b2, off0, off1, off2]} :
     * <ul>
     *   <li>tag 0 — XY-Z, 7 valeurs : {@code [b2, off×3, 0, a, c]} → diag(a, a, c) ;</li>
     *   <li>tag 1 — DIAG, 8 valeurs : {@code [b2, off×3, 1, d0, d1, d2]} → diag(d0, d1, d2) ;</li>
     *   <li>tag 2 — TWO_PARAM, 7 valeurs : {@code [b2, off×3, 2, a, b]} → a·I + b·A
     *       (A = hors-diagonale de 1), moyennes a et b stockées ;</li>
     *   <li>tag 3 — KITAEV, 10 valeurs : {@code [b2, off×3, 3, d0, d1, d2, k, v]} →
     *       diag + exactement une paire hors-diagonale k de valeur v (k ∈ {0,1,2} =
     *       paires 01/02/12), les trois paires jumelles étant <b>exactement égales</b> ;</li>
     *   <li>tag 4 — FULL, 14 valeurs : {@code [b2, off×3, 4, m11..m33]} → 3x3 générale,
     *       repli.</li>
     * </ul>
     * Classification, en ordre strict (il préserve le zéro-dérive des formes exactes) :
     * (1) hors-diagonaux exactement nuls ({@code == 0.}, signed-zero accepté) → formes
     * exactes sans aucune dérive ; (2) forme 2-PARAM testée <b>avec tolérance 1e-12</b>
     * (décision utilisateur explicite) : les moyennes a et b sont stockées à la place des
     * 9 valeurs. La dérive est <b>assumée</b> (≤ 1e-12 par entrée ; mesurée sur le fichier
     * réel : ~6e-17), y compris la perte d'une asymétrie sous-tolérance (la matrice est
     * alors stockée symétrique) ; l'ordre de sommation des moyennes est fixe. Sinon repli :
     * KITAEV (jumeaux exactement égaux, une seule paire non nulle — le contrôle
     * {@code jumeauxEgaux} est impératif : une paire asymétrique rangée symétrique
     * fausserait l'hamiltonien silencieusement) puis FULL, généraux et sans approximation.
     */
    void addInteraction(int b1, int b2, double[][] M, byte[] offset) {
        if (b1 == b2 && offset[0] == 0 && offset[1] == 0 && offset[2] == 0) {
            throw new RuntimeException("L'interaction ne peut pas être locale (b1 == b2 et offset nul).");
        }

        if (interactions.size() <= b1) {
            interactions.add(new ArrayList<>());
        }

        double[][] Mg = MathOps.mult2D(M, getInteractionGtensor(b1));
        final double[][] gtMg = MathOps.mult2D(MathOps.transpose(getInteractionGtensor(b1)), Mg);

        // Prédicats de forme, calculés avant la classification (lisibilité).
        final boolean offNuls = gtMg[0][1] == 0. && gtMg[0][2] == 0. && gtMg[1][0] == 0.
                && gtMg[1][2] == 0. && gtMg[2][0] == 0. && gtMg[2][1] == 0.;
        final boolean jumeauxEgaux = gtMg[0][1] == gtMg[1][0] && gtMg[0][2] == gtMg[2][0]
                && gtMg[1][2] == gtMg[2][1];
        final int pairesNonNulles = (gtMg[0][1] != 0. ? 1 : 0) + (gtMg[0][2] != 0. ? 1 : 0)
                + (gtMg[1][2] != 0. ? 1 : 0);
        final boolean xyZ = offNuls && gtMg[0][0] == gtMg[1][1];
        final boolean diagonaleGenerale = offNuls;
        final double a = (gtMg[0][0] + gtMg[1][1] + gtMg[2][2]) / 3.0;
        final double b = (gtMg[0][1] + gtMg[0][2] + gtMg[1][0]
                + gtMg[1][2] + gtMg[2][0] + gtMg[2][1]) / 6.0;
        final boolean deuxParam = deuxParamATolerance(gtMg, a, b, EPS_DEUX_PARAM);
        final boolean kitaev = jumeauxEgaux && pairesNonNulles == 1;

        // Classification, de la forme la plus contrainte au repli général.
        if (xyZ) {
            // XY-Z : diag(a, a, c) — tag 0, 7 valeurs
            interactions.get(b1).add(new double[] { b2, offset[0], offset[1], offset[2],
                    0, gtMg[0][0], gtMg[2][2] });
        } else if (diagonaleGenerale) {
            // Diagonale générale — tag 1, 8 valeurs
            interactions.get(b1).add(new double[] { b2, offset[0], offset[1], offset[2],
                    1, gtMg[0][0], gtMg[1][1], gtMg[2][2] });
        } else if (deuxParam) {
            // 2-PARAM a.I + b.A (tolérance 1e-12, dérive assumée par l'utilisateur) — tag 2, 7 valeurs
            interactions.get(b1).add(new double[] { b2, offset[0], offset[1], offset[2],
                    2, a, b });
        } else if (kitaev) {
            // KITAEV : diag + exactement une paire hors-diagonale — tag 3, 10 valeurs.
            // Le contrôle jumeauxEgaux est impératif : une paire asymétrique rangée
            // symétrique fausserait l'hamiltonien silencieusement.
            final int k = gtMg[0][1] != 0. ? 0 : (gtMg[0][2] != 0. ? 1 : 2);
            final double v = gtMg[0][1] != 0. ? gtMg[0][1]
                    : (gtMg[0][2] != 0. ? gtMg[0][2] : gtMg[1][2]);
            interactions.get(b1).add(new double[] { b2, offset[0], offset[1], offset[2],
                    3, gtMg[0][0], gtMg[1][1], gtMg[2][2], k, v });
        } else {
            // PLEINE : repli général — tag 4, 14 valeurs
            interactions.get(b1).add(new double[] { b2, offset[0], offset[1], offset[2],
                    4, gtMg[0][0], gtMg[0][1], gtMg[0][2],
                    gtMg[1][0], gtMg[1][1], gtMg[1][2],
                    gtMg[2][0], gtMg[2][1], gtMg[2][2] });
        }
    }

    /** Tolérance 2-PARAM : les 9 coefficients sont chacun à moins de eps de leur moyenne (a, b). */
    private static boolean deuxParamATolerance(double[][] m, double a, double b, double eps) {
        return Math.abs(m[0][0] - a) < eps && Math.abs(m[1][1] - a) < eps && Math.abs(m[2][2] - a) < eps
                && Math.abs(m[0][1] - b) < eps && Math.abs(m[0][2] - b) < eps
                && Math.abs(m[1][0] - b) < eps && Math.abs(m[1][2] - b) < eps
                && Math.abs(m[2][0] - b) < eps && Math.abs(m[2][1] - b) < eps;
    }

    /** Declare un g-tenseur 3x3 du sous-reseau b (le premier declare sert dans addInteraction). */
    void setInteractionGtensor(int b, double[][] M) {
        if (M.length != 3 || M[0].length != 3) {
            throw new RuntimeException("La matrice d'interaction doit être de taille 3x3.");
        }

        // Initialiser le tableau interactionsGtensor pour l'indice b s'il n'est pas
        // déjà initialisé
        if (interactionsGtensor.size() <= b) {
            interactionsGtensor.add(new ArrayList<double[][]>());
        }
        this.interactionsGtensor.get(b).add(M);
    }

    /** g-tenseur du sous-réseau b (le premier déclaré). */
    public double[][] getInteractionGtensor(int b) {
        return interactionsGtensor.get(b).get(0);
    }

    /** Ajoute un site de base et retourne son indice (le premier site vaut 0). */
    int addBasisSite(double[] position) {
        basis.add(position);
        return basis.size() - 1;
    }
}
