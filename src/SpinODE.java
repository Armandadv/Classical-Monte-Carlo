import java.util.Arrays;

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

/**
 * Second membre de la dynamique Landau-Lifshitz pure : dS/dt = -S x H_eff (precession sans
 * amortissement, hbar = 1), H_eff regroupant echange, champ Zeeman et contribution quartique.
 * Etat plat : {@code [site*3 + (Sx, Sy, Sz)]}. L'integration RK4 a pas fixe est dans
 * {@link Rk4PasFixe} : aucun schema d'integration ici.
 */
public class SpinODE implements FirstOrderDifferentialEquations {

    private final Lattice lattice;
    private final int dimension;
    private final int[][] nbBySite;        // [site] -> voisins (référence, pas copie)
    private final double[][][] tmplBySite; // [site] -> gabarits du sous-réseau du site

    public SpinODE(Lattice lattice) {
        this.lattice = lattice;
        this.dimension = lattice.size * 3;
        this.nbBySite = new int[lattice.size][];
        this.tmplBySite = new double[lattice.size][][];
        for (int k = 0; k < lattice.size; k++) {
            this.nbBySite[k] = lattice.getNeighborSites(k);
            this.tmplBySite[k] = lattice.getInteractionMatrices(k);
        }
    }

    /** Dimension du systeme : 3 composantes par site. */
    @Override
    public int getDimension() {
        return dimension;
    }

    /** Remplit yDot avec dS/dt sur l'etat plat y ({@code t} ignore : systeme autonome). */
    @Override
    public void computeDerivatives(double t, double[] y, double[] yDot) {
        eomSpin(lattice, y, yDot);
    }

    /**
     * Dérivée par site. Les caches (nbBySite, tmplBySite) sont construits pour
     * {@code this.lattice} : passer une autre instance lèverait une exception
     * (comportement voulu — échec explicite plutôt que résultat faux).
     */
    public void eomSpin(Lattice lattice, double[] state, double[] yDot) {
        if (lattice != this.lattice) {
            throw new IllegalArgumentException("eomSpin : passer this.lattice (caches construits dessus)");
        }
        final double[] B = lattice.getInteractionField();
        final boolean quartic = lattice.hasQuartic(); // lu UNE fois par appel — hoisté par le JIT
        final double quarticB = lattice.getQuarticB();
        final int n = lattice.size;

        for (int k = 0; k < n; k++) {
            final int[] nb = nbBySite[k];
            final double[][] M = tmplBySite[k];
            final int k3 = 3 * k;
            double dSxdt = 0.0, dSydt = 0.0, dSzdt = 0.0;
            final double Sxk = state[k3], Syk = state[k3 + 1], Szk = state[k3 + 2];

            for (int j = 0; j < nb.length; j++) {
                final int index = nb[j] * 3;
                final double[] mj = M[j];
                double JkjxSj0, JkjxSj1, JkjxSj2;

                if (mj.length == 3) {
                    JkjxSj0 = mj[0] * state[index];
                    JkjxSj1 = mj[1] * state[index + 1];
                    JkjxSj2 = mj[2] * state[index + 2];
                } else if (mj.length == 2) {
                    // 2-PARAM a.I + b.A
                    final double s0 = state[index], s1 = state[index + 1], s2 = state[index + 2];
                    JkjxSj0 = mj[0] * s0 + mj[1] * (s1 + s2);
                    JkjxSj1 = mj[1] * s0 + mj[0] * s1 + mj[1] * s2;
                    JkjxSj2 = mj[1] * (s0 + s1) + mj[0] * s2;
                } else if (mj.length == 5) {
                    // Kitaev : diag(d0, d1, d2) + exactement une paire k de valeur v
                    final double s0 = state[index], s1 = state[index + 1], s2 = state[index + 2];
                    final int kp = (int) mj[3];
                    final double v = mj[4];
                    if (kp == 0) {
                        JkjxSj0 = mj[0] * s0 + v * s1;
                        JkjxSj1 = v * s0 + mj[1] * s1;
                        JkjxSj2 = mj[2] * s2;
                    } else if (kp == 1) {
                        JkjxSj0 = mj[0] * s0 + v * s2;
                        JkjxSj1 = mj[1] * s1;
                        JkjxSj2 = v * s0 + mj[2] * s2;
                    } else {
                        JkjxSj0 = mj[0] * s0;
                        JkjxSj1 = mj[1] * s1 + v * s2;
                        JkjxSj2 = v * s1 + mj[2] * s2;
                    }
                } else {
                    // Pleine symetrique : convention triangle superieur de l'ancien code
                    // (m01 = mj[1], m11 = mj[4], m12 = mj[5], m22 = mj[8]) ; mj[3], mj[6],
                    // mj[7] non lus, comme Spin.exchangeEnergySym et les noyaux Checkerboard.
                    JkjxSj0 = mj[0] * state[index] + mj[1] * state[index + 1] + mj[2] * state[index + 2];
                    JkjxSj1 = mj[1] * state[index] + mj[4] * state[index + 1] + mj[5] * state[index + 2];
                    JkjxSj2 = mj[2] * state[index] + mj[5] * state[index + 1] + mj[8] * state[index + 2];
                }

                dSxdt += JkjxSj0;
                dSydt += JkjxSj1;
                dSzdt += JkjxSj2;

                if (quartic) {
                    final int g = lattice.getQuarticGamma(k, j);
                    final double sjg = state[index + g];
                    final double contrib = 2.0 * quarticB * state[k3 + g] * sjg * sjg;
                    if (g == 0) dSxdt += contrib;
                    else if (g == 1) dSydt += contrib;
                    else dSzdt += contrib;
                }
            }

            dSxdt += B[0];
            dSydt += B[1];
            dSzdt += B[2];

            yDot[k3]     = -(Syk * dSzdt - Szk * dSydt);
            yDot[k3 + 1] = -(Szk * dSxdt - Sxk * dSzdt);
            yDot[k3 + 2] = -(Sxk * dSydt - Syk * dSxdt);
        }
    }
}
