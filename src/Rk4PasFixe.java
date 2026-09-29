/**
 * Intégrateur RK4 à pas fixe, écrit à la main : 4 évaluations de dérivée par pas,
 * sans gestionnaire d'événements ni interpolateur (le framework commons-math coûte
 * ~30 % du temps sur une dérivée de ce coût). L'ordre des opérations par composante
 * réplique commons-math ClassicalRungeKutta (mêmes coefficients, mêmes sommes
 * séquentielles) → résultats bit-exacts.
 *
 * Précondition documentée : bit-exact si (tf - t0) est un multiple exact du pas dt
 * de l'appelant (cas de tous les appelants actuels). Sinon, la répartition des pas
 * diffère de commons-math (ceil + dernier pas clampé vs pas uniforme) → écarts à
 * ~1e-16 sur l'état final.
 */
public class Rk4PasFixe {

    /** dS/dt = f(t, S). */
    @FunctionalInterface
    public interface Derivees { void compute(double t, double[] s, double[] dS); }

    /** Sondes appelées après chaque pas ; retourner true pour arrêter l'intégration. */
    @FunctionalInterface
    public interface Sonde { boolean pas(double t, double[] s, int stepIdx); }

    private final Derivees derivees;

    public Rk4PasFixe(Derivees derivees) { this.derivees = derivees; }

    /**
     * Intègre y de t0 à tf au plus nSteps pas fixes (h = (tf - t0) / nSteps).
     * Modifie y en place (l'état final est écrit dans y). Retourne le nombre de pas
     * effectués (≤ nSteps si la sonde a demandé l'arrêt).
     */
    public int integre(double t0, double tf, int nSteps, double[] y, Sonde sonde) {
        final double h = (tf - t0) / nSteps;
        final int n = y.length;
        final double[] k1 = new double[n], k2 = new double[n], k3 = new double[n],
                k4 = new double[n], tmp = new double[n];
        double t = t0;
        int nPasFait = 0;
        for (int step = 0; step < nSteps; step++) {
            derivees.compute(t, y, k1);
            for (int i = 0; i < n; i++) tmp[i] = y[i] + h * 0.5 * k1[i];
            derivees.compute(t + 0.5 * h, tmp, k2);
            for (int i = 0; i < n; i++) tmp[i] = y[i] + h * 0.5 * k2[i];
            derivees.compute(t + 0.5 * h, tmp, k3);
            for (int i = 0; i < n; i++) tmp[i] = y[i] + h * k3[i];
            derivees.compute(t + h, tmp, k4);
            for (int i = 0; i < n; i++) {
                double somme = (1.0 / 6.0) * k1[i];
                somme += (1.0 / 3.0) * k2[i];
                somme += (1.0 / 3.0) * k3[i];
                somme += (1.0 / 6.0) * k4[i];
                y[i] += h * somme;
            }
            t += h;
            nPasFait = step + 1;
            if (sonde != null && sonde.pas(t, y, step)) break;
        }
        return nPasFait;
    }
}
