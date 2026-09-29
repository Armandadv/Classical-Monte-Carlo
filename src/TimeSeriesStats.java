import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

/**
 * Accumulateur de production a taille dynamique pour le runner adaptatif : avec une
 * thermalisation et une production de longueurs variables (voir
 * {@link AdaptiveThermalization}), le nombre de mesures n'est pas connu a l'avance et
 * un tableau pre-dimensionne ne convient pas.
 *
 * <p>{@link #getMeasurements(int)} renvoie les 9 observables suivantes :</p>
 * <pre>
 *   { &lt;E&gt;/N , 3 sigma_E , Var_pop(E) , &lt;Mx&gt; , 3 sigma_Mx ,
 *     &lt;My&gt; , 3 sigma_My , &lt;Mz&gt; , 3 sigma_Mz }
 * </pre>
 * <p>ou sigma est l'ecart-type d'echantillon (denominateur n-1, {@code getStandardDeviation()}
 * de commons-math) et Var_pop la variance de population ({@code getPopulationVariance()}) de
 * l'energie <b>totale</b> (c'est elle qui alimente {@link AdaptiveCoolingSchedule}).</p>
 *
 * <p>S'y ajoutent des estimateurs corriges de l'auto-correlation
 * ({@link #energyStandardError()}, {@link #energyG()}), qui sont les barres d'erreur
 * physiquement correctes pour une serie Monte Carlo (Wolff 2004) : les "3 sigma" ci-dessus
 * sont des ecarts-types de la <i>distribution</i>, pas des erreurs sur la moyenne, et
 * ignorent la correlation temporelle.</p>
 */
public final class TimeSeriesStats {

    private double[] energy;
    private double[] mx;
    private double[] my;
    private double[] mz;
    private int n;

    /** Constructeur par defaut : capacite initiale 1024. */
    public TimeSeriesStats() { this(1024); }

    /**
     * Construit les 4 series (E, Mx, My, Mz) a la capacite initiale demandee.
     *
     * @param initialCapacity capacite commune aux 4 series, plancher 16 (croissance
     *        par doublement ensuite).
     */
    public TimeSeriesStats(int initialCapacity) {
        int cap = Math.max(16, initialCapacity);
        energy = new double[cap];
        mx = new double[cap];
        my = new double[cap];
        mz = new double[cap];
        n = 0;
    }

    /**
     * Ajoute une mesure.
     *
     * @param magnetization vecteur d'aimantation (au moins 3 composantes).
     */
    public void add(double energy_, double[] magnetization) {
        if (magnetization == null || magnetization.length < 3) {
            throw new IllegalArgumentException("magnetization doit avoir au moins 3 composantes");
        }
        ensureCapacity(n + 1);
        energy[n] = energy_;
        mx[n] = magnetization[0];
        my[n] = magnetization[1];
        mz[n] = magnetization[2];
        n++;
    }

    /** Double la capacite des 4 series tant que {@code need} depasse leur taille courante. */
    private void ensureCapacity(int need) {
        if (need <= energy.length) return;
        int cap = energy.length;
        while (cap < need) cap <<= 1;
        energy = Arrays.copyOf(energy, cap);
        mx = Arrays.copyOf(mx, cap);
        my = Arrays.copyOf(my, cap);
        mz = Arrays.copyOf(mz, cap);
    }

    /** Nombre de mesures accumulees. */
    public int size() { return n; }

    /** Vide l'accumulateur (la capacite est conservee). */
    public void clear() { n = 0; }

    /**
     * Les 9 observables, dans l'ordre decrit en tete de classe.
     *
     * @param latticeSize nombre de sites (pour l'energie par site).
     */
    public double[] getMeasurements(int latticeSize) {
        DescriptiveStatistics statE = new DescriptiveStatistics();
        DescriptiveStatistics statMx = new DescriptiveStatistics();
        DescriptiveStatistics statMy = new DescriptiveStatistics();
        DescriptiveStatistics statMz = new DescriptiveStatistics();

        for (int i = 0; i < n; i++) {
            statE.addValue(energy[i]);
            statMx.addValue(mx[i]);
            statMy.addValue(my[i]);
            statMz.addValue(mz[i]);
        }

        double stdDevE = 3 * statE.getStandardDeviation();
        double meanE = statE.getMean() / latticeSize;
        double varE = statE.getPopulationVariance();

        double meanMx = statMx.getMean();
        double meanMy = statMy.getMean();
        double meanMz = statMz.getMean();
        double stdDevMx = 3 * statMx.getStandardDeviation();
        double stdDevMy = 3 * statMy.getStandardDeviation();
        double stdDevMz = 3 * statMz.getStandardDeviation();

        return new double[] { meanE, stdDevE, varE, meanMx, stdDevMx, meanMy, stdDevMy, meanMz, stdDevMz };
    }

    /**
     * Erreur standard sur &lt;E&gt; corrigee de l'auto-correlation :
     * sqrt(g Var_pop(E) / n) avec g estime par la methode Gamma de Wolff (2004).
     */
    public double energyStandardError() {
        return EquilibrationDetector.standardError(energy, 0, n);
    }

    /** Inefficacite statistique g = 1 + 2 tau (convention de Chodera) de la serie d'energie, methode de Wolff. */
    public double energyG() {
        return EquilibrationDetector.statisticalInefficiencyWolff(energy, 0, n);
    }

    /** Copie de la serie d'energie totale. */
    public double[] energySeries() { return Arrays.copyOf(energy, n); }

    /** Copie de la serie de la composante {@code k} (0 = x, 1 = y, 2 = z) de l'aimantation. */
    public double[] magnetizationSeries(int k) {
        double[] src = switch (k) {
            case 0 -> mx;
            case 1 -> my;
            case 2 -> mz;
            default -> throw new IllegalArgumentException("composante " + k);
        };
        return Arrays.copyOf(src, n);
    }

    /**
     * Variance de population de l'energie <b>totale</b> : c'est la quantite qui alimente
     * {@link AdaptiveCoolingSchedule#next(double, double, int)} (c_v = Var(E)/(N T^2)).
     * Identique a l'element d'indice 2 de {@link #getMeasurements(int)}.
     */
    public double varianceEnergyTotal() {
        return EquilibrationDetector.populationVariance(energy, 0, n);
    }

    /** Diagnostic lisible de la serie d'energie (t0, g, tau, N_eff, Geweke, ...). */
    public String energyDiagnostics() {
        return EquilibrationDetector.diagnostics(energy, 0, n);
    }
}
