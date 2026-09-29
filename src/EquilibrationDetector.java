/**
 * Diagnostics d'equilibration et d'auto-correlation pour series temporelles Monte Carlo.
 *
 * <p>Toutes les methodes sont statiques et operent sur un intervalle semi-ouvert
 * {@code [from, to)} d'un tableau primitif {@code double[]}. Aucune structure
 * {@code List<Double>} n'est utilisee : les boucles internes travaillent directement sur
 * le tableau et n'allouent (au plus) que des temporaires de taille O(W) ou O(n/5) (MSER-5).
 * Le cout des estimateurs est O(n.W) ou W est la fenetre d'auto-correlation retenue
 * (typiquement W &lt;&lt; n des que la serie est stationnaire).</p>
 *
 * <h2>References</h2>
 * <ul>
 *   <li><b>Chodera, J. D.</b> (2016). <i>A simple method for automated equilibration
 *       detection in molecular simulations.</i> J. Chem. Theory Comput. <b>12</b>, 1799-1805.
 *       DOI 10.1021/acs.jctc.5b00784. &mdash; On choisit l'instant de troncature t0 qui
 *       maximise le nombre d'echantillons effectifs N_eff(t0) = (T - t0) / g(t0), ou
 *       g = 1 + 2 tau est l'inefficacite statistique de la serie a partir de t0.
 *       L'estimateur rapide reproduit {@code pymbar.timeseries.statistical_inefficiency(fast=True, mintime=3)}.</li>
 *   <li><b>Wolff, U.</b> (2004). <i>Monte Carlo errors with less errors.</i>
 *       Comput. Phys. Commun. <b>156</b>, 143-153 (erratum <b>176</b>, 383). &mdash;
 *       Methode Gamma avec fenetrage automatique : on somme l'auto-correlation jusqu'a
 *       une fenetre W* determinee par une fonction de cout equilibrant biais et variance.</li>
 *   <li><b>Sokal, A. D.</b> (1997). <i>Monte Carlo methods in statistical mechanics:
 *       foundations and new algorithms.</i> Cours de Troisieme Cycle de la Physique en
 *       Suisse Romande. &mdash; Fenetrage plus ancien et equivalent : W &gt;= c tau_int(W), c ~ 6.</li>
 *   <li><b>Geweke, J.</b> (1992). <i>Evaluating the accuracy of sampling-based approaches
 *       to the calculation of posterior moments.</i> Bayesian Statistics 4, 169-193. &mdash;
 *       Test de convergence z comparant le debut et la fin de la chaine.</li>
 *   <li><b>White, K. P.</b> (1997). <i>An effective truncation heuristic for bias reduction
 *       in simulation output.</i> Simulation <b>69</b>, 323-334 (MSER / MSER-5). Voir aussi
 *       <b>Wang &amp; Glynn</b> (2016), Oper. Res. Lett. <b>44</b>, 447-451, qui montrent que
 *       MSER ne converge pas (pas de concentration) : on l'utilise donc en diagnostic
 *       <i>secondaire</i> uniquement.</li>
 *   <li><b>Gelman, A. &amp; Rubin, D. B.</b> (1992). <i>Inference from iterative simulation
 *       using multiple sequences.</i> Statist. Sci. <b>7</b>, 457-472 (statistique R-chapeau).</li>
 *   <li><b>Katzgraber, H. G.</b> (2009). <i>Introduction to Monte Carlo methods.</i>
 *       arXiv:0905.1629 &mdash; binning logarithmique : les moyennes sur les fenetres
 *       [2^k, 2^(k+1)) doivent coincider aux barres d'erreur pres lorsque la serie est
 *       equilibree (diagnostic qualitatif).</li>
 * </ul>
 *
 * <p>Motivation physique : pres d'une transition de phase le temps de correlation diverge
 * (ralentissement critique, tau ~ xi^z avec z ~ 2 pour une dynamique locale de type
 * Metropolis). Une longueur de thermalisation <i>fixe</i> est donc soit trop courte pres
 * de Tc, soit inutilement longue loin de Tc. Les criteres ci-dessous permettent d'allouer
 * le temps de calcul la ou il est reellement necessaire.</p>
 *
 * <h2>Convention de tau (importante)</h2>
 * <p>Deux definitions du temps d'auto-correlation integre coexistent dans la litterature et
 * different de exactement 1/2 :</p>
 * <pre>
 *   Chodera : tau       = somme_{t &gt;= 1} rho(t)              g = 1 + 2 tau
 *   Wolff   : tau_int^W = 1/2 + somme_{t &gt;= 1} rho(t)        g = 2 tau_int^W
 * </pre>
 * <p>Les deux donnent le <b>meme</b> g (donc la meme barre d'erreur) : seul le nom differe.
 * Cette classe expose partout la convention de <b>Chodera</b> :
 * {@link Result#tauInt()} vaut (g - 1)/2 = tau^Chodera, et
 * {@code AdaptiveThermalization.tauIntSweeps()} en derive. La variable interne de
 * {@link #statisticalInefficiencyWolff} est, elle, un tau_int^Wolff : c'est la seule
 * exception, et elle est convertie en g avant d'etre rendue publique.</p>
 *
 * <h2>Valeurs non finies : alarme, pas silence</h2>
 * <p>Une fenetre contenant un NaN ou un +-Inf donne une moyenne et une variance NaN. Le
 * test de degenerescence {@code !(var &gt; seuil)} etant volontairement NaN-safe, une telle
 * fenetre serait autrement declaree <i>constante</i> et passerait silencieusement tous les
 * criteres d'equilibration &mdash; exactement le contraire du comportement voulu. Chaque point
 * d'entree verifie donc la finitude <b>avant</b> tout calcul de moyenne/variance et signale
 * le probleme de facon visible :</p>
 * <ul>
 *   <li>{@link #detectEquilibration} : {@code Result(from, +Inf, +Inf, 0, false)}
 *       (N_eff = 0 -&gt; echoue a coup sur le critere minNeff, et {@code constant == false}) ;</li>
 *   <li>{@link #statisticalInefficiencyFast} / {@link #statisticalInefficiencyWolff} : +Inf ;</li>
 *   <li>{@link #gewekeZ} : NaN (donc |z| &lt;= zmax est faux) ;</li>
 *   <li>{@link #standardError} : NaN ;</li>
 *   <li>{@link #isConstant} : false.</li>
 * </ul>
 * <p>En amont, {@code AdaptiveThermalization.observe} leve directement une
 * {@code IllegalStateException} : une observable non finie est un bug de la simulation
 * (energie divergente, sigma explose), pas un regime physique a diagnostiquer.</p>
 */
public final class EquilibrationDetector {

    /** Parametre S_tau du fenetrage automatique de Wolff (2004). Valeur usuelle : 1.5. */
    public static final double DEFAULT_S_TAU = 1.5;

    /** Longueur minimale de la queue [t0, to) exploree par {@link #detectEquilibration}. */
    private static final int MIN_TAIL_SAMPLES = 50;

    /** {@code mintime} de pymbar : on n'autorise l'arret sur C_t &lt;= 0 que pour t &gt; 3. */
    private static final int MIN_TIME_FAST = 3;

    /** Valeur de repli pour tau lorsque tau_int^Wolff &lt;= 1/2 (arret immediat du fenetrage). */
    private static final double TINY_TAU = 1e-6;

    private EquilibrationDetector() { }

    /**
     * Resultat d'une detection d'equilibration.
     *
     * @param t0       indice de troncature dans l'indexation <b>originale</b> de {@code x}
     *                 (c'est-a-dire {@code from <= t0 < to}).
     * @param g        inefficacite statistique g = 1 + 2 tau de la queue [t0, to),
     *                 estimee par la methode Gamma de Wolff. Toujours &gt;= 1 (ou +Inf si la
     *                 fenetre contient une valeur non finie).
     * @param tauInt   temps d'auto-correlation integre <b>convention de Chodera</b> :
     *                 tau = somme_{t &gt;= 1} rho(t) = (g - 1) / 2, en unites
     *                 d'echantillons de la serie. Le nom historique du composant est conserve
     *                 pour ne pas casser les appelants ; ce n'est <i>pas</i> le tau_int de
     *                 Wolff, qui vaut tau + 1/2 (voir la section "Convention de tau" de la
     *                 javadoc de classe). Les deux donnent le meme g, donc la meme barre
     *                 d'erreur : seul le nom differe.
     * @param neff     nombre d'echantillons effectifs (to - t0) / g.
     * @param constant vrai si la serie est (numeriquement) constante : aucune information.
     *                 Toujours false pour une fenetre contenant un NaN ou un Inf.
     */
    public record Result(int t0, double g, double tauInt, double neff, boolean constant) { }

    // ------------------------------------------------------------------
    // Utilitaires numeriques
    // ------------------------------------------------------------------

    /** Moyenne arithmetique de {@code x[from, to)}; 0 si l'intervalle est vide. */
    public static double mean(double[] x, int from, int to) {
        int n = to - from;
        if (n <= 0) return 0.0;
        double s = 0.0;
        for (int i = from; i < to; i++) s += x[i];
        return s / n;
    }

    /** Variance de population (denominateur n) de {@code x[from, to)}; 0 si n &lt; 1. */
    public static double populationVariance(double[] x, int from, int to) {
        int n = to - from;
        if (n <= 0) return 0.0;
        return populationVariance(x, from, to, mean(x, from, to));
    }

    private static double populationVariance(double[] x, int from, int to, double mu) {
        int n = to - from;
        if (n <= 0) return 0.0;
        double s = 0.0;
        for (int i = from; i < to; i++) {
            double d = x[i] - mu;
            s += d * d;
        }
        return s / n;
    }

    /**
     * Test de degenerescence numerique d'une fenetre deja connue finie :
     * variance de population &lt;= 1e-300 * max(1, moyenne^2). Une fenetre degeneree ne
     * porte aucune information temporelle (g = 1, tau = 0).
     *
     * <p>La formulation negative {@code !(var > seuil)} est volontaire : elle est NaN-safe.
     * Elle n'est <b>pas</b> un test de NaN pour autant &mdash; c'est {@link #hasNonFinite},
     * appele en amont par tous les points d'entree publics, qui separe le cas
     * "constante" (aucune information) du cas "non finie" (bug), voir la javadoc de classe
     * ("NaN -&gt; alarme, pas silence").</p>
     */
    private static boolean isDegenerate(double var, double mu) {
        return !(var > 1e-300 * Math.max(1.0, mu * mu));
    }

    /**
     * Vrai si {@code x[from, to)} contient au moins une valeur non finie (NaN, +Inf, -Inf).
     *
     * <p>Cout O(n) et une seule passe sans branchement couteux : negligeable devant les
     * O(n.W) des estimateurs d'auto-correlation qui suivent.</p>
     */
    private static boolean hasNonFinite(double[] x, int from, int to) {
        for (int i = from; i < to; i++) {
            if (!Double.isFinite(x[i])) return true;
        }
        return false;
    }

    /**
     * Test de serie (numeriquement) constante : variance de population
     * &lt;= 1e-300 * max(1, moyenne^2).
     *
     * @return false des que la fenetre contient une valeur non finie : un NaN n'est pas une
     *         constante, et le declarer tel ferait passer silencieusement tous les criteres
     *         d'equilibration ("NaN -&gt; alarme, pas silence").
     */
    public static boolean isConstant(double[] x, int from, int to) {
        int n = to - from;
        if (n <= 1) return true;
        if (hasNonFinite(x, from, to)) return false;
        double mu = mean(x, from, to);
        return isDegenerate(populationVariance(x, from, to, mu), mu);
    }

    // ------------------------------------------------------------------
    // Estimateur rapide de Chodera (pymbar)
    // ------------------------------------------------------------------

    /**
     * Inefficacite statistique g = 1 + 2 tau (tau de Chodera) par l'estimateur "rapide" de
     * Chodera (2016), identique a {@code pymbar.timeseries.statistical_inefficiency}
     * avec {@code fast=True, mintime=3}.
     *
     * <p>Algorithme : on note dA = x - moyenne et sigma2 la variance de population.
     * Pour t = 1, 2, 4, 7, 11, ... (increment croissant de 1 a chaque pas) :</p>
     * <pre>
     *   C_t = moyenne(dA[0..N-t) * dA[t..N)) / sigma2
     *   si (C_t &lt;= 0 et t &gt; mintime) -&gt; arret
     *   g += 2 * C_t * (1 - t/N) * increment ;  t += increment ;  increment += 1
     * </pre>
     * <p>Le pas croissant rend le cout O(n sqrt(t_max)) au lieu de O(n t_max) : c'est
     * cet estimateur qui est utilise dans la boucle de recherche de t0, ou il est
     * appele O(n/nskip) fois.</p>
     *
     * <p>g est ici l'inefficacite de la convention de Chodera : g = 1 + 2 tau avec
     * tau = somme_{t &gt;= 1} rho(t).</p>
     *
     * @return g &gt;= 1 ; vaut exactement 1 pour une serie constante ou n &lt; 4 ;
     *         +Inf si la fenetre contient une valeur non finie.
     */
    public static double statisticalInefficiencyFast(double[] x, int from, int to) {
        if (hasNonFinite(x, from, to)) return Double.POSITIVE_INFINITY;
        int n = to - from;
        if (n < 4) return 1.0;
        double mu = mean(x, from, to);
        double sigma2 = populationVariance(x, from, to, mu);
        if (isDegenerate(sigma2, mu)) return 1.0;

        double g = 1.0;
        int t = 1;
        int increment = 1;
        while (t < n - 1) {
            int m = n - t;
            double c = 0.0;
            for (int i = 0; i < m; i++) {
                c += (x[from + i] - mu) * (x[from + i + t] - mu);
            }
            c /= ((double) m * sigma2);
            if (c <= 0.0 && t > MIN_TIME_FAST) break;
            g += 2.0 * c * (1.0 - (double) t / n) * increment;
            t += increment;
            increment += 1;
        }
        if (!Double.isFinite(g) || g < 1.0) return 1.0;
        return g;
    }

    // ------------------------------------------------------------------
    // Methode Gamma de Wolff (fenetrage automatique)
    // ------------------------------------------------------------------

    /** {@link #statisticalInefficiencyWolff(double[], int, int, double)} avec S_tau = 1.5. */
    public static double statisticalInefficiencyWolff(double[] x, int from, int to) {
        return statisticalInefficiencyWolff(x, from, to, DEFAULT_S_TAU);
    }

    /**
     * Inefficacite statistique g par la methode Gamma de Wolff (2004), avec fenetrage
     * automatique.
     *
     * <p><b>Convention.</b> La variable interne {@code tauIntW} de cet algorithme est le
     * tau_int de <i>Wolff</i>, c'est-a-dire 1/2 + somme rho(t). La valeur rendue est donc</p>
     * <pre>
     *   g = 2 tau_int^Wolff(W*) = 1 + 2 tau^Chodera
     * </pre>
     * <p>ce qui est exactement la meme inefficacite (donc la meme barre d'erreur) que celle
     * de l'estimateur rapide : les deux conventions ne different que par le 1/2 porte par le
     * nom de tau, jamais par g. Tout ce qui sort de cette classe (notamment
     * {@link Result#tauInt()}) est exprime en tau de Chodera.</p>
     *
     * <p>Soit Gamma(t) l'auto-covariance (normalisation 1/(n-t)) et
     * tau_int^W(W) = 1/2 + somme_{t=1..W} Gamma(t)/Gamma(0). On definit</p>
     * <pre>
     *   tau(W)   = S_tau / log( (2 tau_int^W(W) + 1) / (2 tau_int^W(W) - 1) )
     *   g_W      = exp(-W / tau(W)) - tau(W) / sqrt(W n)
     * </pre>
     * <p>Le premier terme majore le biais de troncature, le second l'erreur statistique.
     * On retient la premiere fenetre W telle que g_W &lt; 0 (plafond W = n/2). Lorsque
     * tau_int^W(W) &lt;= 1/2 (serie non correlee ou anti-correlee) le logarithme n'est pas
     * defini : on prend tau minuscule, ce qui declenche l'arret immediat.</p>
     *
     * <p>Cet estimateur est nettement moins biaise que l'estimateur rapide (pas de pas
     * croissant) : il est utilise pour la valeur de g <i>rapportee</i>, l'erreur standard
     * et les tests de Geweke.</p>
     *
     * @param sTau parametre de fenetrage (1.5 par defaut ; valeur non finie ou &lt;= 0 -&gt; 1.5).
     * @return g = 2 tau_int^Wolff(W*) = 1 + 2 tau^Chodera, borne inferieurement a 1 ;
     *         +Inf si la fenetre contient une valeur non finie.
     */
    public static double statisticalInefficiencyWolff(double[] x, int from, int to, double sTau) {
        if (hasNonFinite(x, from, to)) return Double.POSITIVE_INFINITY;
        int n = to - from;
        if (n < 4) return 1.0;
        double mu = mean(x, from, to);
        double gamma0 = populationVariance(x, from, to, mu);
        if (isDegenerate(gamma0, mu)) return 1.0;

        double s = (Double.isFinite(sTau) && sTau > 0.0) ? sTau : DEFAULT_S_TAU;
        int wMax = Math.max(1, n / 2);
        double[] gamma = new double[wMax + 1];     // Gamma(0..W) conservees pour la correction de biais
        gamma[0] = gamma0;
        double sumRho = 0.0;
        double tauIntW = 0.5;                       // convention de Wolff : 1/2 + somme rho
        int wOpt = 0;

        for (int w = 1; w <= wMax; w++) {
            int m = n - w;
            double c = 0.0;
            for (int i = 0; i < m; i++) {
                c += (x[from + i] - mu) * (x[from + i + w] - mu);
            }
            double gammaW = c / m;
            gamma[w] = gammaW;
            sumRho += gammaW / gamma0;
            tauIntW = 0.5 + sumRho;
            wOpt = w;
            if (!Double.isFinite(tauIntW)) { tauIntW = 0.5; wOpt = 0; break; }

            double tau;
            if (!(tauIntW > 0.5)) {
                tau = TINY_TAU;
            } else {
                double lg = Math.log((2.0 * tauIntW + 1.0) / (2.0 * tauIntW - 1.0));
                tau = (Double.isFinite(lg) && lg > 0.0) ? s / lg : TINY_TAU;
            }
            double gw = Math.exp(-w / tau) - tau / Math.sqrt((double) w * n);
            if (!(gw >= 0.0)) break;   // NaN-safe : on s'arrete aussi sur NaN
        }

        // Correction de biais a N fini (Wolff 2004, eq. 49 ; UWerr) : chaque Gamma(t) estime avec
        // la moyenne empirique est biaise vers le bas d'environ C_F/N, et le deficit s'accumule
        // sur les W lags sommes. Sans elle g est sous-estime de ~25 % pour n = 100, tau = 4.
        if (wOpt > 0) {
            double cf = gamma[0];
            for (int t = 1; t <= wOpt; t++) cf += 2.0 * gamma[t];
            cf *= 1.0 + (2.0 * wOpt + 1.0) / n;
            double corr = cf / n;
            double g0c = gamma[0] + corr;
            if (Double.isFinite(g0c) && g0c > 0.0 && Double.isFinite(corr)) {
                double sum = 0.0;
                for (int t = 1; t <= wOpt; t++) sum += (gamma[t] + corr) / g0c;
                double corrected = 0.5 + sum;
                if (Double.isFinite(corrected)) tauIntW = corrected;
            }
        }

        double g = 2.0 * tauIntW;                   // = 1 + 2 tau^Chodera
        if (!Double.isFinite(g) || g < 1.0) return 1.0;
        return g;
    }

    // ------------------------------------------------------------------
    // Detection d'equilibration (Chodera 2016)
    // ------------------------------------------------------------------

    /** {@link #detectEquilibration(double[], int, int, int)} avec nskip = max(1, n/200). */
    public static Result detectEquilibration(double[] x, int from, int to) {
        return detectEquilibration(x, from, to, Math.max(1, (to - from) / 200));
    }

    /**
     * Detection automatique de l'instant d'equilibration (Chodera 2016).
     *
     * <p>On balaie les candidats t0 = from, from + nskip, from + 2 nskip, ... tant que la
     * queue [t0, to) contient au moins {@value #MIN_TAIL_SAMPLES} echantillons ; pour chacun on estime g
     * par l'estimateur rapide et on retient le t0 qui <b>maximise</b>
     * N_eff(t0) = (to - t0)/g(t0). Ce critere realise le compromis entre le biais residuel
     * du regime transitoire (qui pousse t0 vers la droite) et la perte d'echantillons
     * (qui pousse t0 vers la gauche).</p>
     *
     * <p>Le g finalement rapporte dans le {@link Result} est <b>recalcule</b> au t0 retenu
     * par la methode Gamma de Wolff, moins biaisee que l'estimateur rapide.</p>
     *
     * @param nskip pas entre deux candidats t0 (&gt;= 1). Le cout total est
     *              O((n/nskip) . n . sqrt(tau)).
     * @return pour une fenetre contenant une valeur non finie :
     *         {@code Result(from, +Inf, +Inf, 0, false)}. N_eff = 0 echoue a coup sur le
     *         critere minNeff de l'appelant, et {@code constant} vaut false : un NaN doit
     *         declencher une alarme, pas passer pour une serie constante.
     */
    public static Result detectEquilibration(double[] x, int from, int to, int nskip) {
        int n = to - from;
        if (n <= 0) return new Result(from, 1.0, 0.0, 0.0, true);
        if (hasNonFinite(x, from, to)) {
            return new Result(from, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, 0.0, false);
        }

        double mu = mean(x, from, to);
        double var = populationVariance(x, from, to, mu);
        boolean constant = isDegenerate(var, mu);
        if (constant || n < 4) {
            return new Result(from, 1.0, 0.0, n, constant);
        }

        int step = Math.max(1, nskip);
        int minTail = MIN_TAIL_SAMPLES;

        int bestT0 = from;
        double bestNeff = -1.0;
        for (int t0 = from; to - t0 >= minTail; t0 += step) {
            double g = statisticalInefficiencyFast(x, t0, to);
            double neff = (to - t0) / g;
            if (neff > bestNeff) {
                bestNeff = neff;
                bestT0 = t0;
            }
        }

        double gW = statisticalInefficiencyWolff(x, bestT0, to);
        double neff = (to - bestT0) / gW;
        return new Result(bestT0, gW, (gW - 1.0) / 2.0, neff, false);
    }

    // ------------------------------------------------------------------
    // Test de Geweke
    // ------------------------------------------------------------------

    /** {@link #gewekeZ(double[], int, int, double, double)} avec 10 % / 50 %. */
    public static double gewekeZ(double[] x, int from, int to) {
        return gewekeZ(x, from, to, 0.1, 0.5);
    }

    /**
     * Score z de Geweke (1992) comparant le debut et la fin de la serie {@code x[from, to)}
     * (typiquement {@code from = t0} renvoye par {@link #detectEquilibration}).
     *
     * <pre>
     *   z = (moyenne_A - moyenne_B) / sqrt( var_A g_A / n_A + var_B g_B / n_B )
     * </pre>
     * <p>ou A est la fraction initiale {@code firstFrac} et B la fraction finale
     * {@code lastFrac}. Les variances sont des variances de population et g est estime
     * independamment sur chaque fenetre par la methode Gamma de Wolff (c'est la correction
     * d'auto-correlation qui rend le test valide pour une chaine Markovienne). Sous
     * l'hypothese de stationnarite z suit approximativement N(0,1) ; on rejette pour
     * |z| &gt;= 2.</p>
     *
     * @return z, ou 0 si une des deux fenetres a une variance nulle ou si le calcul
     *         produit une valeur non finie (comportement "sans-alarme") ;
     *         <b>NaN</b> si la fenetre contient une valeur non finie, de sorte que le test
     *         {@code |z| <= zmax} de l'appelant echoue ("NaN -&gt; alarme, pas silence").
     */
    public static double gewekeZ(double[] x, int from, int to, double firstFrac, double lastFrac) {
        int n = to - from;
        if (n > 0 && hasNonFinite(x, from, to)) return Double.NaN;   // NaN -> alarme, avant tout raccourci
        if (n < 8) return 0.0;
        double fA = (Double.isFinite(firstFrac) && firstFrac > 0.0) ? firstFrac : 0.1;
        double fB = (Double.isFinite(lastFrac) && lastFrac > 0.0) ? lastFrac : 0.5;

        int nA = (int) Math.max(2, Math.floor(fA * n));
        int nB = (int) Math.max(2, Math.floor(fB * n));
        nA = Math.min(nA, n);
        nB = Math.min(nB, n);

        int aFrom = from, aTo = from + nA;
        int bFrom = to - nB, bTo = to;

        double mA = mean(x, aFrom, aTo);
        double mB = mean(x, bFrom, bTo);
        double vA = populationVariance(x, aFrom, aTo, mA);
        double vB = populationVariance(x, bFrom, bTo, mB);
        if (!(vA > 0.0) || !(vB > 0.0)) return 0.0;

        double gA = statisticalInefficiencyWolff(x, aFrom, aTo);
        double gB = statisticalInefficiencyWolff(x, bFrom, bTo);

        double denom = Math.sqrt(vA * gA / nA + vB * gB / nB);
        if (!(denom > 0.0)) return 0.0;
        double z = (mA - mB) / denom;
        return Double.isFinite(z) ? z : 0.0;
    }

    // ------------------------------------------------------------------
    // MSER-5 (diagnostic secondaire)
    // ------------------------------------------------------------------

    /**
     * Heuristique MSER-5 de White (1997) : moyennes par lots de 5, puis
     * <pre>
     *   d* = argmin_d  [ somme_{i &gt;= d} (Y_i - Ybar(d))^2 ] / (n_b - d)^2
     * </pre>
     * <p>ou Y_i sont les moyennes de lots et n_b leur nombre ; on impose au moins 5 lots
     * restants. Implementation en O(n) (sommes suffixes accumulees a la volee).</p>
     *
     * <p><b>Centrage prealable.</b> Les moyennes de lots sont recentrees sur leur moyenne
     * globale avant la boucle suffixe. Le critere ne depend que des ecarts (Y_i - Ybar(d)) :
     * il est donc invariant par translation, et le centrage est mathematiquement neutre.
     * Numeriquement il est indispensable : la serie visee est une energie <b>totale</b> de
     * l'ordre de -10^3 a -10^6 avec des fluctuations O(1), si bien que la forme brute
     * {@code sq - sum^2/m} soustrait deux nombres de l'ordre de 10^12 pour en extraire un
     * resultat de l'ordre de 1 &mdash; soit une perte de ~12 chiffres significatifs sur 16, et
     * regulierement un {@code ss} negatif qu'il fallait rabattre a 0. Une fois centre,
     * {@code sq} et {@code sum^2/m} sont du meme ordre que le resultat : l'annulation
     * catastrophique et le garde-fou {@code ss < 0} disparaissent tous les deux.</p>
     *
     * <p><b>Diagnostic secondaire uniquement</b> : Wang &amp; Glynn (2016) montrent que MSER
     * ne se concentre pas autour du vrai point de troncature ; on l'utilise comme
     * verification croisee de {@link #detectEquilibration}, pas comme critere d'arret.</p>
     *
     * @return indice de troncature dans l'indexation originale, multiple de 5 a partir de
     *         {@code from} ; vaut {@code from} si la serie est trop courte (n &lt; 25).
     */
    public static int mser5(double[] x, int from, int to) {
        int n = to - from;
        if (n < 25) return from;
        int nb = n / 5;
        if (nb < 5) return from;

        double[] y = new double[nb];
        double gsum = 0.0;
        for (int b = 0; b < nb; b++) {
            int base = from + 5 * b;
            double s = x[base] + x[base + 1] + x[base + 2] + x[base + 3] + x[base + 4];
            y[b] = s / 5.0;
            gsum += y[b];
        }
        // Recentrage sur la moyenne globale des lots : neutre pour le critere (invariant par
        // translation), decisif pour la precision de sq - sum^2/m (voir javadoc).
        double gmean = gsum / nb;
        for (int b = 0; b < nb; b++) y[b] -= gmean;

        double sum = 0.0, sq = 0.0;
        int bestD = 0;
        double best = Double.POSITIVE_INFINITY;
        for (int d = nb - 1; d >= 0; d--) {
            sum += y[d];
            sq += y[d] * y[d];
            if (d > nb - 5) continue;          // il faut au moins 5 lots restants
            int m = nb - d;
            double ss = sq - sum * sum / m;    // >= 0 : les y sont centres, plus d'annulation
            double crit = ss / ((double) m * m);
            if (crit < best) {
                best = crit;
                bestD = d;
            }
        }
        return from + 5 * bestD;
    }

    // ------------------------------------------------------------------
    // Gelman-Rubin
    // ------------------------------------------------------------------

    /**
     * Statistique R-chapeau de Gelman &amp; Rubin (1992) pour M chaines de meme longueur N :
     * <pre>
     *   B = N/(M-1) somme_m (moyenne_m - moyenne_globale)^2
     *   W = moyenne des variances intra-chaine (denominateur N-1)
     *   V = (N-1)/N W + B/N
     *   R = sqrt(V/W)
     * </pre>
     * <p>R -&gt; 1 quand les chaines ont oublie leur condition initiale ; R &gt; 1.05-1.1
     * signale que les chaines n'ont pas melange (utile pour des recuits paralleles ou des
     * replicas partant de configurations differentes).</p>
     *
     * @return R, ou NaN si moins de 2 chaines, longueurs inegales/insuffisantes ou W = 0.
     */
    public static double gelmanRubin(double[][] chains) {
        if (chains == null || chains.length < 2) return Double.NaN;
        int m = chains.length;
        if (chains[0] == null) return Double.NaN;
        int n = chains[0].length;
        if (n < 2) return Double.NaN;
        for (double[] c : chains) {
            if (c == null || c.length != n) return Double.NaN;
        }

        double[] means = new double[m];
        double w = 0.0, grand = 0.0;
        for (int j = 0; j < m; j++) {
            double mu = mean(chains[j], 0, n);
            means[j] = mu;
            grand += mu;
            double sv = 0.0;
            for (int i = 0; i < n; i++) {
                double d = chains[j][i] - mu;
                sv += d * d;
            }
            w += sv / (n - 1);
        }
        grand /= m;
        w /= m;
        if (!(w > 0.0)) return Double.NaN;

        double b = 0.0;
        for (int j = 0; j < m; j++) {
            double d = means[j] - grand;
            b += d * d;
        }
        b *= (double) n / (m - 1);

        double v = (n - 1.0) / n * w + b / n;
        double r = Math.sqrt(v / w);
        return Double.isFinite(r) ? r : Double.NaN;
    }

    // ------------------------------------------------------------------
    // Binning logarithmique (Katzgraber)
    // ------------------------------------------------------------------

    /**
     * Moyennes sur les fenetres logarithmiques [2^k, 2^(k+1)) <b>relatives a {@code from}</b>
     * (derniere fenetre tronquee a {@code to}).
     *
     * <p>Note : l'echantillon d'offset 0 n'appartient a aucune fenetre (la premiere fenetre
     * est [1, 2)), conformement au schema [2^k, 2^(k+1)) ; c'est sans consequence pour un
     * diagnostic.</p>
     */
    public static double[] logBinMeans(double[] x, int from, int to) {
        int n = to - from;
        if (n <= 1) return new double[0];
        int nb = 0;
        for (int lo = 1; lo < n; lo <<= 1) nb++;
        double[] out = new double[nb];
        int k = 0;
        for (int lo = 1; lo < n; lo <<= 1) {
            int hi = Math.min(n, lo << 1);
            out[k++] = mean(x, from + lo, from + hi);
        }
        return out;
    }

    /**
     * {@link #logBinsAgree(double[], int, int, int, double)} avec lastBins = 3 et nSigma = 3.
     *
     * <p>On prend 3 sigma (et non 2) pour le raccourci : comparer 3 cases deux a deux fait
     * 3 tests simultanes, et un seuil a 2 sigma donnerait ~13 % de fausses alarmes sur une
     * serie pourtant equilibree. Utiliser la surcharge complete pour un seuil explicite.</p>
     */
    public static boolean logBinsAgree(double[] x, int from, int to) {
        return logBinsAgree(x, from, to, 3, 3.0);
    }

    /**
     * Critere de binning logarithmique (Katzgraber) : les {@code lastBins} dernieres
     * fenetres [2^k, 2^(k+1)) doivent avoir des moyennes compatibles a {@code nSigma}
     * pres, l'erreur de chaque fenetre etant l'erreur standard corrigee de
     * l'auto-correlation ({@link #standardError}).
     *
     * <p>Si une derive residuelle subsiste, les fenetres successives (qui doublent de
     * longueur, donc balaient des epoques de plus en plus tardives) donnent des moyennes
     * incompatibles.</p>
     *
     * @return true si toutes les paires sont compatibles (ou si le nombre de fenetres est
     *         insuffisant pour conclure).
     */
    public static boolean logBinsAgree(double[] x, int from, int to, int lastBins, double nSigma) {
        int n = to - from;
        if (n <= 1 || lastBins < 2) return true;
        int nb = 0;
        for (int lo = 1; lo < n; lo <<= 1) nb++;
        if (nb < lastBins) return true;

        int firstBin = nb - lastBins;
        double[] mu = new double[lastBins];
        double[] se = new double[lastBins];
        int k = 0, idx = 0;
        for (int lo = 1; lo < n; lo <<= 1, idx++) {
            if (idx < firstBin) continue;
            int hi = Math.min(n, lo << 1);
            int a = from + lo, b = from + hi;
            mu[k] = mean(x, a, b);
            se[k] = standardError(x, a, b);
            k++;
        }

        double ns = (Double.isFinite(nSigma) && nSigma > 0.0) ? nSigma : 2.0;
        for (int i = 0; i < lastBins; i++) {
            for (int j = i + 1; j < lastBins; j++) {
                double d = Math.abs(mu[i] - mu[j]);
                double sA = Double.isFinite(se[i]) ? se[i] : 0.0;
                double sB = Double.isFinite(se[j]) ? se[j] : 0.0;
                double s = Math.sqrt(sA * sA + sB * sB);
                if (!(s > 0.0)) {
                    if (d > 0.0) return false;
                    continue;
                }
                if (d > ns * s) return false;
            }
        }
        return true;
    }

    /**
     * Erreur standard de la moyenne corrigee de l'auto-correlation :
     * sqrt(g * var / n), avec var la variance de population et g l'inefficacite
     * statistique estimee par la methode Gamma de Wolff. C'est l'erreur "physique"
     * a rapporter pour une observable Monte Carlo (Wolff 2004).
     *
     * @return 0 pour une serie constante ou n &lt; 2 ; NaN si l'intervalle est vide
     *         <b>ou</b> s'il contient une valeur non finie ("NaN -&gt; alarme, pas silence").
     */
    public static double standardError(double[] x, int from, int to) {
        int n = to - from;
        if (n <= 0) return Double.NaN;
        if (hasNonFinite(x, from, to)) return Double.NaN;
        if (n < 2) return 0.0;
        double mu = mean(x, from, to);
        double var = populationVariance(x, from, to, mu);
        if (!(var > 0.0)) return 0.0;
        double g = statisticalInefficiencyWolff(x, from, to);
        double se = Math.sqrt(g * var / n);
        return Double.isFinite(se) ? se : 0.0;
    }

    /**
     * Chaine de diagnostic lisible resumant l'etat d'une serie : t0, g, tau (Chodera), N_eff,
     * z de Geweke, MSER-5 et accord des fenetres logarithmiques.
     */
    public static String diagnostics(double[] x, int from, int to) {
        Result r = detectEquilibration(x, from, to);
        double z = gewekeZ(x, r.t0(), to);
        int mser = mser5(x, from, to);
        boolean bins = logBinsAgree(x, r.t0(), to);
        return String.format(
                "n=%d t0=%d g=%.3f tau=%.3f neff=%.1f geweke=%.2f mser5=%d logbins=%s%s",
                to - from, r.t0(), r.g(), r.tauInt(), r.neff(), z, mser,
                bins ? "ok" : "DESACCORD", r.constant() ? " [constante]" : "");
    }
}
