/**
 * Detecteur de stagnation ("je ne peux plus ameliorer l'etat courant") pour une chaine de
 * Monte Carlo Metropolis, destine a remplacer les plafonds fixes
 * ({@code maxSweeps = 200 000}, {@code maxProductionSweeps = 400 000}) par un arret motive
 * quand poursuivre ne peut statistiquement plus rien apporter. Voir spec_stagnation.md,
 * section "Criteres d'arret <i>je ne peux plus ameliorer</i>".
 *
 * <p>Trois criteres independants, chacun expose comme primitive statique testable
 * isolement :</p>
 * <ul>
 *   <li><b>Derive invisible</b> ({@link #drift}) : regression lineaire OLS de l'observable
 *       sur l'indice, sur la fenetre d'analyse. L'erreur de la pente est gonflee par
 *       l'inefficacite statistique g de Wolff ({@link EquilibrationDetector#statisticalInefficiencyWolff})
 *       parce qu'un estimateur lineaire sur une serie autocorrelee a une variance plus
 *       grande que dans le cas i.i.d. d'un facteur g : g est le nombre d'echantillons
 *       correles par echantillon independant, donc l'erreur (ecart-type) d'un estimateur
 *       lineaire (combinaison lineaire des x_i, tout comme une moyenne empirique) est
 *       gonflee de sqrt(g) par rapport au cas i.i.d. &mdash; c'est exactement le meme
 *       raisonnement que celui de {@link EquilibrationDetector#standardError}.</li>
 *   <li><b>Derive negligeable en pratique</b> ({@link #negligibleDrift}) : meme en
 *       depensant tout le budget de sweeps restant, l'observable ne bougerait pas de plus
 *       d'une erreur standard &mdash; continuer ne peut donc pas ameliorer le resultat de
 *       facon detectable.</li>
 *   <li><b>Saturation de N_eff</b> ({@link #neffSaturated}) : si l'inefficacite statistique
 *       g croit proportionnellement au nombre d'echantillons (regime critique, chaine
 *       quasi-gelee), N_eff = n/g plafonne et etendre la chaine n'ameliore plus la
 *       precision statistique.</li>
 * </ul>
 *
 * <p>A cela s'ajoute la memoire de configuration : {@link SpinOverlap#analyse} indique si
 * le recouvrement q(t) avec la configuration heritee a atteint son plateau, et si ce
 * plateau est proche de 1 (chaine gelee, Ogielski 1985 &mdash; voir la javadoc de
 * {@link SpinOverlap}).</p>
 *
 * <p><b>Hysteresis.</b> {@link #update} ne declare un verdict {@code STAGNANT_*} qu'apres
 * {@link Config#minStagnantChecks} verifications consecutives ou le candidat est stagnant :
 * un palier transitoire (fluctuation qui ralentit temporairement la derive mesuree, ou
 * franchissement fortuit du seuil de plateau) ne doit pas etre confondu avec un vrai
 * blocage. Le risque physique dominant est d'arreter une thermalisation qui progressait
 * encore, jamais d'en prolonger une qui a fini : en cas de doute, le detecteur doit donc
 * pencher du cote "continuer".</p>
 *
 * <p><b>Valeurs non finies.</b> {@link #update} leve une {@link IllegalStateException} des
 * qu'un NaN ou un +-Inf apparait dans l'une des series fournies, exactement comme
 * {@code AdaptiveThermalization.observe} : une observable non finie signale un bug de la
 * simulation (energie divergente), pas un regime de stagnation.</p>
 */
public final class StagnationDetector {

    /**
     * Parametres du detecteur. POJO mutable avec setters fluides (meme style que
     * {@code AdaptiveThermalization.Config}).
     */
    public static final class Config {
        /** |z| en dessous duquel la derive OLS est statistiquement invisible. */
        public double driftZmax = 1.0;
        /** Derive extrapolee sur le budget restant, en unites d'erreur standard. */
        public double driftTolerance = 1.0;
        /** Echantillons minimum dans la fenetre d'analyse ; en dessous : INSUFFICIENT_DATA. */
        public int minWindow = 40;
        /** Fraction de la serie (la plus recente) retenue comme fenetre d'analyse. */
        public double windowFraction = 0.5;
        /** q_inf au-dela duquel l'etat est declare gele (STAGNANT_FROZEN). */
        public double frozenOverlap = 0.90;
        /** Fraction de plateau passee a {@link SpinOverlap#analyse}. */
        public double plateauFraction = 0.5;
        /** Gain relatif de N_eff en dessous duquel les extensions de chaine saturent. */
        public double neffGainFraction = 0.15;
        /** Nombre de verdicts stagnants consecutifs requis avant STAGNANT_* (hysteresis). */
        public int minStagnantChecks = 2;

        public Config driftZmax(double v) { this.driftZmax = v; return this; }
        public Config driftTolerance(double v) { this.driftTolerance = v; return this; }
        public Config minWindow(int v) { this.minWindow = v; return this; }
        public Config windowFraction(double v) { this.windowFraction = v; return this; }
        public Config frozenOverlap(double v) { this.frozenOverlap = v; return this; }
        public Config plateauFraction(double v) { this.plateauFraction = v; return this; }
        public Config neffGainFraction(double v) { this.neffGainFraction = v; return this; }
        public Config minStagnantChecks(int v) { this.minStagnantChecks = v; return this; }

        /** Copie profonde (le POJO ne contient que des scalaires). */
        public Config copy() {
            Config c = new Config();
            c.driftZmax = driftZmax;
            c.driftTolerance = driftTolerance;
            c.minWindow = minWindow;
            c.windowFraction = windowFraction;
            c.frozenOverlap = frozenOverlap;
            c.plateauFraction = plateauFraction;
            c.neffGainFraction = neffGainFraction;
            c.minStagnantChecks = minStagnantChecks;
            return c;
        }
    }

    /** Verdict de {@link #update}. */
    public enum Verdict { INSUFFICIENT_DATA, PROGRESSING, STAGNANT_PLATEAU, STAGNANT_FROZEN }

    /**
     * Resultat d'une regression lineaire OLS de {@code x} sur l'indice, sur {@code [from, to)}.
     *
     * @param slopePerSample pente b (unite de x par echantillon).
     * @param totalChange    variation extrapolee sur toute la fenetre : {@code b * (to - from - 1)}.
     * @param stdErr         erreur standard de {@code totalChange}, inefficacite statistique
     *                       de Wolff incluse.
     * @param z              {@code totalChange / stdErr} (0 si les deux sont nuls,
     *                       +-infini si {@code stdErr} est nul mais pas {@code totalChange}).
     * @param flat           {@code |z| <= driftZmax} pour le {@code driftZmax} utilise.
     * @param from           borne inferieure (incluse) de la fenetre analysee.
     * @param to             borne superieure (exclue) de la fenetre analysee.
     */
    public record Drift(double slopePerSample, double totalChange, double stdErr, double z,
                         boolean flat, int from, int to) { }

    /** Prevision de sweeps necessaires pour atteindre un N_eff cible, a g constant. */
    public record Forecast(double sweepsToTarget, boolean reachable) { }

    // ==================================================================
    // Primitives statiques, sans etat
    // ==================================================================

    /** {@link #drift(double[], int, int, double)} avec {@code driftZmax = 1.0}. */
    public static Drift drift(double[] x, int from, int to) {
        return drift(x, from, to, 1.0);
    }

    /**
     * Regression lineaire OLS de {@code x[from, to)} sur l'indice relatif
     * {@code i = 0, ..., (to-from-1)}, avec erreur de pente gonflee par l'inefficacite
     * statistique de Wolff (voir la javadoc de classe pour la justification de sqrt(g)).
     *
     * <pre>
     *   b        = Sxy / Sxx                                  (pente OLS standard)
     *   sigma^2  = SS_res / (n - 2)                            (variance residuelle)
     *   Var(b)   = g . sigma^2 / Sxx                            (g = inefficacite de Wolff)
     *   totalChange = b . (n - 1)
     *   stdErr      = sqrt(Var(b)) . (n - 1)
     *   z           = totalChange / stdErr
     * </pre>
     *
     * @return {@code Drift(0,0,0,0,true,from,to)} si {@code to - from < 3} ou si la serie
     *         est (numeriquement) constante sur la fenetre : aucune pente n'y est
     *         estimable, et l'absence d'information doit se lire comme "pas de derive
     *         detectable", pas comme une erreur.
     */
    public static Drift drift(double[] x, int from, int to, double driftZmax) {
        int n = to - from;
        if (n < 3 || EquilibrationDetector.isConstant(x, from, to)) {
            return new Drift(0.0, 0.0, 0.0, 0.0, true, from, to);
        }

        double meanX = EquilibrationDetector.mean(x, from, to);
        double meanI = (n - 1) / 2.0;
        double sxx = 0.0, sxy = 0.0;
        for (int i = 0; i < n; i++) {
            double di = i - meanI;
            sxx += di * di;
            sxy += di * (x[from + i] - meanX);
        }
        double b = sxy / sxx;
        double a = meanX - b * meanI;

        // g doit etre estime sur les RESIDUS de la regression, pas sur la serie brute x.
        // x contient la tendance (c'est justement ce qu'on cherche a detecter) : une
        // tendance lineaire, vue par la methode Gamma de Wolff comme une fonction
        // d'auto-covariance qui decroit tres lentement (elle ne s'annule qu'au bout de
        // l'ordre de n), est interpretee comme une auto-correlation enorme -- tau_int de
        // l'ordre de n, donc g qui explose avec n au lieu de converger vers une valeur finie
        // caracteristique du bruit. stdErr (qui contient sqrt(g)) explose alors avec elle, et
        // |z| = |totalChange| / stdErr s'effondre : plus la derive est nette et longue, plus
        // "flat" a des chances de devenir vrai a tort. C'est le sens dangereux de l'erreur
        // (faux positif de stagnation sur une chaine qui thermalise encore). La quantite
        // statistiquement correcte est l'auto-correlation du bruit residuel r_i =
        // x[from+i] - (a + b.i), *apres* avoir retire la tendance estimee -- c'est la
        // correction standard pour l'erreur d'un estimateur lineaire sur une serie
        // autocorrelee (cf. Newey-West ; ici on reutilise directement la methode Gamma de
        // Wolff (2004), mais appliquee aux residus plutot qu'a la serie brute).
        double[] resid = new double[n];
        double ssRes = 0.0;
        for (int i = 0; i < n; i++) {
            double r = x[from + i] - (a + b * i);
            resid[i] = r;
            ssRes += r * r;
        }
        double sigma2 = ssRes / (n - 2);
        double g = EquilibrationDetector.statisticalInefficiencyWolff(resid, 0, n);
        double varB = g * sigma2 / sxx;

        double lengthSpan = n - 1;
        double totalChange = b * lengthSpan;
        double stdErr = Math.sqrt(varB) * Math.abs(lengthSpan);

        double z;
        if (Double.isFinite(stdErr) && stdErr > 0.0) {
            z = totalChange / stdErr;
        } else if (totalChange == 0.0) {
            z = 0.0;
        } else {
            // stdErr nul (ou non fini) mais un changement reel est mesure : une derive
            // exacte, sans aucun bruit residuel estimable, est par definition infiniment
            // significative -- jamais "invisible".
            z = (totalChange > 0.0) ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
        }
        boolean flat = Double.isFinite(z) && Math.abs(z) <= driftZmax;
        return new Drift(b, totalChange, stdErr, z, flat, from, to);
    }

    /**
     * Vrai si, meme en depensant tout le {@code budgetSamples} restant au rythme de derive
     * mesure, l'observable bougerait de moins de {@code tolerance} erreurs standard :
     * <pre>  |slopePerSample| . budgetSamples &lt;= tolerance . standardError(x, from, to)  </pre>
     * Sens physique : continuer ne peut alors pas ameliorer le resultat de facon
     * detectable, meme sans regarder si la derive est deja "invisible" sur la fenetre
     * ecoulee ({@link #drift}) &mdash; c'est un critere complementaire, pas equivalent : une
     * derive lente mais tres significative statistiquement (grande fenetre, faible bruit)
     * peut rester en dessous du seuil pratique si le budget restant est court.
     */
    public static boolean negligibleDrift(double[] x, int from, int to, double budgetSamples, double tolerance) {
        double slope = drift(x, from, to).slopePerSample();
        double se = EquilibrationDetector.standardError(x, from, to);
        double projected = Math.abs(slope) * budgetSamples;
        double allowance = tolerance * se;
        return Double.isFinite(projected) && Double.isFinite(allowance) && projected <= allowance;
    }

    /**
     * Vrai si l'extension de chaine qui a fait passer N_eff de {@code prev} a {@code cur}
     * n'a apporte qu'un gain relatif inferieur a {@code gainFraction} :
     * {@code (cur - prev) <= gainFraction . prev}.
     *
     * @return false si {@code prev <= 0} (pas de base de comparaison valable).
     */
    public static boolean neffSaturated(double neffPrevious, double neffCurrent, double gainFraction) {
        if (!(neffPrevious > 0.0)) return false;
        return (neffCurrent - neffPrevious) <= gainFraction * neffPrevious;
    }

    /**
     * Extrapole le nombre de sweeps necessaires pour atteindre {@code target} echantillons
     * effectifs, en supposant que N_eff croit lineairement avec le nombre de sweeps a
     * inefficacite statistique g fixee : {@code sweepsToTarget = sweepsDone . target / neffNow}.
     *
     * @return {@code sweepsToTarget = +Inf, reachable = false} si {@code neffNow <= 0}
     *         (aucune extrapolation possible sans N_eff de reference).
     */
    public static Forecast forecastNeff(double neffNow, long sweepsDone, double target, long budgetSweeps) {
        if (!(neffNow > 0.0)) return new Forecast(Double.POSITIVE_INFINITY, false);
        double sweepsToTarget = sweepsDone * target / neffNow;
        boolean reachable = sweepsToTarget <= budgetSweeps;
        return new Forecast(sweepsToTarget, reachable);
    }

    // ==================================================================
    // Detecteur avec hysteresis (etat de l'instance)
    // ==================================================================

    private final Config cfg;

    private Verdict verdict = Verdict.INSUFFICIENT_DATA;
    private Drift lastDrift = new Drift(0.0, 0.0, 0.0, 0.0, true, 0, 0);
    private SpinOverlap.Decay lastDecay;
    private SpinOverlap.Decay lastFullDecay;
    private int consecutiveStagnant;

    /**
     * @param cfg configuration ; {@code null} -&gt; {@code new Config()}. Une copie
     *            defensive est conservee.
     * @throws IllegalArgumentException si un parametre n'a pas de sens (valeurs qui
     *         doivent etre strictement positives fournies &lt;= 0, {@code minWindow < 1},
     *         {@code minStagnantChecks < 1}, {@code frozenOverlap > 1} ou
     *         {@code plateauFraction > 1}).
     */
    public StagnationDetector(Config cfg) {
        Config c = (cfg == null) ? new Config() : cfg.copy();
        requirePositive(c.driftZmax, "driftZmax");
        requirePositive(c.driftTolerance, "driftTolerance");
        if (c.minWindow < 1) {
            throw new IllegalArgumentException("minWindow doit etre >= 1, recu " + c.minWindow);
        }
        requirePositive(c.windowFraction, "windowFraction");
        // q est un recouvrement de spins unitaires : q_inf <= 1 toujours. Un frozenOverlap
        // > 1 rendrait STAGNANT_FROZEN silencieusement inatteignable (decay.qInf() >= seuil
        // ne serait jamais vrai) -- mieux vaut le refuser explicitement qu'echouer en silence.
        requireInUnitInterval(c.frozenOverlap, "frozenOverlap");
        // plateauFraction est passe tel quel a SpinOverlap.analyse ; on le valide ici plutot
        // que de le corriger silencieusement, contrairement a SpinOverlap.analyse (primitive
        // statique) qui, elle, retombe sur 0.5 par defaut sans lever d'exception.
        requireInUnitInterval(c.plateauFraction, "plateauFraction");
        requirePositive(c.neffGainFraction, "neffGainFraction");
        if (c.minStagnantChecks < 1) {
            throw new IllegalArgumentException("minStagnantChecks doit etre >= 1, recu " + c.minStagnantChecks);
        }
        this.cfg = c;
    }

    private static void requirePositive(double v, String name) {
        if (!(v > 0.0)) {
            throw new IllegalArgumentException(name + " doit etre > 0, recu " + v);
        }
    }

    /** Valide {@code v} dans (0, 1] : refuse aussi bien les valeurs <= 0 que > 1. */
    private static void requireInUnitInterval(double v, String name) {
        if (!(v > 0.0) || v > 1.0) {
            throw new IllegalArgumentException(name + " doit etre dans (0, 1], recu " + v);
        }
    }

    /** Remet le detecteur a l'etat initial (verdict INSUFFICIENT_DATA, compteur a 0). */
    public void reset() {
        verdict = Verdict.INSUFFICIENT_DATA;
        lastDrift = new Drift(0.0, 0.0, 0.0, 0.0, true, 0, 0);
        lastDecay = null;
        lastFullDecay = null;
        consecutiveStagnant = 0;
    }

    /**
     * Analyse la fenetre la plus recente des series fournies et met a jour le verdict.
     *
     * <p>La fenetre d'analyse est constituee des
     * {@code max(minWindow, round(windowFraction * n))} derniers echantillons d'energie
     * (bornee par {@code n = nEnergy}). En dessous de {@code minWindow} echantillons,
     * verdict = {@code INSUFFICIENT_DATA} et le compteur d'hysteresis est remis a zero.</p>
     *
     * <p>Si {@code nOverlap > 0}, deux analyses distinctes de {@link SpinOverlap#analyse}
     * sont calculees et exposees separement, car elles repondent a deux questions
     * differentes :</p>
     * <ul>
     *   <li>{@link #lastDecay()} : la meme construction de fenetre (bornee par
     *       {@code nOverlap}) qu'utilisee pour la derive d'energie, c'est-a-dire seulement
     *       les echantillons les plus recents. Question posee : "q bouge-t-il encore
     *       <b>maintenant</b> ?" -- c'est elle qui alimente le verdict (candidat stagnant,
     *       {@code STAGNANT_FROZEN}) ; juger la memoire de configuration sur des echantillons
     *       trop anciens n'aurait pas de sens pour cette decision.</li>
     *   <li>{@link #lastFullDecay()} : {@code SpinOverlap.analyse(overlap, 0, nOverlap, ...)},
     *       sur la serie <b>complete</b> depuis le debut de la temperature. Question posee :
     *       "combien de sweeps a-t-il fallu pour oublier l'etat herite ?" -- c'est tau_q
     *       mesure depuis q(0) = 1 qui donne le critere externe "deux temperatures
     *       decorrelees" (thermSweeps &gt;= K . tau_q, voir spec_stagnation.md). La fenetre
     *       tardive de {@code lastDecay()} ne contient plus la decroissance depuis 1 : y
     *       chercher tau_q y vaudrait ~0 par construction, ce qui est inutilisable pour ce
     *       cablage externe. Le cout supplementaire est O(nOverlap), negligeable devant la
     *       recherche de t0 de Chodera sur le canal energie.</li>
     * </ul>
     *
     * <p>Un candidat est <b>stagnant</b> si la derive d'energie est invisible ou
     * pratiquement negligeable sur le budget restant, <b>et</b> (canal overlap absent, ou
     * le recouvrement a atteint son plateau). Le compteur d'hysteresis n'avance que sur des
     * candidats stagnants <b>consecutifs</b> ; le verdict ne devient {@code STAGNANT_*} que
     * lorsque ce compteur atteint {@link Config#minStagnantChecks}, auquel cas
     * {@code STAGNANT_FROZEN} est rendu si un canal overlap est present et que son plateau
     * {@code q_inf >= frozenOverlap}, {@code STAGNANT_PLATEAU} sinon.</p>
     *
     * @param energy      serie d'energie (ou toute observable scalaire pertinente).
     * @param nEnergy     nombre d'echantillons valides dans {@code energy} (&gt;= 0).
     * @param overlap     serie de recouvrement q(t) ; ignoree si {@code nOverlap <= 0}
     *                    (peut alors etre {@code null}).
     * @param nOverlap    nombre d'echantillons valides dans {@code overlap} ; {@code <= 0}
     *                    signifie "canal absent".
     * @param budgetSamples nombre d'echantillons d'energie restants dans le budget alloue,
     *                    passe a {@link #negligibleDrift}.
     * @return le nouveau {@link #verdict()}.
     * @throws IllegalArgumentException si {@code nEnergy < 0}, si {@code energy} est trop
     *         court pour {@code nEnergy}, ou si {@code nOverlap > 0} et {@code overlap} est
     *         trop court pour {@code nOverlap}.
     * @throws IllegalStateException si une valeur non finie (NaN ou +-Inf) est presente
     *         dans les {@code nEnergy} (resp. {@code nOverlap}) premiers echantillons.
     */
    public Verdict update(double[] energy, int nEnergy, double[] overlap, int nOverlap, double budgetSamples) {
        if (nEnergy < 0) {
            throw new IllegalArgumentException("nEnergy doit etre >= 0, recu " + nEnergy);
        }
        if (nEnergy > 0 && (energy == null || energy.length < nEnergy)) {
            throw new IllegalArgumentException("energy trop court pour nEnergy=" + nEnergy);
        }
        boolean hasOverlap = nOverlap > 0;
        if (hasOverlap && (overlap == null || overlap.length < nOverlap)) {
            throw new IllegalArgumentException("overlap trop court pour nOverlap=" + nOverlap);
        }
        requireFinite(energy, nEnergy, "energy");
        if (hasOverlap) requireFinite(overlap, nOverlap, "overlap");

        int n = nEnergy;
        if (n < cfg.minWindow) {
            verdict = Verdict.INSUFFICIENT_DATA;
            consecutiveStagnant = 0;
            return verdict;
        }

        int[] win = analysisWindow(n);
        int from = win[0], to = win[1];
        Drift d = drift(energy, from, to, cfg.driftZmax);
        lastDrift = d;

        boolean energyFlat = d.flat() || negligibleDrift(energy, from, to, budgetSamples, cfg.driftTolerance);

        SpinOverlap.Decay decay = null;
        boolean overlapOk = true;
        if (hasOverlap) {
            int[] owin = analysisWindow(nOverlap);
            decay = SpinOverlap.analyse(overlap, owin[0], owin[1], cfg.plateauFraction);
            overlapOk = decay.plateaued();
            // Serie complete depuis le debut de la temperature (q(0) = 1) : c'est elle qui
            // donne tau_q pour le critere externe "deux temperatures decorrelees", voir la
            // javadoc de cette methode. Cout O(nOverlap), sans lien avec le verdict rendu ici.
            lastFullDecay = SpinOverlap.analyse(overlap, 0, nOverlap, cfg.plateauFraction);
        } else {
            lastFullDecay = null;
        }
        lastDecay = decay;

        boolean candidate = energyFlat && overlapOk;
        consecutiveStagnant = candidate ? consecutiveStagnant + 1 : 0;

        if (consecutiveStagnant >= cfg.minStagnantChecks) {
            boolean frozen = hasOverlap && decay.qInf() >= cfg.frozenOverlap;
            verdict = frozen ? Verdict.STAGNANT_FROZEN : Verdict.STAGNANT_PLATEAU;
        } else {
            verdict = Verdict.PROGRESSING;
        }
        return verdict;
    }

    /** Fenetre {@code [from, to)} des {@code max(minWindow, round(windowFraction*n))} derniers points. */
    private int[] analysisWindow(int n) {
        int len = Math.max(cfg.minWindow, (int) Math.round(cfg.windowFraction * n));
        if (len > n) len = n;
        if (len < 1) len = Math.min(1, n);
        return new int[] { n - len, n };
    }

    private static void requireFinite(double[] x, int n, String name) {
        for (int i = 0; i < n; i++) {
            if (!Double.isFinite(x[i])) {
                throw new IllegalStateException(
                        "valeur non finie dans " + name + "[" + i + "] = " + x[i]);
            }
        }
    }

    /** Dernier verdict calcule par {@link #update} (INSUFFICIENT_DATA avant le premier appel). */
    public Verdict verdict() { return verdict; }

    /** Vrai si {@link #verdict()} est {@code STAGNANT_PLATEAU} ou {@code STAGNANT_FROZEN}. */
    public boolean stagnant() {
        return verdict == Verdict.STAGNANT_PLATEAU || verdict == Verdict.STAGNANT_FROZEN;
    }

    /** Derniere regression de derive calculee (fenetre [0,0), flat=true, avant le premier appel utile). */
    public Drift lastDrift() { return lastDrift; }

    /**
     * Dernier resultat de {@link SpinOverlap#analyse} sur la fenetre d'analyse recente
     * (memes bornes que {@link #lastDrift()}) : "q bouge-t-il encore maintenant ?" -- c'est
     * celui qui alimente le verdict. {@code null} si aucun canal overlap n'a ete fourni.
     */
    public SpinOverlap.Decay lastDecay() { return lastDecay; }

    /**
     * Dernier resultat de {@link SpinOverlap#analyse} sur la serie d'overlap <b>complete</b>
     * depuis le debut de la temperature (q(0) = 1) : "combien de sweeps a-t-il fallu pour
     * oublier l'etat herite ?" -- c'est cette valeur (pas {@link #lastDecay()}) qui donne
     * tau_q pour le critere externe de decorrelation entre deux temperatures. {@code null}
     * si aucun canal overlap n'a ete fourni.
     */
    public SpinOverlap.Decay lastFullDecay() { return lastFullDecay; }

    /** Nombre de verdicts candidats-stagnants consecutifs accumules (hysteresis). */
    public int consecutiveStagnant() { return consecutiveStagnant; }

    /** Copie defensive de la configuration effective. */
    public Config config() { return cfg.copy(); }

    /** Resume lisible sur une ligne, destine aux journaux de simulation. */
    public String report() {
        StringBuilder sb = new StringBuilder(160);
        sb.append("StagnationDetector[verdict=").append(verdict)
          .append(" consecutif=").append(consecutiveStagnant).append('/').append(cfg.minStagnantChecks)
          .append(String.format(" drift(z=%.2f flat=%b)", lastDrift.z(), lastDrift.flat()));
        if (lastDecay != null) {
            sb.append(String.format(" overlap(qInf=%.3f tau=%.1f plateau=%b)",
                    lastDecay.qInf(), lastDecay.tauSamples(), lastDecay.plateaued()));
        } else {
            sb.append(" overlap=absent");
        }
        sb.append(']');
        return sb.toString();
    }
}
