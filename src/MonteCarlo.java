
import java.util.random.RandomGenerator;
import java.util.random.RandomGeneratorFactory;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

/**
 * Orchestrateur du recuit adaptatif : par temperature, thermalisation pilotee
 * ({@link AdaptiveThermalization}) puis production prolongee jusqu'a un N_eff cible ; bloc
 * S(Q,w) a la derniere temperature et precession de fin. Voir README_PARALLEL.md §5-7
 * (§10 pour le diagnostic de stagnation).
 */
public class MonteCarlo {

    // Historiques du sequentiel : thermalizationSweeps et thermalizationRate ne servent plus
    // qu'au calcul de nSQW (le pilotage de la thermalisation est desormais dans
    // AdaptiveThermalization) ; measurementRateSQW reste la cadence des instantanes du bloc S(Q,w).
    final int thermalizationSweeps;
    final int measurementRateSQW;
    final int thermalizationRate;
    private int measurementIdxSQW;
    private int MCIdx;
    final private double initTemp;
    final private double endTemp;
    /** Lignes d'observables, une par temperature : {@link #N_COLUMNS} colonnes (cf. COL_*). */
    List<double[]> observablesData = new ArrayList<>();
    private int nSQW;
    private float dTemp;
    private final int H,K;
    private DynamicStructureFactor sqw;

    private final RandomGenerator random = RandomGenerator.of("Xoshiro256PlusPlus");
    /** Garde du tirage unique de la configuration initiale (cf. runParallel). */
    boolean isInitialized;

    /**
     * Valeurs historiques du sequentiel : recuit geometrique de T = 10 a 0.001 (x0.95),
     * instantane S(Q,w) tous les {@code measurementRateSQW} = 1000 balayages (d'ou
     * {@code nSQW = 40}), dimensions honeycomb H = K = 3 des facteurs de structure.
     * {@link #runParallel} outrepasse ces bornes via {@link ParallelOptions}.
     */
    public MonteCarlo() {
        this.measurementRateSQW = 1000;
        this.thermalizationSweeps = 40_000;
        this.thermalizationRate = 80_000;
        this.measurementIdxSQW = 0;
        this.MCIdx = 0;
        this.initTemp = 10.d;
        this.endTemp = 0.001d;
        this.nSQW = (int) (thermalizationRate - thermalizationSweeps) / measurementRateSQW ;
        this.dTemp = 0.95f;
        this.isInitialized = false;
        this.H = 3;
        this.K= 3;
    }

    // ------------------------------------------------------------------ colonnes de sortie

    /*
     * Indices des colonnes de observablesData. Les 16 premieres constituent le layout
     * historique des seize colonnes d'observables (fige) ; les suivantes sont les diagnostics
     * du runner, reprises en clair dans EquilibrationLog.
     */

    /** Colonne 0 : index 1-based de la temperature. */
    public static final int COL_MCIDX = 0;
    /** Colonne 1 : taux d'acceptation moyen. */
    public static final int COL_ACC = 1;
    /** Colonne 2 : sigma du mouvement gaussien en fin de temperature. */
    public static final int COL_SIGMA = 2;
    /** Colonne 3 : temperature. */
    public static final int COL_T = 3;
    /** Colonne 4 : H_x (base de Kitaev). */
    public static final int COL_HX = 4;
    /** Colonne 5 : H_y. */
    public static final int COL_HY = 5;
    /** Colonne 6 : H_z. */
    public static final int COL_HZ = 6;
    /** Colonne 7 : &lt;E&gt;/N. */
    public static final int COL_E = 7;
    /** Colonne 8 : 3 sigma_E (ecart-type de distribution, pas une erreur sur la moyenne). */
    public static final int COL_STDE = 8;
    /** Colonne 9 : variance de population de l'energie <b>totale</b>. */
    public static final int COL_VARE = 9;
    /** Colonne 10 : &lt;M_x&gt;. */
    public static final int COL_MX = 10;
    /** Colonne 11 : 3 sigma_Mx. */
    public static final int COL_STDMX = 11;
    /** Colonne 12 : &lt;M_y&gt;. */
    public static final int COL_MY = 12;
    /** Colonne 13 : 3 sigma_My. */
    public static final int COL_STDMY = 13;
    /** Colonne 14 : &lt;M_z&gt;. */
    public static final int COL_MZ = 14;
    /** Colonne 15 : 3 sigma_Mz. */
    public static final int COL_STDMZ = 15;

    /** Colonne 16 : sweeps de thermalisation reellement consommes. */
    public static final int COL_THERM = 16;
    /** Colonne 17 : t0, debut de la zone equilibree (Chodera). */
    public static final int COL_T0 = 17;
    /** Colonne 18 : tau_int de la serie de <b>thermalisation</b>, en sweeps. */
    public static final int COL_TAU = 18;
    /** Colonne 19 : g = 1 + 2 tau de la serie de thermalisation. */
    public static final int COL_G = 19;
    /** Colonne 20 : N_eff de la queue de la serie de thermalisation. */
    public static final int COL_NEFF = 20;
    /** Colonne 21 : sweeps de production totaux (bloc planifie + extensions). */
    public static final int COL_PROD = 21;
    /** Colonne 22 : 1 si l'equilibration a ete detectee, 0 si coupee par {@code maxSweeps}. */
    public static final int COL_CONVERGED = 22;
    /** Colonne 23 : secondes passees a cette temperature. */
    public static final int COL_SECONDS = 23;
    /** Colonne 24 : cadence de mesure de production, en sweeps. */
    public static final int COL_MEASURE_EVERY = 24;
    /** Colonne 25 : nombre de mesures de production accumulees. */
    public static final int COL_NMEAS = 25;
    /** Colonne 26 : sigma a la fin de la thermalisation. */
    public static final int COL_SIGMA_THERM = 26;
    /** Colonne 27 : N_eff <b>reellement obtenu</b> sur la serie de production. */
    public static final int COL_NEFF_ACHIEVED = 27;
    /** Colonne 28 : g de la serie de production (Wolff). */
    public static final int COL_G_PROD = 28;

    /** Nombre total de colonnes de {@code observablesData} en mode adaptatif. */
    public static final int N_COLUMNS = 29;

    /** Nombre d'observables de {@code TimeSeriesStats.getMeasurements(int)} : colonnes COL_E a COL_STDMZ. */
    public static final int N_OBS = 9;

    /**
     * Fenetre d'integration S(Q,w) historique (hbar/meV) : ndt = W/(dt*MR) = 100.
     * Avec --wmax X &gt; 0 : MR = round(pi/(X*dt)), ndt = round(W/Delta) ; --wmax 0 : ndt = 1.
     */
    public static final int SQW_INTEGRATION_WINDOW = 20;

    // ------------------------------------------------------------------ accesseurs du recuit

    /** Temperature initiale du recuit, exposee pour les tests. */
    public double initTemp() { return initTemp; }

    /** Temperature finale (exclue : la boucle s'arrete des que {@code T <= endTemp}). */
    public double endTemp() { return endTemp; }

    /** Ratio geometrique {@code T_{k+1} / T_k} du refroidissement historique. */
    public double dTemp() { return dTemp; }

    /** Instantanes S(Q,w) au-dela duquel {@code processSQW} ecrit son fichier ; sqwBlock en prend {@code nSQW + 1}. */
    public int nSQW() { return nSQW; }

    /**
     * Options du runner adaptatif {@link MonteCarlo#runParallel}. Defauts : refroidissement
     * geometrique, precession en fin de temperature, S(Q,w) a la derniere temperature,
     * sigma fige en production ; sur-relaxation desactivee ici (AppParallel l'active).
     */
    public static final class ParallelOptions {
        /** Parametres du controleur d'equilibration ; {@code channels} force a 2 sur une copie, l'objet fourni n'est jamais modifie. */
        public AdaptiveThermalization.Config thermalization = new AdaptiveThermalization.Config();
        /** {@code null} : refroidissement geometrique {@code T *= dTemp} (schema historique). */
        public AdaptiveCoolingSchedule cooling = null;
        /** {@code null} : pas de journal CSV. Le journal n'est pas ferme par le runner. */
        public EquilibrationLog log = null;
        /** Repertoire racine des fichiers S(Q,w) ; {@code null} = chemin code en dur de DynamicStructureFactor. */
        public String sqwOutputDir = null;
        /**
         * {@code true} (defaut) : sigma fige en production — l'adapter rendrait la chaine
         * non markovienne (proposition dependant de l'historique des taux de rejet).
         */
        public boolean freezeSigmaInProduction = true;
        /** Balayages de sur-relaxation microcanonique (Creutz 1987) apres chaque balayage Metropolis ; 0 = desactive. */
        public int overRelaxationPerSweep = 0;
        /** {@code false} : aucun instantane S(Q,w) (indispensable pour les tests). */
        public boolean computeSQW = true;
        /** {@code true} : {@code spin.Precession1D()} en fin de temperature. */
        public boolean precessionAtEnd = true;
        /** NaN (defaut) : parametres S(Q,w) historiques. Sinon w_max vise, en unites de Spin. */
        public double wMax = Double.NaN;
        /** Temperature initiale du recuit si fournie (sinon NaN = historique, 10) ; validee &gt; 0 et &gt; endTemp. */
        public double initTemp = Double.NaN;
        /** Temperature finale si fournie (sinon NaN = historique, 0.001) : borne d'arret et temperature du bloc S(Q,w). */
        public double endTemp = Double.NaN;
        /** Fenetre d'integration S(Q,w) si fournie (sinon NaN = {@link #SQW_INTEGRATION_WINDOW}) ; fixe dw = 2 pi / W. */
        public double window = Double.NaN;
        /**
         * Nombre d'instantanes du bloc S(Q,w) si fourni ({@code null} = defaut : le champ
         * historique {@code nSQW} = 40) : le bloc prend N+1 mesures spectrales moyennees et
         * DynamicStructureFactor n'ecrit son Avro qu'a la (N+1)-ieme (son compteur Nmesures
         * doit depasser la valeur passee a {@code processSQW}). Option CLI --nSQW ; valide
         * &gt;= 2 (et &lt;= 1 000 000, garde-fou d'overflow) cote parseur AppParallel.
         */
        public Integer nSQW = null;
        /**
         * Composante z du vecteur de diffusion, en r.l.u. (option --qz d'AppParallel ; defaut 0
         * = honeycomb 2D). Enregistree comme propriete "Qzz" dans tout Avro de structure ecrit :
         * structure_dipolaire toujours ; structure_couleur/dimere seulement si QUARTIC actif.
         */
        public double qz = 0.0;
        /** Non null : appele en fin de production de chaque temperature (config equilibree, avant la precession) — instantane des spins. */
        public java.util.function.BiConsumer<Double, Lattice> spinSnapshotListener = null;
        /** {@code true} : barre de progression et ligne de resume par temperature sur stdout. */
        public boolean verbose = true;
        /**
         * Graine de la configuration initiale. {@code null} : RNG de l'instance (graine
         * aleatoire). Fixee, avec celle du sweeper : observables identiques bit a bit.
         */
        public Long initialConfigSeed = null;
        /**
         * Crochet de test : appele A LA PLACE de {@code solveODE1D} + {@code processSQW} pour
         * chaque instantane S(Q,w), sans entree-sortie (processSQW ecrit dans un repertoire
         * code en dur). {@code null} en production.
         */
        java.util.function.IntConsumer sqwSnapshotListener = null;

        // Diagnostic de stagnation : champs additifs, defauts = comportement historique bit a bit.

        /**
         * {@code null} (defaut) : aucun diagnostic. Non null : detecteur actif par temperature,
         * cadence calquee sur le controleur ; q(t) est tenu localement et jamais ajoute aux
         * canaux du controleur (decroissance monotone, pas un plateau stationnaire).
         */
        public StagnationDetector.Config stagnation = null;
        /** {@code true} : arrete la thermalisation sur verdict STAGNANT_* du detecteur. Sans effet si {@code stagnation == null}. */
        public boolean stopThermalizationOnStagnation = false;
        /**
         * {@code true} : arrete les extensions de production quand le gain de N_eff sature
         * ({@code neffGainFraction}) ; la premiere extension n'est jamais coupee. Defaut de
         * StagnationDetector.Config (0.15) si {@code stagnation == null}.
         */
        public boolean stopProductionOnSaturation = false;
        /** Journal de decorrelation (une ligne par temperature) ; aucune ecriture sans detecteur actif. Non ferme par le runner. */
        public StagnationLog stagnationLog = null;
        /**
         * K_q du critere {@code thermSweeps >= K_q * tauOverlapSweeps} ; defaut 3 ~ memoire
         * residuelle exp(-3) ~ 5 %. Voir Ogielski, Phys. Rev. B 32, 7384 (1985).
         */
        public double memoryMultiplier = 3.0;
    }

    /**
     * Runner adaptatif et parallele : balayages delegues a {@link CheckerboardMetropolis},
     * thermalisation et production dimensionnees en ligne ({@link AdaptiveThermalization}),
     * production prolongee jusqu'a la cible de N_eff — mesuree sur l'energie seule (la
     * direction de |M| diffuse librement, mode de Goldstone : exiger le meme N_eff ferait
     * exploser le cout). Par temperature : thermalisation, production (+ extensions),
     * temperature suivante, S(Q,w) si {@code nextT <= tEnd} (decide par nextT : une fois
     * et une seule, a la derniere temperature), precession, stockage.
     *
     * @param sweeper construit sur le <b>meme</b> objet {@code lattice} ; non ferme par cette methode.
     * @param opt     options ; {@code null} = {@code new ParallelOptions()}.
     * @throws IllegalArgumentException parametre nul ou incoherent, sweeper sur un autre
     *         reseau, ou computeSQW actif sans {@code uc}/{@code output}.
     */
    public void runParallel(Lattice lattice, Spin spin, UnitCell uc, String output,
                            double[] H_field, double B,
                            CheckerboardMetropolis sweeper, ParallelOptions opt) {
        final ParallelOptions options = (opt == null) ? new ParallelOptions() : opt;

        // Bornes du recuit : options si fournies, sinon valeurs historiques du constructeur.
        final double tInit = Double.isNaN(options.initTemp) ? this.initTemp : options.initTemp;
        final double tEnd = Double.isNaN(options.endTemp) ? this.endTemp : options.endTemp;
        if (tInit <= tEnd) {
            throw new IllegalArgumentException(
                    "initTemp (" + tInit + ") doit etre > endTemp (" + tEnd + ")");
        }

        // Validation immediate : ces erreurs sinon plantent apres des heures de calcul.
        if (lattice == null || spin == null || sweeper == null) {
            throw new IllegalArgumentException("lattice, spin et sweeper sont obligatoires");
        }
        if (H_field == null || H_field.length < 3) {
            // Les colonnes 4..6 lisent H_field[0..2] a chaque temperature.
            throw new IllegalArgumentException("H_field doit avoir au moins 3 composantes"
                    + (H_field == null ? " (null)" : " (longueur " + H_field.length + ")"));
        }
        if (sweeper.lattice() != lattice) {
            // Un sweeper sur un autre Lattice thermaliserait une copie des spins.
            throw new IllegalArgumentException(
                    "le sweeper ne porte pas sur ce reseau (CheckerboardMetropolis.lattice() != lattice)");
        }
        if (options.computeSQW) {
            // uc/output indispensables des le premier instantane, a la derniere temperature.
            if (uc == null) {
                throw new IllegalArgumentException(
                        "computeSQW = true exige une UnitCell (DynamicStructureFactor.processSQW)");
            }
            if (output == null) {
                throw new IllegalArgumentException(
                        "computeSQW = true exige un nom de sortie (DynamicStructureFactor.processSQW)");
            }
        }

        final int latticeSize = lattice.size;

        // Initialisation du reseau (au generateur pres, cf. initialConfigSeed).
        if (!isInitialized) {
            final RandomGenerator init = (options.initialConfigSeed == null)
                    ? this.random
                    : RandomGeneratorFactory.of("Xoshiro256PlusPlus").create(options.initialConfigSeed);
            for (int i = 0; i < latticeSize; i++) {
                lattice.setSpin1D(i, spin.SRG7(init));
            }
            isInitialized = true;
        }
        // DynamicStructureFactor alloue seulement dans sqwBlock (~280 Mo par replique en 30x30).

        // Controleur d'equilibration : cree une fois, remis a zero a chaque temperature.
        final AdaptiveThermalization.Config cfg =
                (options.thermalization == null ? new AdaptiveThermalization.Config()
                                                : options.thermalization).copy();
        cfg.channels = 2; // energie et |M|
        final AdaptiveThermalization ctrl = new AdaptiveThermalization(cfg);

        // Stagnation : cree une fois ; null tant qu'aucun diagnostic n'est demande.
        final boolean stagnationActive = options.stagnation != null;
        final StagnationDetector stag = stagnationActive
                ? new StagnationDetector(options.stagnation) : null;
        final SpinOverlap overlapRef = stagnationActive ? new SpinOverlap(latticeSize) : null;
        // Le journal de decorrelation exige un detecteur actif.
        final boolean stagnationLogActive = stagnationActive && options.stagnationLog != null;
        // Defaut de StagnationDetector.Config si stagnation == null.
        final double neffGainFraction = stagnationActive
                ? options.stagnation.neffGainFraction
                : new StagnationDetector.Config().neffGainFraction;
        // Serie q(t) : tableau primitif a croissance geometrique, pas de ArrayList<Double> dans la boucle chaude.
        double[] qSeries = stagnationActive ? new double[1024] : null;
        int qN = 0;

        final double logInit = Math.log(tInit);
        final double logSpan = logInit - Math.log(tEnd);

        // Meme regle de sortie que « for (Temp = initTemp; Temp > endTemp; Temp *= dTemp) ».
        double Temp = tInit;
        while (Temp > tEnd) {

            ctrl.reset();
            double sig = 60.d;
            long thermSweeps = 0L;
            final long tStartNanos = System.nanoTime();

            // Reference q(t) : configuration heritee de la temperature precedente (au premier
            // T : configuration initiale). Remise a zero ici, pas a la construction.
            if (stagnationActive) {
                overlapRef.capture(lattice);
                qN = 0;
                stag.reset();
            }
            int lastStagnationDetections = 0;
            boolean stagnationStop = false;

            // ------------------------------------------------------------- thermalisation
            // Garde-fou : au-dela de maxSweeps + sampleEvery, le controleur a forcement statue.
            final long thermCap = (long) cfg.maxSweeps + cfg.sampleEvery;
            // stagnationStop ne devient vrai qu'apres au moins un sweep : la garde reste valable.
            while (!ctrl.isEquilibrated() && !stagnationStop) {
                if (thermSweeps >= thermCap) {
                    throw new IllegalStateException("thermalisation sans decision apres "
                            + thermSweeps + " sweeps (maxSweeps=" + cfg.maxSweeps + ")");
                }
                CheckerboardMetropolis.SweepStats st = sweeper.sweep(Temp, sig);
                for (int k = 0; k < options.overRelaxationPerSweep; k++) {
                    sweeper.overRelaxationSweep();
                }
                thermSweeps++;
                sig = CheckerboardMetropolis.nextSigma(sig, st.rejectionRate());
                if (ctrl.wantsSample(thermSweeps)) {
                    ctrl.observe(thermSweeps, sweeper.energy(), sweeper.magnetizationNorm());
                    if (stagnationActive) {
                        // q(t) echantillonne au meme instant que ctrl.observe() : indices alignes.
                        qSeries = growSeries(qSeries, qN + 1);
                        qSeries[qN++] = overlapRef.overlap(lattice);

                        // Verification cadencee sur detectionsRun() : reutilise la cadence du controleur.
                        final int detections = ctrl.detectionsRun();
                        if (detections != lastStagnationDetections) {
                            lastStagnationDetections = detections;
                            final double[] energySeries = ctrl.channelSeries(0);
                            // Budget en ECHANTILLONS (pas en sweeps) restant avant maxSweeps.
                            final long budgetSamples = Math.max(0L,
                                    ((long) cfg.maxSweeps - thermSweeps) / cfg.sampleEvery);
                            stag.update(energySeries, energySeries.length, qSeries, qN,
                                        (double) budgetSamples);
                            stagnationStop = options.stopThermalizationOnStagnation
                                    && stag.stagnant();
                        }
                    }
                }
            }
            final double sigmaEndOfThermalization = sig;

            // Raison d'arret : pour le seul journal de decorrelation.
            final String stopReason;
            if (stagnationStop) {
                stopReason = stag.verdict().name(); // "STAGNANT_PLATEAU" ou "STAGNANT_FROZEN"
            } else if (ctrl.converged()) {
                stopReason = "EQUILIBRATED";
            } else {
                stopReason = "MAX_SWEEPS";
            }

            // ------------------------------------------------------------- production
            // Cadence de mesure adaptative (~2 tau) : productionSweeps() en tient compte.
            final int measureEvery = Math.max(1, ctrl.productionMeasureEvery());
            final int plannedProd = ctrl.productionSweeps();

            final TimeSeriesStats stats = new TimeSeriesStats();
            final ProductionAccumulator acc = new ProductionAccumulator();
            sig = runProductionBlock(sweeper, options, Temp, sig, 0L, plannedProd, measureEvery,
                                     stats, acc);
            long totalProd = plannedProd;

            // Filet : au moins une mesure, sinon getMeasurements() renverrait des NaN.
            if (stats.size() == 0) {
                stats.add(sweeper.energy(), sweeper.magnetization());
            }

            // Extension : livrer le N_eff promis, pas seulement l'avoir prevu.
            double gProduction = productionG(stats);
            double neffAchieved = stats.size() / gProduction;
            // firstExtension : la premiere extension n'a pas de precedente a laquelle se comparer.
            boolean neffSaturated = false;
            boolean firstExtension = true;
            while (neffAchieved < cfg.targetNeffProduction && totalProd < cfg.maxProductionSweeps) {
                // Manque converti en sweeps ; plancher a 25 % du bloc planifie contre les
                // extensions minuscules (g bruitee, recalcul O(n log n)).
                final double needed = (cfg.targetNeffProduction - neffAchieved)
                        * gProduction * measureEvery;
                long extra = (long) Math.max(Math.ceil(needed), Math.ceil(0.25d * plannedProd));
                extra = Math.min(extra, cfg.maxProductionSweeps - totalProd);
                if (extra <= 0L) break;
                final double neffBeforeExtension = neffAchieved;
                sig = runProductionBlock(sweeper, options, Temp, sig, totalProd, (int) extra,
                                         measureEvery, stats, acc);
                totalProd += extra;
                gProduction = productionG(stats);
                neffAchieved = stats.size() / gProduction;
                if (options.stopProductionOnSaturation && !firstExtension
                        && StagnationDetector.neffSaturated(neffBeforeExtension, neffAchieved,
                                                            neffGainFraction)) {
                    neffSaturated = true;
                    break;
                }
                firstExtension = false;
            }

            // ------------------------------------------------------------- temperature suivante
            final double nextT = (options.cooling != null)
                    ? options.cooling.next(Temp, stats.varianceEnergyTotal(), latticeSize)
                    : Temp * this.dTemp;

            // ------------------------------------------------------------- S(Q,w)
            // « Derniere temperature » n'est plus une estimation : nextT est connu.
            if (options.computeSQW && nextT <= tEnd) {
                sqwBlock(lattice, spin, uc, output, B, sweeper, options, Temp, sig);
            }

            // Instantane des spins (--spins-per-T) : config equilibree, avant la precession.
            if (options.spinSnapshotListener != null) {
                options.spinSnapshotListener.accept(Temp, lattice);
            }

            if (options.precessionAtEnd) {
                spin.Precession1D();
            }

            // ------------------------------------------------------------- stockage
            MCIdx++;
            final double[] measurements = stats.getMeasurements(latticeSize);
            if (measurements.length != N_OBS) {
                // Sinon l'arraycopy ecraserait les colonnes de diagnostics.
                throw new IllegalStateException("TimeSeriesStats.getMeasurements() a renvoye "
                        + measurements.length + " valeurs, " + N_OBS + " attendues");
            }
            final double acceptance = acc.meanAcceptance();
            final double seconds = (System.nanoTime() - tStartNanos) * 1e-9;
            final double ePerSite = measurements[0];
            final double eErrPerSite = stats.energyStandardError() / latticeSize;
            final double mNorm = Math.sqrt(measurements[3] * measurements[3]
                    + measurements[5] * measurements[5]
                    + measurements[7] * measurements[7]);

            final TemperaturePoint point = new TemperaturePoint(
                    MCIdx, Temp, thermSweeps, ctrl.t0Sweep(), ctrl.tauIntSweeps(),
                    ctrl.statisticalInefficiency(), ctrl.neff(), ctrl.gewekeZ(0), ctrl.gewekeZ(1),
                    (int) Math.min(Integer.MAX_VALUE, totalProd), ctrl.converged(), acceptance,
                    sig, sigmaEndOfThermalization, ePerSite, eErrPerSite, mNorm, seconds,
                    measureEvery, stats.size(), neffAchieved, gProduction,
                    // Observables pour le CSV d'equilibration (fusion) : champ + serie de production.
                    H_field[0], H_field[1], H_field[2],
                    measurements[1], measurements[2],
                    measurements[3], measurements[4], measurements[5],
                    measurements[6], measurements[7], measurements[8]);

            final double[] combinedData = new double[N_COLUMNS];
            combinedData[COL_MCIDX] = MCIdx;
            combinedData[COL_ACC] = acceptance;
            combinedData[COL_SIGMA] = sig;
            combinedData[COL_T] = Temp;
            combinedData[COL_HX] = H_field[0];
            combinedData[COL_HY] = H_field[1];
            combinedData[COL_HZ] = H_field[2];
            System.arraycopy(measurements, 0, combinedData, COL_E, N_OBS);
            final double[] extras = point.extras();
            System.arraycopy(extras, 0, combinedData, COL_THERM, extras.length);
            observablesData.add(combinedData);

            if (options.log != null) {
                options.log.write(point);
            }

            // ------------------------------------------------------------- diagnostic de stagnation
            // Journal separe, ecrit seulement si detecteur ET journal actifs.
            if (stagnationLogActive) {
                final StagnationDetector.Drift lastDrift = stag.lastDrift();
                final SpinOverlap.Decay lastFullDecay = stag.lastFullDecay();

                // tauOverlapSweeps = MULTIPLICATION simple (tau_q mesure en echantillons sur q(t))
                // — NE PAS utiliser l'inversion exponentielle de tauSweepsFromSamples (facteur
                // logarithmique d'erreur).
                final double tauOverlapSweeps = (lastFullDecay != null)
                        ? lastFullDecay.tauSamples() * cfg.sampleEvery
                        : Double.NaN;
                // tau_q = 0 => +Inf : memoire nulle, infiniment mieux oubliee que le seuil exige.
                final double memoryRatio = thermSweeps / tauOverlapSweeps;
                final boolean temperaturesDecorrelated = lastFullDecay != null
                        && lastFullDecay.plateaued()
                        && memoryRatio >= options.memoryMultiplier;

                // Conversion g -> tau de la serie de production : tau_sweeps = (g-1)/2 * Delta.
                final double tauProductionSweeps = (gProduction - 1.0d) / 2.0d * measureEvery;
                // g_Delta <= 2 <=> measureEvery >= 2 tau_prod <=> measurementsPerTau >= 1.
                final double measurementsPerTau = measureEvery / (2.0d * tauProductionSweeps);
                final boolean measurementsDecorrelated = measurementsPerTau >= 1.0d;

                final StagnationPoint stagPoint = new StagnationPoint(
                        MCIdx, Temp, thermSweeps, totalProd,
                        lastFullDecay != null ? lastFullDecay.qInf() : Double.NaN,
                        tauOverlapSweeps, memoryRatio, temperaturesDecorrelated,
                        lastDrift.z(), lastDrift.slopePerSample() / cfg.sampleEvery,
                        stopReason, ctrl.converged(),
                        measureEvery, tauProductionSweeps, measurementsPerTau,
                        measurementsDecorrelated, neffAchieved, neffSaturated, seconds,
                        H_field[0], H_field[1], H_field[2]);
                options.stagnationLog.write(stagPoint);
            }

            if (options.verbose) {
                printTemperatureSummary(point, logInit, logSpan, tEnd);
            }

            Temp = nextT;
        }
    }

    /** Taux d'acceptation moyen sur les sweeps reellement executes (blocs d'extension compris). */
    private static final class ProductionAccumulator {
        private double sum;
        private long sweeps;

        void add(double acceptanceRate) { sum += acceptanceRate; sweeps++; }

        double meanAcceptance() { return sweeps == 0L ? 0.d : sum / sweeps; }
    }

    /**
     * Un bloc de production : {@code sweeps} balayages, une mesure tous les {@code measureEvery},
     * acceptation cumulee. Grille comptee depuis le debut de la production : g (Wolff) garde un sens.
     *
     * @return sigma en fin de bloc (inchange si {@code freezeSigmaInProduction}).
     */
    private double runProductionBlock(CheckerboardMetropolis sweeper, ParallelOptions options,
                                      double Temp, double sigma, long sweepOffset, int sweeps,
                                      int measureEvery, TimeSeriesStats stats,
                                      ProductionAccumulator acc) {
        double sig = sigma;
        for (int p = 1; p <= sweeps; p++) {
            final long sweepIndex = sweepOffset + p;
            CheckerboardMetropolis.SweepStats st = sweeper.sweep(Temp, sig);
            for (int k = 0; k < options.overRelaxationPerSweep; k++) {
                sweeper.overRelaxationSweep();
            }
            acc.add(st.acceptanceRate());
            if (!options.freezeSigmaInProduction) {
                sig = CheckerboardMetropolis.nextSigma(sig, st.rejectionRate());
            }
            if (sweepIndex % measureEvery == 0) {
                stats.add(sweeper.energy(), sweeper.magnetization());
            }
        }
        return sig;
    }

    /**
     * g de la serie de production (Gamma de Wolff), ramene a 1 si non exploitable
     * (n &lt;= 1, serie gelee, estimation negative) : {@code N_eff = n}.
     */
    private static double productionG(TimeSeriesStats stats) {
        final double g = stats.energyG();
        return (Double.isFinite(g) && g >= 1.d) ? g : 1.d;
    }

    /** Double la capacite si {@code need} depasse la longueur (reserve a qSeries : pas d'allocation par sweep). */
    private static double[] growSeries(double[] arr, int need) {
        if (need <= arr.length) return arr;
        int cap = arr.length;
        while (cap < need) cap <<= 1;
        return Arrays.copyOf(arr, cap);
    }

    /**
     * Bloc S(Q,w) de la derniere temperature : exactement {@code (N + 1) * measurementRateSQW}
     * balayages, un instantane tous les {@code measurementRateSQW} (contrainte de processSQW,
     * qui n'ecrit qu'au-dela de N), ou N est le nombre d'instantanes : {@code opt.nSQW}
     * (CLI --nSQW) s'il est fourni, sinon la valeur historique du champ {@code nSQW} (40).
     * DynamicStructureFactor est reconstruit A CHAQUE appel :
     * son compteur Nmesures ne se remet pas a zero, une instance reutilisee reecrirait le
     * fichier a chaque instantane. Aucune mesure thermodynamique ici.
     */
    private void sqwBlock(Lattice lattice, Spin spin, UnitCell uc, String output, double B,
                          CheckerboardMetropolis sweeper, ParallelOptions opt,
                          double Temp, double sigma) {

        // Sampling S(Q,w) : wMax > 0 => Delta = pi/wMax (Nyquist) ; wMax = 0 => statique, ndt = 1.
        final double dt = spin.getDt();
        // Fenetre : --window si fournie, sinon 20 ; fixe dw = 2 pi / fenetre.
        final double fenetre = Double.isNaN(opt.window) ? SQW_INTEGRATION_WINDOW : opt.window;
        final int measurementRate;
        final int ndt;
        final double delta;
        if (!Double.isNaN(opt.wMax) && opt.wMax > 0) {
            final double deltaCible = Math.PI / opt.wMax;   // convention standard : Nyquist = pi/Delta
            measurementRate = Math.max(1, (int) Math.round(deltaCible / dt));
            delta = dt * measurementRate;
            ndt   = (int) Math.round(fenetre / delta);
        } else if (!Double.isNaN(opt.wMax)) {
            measurementRate = Math.max(1, (int) (fenetre / dt));
            delta = fenetre;
            ndt   = 1;
        } else {
            // Historique : les 100 instantanes sont tous FFTes (resolution x2 vs l'ancien code
            // => spectres non comparables aux fichiers d'avant le 12/09/2026).
            measurementRate = 20;
            delta = dt * measurementRate;
            ndt   = (int) (fenetre / dt) / measurementRate;  // = 100
        }
        // Garde-fou : une fenetre derivee trop courte (cas atteignable : --domega enorme sans
        // --wmax) ne doit jamais construire un DynamicStructureFactor avec time < 2. La
        // branche statique (wMax = 0 => ndt = 1, coupe unique a w = 0) est exempte : son
        // ndt = 1 est la convention, pas un defaut de reglage.
        if (ndt < 2 && (Double.isNaN(opt.wMax) || opt.wMax > 0)) {
            throw new IllegalStateException("S(Q,w) : fenetre " + fenetre + " / delta " + delta
                    + " => ndt = " + ndt + " (< 2) : au moins 2 points spectraux sont requis"
                    + " (augmenter la fenetre W, ou baisser --domega)");
        }
        // Nombre d'instantanes du bloc : override ParallelOptions (CLI --nSQW), sinon la valeur
        // historique du champ. Resolu UNE fois ici : les QUATRE sites ci-dessous (budget de
        // sweeps, les deux processSQW, la barre de progression) doivent voir la MEME valeur,
        // sinon DynamicStructureFactor n'ecrit jamais son Avro (il n'ecrit qu'a la (N+1)-ieme
        // mesure, voir son compteur Nmesures).
        final int nSQWEffectif = (opt.nSQW != null) ? opt.nSQW : nSQW;
        this.sqw = new DynamicStructureFactor(H, K, ndt, delta, lattice);
        if (opt.sqwOutputDir != null) this.sqw.setOutputDirPath(opt.sqwOutputDir);
        // Canaux couleur/dimere : seulement si le terme quartique est actif (ligne QUARTIC de
        // l'input, b != 0). Le canal dimere mesure tau_b = (S_i^g S_j^g)^2, la variable de
        // liaison du terme quartique lui-meme ; sans QUARTIC, ni calcules ni ecrits.
        final boolean quartic = lattice.hasQuartic();
        // Canal couleur : n_i^a = (S_i^a)^2, statique, ecriture Avro apres la boucle.
        final ColorStructureFactor colorSF = quartic ? new ColorStructureFactor(H, K, lattice) : null;
        if (quartic && opt.sqwOutputDir != null) colorSF.setOutputDirPath(opt.sqwOutputDir);
        // Canal dimere : tau_b = (S_i^g S_j^g)^2 aux milieux de liaisons (pinch point).
        final DimerStructureFactor dimerSF = quartic ? new DimerStructureFactor(H, K, lattice) : null;
        if (quartic && opt.sqwOutputDir != null) dimerSF.setOutputDirPath(opt.sqwOutputDir);
        // Traceabilite --qz : tout Avro de structure ecrit porte la propriete "Qzz" (dipolaire
        // toujours ; couleur/dimere quand ils existent, i.e. QUARTIC actif).
        sqw.setQzz(opt.qz);
        if (quartic) {
            colorSF.setQzz(opt.qz);
            dimerSF.setQzz(opt.qz);
        }
        final int sweeps = (nSQWEffectif + 1) * measurementRateSQW;
        double sig = sigma;
        measurementIdxSQW = 0;
        for (int p = 1; p <= sweeps; p++) {
            CheckerboardMetropolis.SweepStats st = sweeper.sweep(Temp, sig);
            for (int k = 0; k < opt.overRelaxationPerSweep; k++) {
                sweeper.overRelaxationSweep();
            }
            if (!opt.freezeSigmaInProduction) {
                sig = CheckerboardMetropolis.nextSigma(sig, st.rejectionRate());
            }
            if (p % measurementRateSQW == 0) {
                measurementIdxSQW += 1;
                if (opt.sqwSnapshotListener != null) {
                    // Mode test : ni RK4 ni ecriture de fichier.
                    opt.sqwSnapshotListener.accept(measurementIdxSQW);
                } else if (ndt == 1) {
                    // Statique : snapshot normalise direct, sans RK4.
                    final double[] src = lattice.getSpin1Dlattice();
                    final double[] S_t = new double[src.length];
                    for (int i = 0; i < lattice.size; i++) {
                        final double n = MathOps.norm(src, i * 3);
                        S_t[3 * i] = src[3 * i] / n;
                        S_t[3 * i + 1] = src[3 * i + 1] / n;
                        S_t[3 * i + 2] = src[3 * i + 2] / n;
                    }
                    sqw.processSQW(uc, S_t, output, nSQWEffectif, B);
                    if (quartic) {
                        colorSF.processConfig(uc, S_t);
                        dimerSF.processConfig(uc, S_t);
                    }
                } else {
                    final double[] S_t = spin.solveODE1D(0, fenetre,
                            measurementRate, ndt);
                    sqw.processSQW(uc, S_t, output, nSQWEffectif, B);
                    if (quartic) {
                        colorSF.processConfig(uc, S_t);
                        dimerSF.processConfig(uc, S_t);
                    }
                }
                if (opt.verbose) {
                    System.out.print("\r" + measurementIdxSQW + "/" + nSQWEffectif);
                }
            }
        }
        // Ecriture Avro couleur/dimere apres les nSQW+1 mesures ; pas en mode test.
        // Noms deterministes : le repertoire (opt.sqwOutputDir) porte l'identite de la replique.
        // Sans terme quartique, ces canaux n'ont pas ete alloues : rien a ecrire.
        if (quartic && opt.sqwSnapshotListener == null) {
            colorSF.writeAvro("structure_couleur_h" + String.format(java.util.Locale.US, "%.6f", B) + ".avro");
            dimerSF.writeAvro("structure_dimere_h" + String.format(java.util.Locale.US, "%.6f", B) + ".avro");
        }
    }

    /** Ligne de resume et barre de progression d'une temperature (mode verbeux). */
    private void printTemperatureSummary(TemperaturePoint p, double logInit, double logSpan,
                                         double tEnd) {
        System.out.print("\r");
        System.out.println(String.format(java.util.Locale.US,
                "%s T=%9.5f therm=%7d t0=%7d tau=%8.2f g=%8.2f neff=%8.1f %s "
                + "prod=%7d dt=%4d nmes=%5d neffObt=%8.1f gProd=%7.2f acc=%.4f sigma=%9.5f "
                + "E/N=%12.6f +- %.6f |M|=%10.4f %7.2fs",
                Thread.currentThread().getName(), p.T(), p.thermSweeps(), p.t0Sweep(),
                p.tauSweeps(), p.g(), p.neff(), p.converged() ? "conv" : "CUT ", p.prodSweeps(),
                p.measureEvery(), p.nMeasurements(), p.neffAchieved(), p.gProduction(),
                p.acceptance(), p.sigma(), p.ePerSite(), p.eErr(), p.mNorm(), p.seconds()));

        // Progression en log(T) : le nombre de temperatures n'est pas connu a l'avance.
        final int barLength1 = 50;
        final double frac = logSpan > 0.d ? (logInit - Math.log(p.T())) / logSpan : 1.d;
        final int progress1 = (int) (Math.max(0.d, Math.min(1.d, frac)) * barLength1);
        System.out.print("\r" + Thread.currentThread().getName() + ": [");
        for (int i = 0; i < barLength1; i++) {
            if (i < progress1) {
                System.out.print("=");
            } else if (i == progress1) {
                System.out.print(">");
            } else {
                System.out.print(" ");
            }
        }
        System.out.print(String.format(java.util.Locale.US, "] T = %.5f / %.5f ",
                p.T(), tEnd));
    }

}
