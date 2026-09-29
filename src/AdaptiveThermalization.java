import java.util.Arrays;

/**
 * Controleur d'equilibration <i>en ligne</i> pour une temperature donnee d'un recuit
 * Monte Carlo : il remplace la longueur de thermalisation fixe
 * ({@code thermalizationSweeps = 40 000}) par un critere automatique fonde sur la
 * litterature (voir {@link EquilibrationDetector} pour les references completes).
 *
 * <p>Principe : la simulation echantillonne une ou plusieurs observables (typiquement
 * l'energie et |M|) tous les {@code sampleEvery} sweeps. A une cadence <b>geometrique</b>
 * (seuils {@code checkEvery}, puis {@code ceil(n . checkGrowth)} avec checkGrowth = 1.2 par
 * defaut ; voir {@link Config#checkGrowth}), et une fois {@code minSweeps} <i>ecoules</i>
 * atteints, le controleur lance <b>pour chaque canal</b> :</p>
 * <ol>
 *   <li>la detection d'equilibration de Chodera (2016) : t0 maximisant
 *       N_eff = (n - t0)/g, g estime par l'estimateur rapide de pymbar puis raffine par la
 *       methode Gamma de Wolff (2004) ;</li>
 *   <li>trois criteres d'acceptation :
 *     <ul>
 *       <li><b>longueur</b> : (n - t0) . sampleEvery &gt;= K . tau_sweeps, avec
 *           tau_sweeps = {@link #tauSweepsFromSamples} (conversion exponentielle, et
 *           <i>non</i> une simple multiplication par sampleEvery) et
 *           K = {@code tauMultiplier} (defaut 20) ; on exige donc une queue longue d'au
 *           moins ~20 temps d'auto-correlation ;</li>
 *       <li><b>stationnarite</b> : |z| &lt;= {@code gewekeZmax} (Geweke 1992) evalue sur la
 *           <b>seconde moitie</b> de la queue, [t0 + f (n - t0), n) avec
 *           f = {@code gewekeSkipFraction} = 0.5 par defaut (voir
 *           {@link Config#gewekeSkipFraction} : sur [t0, n) le test serait structurellement
 *           incompatible avec le t0 de Chodera). En dessous de
 *           {@code minGewekeSamples} points cette fenetre est trop courte pour que z soit
 *           interpretable : le critere est declare <i>non evaluable</i> et le canal ne passe
 *           pas (voir {@link Config#minGewekeSamples}) ;</li>
 *       <li><b>statistique</b> : N_eff &gt;= {@code minNeff}.</li>
 *     </ul>
 *   </li>
 * </ol>
 * <p>L'equilibration est declaree lorsque <b>tous</b> les canaux passent les trois criteres.
 * Lorsque {@code maxSweeps} est atteint, {@link #isEquilibrated()} passe a true dans tous les
 * cas ; une derniere detection est alors lancee si la plus recente est perimee, et
 * {@link #converged()} porte <b>son</b> verdict. Cette detection de derniere minute peut tres
 * bien reussir &mdash; la cadence geometrique espace les verifications et rien ne garantit que
 * l'une d'elles tombe sur le dernier echantillon &mdash; auquel cas la thermalisation est bel et
 * bien convergee et ne doit pas etre signalee comme coupee. Si elle echoue (ou si elle avait
 * deja echoue sur ce meme nombre d'echantillons), {@code converged()} vaut false : l'appelant
 * decide (poursuivre, marquer le point de temperature comme suspect, etc.). Les estimations de
 * t0 / tau les plus recentes restent disponibles.</p>
 *
 * <p>Interet physique : pres de Tc le temps de correlation diverge (tau ~ xi^z, z ~ 2 pour
 * une dynamique locale) et une thermalisation fixe est trop courte ; loin de Tc elle est
 * inutilement longue. Ce controleur concentre automatiquement le temps de calcul la ou il
 * est necessaire.</p>
 *
 * <p><b>Canaux constants.</b> Un canal numeriquement constant (systeme gele a tres basse
 * temperature, ou observable identiquement nulle par symetrie) ne porte aucune information
 * temporelle : ni tau, ni le z de Geweke, ni N_eff n'y sont definis. Un tel canal
 * <b>passe</b> donc les trois criteres (avec t0 = 0, tau = 0, g = 1, z = 0, N_eff = n)
 * plutot que de bloquer indefiniment la thermalisation jusqu'a maxSweeps. Attention : un
 * systeme peut etre "gele mais pas equilibre" au sens vitreux ; aucun critere sur une seule
 * trajectoire ne peut le detecter (une thermalisation de longueur fixe ne le detecte pas
 * davantage). Pour ce cas de figure, utiliser plusieurs repliques partant de configurations
 * differentes et {@link EquilibrationDetector#gelmanRubin(double[][])}.</p>
 *
 * <p><b>Observables non finies.</b> Un NaN ou un Inf n'est <i>pas</i> un canal constant :
 * c'est un bug de la simulation (energie divergente, sigma qui explose, division par T = 0).
 * {@link #observe} leve donc immediatement une {@link IllegalStateException} plutot que de
 * laisser la serie passer silencieusement les criteres. Voir aussi la section
 * "Valeurs non finies" de {@link EquilibrationDetector}.</p>
 *
 * <p><b>Origine des sweeps.</b> Les indices passes a {@link #observe} n'ont besoin d'etre
 * que <i>croissants</i> : le controleur les rapporte a une origine
 * {@code sweepOrigin = premierSweep - sampleEvery} capturee au premier {@code observe()}
 * suivant {@link #reset()}. Les seuils {@code minSweeps} / {@code maxSweeps} et
 * {@link #sweepsObserved()} portent sur ce temps <i>ecoule</i>, de sorte qu'un appelant qui
 * compte a partir de 1 avec des mesures a 10, 20, ... obtient exactement
 * {@code ecoule == sweep}, et qu'un appelant dont l'horloge demarre a 1 000 000 obtient
 * exactement le meme comportement.</p>
 *
 * <p>Cette classe ne depend d'aucune classe de reseau : elle ne manipule que des series
 * temporelles (tableaux primitifs a croissance geometrique).</p>
 */
public final class AdaptiveThermalization {

    /**
     * Parametres du controleur. POJO mutable avec setters fluides ; valeurs par defaut
     * raisonnables pour un Metropolis local sur quelques milliers de sites.
     */
    public static final class Config {
        /** Nombre minimal de sweeps avant toute tentative de detection. */
        public int minSweeps = 1_000;
        /** Plafond dur de sweeps de thermalisation (securite). */
        public int maxSweeps = 200_000;
        /** Nombre de sweeps entre deux echantillons de la serie temporelle. */
        public int sampleEvery = 10;
        /** Nombre d'echantillons avant la premiere execution de la detection. */
        public int checkEvery = 50;
        /**
         * Croissance geometrique de la cadence de verification.
         *
         * <p>Les seuils de declenchement sont {@code checkEvery}, puis
         * {@code ceil(precedent * checkGrowth)}, strictement croissants. Une cadence
         * <i>fixe</i> (une verification tous les {@code checkEvery} echantillons) fait
         * n/checkEvery verifications, soit 400 pour n = 20 000 ; comme chaque verification
         * coute O(maxCandidates . n . sqrt(tau)) et que tau peut croitre comme n sur une
         * serie qui derive, la somme est dominee par les verifications tardives et atteint
         * ~3.10^10 operations (~28 s mesurees). Avec une croissance g le nombre de
         * verifications tombe a log(n/checkEvery)/log(g), soit ~33 pour g = 1.2 et
         * n = 20 000.</p>
         *
         * <p>Prix a payer : le depassement de latence de detection vaut au plus un facteur g
         * (20 % pour g = 1.2) sur la longueur de thermalisation, puisque l'equilibration est
         * constatee au plus tard au seuil suivant. C'est negligeable devant le cout d'une
         * cadence fixe en O(n^1.5) par verification.</p>
         *
         * <p>Valeur &lt;= 1 : on retombe sur la cadence fixe historique
         * ({@code n % checkEvery == 0}).</p>
         */
        public double checkGrowth = 1.2;
        /** K : la queue doit couvrir au moins K temps d'auto-correlation. */
        public double tauMultiplier = 20.0;
        /** Seuil |z| du test de Geweke. */
        public double gewekeZmax = 2.0;
        /**
         * Fraction initiale de la queue [t0, n) exclue du test de Geweke.
         *
         * <p>Le t0 de Chodera <b>maximise N_eff</b> : il accepte donc deliberement un biais
         * residuel de l'ordre d'une erreur standard en t0 (c'est tout l'interet du critere,
         * qui echange un peu de biais contre beaucoup d'echantillons). Appliquer Geweke a
         * partir de t0 exactement met la fenetre A pile sur ce biais accepte ; comme le biais
         * est fixe par t0 alors que l'erreur standard decroit en 1/sqrt(n), |z| <b>croit</b>
         * avec n et finit par rejeter toute troncature finie d'un transitoire exponentiel.
         * Les deux criteres sont alors structurellement incompatibles.</p>
         *
         * <p>Le residu decroit sur un temps de relaxation ; des lors que la queue satisfait le
         * critere de longueur (L &gt;= K tau), sa <b>seconde moitie</b> en est donc exempte.
         * On evalue Geweke sur {@code [t0 + floor(f (n - t0)), n)} avec f = 0.5 : le biais
         * accepte n'y contribue plus, tandis qu'une vraie derive lente (mode de relaxation
         * bien plus lent que tau_int : grossissement de domaines, systeme vitreux) continue de
         * se manifester comme une difference entre les deux fenetres tardives.</p>
         *
         * <p>f = 0 restaure le comportement historique (Geweke a partir de t0).</p>
         */
        public double gewekeSkipFraction = 0.5;
        /** Nombre minimal d'echantillons effectifs dans la queue. */
        public double minNeff = 50.0;
        /**
         * Taille minimale de la fenetre de Geweke (apres application de
         * {@link #gewekeSkipFraction}) en deca de laquelle le test est declare
         * <b>non evaluable</b>.
         *
         * <p>Pourquoi. Le z de Geweke compare la fraction initiale A (10 % de la fenetre) a la
         * fraction finale B (50 %). Sur une fenetre de 25 echantillons, A ne contient que
         * {@code max(2, floor(0.1 * 25)) = 2} points : sa variance est estimee sur 2 valeurs et
         * l'inefficacite statistique n'y est evidemment pas estimable. Le z resultant n'est plus
         * du tout N(0,1) et rejette ~24 % des series i.i.d. parfaitement stationnaires. Or ces
         * fenetres minuscules sont exactement celles des <i>premieres</i> verifications, ou une
         * fausse alarme coute cher : elle repousse le verrou au seuil geometrique suivant.</p>
         *
         * <p>Comportement. Si {@code (n - t0 - skip) < minGewekeSamples}, le canal <b>ne passe
         * pas</b> cette verification (z est rapporte NaN) : la thermalisation continue
         * simplement jusqu'a disposer d'assez d'echantillons. C'est le choix conservateur, et
         * il ne peut pas bloquer une serie reellement equilibree : le critere
         * {@code N_eff >= minNeff} (50 par defaut) implique deja une queue d'au moins 50
         * echantillons, donc une fenetre de Geweke d'au moins 25 avec f = 0.5 &mdash; le seuil
         * de 20 est atteint bien avant que les autres criteres ne le soient.</p>
         *
         * <p>Valeur &lt;= 0 : le test est toujours evalue (comportement historique).</p>
         */
        public int minGewekeSamples = 20;
        /** Pas entre candidats t0 dans la recherche de Chodera. */
        public int nskip = 5;
        /**
         * Plafond du nombre de candidats t0 examines a chaque detection.
         *
         * <p>La recherche de Chodera coute O((n/nskip) . n . sqrt(tau)) par canal et par
         * verification. Avec un {@code nskip} fixe ce cout croit en O(n^2 sqrt(tau)) : sur une
         * temperature qui ne converge pas, n atteint maxSweeps/sampleEvery (20 000 echantillons
         * pour les valeurs par defaut) et une seule verification deviendrait aussi chere que
         * les sweeps Monte Carlo qu'elle est censee surveiller.</p>
         *
         * <p>On utilise donc un pas effectif
         * {@code nskipEff = max(nskip, n / maxCandidates)} : le nombre de candidats est borne
         * par ~{@code maxCandidates} et chaque verification reste en
         * O(maxCandidates . n . sqrt(tau)). La resolution sur t0 se degrade proportionnellement
         * a n (elle vaut n/maxCandidates echantillons), ce qui est sans consequence : t0 n'a
         * pas besoin d'etre connu a mieux que ~0.5 % de la serie.</p>
         */
        public int maxCandidates = 200;
        /** Nombre de canaux observes (p.ex. 2 : energie et |M|). Tous doivent converger. */
        public int channels = 2;
        /** Cible de N_eff pour la phase de production. */
        public double targetNeffProduction = 400.0;
        /** Plancher de sweeps de production. */
        public int minProductionSweeps = 1_000;
        /** Plafond de sweeps de production. */
        public int maxProductionSweeps = 400_000;
        /**
         * Cadence de mesure (en sweeps) pendant la production ; sert de <b>plafond</b> lorsque
         * {@link #adaptiveMeasureInterval} est actif.
         */
        public int productionMeasureEvery = 100;
        /**
         * Adapte la cadence de mesure au temps d'auto-correlation mesure.
         *
         * <p>Mesurer tous les ~2 tau rend les mesures successives quasi independantes
         * (g_Delta = 1 + 2 tau/Delta = 2), si bien que N_eff = 400 coute ~400 . 4 tau sweeps.
         * Loin de Tc, tau vaut 1 a 2 sweeps : la production tombe a ~1 a 3 k sweeps au lieu
         * des 40 000 fixes de la version historique. Pres de Tc, tau diverge et la production
         * depasse automatiquement 40 000. C'est exactement la concentration du temps de calcul
         * recherchee.</p>
         *
         * <p>Le plancher {@link #minMeasureEvery} existe parce qu'une mesure n'est pas
         * gratuite : energie + aimantation coutent environ un demi-sweep, donc mesurer a
         * chaque sweep augmenterait le cout total de ~50 % sans gagner d'information.</p>
         */
        public boolean adaptiveMeasureInterval = true;
        /** Cadence de mesure minimale (en sweeps) : une mesure coute ~0.5 sweep. */
        public int minMeasureEvery = 2;
        /**
         * Plancher <b>prudent</b> de tau utilise par le planificateur de production, exprime
         * en fraction de {@code sampleEvery} : le planificateur travaille avec
         * {@code tauPlan = max(tauSweeps, tauPlanFloorFraction * sampleEvery)}.
         *
         * <p>Raison. La serie de thermalisation est echantillonnee tous les Delta_s sweeps,
         * donc l'auto-correlation n'y est resolue que pour des retards multiples de Delta_s.
         * Un tau_echantillons mesure a ~0 ne signifie <b>pas</b> "tau = 0 sweep" mais
         * seulement "tau non resolu, quelque part sous ~Delta_s/3" : le premier point
         * mesurable, rho(Delta_s), est deja tombe dans le bruit. Planifier la production sur
         * tau = 0 reviendrait a affirmer que des mesures espacees de 2 sweeps sont
         * independantes, ce que rien dans les donnees ne soutient.</p>
         *
         * <p>Avec la valeur par defaut 0.25 et Delta_s = 10, le planificateur suppose donc
         * au moins tau = 2.5 sweeps : Delta = ceil(2 tau) = 5 et
         * production = 400 . (5 + 5) = 4 000 sweeps pour N_eff = 400. C'est encore 20 fois
         * moins que les 40 000 sweeps fixes de la version historique, tout en refusant le
         * plancher absurde qu'aurait donne tau = 0 (Delta = 2, production = 1 000).</p>
         *
         * <p>Valeur &lt;= 0 : plus de plancher (le planificateur croit tau sur parole).</p>
         */
        public double tauPlanFloorFraction = 0.25;

        public Config minSweeps(int v) { this.minSweeps = v; return this; }
        public Config maxSweeps(int v) { this.maxSweeps = v; return this; }
        public Config sampleEvery(int v) { this.sampleEvery = v; return this; }
        public Config checkEvery(int v) { this.checkEvery = v; return this; }
        public Config checkGrowth(double v) { this.checkGrowth = v; return this; }
        public Config tauMultiplier(double v) { this.tauMultiplier = v; return this; }
        public Config gewekeZmax(double v) { this.gewekeZmax = v; return this; }
        public Config gewekeSkipFraction(double v) { this.gewekeSkipFraction = v; return this; }
        public Config minNeff(double v) { this.minNeff = v; return this; }
        public Config minGewekeSamples(int v) { this.minGewekeSamples = v; return this; }
        public Config nskip(int v) { this.nskip = v; return this; }
        public Config maxCandidates(int v) { this.maxCandidates = v; return this; }
        public Config channels(int v) { this.channels = v; return this; }
        public Config targetNeffProduction(double v) { this.targetNeffProduction = v; return this; }
        public Config minProductionSweeps(int v) { this.minProductionSweeps = v; return this; }
        public Config maxProductionSweeps(int v) { this.maxProductionSweeps = v; return this; }
        public Config productionMeasureEvery(int v) { this.productionMeasureEvery = v; return this; }
        public Config adaptiveMeasureInterval(boolean v) { this.adaptiveMeasureInterval = v; return this; }
        public Config minMeasureEvery(int v) { this.minMeasureEvery = v; return this; }
        public Config tauPlanFloorFraction(double v) { this.tauPlanFloorFraction = v; return this; }

        /** Copie profonde (le POJO ne contient que des scalaires). */
        public Config copy() {
            Config c = new Config();
            c.minSweeps = minSweeps;
            c.maxSweeps = maxSweeps;
            c.sampleEvery = sampleEvery;
            c.checkEvery = checkEvery;
            c.checkGrowth = checkGrowth;
            c.tauMultiplier = tauMultiplier;
            c.gewekeZmax = gewekeZmax;
            c.gewekeSkipFraction = gewekeSkipFraction;
            c.minNeff = minNeff;
            c.minGewekeSamples = minGewekeSamples;
            c.nskip = nskip;
            c.maxCandidates = maxCandidates;
            c.channels = channels;
            c.targetNeffProduction = targetNeffProduction;
            c.minProductionSweeps = minProductionSweeps;
            c.maxProductionSweeps = maxProductionSweeps;
            c.productionMeasureEvery = productionMeasureEvery;
            c.adaptiveMeasureInterval = adaptiveMeasureInterval;
            c.minMeasureEvery = minMeasureEvery;
            c.tauPlanFloorFraction = tauPlanFloorFraction;
            return c;
        }
    }

    /**
     * Instantane de l'etat du controleur.
     *
     * @param sweepsObserved nombre de sweeps <b>ecoules</b> (voir {@link #sweepsObserved()}).
     * @param gewekeZ copie defensive des scores z par canal.
     */
    public record Status(long sweepsObserved, int samples, boolean equilibrated, boolean converged,
                         long t0Sweep, double tauIntSweeps, double statisticalInefficiency,
                         double neff, double[] gewekeZ, int productionMeasureEvery,
                         int productionSweeps) { }

    private final Config cfg;

    private double[][] series;
    private int n;
    private long firstSweep;
    private long lastSweep;
    /**
     * Origine du temps : {@code firstSweep - sampleEvery}, capturee au premier
     * {@code observe()} suivant {@link #reset()}. Voir la section "Origine des sweeps" de la
     * javadoc de classe.
     */
    private long sweepOrigin;

    private final int[] t0Index;
    private final double[] gCh;
    private final double[] tauCh;
    private final double[] neffCh;
    private final double[] zCh;

    private boolean equilibrated;
    private boolean converged;
    private int lastDetectionSamples;
    private int nextThreshold;
    private int detectionsRun;

    public AdaptiveThermalization(Config config) {
        if (config == null) throw new IllegalArgumentException("config null");
        this.cfg = config.copy();
        if (this.cfg.channels < 1) this.cfg.channels = 1;
        if (this.cfg.sampleEvery < 1) this.cfg.sampleEvery = 1;
        if (this.cfg.checkEvery < 1) this.cfg.checkEvery = 1;
        if (this.cfg.nskip < 1) this.cfg.nskip = 1;
        // Ces deux incoherences rendraient la boucle de thermalisation de l'appelant infinie :
        // isEquilibrated() ne peut devenir vrai que dans observe(), qui n'est atteint qu'aux
        // multiples de sampleEvery. On refuse (comme AdaptiveCoolingSchedule) plutot que de borner.
        if (this.cfg.sampleEvery > this.cfg.maxSweeps) {
            throw new IllegalArgumentException("sampleEvery (" + this.cfg.sampleEvery
                    + ") > maxSweeps (" + this.cfg.maxSweeps + ") : aucun echantillon ne serait pris");
        }
        if (this.cfg.minSweeps > this.cfg.maxSweeps) {
            throw new IllegalArgumentException("minSweeps (" + this.cfg.minSweeps
                    + ") > maxSweeps (" + this.cfg.maxSweeps + ")");
        }
        int ch = this.cfg.channels;
        this.series = new double[ch][256];
        this.t0Index = new int[ch];
        this.gCh = new double[ch];
        this.tauCh = new double[ch];
        this.neffCh = new double[ch];
        this.zCh = new double[ch];
        reset();
    }

    /** Constructeur avec configuration par defaut. */
    public AdaptiveThermalization() { this(new Config()); }

    /** Remet le controleur a zero (les tableaux gardent leur capacite). */
    public void reset() {
        n = 0;
        firstSweep = 0L;
        lastSweep = 0L;
        sweepOrigin = 0L;
        equilibrated = false;
        converged = false;
        lastDetectionSamples = -1;
        nextThreshold = Math.max(1, cfg.checkEvery);
        detectionsRun = 0;
        Arrays.fill(t0Index, 0);
        Arrays.fill(gCh, 1.0);
        Arrays.fill(tauCh, 0.0);
        Arrays.fill(neffCh, 0.0);
        Arrays.fill(zCh, 0.0);
    }

    /** Configuration effective (copie interne : la modifier n'a pas d'effet retroactif). */
    public Config config() { return cfg.copy(); }

    /** Vrai si ce sweep doit produire un echantillon (sweep % sampleEvery == 0). */
    public boolean wantsSample(long sweep) {
        return sweep % cfg.sampleEvery == 0L;
    }

    /**
     * Enregistre un echantillon et, le cas echeant, relance la detection.
     *
     * <p>Les indices {@code sweep} n'ont besoin d'etre que <b>croissants</b> : leur origine
     * est arbitraire et le controleur la neutralise (voir la section "Origine des sweeps"
     * de la javadoc de classe). Les seuils {@code minSweeps} / {@code maxSweeps} portent
     * donc sur {@code sweep - sweepOrigin}, pas sur {@code sweep} brut.</p>
     *
     * @param values une valeur par canal ; {@code values.length} doit valoir
     *               {@code config().channels}.
     * @throws IllegalArgumentException si le nombre de canaux est incorrect.
     * @throws IllegalStateException    si une valeur n'est pas finie (NaN ou +-Inf). Une
     *                                  observable non finie est un bug de la simulation, pas
     *                                  un canal constant : elle doit declencher une alarme et
     *                                  non passer silencieusement les criteres.
     */
    public void observe(long sweep, double... values) {
        if (values == null || values.length != cfg.channels) {
            throw new IllegalArgumentException(
                    "observe() attend " + cfg.channels + " valeurs, recu "
                            + (values == null ? "null" : String.valueOf(values.length)));
        }
        for (int c = 0; c < cfg.channels; c++) {
            if (!Double.isFinite(values[c])) {
                throw new IllegalStateException("observable non finie (canal " + c
                        + ", sweep " + sweep + ") : " + values[c]);
            }
        }
        if (n > 0 && sweep - lastSweep != cfg.sampleEvery) {
            // Tout le controleur (t0Sweep, tauSweeps, criteres de longueur) suppose une grille
            // reguliere de pas sampleEvery : on l'impose au lieu de le supposer.
            throw new IllegalArgumentException("observe() attend des sweeps espaces de sampleEvery="
                    + cfg.sampleEvery + " (utiliser wantsSample) : " + lastSweep + " -> " + sweep);
        }
        ensureCapacity(n + 1);
        for (int c = 0; c < cfg.channels; c++) series[c][n] = values[c];
        if (n == 0) {
            firstSweep = sweep;
            // Le premier echantillon clot le premier intervalle de sampleEvery sweeps :
            // l'origine est donc un sampleEvery avant lui, et un appelant qui compte a
            // partir de 1 avec des mesures a 10, 20, ... obtient ecoule == sweep.
            sweepOrigin = sweep - cfg.sampleEvery;
        }
        n++;
        lastSweep = sweep;

        if (equilibrated) return;

        long elapsed = sweep - sweepOrigin;

        if (maybeRunDetection(elapsed)) {
            equilibrated = true;
            converged = true;
            return;
        }

        if (elapsed >= cfg.maxSweeps) {
            // La detection lancee ici pour rafraichir les estimations peut parfaitement
            // reussir : la cadence geometrique espace les verifications, et rien ne garantit
            // que la derniere tombe sur un echantillon de controle. Forcer converged = false
            // sans regarder son verdict signalait alors a tort une thermalisation coupee.
            // Si une detection a deja tourne sur ce meme n (via maybeRunDetection), elle a
            // echoue -- sinon on serait deja sorti plus haut -- et converged reste false.
            converged = (lastDetectionSamples != n) && runDetection();
            equilibrated = true;
        }
    }

    /**
     * Lance une detection si le budget de sweeps <b>et</b> la cadence l'autorisent ;
     * retourne true si tous les criteres passent.
     *
     * <p>La garde {@code elapsed < minSweeps} est un retour anticipe explicite et non un
     * membre gauche de {@code &&} : {@link #takeCheckSlot()} <i>consomme</i> un creneau de
     * verification (il avance {@code nextThreshold}), et faire dependre cette mutation de
     * l'ordre d'evaluation d'un court-circuit est exactement le genre de couplage qui casse
     * a la premiere reecriture de la condition.</p>
     */
    private boolean maybeRunDetection(long elapsed) {
        if (elapsed < cfg.minSweeps) return false;
        if (!takeCheckSlot()) return false;
        return runDetection();
    }

    /**
     * Nombre d'echantillons de la queue exclus du test de Geweke (voir
     * {@link Config#gewekeSkipFraction}). Borne de facon a laisser au moins 8 echantillons.
     */
    private int gewekeSkip(int tailLen) {
        double f = cfg.gewekeSkipFraction;
        if (!Double.isFinite(f) || f <= 0.0) return 0;
        if (f > 0.9) f = 0.9;
        int skip = (int) Math.floor(f * tailLen);
        return Math.max(0, Math.min(skip, tailLen - 8));
    }

    /**
     * <b>Consomme</b> un creneau de verification si le nombre d'echantillons courant en
     * declenche un, en avancant le seuil ; retourne false sinon (et ne mute alors rien).
     *
     * <p>Le nom porte la mutation : cette methode n'est pas un predicat pur, et ne doit
     * jamais etre appelee derriere un court-circuit (voir {@link #maybeRunDetection}).</p>
     *
     * <p>Cadence geometrique : seuils {@code checkEvery}, puis {@code ceil(n * checkGrowth)}
     * calcule a partir du n <b>reel</b> de la verification qui vient de se produire. Ainsi,
     * si {@code minSweeps} n'est atteint qu'apres le premier seuil, la verification a lieu au
     * premier echantillon suivant {@code minSweeps} et la suite geometrique repart de la.</p>
     */
    private boolean takeCheckSlot() {
        if (!(cfg.checkGrowth > 1.0)) {
            return n % cfg.checkEvery == 0;          // cadence fixe historique
        }
        if (n < nextThreshold) return false;
        int next = (int) Math.ceil((double) n * cfg.checkGrowth);
        nextThreshold = Math.max(n + 1, next);       // strictement croissant
        return true;
    }

    /**
     * Execute la detection sur tous les canaux ; retourne true si tous les criteres passent.
     *
     * <p>Le pas entre candidats t0 est plafonne (voir {@link Config#maxCandidates}) afin que
     * le cout d'une verification ne croisse pas quadratiquement avec la longueur de la serie.</p>
     */
    private boolean runDetection() {
        lastDetectionSamples = n;
        detectionsRun++;
        boolean all = true;
        int nskipEff = Math.max(cfg.nskip, n / Math.max(1, cfg.maxCandidates));
        for (int c = 0; c < cfg.channels; c++) {
            double[] s = series[c];
            EquilibrationDetector.Result r = EquilibrationDetector.detectEquilibration(s, 0, n, nskipEff);

            if (r.constant()) {
                // Canal numeriquement constant : il ne porte aucune information temporelle,
                // donc aucun critere ne peut y etre evalue. Il passe (voir javadoc de classe).
                t0Index[c] = 0;
                gCh[c] = 1.0;
                tauCh[c] = 0.0;
                neffCh[c] = n;
                zCh[c] = 0.0;
                continue;
            }

            t0Index[c] = r.t0();
            gCh[c] = r.g();
            tauCh[c] = r.tauInt();
            neffCh[c] = r.neff();

            // Fenetre de Geweke : seconde moitie de la queue, et seulement si elle est assez
            // longue pour que z soit interpretable (voir Config#minGewekeSamples). Sinon on
            // rapporte NaN, ce qui fait echouer |z| <= zmax sans code supplementaire : le
            // canal ne passe pas et la thermalisation continue.
            int tailLen = n - r.t0();
            int skip = gewekeSkip(tailLen);
            zCh[c] = (tailLen - skip >= cfg.minGewekeSamples)
                    // Fenetres 25 % / 50 % (et non 10 % / 50 %) : avec la fenetre minimale de
                    // minGewekeSamples = 20 points, 10 % ne donnerait que nA = 2 points.
                    ? EquilibrationDetector.gewekeZ(s, r.t0() + skip, n, 0.25, 0.5)
                    : Double.NaN;

            double tauSweeps = tauSweepsFromSamples(r.tauInt(), cfg.sampleEvery);
            double tailSweeps = (double) tailLen * cfg.sampleEvery;
            boolean okLength = tailSweeps >= cfg.tauMultiplier * tauSweeps;
            boolean okGeweke = Math.abs(zCh[c]) <= cfg.gewekeZmax;   // NaN -> false
            boolean okNeff = r.neff() >= cfg.minNeff;
            if (!(okLength && okGeweke && okNeff)) all = false;
        }
        return all;
    }

    private void ensureCapacity(int need) {
        if (need <= series[0].length) return;
        int cap = series[0].length;
        while (cap < need) cap <<= 1;
        for (int c = 0; c < cfg.channels; c++) {
            series[c] = Arrays.copyOf(series[c], cap);
        }
    }

    // ------------------------------------------------------------------
    // Accesseurs
    // ------------------------------------------------------------------

    /**
     * Vrai lorsque la thermalisation est declaree terminee : soit tous les criteres sont
     * passes ({@link #converged()} = true), soit {@code maxSweeps} a ete atteint
     * ({@link #converged()} = false).
     */
    public boolean isEquilibrated() { return equilibrated; }

    /** Vrai si l'equilibration a ete <i>detectee</i> et non simplement imposee par maxSweeps. */
    public boolean converged() { return converged; }

    /** Sweep a partir duquel la serie est consideree equilibree (max sur les canaux). */
    public long t0Sweep() {
        int worst = 0;
        for (int c = 0; c < cfg.channels; c++) worst = Math.max(worst, t0Index[c]);
        return firstSweep + (long) worst * cfg.sampleEvery;
    }

    /**
     * Temps d'auto-correlation integre exprime en <b>sweeps</b> (max sur les canaux),
     * convention de Chodera (tau = somme_{t &gt;= 1} rho(t), g = 1 + 2 tau).
     *
     * <p>Voir {@link #tauSweepsFromSamples} : la conversion n'est <i>pas</i> une
     * multiplication par {@code sampleEvery}.</p>
     */
    public double tauIntSweeps() {
        double worst = 0.0;
        for (int c = 0; c < cfg.channels; c++) worst = Math.max(worst, tauCh[c]);
        return tauSweepsFromSamples(worst, cfg.sampleEvery);
    }

    /**
     * Convertit un tau de Chodera exprime en <b>echantillons</b> (espacement
     * Delta_s = {@code sampleEvery} sweeps) en un tau de Chodera exprime en <b>sweeps</b>.
     *
     * <h3>Pourquoi ce n'est pas une multiplication</h3>
     * <p>{@code tauCh[c]} est un tau <i>discret</i> : tau_s = somme_{k &gt;= 1} rho(k Delta_s).
     * Multiplier par Delta_s revient a approximer cette somme par une integrale, ce qui n'est
     * legitime que pour tau &gt;&gt; Delta_s. Des que tau est de l'ordre de Delta_s l'erreur
     * atteint un facteur 2, et le planificateur de production se trompe d'autant sur N_eff.</p>
     *
     * <p>Pour une auto-correlation exponentielle rho(t) = exp(-t / tau_e) (tau_e en sweeps) :</p>
     * <pre>
     *   tau_s = somme_{k &gt;= 1} exp(-k Delta_s / tau_e) = 1 / (exp(Delta_s / tau_e) - 1)
     *   =&gt;  tau_e = Delta_s / log1p(1 / tau_s)
     *   et le tau integre en sweeps vaut  tau_sw = 1 / (exp(1 / tau_e) - 1)  (~ tau_e - 1/2).
     * </pre>
     * <p>Exemples pour Delta_s = 10 : tau_s = 100 -&gt; 1004 sweeps (la multiplication naive
     * donnait 1000, erreur 0.4 %) ; tau_s = 1 -&gt; 13.9 sweeps (naive : 10, erreur 40 %) ;
     * tau_s = 0.1 -&gt; 3.7 sweeps (naive : 1, erreur x3.7). La conversion exacte est donc
     * systematiquement <i>plus prudente</i> que la multiplication, ce qui est le bon sens
     * de l'erreur pour un critere de longueur de queue.</p>
     *
     * <p>Attention : tau_s ~ 0 ne signifie pas "tau nul" mais "tau non resolu sous
     * ~Delta_s/3". Le planificateur de production applique pour cela un plancher separe,
     * {@link Config#tauPlanFloorFraction}.</p>
     *
     * @param tauSamples tau de Chodera en echantillons ; &lt;= 0 (ou NaN) -&gt; 0.
     * @param sampleEvery espacement Delta_s en sweeps (&lt; 1 traite comme 1).
     * @return tau de Chodera en sweeps, &gt;= 0 (+Inf si {@code tauSamples} est +Inf).
     */
    public static double tauSweepsFromSamples(double tauSamples, int sampleEvery) {
        if (!(tauSamples > 0.0)) return 0.0;              // 0, negatif, NaN
        int ds = Math.max(1, sampleEvery);
        if (ds == 1) return tauSamples;                   // deja en sweeps (la formule est exacte,
                                                          // ce raccourci evite deux transcendantes)
        double tauE = ds / Math.log1p(1.0 / tauSamples);  // constante exponentielle, en sweeps
        if (!(tauE > 0.0)) return 0.0;                    // tauSamples denormal -> log1p(Inf) = Inf
        if (!Double.isFinite(tauE)) return (double) ds * tauSamples;   // 1/tauSamples a sous-deborde
        double d = Math.expm1(1.0 / tauE);                // expm1 : pas d'annulation pour 1/tauE -> 0
        if (!(d > 0.0)) return tauE;                      // 1/tauE a sous-deborde : tau_sw ~ tauE
        double tauSw = 1.0 / d;
        return Double.isFinite(tauSw) ? tauSw : tauE;
    }

    /** Inefficacite statistique g = 1 + 2 tau (max sur les canaux). */
    public double statisticalInefficiency() {
        double worst = 1.0;
        for (int c = 0; c < cfg.channels; c++) worst = Math.max(worst, gCh[c]);
        return worst;
    }

    /** Nombre d'echantillons effectifs de la queue (min sur les canaux). */
    public double neff() {
        double worst = Double.POSITIVE_INFINITY;
        for (int c = 0; c < cfg.channels; c++) worst = Math.min(worst, neffCh[c]);
        return Double.isFinite(worst) ? worst : 0.0;
    }

    /** Dernier score z de Geweke du canal {@code ch}. */
    public double gewekeZ(int ch) { return zCh[ch]; }

    /**
     * Nombre de sweeps <b>ecoules</b> depuis le dernier {@link #reset()}, c'est-a-dire
     * {@code dernierSweep - sweepOrigin} (0 si aucun echantillon). C'est la quantite
     * comparee a {@code minSweeps} / {@code maxSweeps}, et elle ne depend pas de l'origine
     * arbitraire choisie par l'appelant.
     */
    public long sweepsObserved() { return n == 0 ? 0L : lastSweep - sweepOrigin; }

    /** Nombre d'echantillons stockes. */
    public int samples() { return n; }

    /**
     * Nombre de detections effectivement executees depuis le dernier {@link #reset()}.
     * Utile pour verifier que la cadence geometrique borne bien le cout de surveillance.
     */
    public int detectionsRun() { return detectionsRun; }

    /** Copie de la serie du canal {@code ch} (longueur {@link #samples()}). */
    public double[] channelSeries(int ch) { return Arrays.copyOf(series[ch], n); }

    /**
     * Nombre de sweeps de production recommande.
     *
     * <p>Derivation : si l'on mesure tous les Delta sweeps, les mesures successives sont
     * correlees avec une inefficacite g_Delta = 1 + 2 tau / Delta (tau exprime en sweeps).
     * Avec nMeas mesures on obtient N_eff = nMeas / g_Delta, donc le nombre de sweeps vaut</p>
     * <pre>
     *   sweeps = nMeas . Delta = N_eff . g_Delta . Delta = N_eff . (Delta + 2 tau)
     * </pre>
     * <p>d'ou {@code ceil(targetNeffProduction * (productionMeasureEvery() + 2 tauPlan))}
     * ou Delta = {@link #productionMeasureEvery()} est la cadence <b>effective</b>,
     * borne a [minProductionSweeps, maxProductionSweeps]. Pres de Tc, tau diverge et le
     * terme 2 tau domine : la production s'allonge automatiquement, alors que loin de Tc
     * on retombe sur N_eff . Delta.</p>
     *
     * <p>tauPlan est le tau <b>plancher</b> du planificateur,
     * {@code max(tauIntSweeps(), tauPlanFloorFraction * sampleEvery)} : voir
     * {@link Config#tauPlanFloorFraction} pour la raison (un tau mesure a ~0 signifie
     * "non resolu sous ~Delta_s/3", pas "nul").</p>
     */
    public int productionSweeps() { return productionSweeps(cfg, tauIntSweeps()); }

    /**
     * Cadence de mesure effective pendant la production, en sweeps :
     * {@code adaptiveMeasureInterval ? clamp(ceil(2 tauPlan), minMeasureEvery, productionMeasureEvery)
     * : productionMeasureEvery}, ou tauPlan applique le plancher
     * {@link Config#tauPlanFloorFraction}.
     *
     * <p>Mesurer tous les 2 tau donne g_Delta = 1 + 2 tau/Delta = 2 : deux mesures
     * consecutives sont quasi independantes, ce qui est le meilleur compromis entre le cout
     * des mesures et le gain d'information. Voir {@link Config#adaptiveMeasureInterval}.</p>
     */
    public int productionMeasureEvery() { return productionMeasureEvery(cfg, tauIntSweeps()); }

    /**
     * tau de planification : {@code max(tauIntSweeps, tauPlanFloorFraction * sampleEvery)},
     * les valeurs non finies ou negatives etant traitees comme 0.
     *
     * <p>Ce plancher est le seul endroit ou le planificateur se protege du fait qu'un tau
     * mesure sur une serie echantillonnee tous les Delta_s sweeps n'a de sens que jusqu'a
     * ~Delta_s/3 : en dessous, "tau = 0" veut dire "non resolu". Voir
     * {@link Config#tauPlanFloorFraction}.</p>
     */
    private static double tauPlan(Config cfg, double tauIntSweeps) {
        double t = (Double.isFinite(tauIntSweeps) && tauIntSweeps > 0.0) ? tauIntSweeps : 0.0;
        double f = cfg.tauPlanFloorFraction;
        if (!Double.isFinite(f) || f <= 0.0) return t;
        return Math.max(t, f * Math.max(1, cfg.sampleEvery));
    }

    /** Version pure de {@link #productionMeasureEvery()} pour un tau hypothetique (planification). */
    public static int productionMeasureEvery(Config cfg, double tauIntSweeps) {
        if (!cfg.adaptiveMeasureInterval) return cfg.productionMeasureEvery;
        double t = tauPlan(cfg, tauIntSweeps);
        long d = (long) Math.ceil(2.0 * t);
        if (d < cfg.minMeasureEvery) d = cfg.minMeasureEvery;
        if (d > cfg.productionMeasureEvery) d = cfg.productionMeasureEvery;
        return (int) d;
    }

    /** Version pure de {@link #productionSweeps()} pour un tau hypothetique (planification). */
    public static int productionSweeps(Config cfg, double tauIntSweeps) {
        double t = tauPlan(cfg, tauIntSweeps);
        double v = cfg.targetNeffProduction * (productionMeasureEvery(cfg, tauIntSweeps) + 2.0 * t);
        if (!Double.isFinite(v) || v < 0.0) v = cfg.minProductionSweeps;
        long s = (long) Math.ceil(v);
        if (s < cfg.minProductionSweeps) s = cfg.minProductionSweeps;
        if (s > cfg.maxProductionSweeps) s = cfg.maxProductionSweeps;
        return (int) s;
    }

    /** Instantane complet de l'etat (copie defensive des z). */
    public Status status() {
        return new Status(sweepsObserved(), n, equilibrated, converged, t0Sweep(), tauIntSweeps(),
                statisticalInefficiency(), neff(), Arrays.copyOf(zCh, zCh.length),
                productionMeasureEvery(), productionSweeps());
    }

    /** Resume sur une ligne, destine aux journaux de simulation. */
    public String report() {
        StringBuilder sb = new StringBuilder(160);
        sb.append("AdaptiveThermalization[sweeps=").append(sweepsObserved())
          .append(" n=").append(n)
          .append(" t0=").append(t0Sweep())
          .append(String.format(" tau=%.1f sw", tauIntSweeps()))
          .append(String.format(" g=%.2f", statisticalInefficiency()))
          .append(String.format(" neff=%.1f", neff()))
          .append(" z=");
        for (int c = 0; c < cfg.channels; c++) {
            if (c > 0) sb.append(',');
            sb.append(String.format("%.2f", zCh[c]));
        }
        sb.append(equilibrated ? " EQUILIBRE" : " en cours")
          .append(converged ? " (converge)" : (equilibrated ? " (NON converge: maxSweeps)" : ""))
          .append(" mesure=").append(productionMeasureEvery()).append("sw")
          .append(" prod=").append(productionSweeps())
          .append(']');
        return sb.toString();
    }
}
