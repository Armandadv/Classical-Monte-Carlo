import java.io.File;
import java.nio.file.Path;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.random.RandomGenerator;
import java.util.random.RandomGeneratorFactory;

/**
 * Point d'entree parallele : construction de replique via {@link ReplicaFactory},
 * echantillonnage par {@link CheckerboardMetropolis} orchestre par {@link MonteCarlo#runParallel}.
 *
 * <h2>Ordonnancement hybride</h2>
 * <p>Deux niveaux de parallelisme coexistent :</p>
 * <ul>
 *   <li><b>entre repliques</b> : sans aucune synchronisation, mais chaque replique porte son
 *       propre reseau et son propre S(Q,w) — c'est ce niveau qui sature le bus memoire ;</li>
 *   <li><b>a l'interieur d'un balayage</b> : {@code threadsPerReplica} threads se partagent une
 *       classe de couleur, avec une barriere par classe — cout de synchronisation, mais un seul
 *       jeu de donnees en cache.</li>
 * </ul>
 * <p>On lance donc {@code concurrentReplicas = totalThreads / threadsPerReplica} repliques a la
 * fois. Les workers de {@link CheckerboardMetropolis} attendent par {@code LockSupport.parkNanos}
 * apres une phase d'attente active : les sur-souscrire degrade fortement le debit. Le produit
 * {@code threadsPerReplica * concurrentReplicas} ne doit donc jamais depasser le nombre de coeurs
 * reellement disponibles. Utiliser {@code --bench} pour choisir {@code threadsPerReplica}.</p>
 *
 * <h2>Usage</h2>
 * <pre>
 * java -cp "bin:lib/*" AppParallel &lt;inputDir&gt; &lt;outputDir&gt; [totalThreads] [threadsPerReplica]
 *        [La] [Lb] [hStart] [hEnd] [hStep] [seed] [options]
 * java -cp "bin:lib/*" AppParallel --bench La Lb [J1 J2 J3]
 * java -cp "bin:lib/*" AppParallel --bench &lt;fichierInput&gt; La Lb
 * </pre>
 * <p>La seconde forme de {@code --bench} mesure le debit sur le <b>modele reellement lu</b> dans
 * {@code fichierInput} (recette {@link ReplicaFactory}, comme le mode production), au lieu du
 * modele honeycomb synthetique J1/J2/J3 de la premiere forme : c'est la coloration et le nombre de
 * voisins par site du fichier d'entree qui determinent le debit, pas ceux d'un honeycomb generique.
 * Les deux formes se distinguent par le premier argument : un entier declenche la forme synthetique,
 * tout ce qui n'en est pas un (un chemin de fichier) declenche la forme fichier.</p>
 *
 * <h2>Codes de retour</h2>
 * <p>{@link #run(String[])} renvoie 0 en cas de succes, 1 si au moins une replique a echoue,
 * 2 pour une erreur d'usage. Elle ne leve rien et n'appelle jamais {@code System.exit} : c'est
 * {@link #main(String[])} qui traduit le code en sortie de processus, ce qui permet aux tests
 * d'appeler {@code run} directement.</p>
 */
public final class AppParallel {

    private AppParallel() {}

    /** Code de retour : tout s'est bien passe. */
    public static final int EXIT_OK = 0;
    /** Code de retour : au moins une replique a echoue (le detail est imprime). */
    public static final int EXIT_FAILURE = 1;
    /** Code de retour : erreur d'usage (arguments invalides). */
    public static final int EXIT_USAGE = 2;

    /**
     * Garde-fou sur le nombre de valeurs de champ. Une faute de frappe sur {@code hStep}
     * ({@code 0.0000125} au lieu de {@code 0.0125}) engendrerait 80 000 fois plus de repliques
     * que voulu, chacune allouant un reseau : mieux vaut une erreur d'usage immediate qu'un
     * {@code OutOfMemoryError} une heure plus tard.
     */
    private static final int MAX_FIELDS = 100_000;

    /** Nombre de tranches par couleur (voir {@code replica}) : reproductibilite independante des threads. */
    static final int CHUNKS_PER_COLOR = 16;

    private static final String USAGE = String.join("\n",
        "Usage :",
        "  java -cp \"bin:lib/*\" AppParallel <inputDir> <outputDir> [totalThreads] [threadsPerReplica]",
        "         [La] [Lb] [hStart] [hEnd] [hStep] [seed] [options]",
        "  java -cp \"bin:lib/*\" AppParallel --bench <La> <Lb> [J1 J2 J3]",
        "  java -cp \"bin:lib/*\" AppParallel --bench <fichierInput> <La> <Lb>",
        "",
        "Positionnels (valeurs par defaut entre parentheses) :",
        "  inputDir           repertoire contenant les fichiers bcaoExplor_*",
        "  outputDir          repertoire de sortie (cree au besoin)",
        "  totalThreads       (Runtime.availableProcessors()) coeurs a utiliser au total",
        "  threadsPerReplica  (4)     threads par replique ; cf. --bench",
        "  La, Lb             (30 30) taille du reseau (x2 sous-reseaux)",
        "  hStart hEnd hStep  (0.0 1.0 0.0125) boucle en champ, h dans [hStart, hEnd[ ; actifs",
        "                     UNIQUEMENT avec --field (voir Options) ; sans --field, ni lus ni",
        "                     valides (1 seule replique a h=0)",
        "  seed               (12345) graine maitresse",
        "",
        "Options (chaque paire --x / --no-x est imprimee avec sa valeur effective par le",
        "resume avant tout calcul) :",
        "  --sqw / --no-sqw   (defaut : --sqw) evalue (ou non) S(Q,w) a la derniere temperature",
        "  --field / --no-field  (defaut : --no-field) RUPTURE de comportement : par defaut,",
        "                     PAS de boucle en champ (1 seule replique a h=0) ; --field retablit",
        "                     le balayage h_k = hStart + k*hStep, k = 0..nH-1 (avant cette",
        "                     option, 80 champs etaient balayes par defaut). Sans --field,",
        "                     hStart/hEnd/hStep sont ignores.",
        "  --annealing / --no-annealing  (defaut : --annealing) execute (ou non) le recuit",
        "                     (thermalisation + production). --no-annealing necessite --ac",
        "                     (susceptibilite AC, voir le bloc dedie plus bas) : sans mesure",
        "                     alternative -> EXIT_USAGE.",
        "  --adaptive-cooling / --no-adaptive-cooling  (defaut : --no-adaptive-cooling, soit le",
        "                     refroidissement geometrique 0.95) ; --cooling est un alias de",
        "                     --adaptive-cooling.",
        "  --freeze-sigma / --no-freeze-sigma  (defaut : --freeze-sigma, sigma fige pendant la",
        "                     production) ; --legacy-sigma est un alias de --no-freeze-sigma",
        "                     (comportement historique, sigma adaptatif en production).",
        "  --overrelax N / --no-overrelax  (defaut : --overrelax 2) balayages de sur-relaxation",
        "                     apres chaque balayage Metropolis ; N=0 equivaut a --no-overrelax.",
        "                     La valeur par defaut divise le cout du recuit par ~35.",
        "  --precession-at-end / --no-precession-at-end  (defaut : --precession-at-end) relaxation classique",
        "                     (spin.Precession1D) en fin de chaque temperature",
        "  --stagnation / --no-stagnation  (defaut : --no-stagnation) active (ou non) le",
        "                     detecteur de stagnation (StagnationDetector.Config), exactement ce",
        "                     que --preset screening fait pour ce seul reglage (les autres",
        "                     reglages du preset ne sont pas touches par cette option seule)",
        "  --stop-therm-on-stagnation / --no-stop-therm-on-stagnation  (defaut : desactive)",
        "                     arrete la thermalisation des que le detecteur conclut a la",
        "                     stagnation",
        "  --stop-prod-on-saturation / --no-stop-prod-on-saturation  (defaut : desactive)",
        "                     arrete l'extension de production des que le gain de N_eff sature",
        "                     Ces deux arrets anticipes exigent un detecteur de stagnation actif",
        "                     (--stagnation ou --preset screening) : sinon EXIT_USAGE (aucune",
        "                     activation implicite du detecteur). Le detecteur DIAGNOSTIC installe",
        "                     par --decorrelation-log seul ne compte pas : il ne declenche jamais",
        "                     d'arret, ces options n'auraient donc aucun effet.",
        "  --field-dir x,y,z  direction du champ dans le tenseur g (defaut -1,1,0), vecteur non nul",
        "  --minTherm N       minSweeps du controleur d'equilibration",
        "  --maxTherm N       maxSweeps du controleur d'equilibration",
        "  --targetNeff N     N_eff vise pendant la production",
        "  --wmax X           w_max vise (unites de Spin) : derive measurementRate/N. Absent = historique.",
        "                     X = 0 : mode statique, coupe unique a w = 0 (S(q) thermique).",
        "  --qz X             composante z du vecteur de diffusion (r.l.u.), defaut 0 ;",
        "                     enregistree dans les fichiers de structure (propriete Qzz).",
        "  --initTemp T       temperature initiale du recuit (defaut historique 10).",
        "  --endTemp T        temperature finale du recuit (defaut historique 0.001) ; le bloc",
        "                     S(Q,w) est mesure a cette derniere temperature. Exige > 0 et < initTemp.",
        "  --window W         duree de la fenetre d'integration par instantane S(Q,w) (defaut 20) :",
        "                     fixe la resolution dw = 2 pi/W et ~ W*wmax/pi points dans la bande.",
        "  --domega X         resolution VISEE en energie du S(Q,w), en unites hbarre = 1 :",
        "                     ecriture alternative de --window, W = 2 pi/X (la largeur des raies",
        "                     de Bragg suit X). Incompatible avec --window (n'en fournir qu'un,",
        "                     sinon EXIT_USAGE). Avec --wmax > 0, exiger X < wmax (au moins 2",
        "                     points dans la bande, sinon EXIT_USAGE).",
        "  --nSQW N           nombre d'instantanes du bloc S(Q,w) (defaut : 40 historique) :",
        "                     N+1 mesures spectrales moyennees, Avro ecrit a la (N+1)-ieme ;",
        "                     2 <= N <= 1 000 000, cout du bloc proportionnel a N.",
        "  --seed-random      graine maitresse tiree au hasard a chaque lancement (precede la graine",
        "                     positionnelle) ; la graine effective est affichee au demarrage : la",
        "                     redonner en argument 10 reproduit exactement ce run.",
        "  --spins-per-T      instantane des spins (avro) apres la production de CHAQUE temperature,",
        "                     nomme avec T ; sans l'option, seul l'etat final est sauvegarde.",
        "  --gewekeSkip F     fraction initiale ecartee par le test de Geweke (defaut 0.5)",
        "  --preset P         screening|standard|strict (defaut standard = comportement actuel).",
        "                     screening : criblage rapide, arret anticipe des que le detecteur de",
        "                       stagnation conclut qu'il ne peut plus ameliorer l'etat courant",
        "                       (targetNeff=100, maxTherm=maxProd=50 000, arrets anticipes actifs).",
        "                     standard  : valeurs actuelles, aucun detecteur, aucun arret anticipe.",
        "                     strict    : chaine longue et decorrelee (targetNeff=1000, minNeff=200,",
        "                       tauMultiplier=50, maxTherm=500 000, maxProd=2 000 000), pas de",
        "                       detecteur de stagnation (le budget est suppose suffisant).",
        "                     Les options explicites (--targetNeff, --maxTherm, --stagnation, ...)",
        "                     placees APRES --preset sur la ligne de commande l'emportent sur le",
        "                     preset ; placees avant, c'est le preset qui l'emporte (les deux se",
        "                     mutent dans l'ordre de lecture de la ligne de commande). Le parseur",
        "                     est lineaire et n'a pas de logique de precedence dediee : la meme",
        "                     regle -- le dernier lu gagne -- s'applique entre --preset,",
        "                     --stagnation/--no-stagnation et --decorrelation-log.",
        "  --decorrelation-log  ecrit un journal decorrelation.csv (StagnationLog) a cote du",
        "                     journal d'equilibration : q_inf, tau_q, decorrelation entre",
        "                     temperatures et entre mesures de production. Si aucun detecteur de",
        "                     stagnation n'est deja actif (--stagnation ou --preset screening), en",
        "                     active un en mode DIAGNOSTIC SEUL (le detecteur mesure mais ne",
        "                     declenche jamais d'arret anticipe) : la trajectoire n'est donc jamais",
        "                     modifiee par cette seule option.",
        "  --memoryMultiplier x  (3.0) K_q : seuil thermSweeps >= K_q * tau_q pour declarer deux",
        "                     temperatures decorrelees (voir README_PARALLEL.md).",
        "",
        "Susceptibilite AC (chi'(T,omega), chi''(T,omega), voir ACSusceptibility/ACPoint) :",
        "  --ac / --no-ac     (defaut : --no-ac) calcule, en plus du recuit (ou a sa place avec",
        "                     --no-annealing), la reponse a un champ oscillant H(t) = hStatic +",
        "                     h0.sin(omega t).acDrive. omega = 2 pi / periodSweeps est exprimee en",
        "                     RADIANS PAR SWEEP (le temps est compte en sweeps Metropolis, aucune",
        "                     conversion en Hz n'est faite) : deux runs a des periodes differentes",
        "                     ne se comparent donc que qualitativement, jamais sur une echelle de",
        "                     frequence physique. Chaque periode de --acPeriods refait une echelle",
        "                     de refroidissement COMPLETE depuis --acTempMax (redondant en temps",
        "                     de calcul, necessaire a la correction physique : voir la javadoc de",
        "                     ACSusceptibility). chi'/chi'' sont normalisees par H0 et par site",
        "                     (grandeur intensive, voir ACPoint). Ecrit un CSV separe",
        "                     ac_susceptibility_*.csv (ACSusceptibilityLog), independant de",
        "                     equilibration.csv. --no-annealing exige --ac (sinon EXIT_USAGE,",
        "                     voir --annealing ci-dessus) ; --ac est compatible avec --field (une",
        "                     mesure AC par valeur de champ statique).",
        "  --acReplicas N     (8) repliques AC independantes par (fichier, champ statique) ; le",
        "                     CSV d'un (fichier, champ) n'est ecrit que si TOUTES ses repliques",
        "                     reussissent (sinon EXIT_FAILURE, comme les repliques de recuit).",
        "  --acH0 x           (0.02) amplitude du champ oscillant ; x > 0.",
        "  --acPeriods p1,p2,...  (64,128,256) periodes d'excitation, en sweeps, entiers >= 2,",
        "                     liste non vide.",
        "  --acTransientPeriods N  (5) periodes de regime transitoire ecartees avant la",
        "                     demodulation ; N >= 0.",
        "  --acMeasurePeriods N    (30) periodes de mesure accumulees dans la demodulation ;",
        "                     N >= 1.",
        "  --acTempMax x / --acTempMin x / --acTempStep x  (2.0 / 0.05 / 0.05) grille de",
        "                     temperature AC (lineaire decroissante), independante de la grille",
        "                     de recuit ; il faut acTempMax > acTempMin > 0 et acTempStep > 0.",
        "  --acFieldDir x,y,z  direction du champ oscillant dans le tenseur g (defaut : celle de",
        "                     --field-dir si non precisee) ; vecteur non nul.",
        "  --quiet            supprime la barre de progression et le resume par temperature",
        "  -h, --help         cette aide");

    // ------------------------------------------------------------------ parametres analyses

    private static final class Args {
        String inputDir;
        String outputDir;
        int totalThreads = Runtime.getRuntime().availableProcessors();
        int threadsPerReplica = 4;
        int La = 30, Lb = 30;
        double hStart = 0.0, hEnd = 1.0, hStep = 0.0125;
        long seed = 12345L;

        // ---- (1) --sqw / --no-sqw (defaut : true) -- S(Q,w) a la derniere temperature.
        boolean computeSQW = true;

        // ---- --wmax (defaut : NaN = non fourni -> parametres de mesure historiques).
        // 0 explicite = mode statique : coupe unique a w = 0, soit le facteur de structure
        // statique S(q) (moyenne thermique de |S_perp(q)|^2 sur les configurations MC).
        // Negatif : rejet a la validation. ----
        double wMax = Double.NaN;

        // ---- --initTemp / --endTemp (defaut : NaN = non fourni -> bornes historiques du
        // recuit, 10 / 0.001). La derniere temperature est aussi celle du bloc S(Q,w).
        // <= 0 ou initTemp <= endTemp : rejet a la validation. ----
        double initTemp = Double.NaN;
        double endTemp = Double.NaN;

        // ---- --window (defaut : NaN = non fourni -> fenetre historique 20). Duree de
        // l'integration LL par instantane S(Q,w) : fixe la resolution dw = 2 pi / W et le
        // nombre de points en omega (~ W * wmax / pi dans la bande positive). ----
        double window = Double.NaN;

        // ---- --domega (defaut : NaN = non fourni). Resolution VISEE en energie du S(Q,w),
        // en unites hbarre = 1 : ecriture alternative de --window, W = 2 pi / domega (voir
        // fenetreDepuisDeltaOmega). Converti en a.window a la fin de parse(), apres
        // validation : toute la suite (ParallelOptions.window, sqwBlock) ne voit que la
        // fenetre. Rejets : non finie ou <= 0, concomitance avec --window, et avec
        // --wmax > 0 : domega >= wMax (moins de 2 points dans la bande). ----
        double domega = Double.NaN;

        // ---- --nSQW (defaut : null = 40 historique). Nombre d'instantanes du bloc S(Q,w) :
        // N+1 mesures spectrales moyennees, l'Avro est ecrit a la (N+1)-ieme ; le bloc
        // planifie (N+1) * measurementRateSQW balayages, d'ou le garde-fou de validation
        // (>= 2 : une moyenne d'un seul instantane n'a pas de sens ; <= 1 000 000 :
        // overflow int impossible). ----
        Integer nSQW = null;

        // ---- --seed-random (defaut : false = graine fixe, reproductible). true : la graine
        // maitresse est tiree au hasard a chaque lancement (apres parse, avant printOptions,
        // donc la graine effective est affichee ; la redonner en argument 10 rend le run
        // reproductible). Precedence sur la graine positionnelle. ----
        boolean seedRandom = false;

        // ---- --spins-per-T (defaut : false). true : un instantane des spins (avro) est ecrit
        // a chaque temperature du recuit, apres la production, nomme avec la temperature. ----
        boolean spinsPerT = false;

        // ---- --qz (defaut : 0.0) composante z du vecteur de diffusion, en r.l.u. Defaut 0 =
        // honeycomb 2D. Non finie (NaN/Infini) : rejet a la validation. Enregistree comme
        // propriete "Qzz" dans les trois Avro de structure (traceabilite ; inertie physique
        // hors DynamicStructureFactor, dont le deltaZ = z_B − z_A vient desormais de
        // l'input — a qzz = 0 le terme de phase z s'annule). ----
        double qz = 0.0;

        // ---- (2) --field / --no-field (defaut : FALSE) : aucune boucle en champ sans --field
        // (une seule replique a h=0) ; --field active le balayage hStart/hEnd/hStep.
        boolean runUnderField = false;

        // ---- (3) --annealing / --no-annealing (defaut : true) : execute le recuit
        // (thermalisation + production) ; --no-annealing exige --ac (voir parse()).
        boolean runAnnealing = true;

        // ---- (4) --adaptive-cooling / --no-adaptive-cooling (defaut : false, geometrique
        // 0.95) ; --cooling reste un alias de --adaptive-cooling.
        boolean adaptiveCooling = false;

        // ---- (5) --freeze-sigma / --no-freeze-sigma (defaut : true) ; --legacy-sigma reste un
        // alias de --no-freeze-sigma. Affecte opt.freezeSigmaInProduction.
        boolean freezeSigma = true;

        // ---- (6) --overrelax N / --no-overrelax (equivalent a --overrelax 0). N>=1 =>
        // useOverRelaxation=true, N=0 => false. Defaut : true avec N=2. Valeur effective
        // (opt.overRelaxationPerSweep) = useOverRelaxation ? overrelaxN : 0.
        /**
         * Deux sur-relaxations micro-canoniques par balayage Metropolis, par defaut.
         *
         * <p>Mesure sur le recuit complet 30x30 : tau passe de ~100-200 sweeps a ~1 sweep aux
         * basses temperatures, et le recuit entier tombe de 10,9 M a 0,31 M balayages Metropolis
         * — soit un facteur 35 — avec un &lt;E&gt;(T) identique a 1,3 sigma pres sur toute la
         * gamme. La sur-relaxation est micro-canonique (elle conserve exactement l'energie) donc
         * elle ne peut pas biaiser l'ensemble ; elle ne fait que decorreler plus vite. Il n'y a
         * aucune raison de ne pas l'activer par defaut. {@code --no-overrelax} (ou
         * {@code --overrelax 0}) la desactive.</p>
         */
        int overrelaxN = 2;
        boolean useOverRelaxation = true;

        boolean quiet = false;
        double[] fieldDir = { -1.d, 1.d, 0.d };
        Integer minTherm, maxTherm;
        Double targetNeff;
        Double gewekeSkip;

        // ---- (7) --precession-at-end / --no-precession-at-end (defaut : true) -- expose
        // opt.precessionAtEnd.
        boolean precessionAtEnd = true;

        // ---- Diagnostic de decorrelation et arret par stagnation (voir --preset, ci-dessous).
        /** Nom du preset effectif, pour {@link #printOptions}. "standard" = comportement actuel. */
        String preset = "standard";
        /**
         * (8) --stagnation / --no-stagnation (defaut : false). Positionne directement par cette
         * option, ou par {@link #applyPreset} (screening) : le parseur est lineaire, donc la
         * derniere des deux lue sur la ligne de commande l'emporte (voir USAGE).
         */
        boolean stagnationDetection = false;
        /** (9) --stop-therm-on-stagnation / --no-stop-therm-on-stagnation (defaut : false). */
        boolean stopThermalizationOnStagnation = false;
        /** (10) --stop-prod-on-saturation / --no-stop-prod-on-saturation (defaut : false). */
        boolean stopProductionOnSaturation = false;
        Integer maxProductionSweeps;
        Double minNeff;
        Double tauMultiplier;
        /** {@code --decorrelation-log}. */
        boolean decorrelationLog = false;
        Double memoryMultiplier;

        // ---- (11) --ac / --no-ac (defaut : false) -- voir le bloc dedie de USAGE.
        /** Calcule la susceptibilite AC (chi', chi'') en plus (ou, avec --no-annealing, a la place) du recuit. */
        boolean computeAC = false;
        /** Repliques AC independantes par (fichier, champ statique) ; {@code --acReplicas}. */
        int acReplicas = 8;
        /** Amplitude du champ oscillant ({@code ACSusceptibility.Config#h0}) ; {@code --acH0}. */
        double acH0 = 0.02;
        /** Periodes d'excitation, en sweeps ({@code ACSusceptibility.Config#periodsSweeps}) ; {@code --acPeriods}. */
        int[] acPeriods = { 64, 128, 256 };
        /** Periodes de regime transitoire ecartees avant la demodulation ; {@code --acTransientPeriods}. */
        int acTransientPeriods = 5;
        /** Periodes de mesure accumulees dans la demodulation ; {@code --acMeasurePeriods}. */
        int acMeasurePeriods = 30;
        /** Temperature de depart (la plus haute) de chaque echelle de refroidissement AC ; {@code --acTempMax}. */
        double acTempMax = 2.0;
        /** Temperature d'arret (la plus basse) de la grille AC ; {@code --acTempMin}. */
        double acTempMin = 0.05;
        /** Pas de la grille lineaire de temperature AC ; {@code --acTempStep}. */
        double acTempStep = 0.05;
        /** Direction du champ oscillant ; {@code null} => retombe sur {@link #fieldDir}. {@code --acFieldDir}. */
        double[] acFieldDir = null;
    }

    /**
     * Applique un preset a {@code a} : "standard" ne modifie rien (comportement historique,
     * y compris si un preset different avait ete choisi plus tot sur la meme ligne de commande
     * -- il remet donc explicitement a leurs defauts les seuls champs qu'un preset peut avoir
     * touches). "screening" et "strict" sont decrits dans {@link #USAGE}.
     *
     * <p>Precedence avec les options explicites ({@code --targetNeff}, {@code --maxTherm}, ...) :
     * {@link #parse} applique chaque option dans l'ordre ou elle apparait sur la ligne de
     * commande, en mutant directement les memes champs de {@code Args}. Une option explicite
     * placee APRES {@code --preset} ecrase donc la valeur du preset ; placee AVANT, c'est le
     * preset qui l'ecrase. Aucune logique de precedence dediee n'est necessaire : c'est le seul
     * comportement qu'un parseur lineaire peut avoir, et c'est celui documente dans
     * {@link #USAGE}.</p>
     *
     * @throws IllegalArgumentException si {@code name} n'est ni "screening", ni "standard", ni
     *         "strict".
     */
    private static void applyPreset(Args a, String name) {
        switch (name) {
            case "standard" -> {
                a.stagnationDetection = false;
                a.stopThermalizationOnStagnation = false;
                a.stopProductionOnSaturation = false;
                a.maxProductionSweeps = null;
                a.minNeff = null;
                a.tauMultiplier = null;
            }
            case "screening" -> {
                // Criblage rapide : le detecteur de stagnation est la strategie principale de
                // limitation du cout, targetNeff/maxTherm ne sont plus qu'un filet de securite.
                a.targetNeff = 100.0;
                a.maxTherm = 50_000;
                a.maxProductionSweeps = 50_000;
                a.stagnationDetection = true;
                a.stopThermalizationOnStagnation = true;
                a.stopProductionOnSaturation = true;
            }
            case "strict" -> {
                // Chaine longue et decorrelee : pas de detecteur de stagnation, le budget est
                // suppose suffisant pour atteindre la cible sans avoir besoin de le constater.
                a.targetNeff = 1000.0;
                a.minNeff = 200.0;
                a.tauMultiplier = 50.0;
                a.maxProductionSweeps = 2_000_000;
                a.maxTherm = 500_000;
                a.stagnationDetection = false;
                a.stopThermalizationOnStagnation = false;
                a.stopProductionOnSaturation = false;
            }
            default -> throw new IllegalArgumentException(
                    "--preset attend screening, standard ou strict, recu " + name);
        }
    }

    /**
     * Point d'entree du processus : delegue a {@link #run(String[])} et ne sort en erreur que
     * si celle-ci renvoie un code non nul.
     */
    public static void main(String[] args) {
        final int code = run(args);
        if (code != EXIT_OK) System.exit(code);
    }

    /**
     * Execute le lot. Ne leve rien : les echecs de repliques sont collectes, resumes sur
     * stdout, et traduits par {@link #EXIT_FAILURE}.
     *
     * @return {@link #EXIT_OK}, {@link #EXIT_FAILURE} ou {@link #EXIT_USAGE}.
     */
    public static int run(String[] args) {
        if (args == null || args.length == 0) {
            System.out.println(USAGE);
            return EXIT_USAGE;
        }
        if ("-h".equals(args[0]) || "--help".equals(args[0])) {
            System.out.println(USAGE);
            return EXIT_OK;
        }
        if ("--bench".equals(args[0])) {
            return bench(Arrays.copyOfRange(args, 1, args.length));
        }

        final Args a;
        try {
            a = parse(args);
        } catch (IllegalArgumentException e) {
            System.out.println("Argument invalide : " + e.getMessage());
            System.out.println();
            System.out.println(USAGE);
            return EXIT_USAGE;
        }
        if (a.seedRandom) {
            // Graine maitresse tiree au hasard UNE fois par lancement, AVANT printOptions :
            // la graine effective est affichee, et la redonner en argument positionnel 10
            // (ou la remettre dans le script) rend le lancement reproductible.
            a.seed = new java.security.SecureRandom().nextLong();
        }
        return execute(a);
    }

    private static Args parse(String[] argv) {
        final Args a = new Args();
        final List<String> pos = new ArrayList<>();
        for (int i = 0; i < argv.length; i++) {
            final String s = argv[i];
            switch (s) {
                case "--sqw"         -> a.computeSQW = true;
                case "--no-sqw"      -> a.computeSQW = false;
                case "--field"       -> a.runUnderField = true;
                case "--no-field"    -> a.runUnderField = false;
                case "--annealing"   -> a.runAnnealing = true;
                case "--no-annealing"-> a.runAnnealing = false;
                case "--adaptive-cooling", "--cooling" -> a.adaptiveCooling = true;
                case "--no-adaptive-cooling" -> a.adaptiveCooling = false;
                case "--freeze-sigma"-> a.freezeSigma = true;
                case "--no-freeze-sigma", "--legacy-sigma" -> a.freezeSigma = false;
                case "--precession-at-end"   -> a.precessionAtEnd = true;
                case "--no-precession-at-end"-> a.precessionAtEnd = false;
                case "--stagnation"      -> a.stagnationDetection = true;
                case "--no-stagnation"   -> a.stagnationDetection = false;
                case "--stop-therm-on-stagnation"    -> a.stopThermalizationOnStagnation = true;
                case "--no-stop-therm-on-stagnation" -> a.stopThermalizationOnStagnation = false;
                case "--stop-prod-on-saturation"      -> a.stopProductionOnSaturation = true;
                case "--no-stop-prod-on-saturation"   -> a.stopProductionOnSaturation = false;
                case "--quiet"       -> a.quiet = true;
                case "--overrelax"   -> {
                    final int n = Integer.parseInt(need(argv, ++i, s));
                    a.overrelaxN = n;
                    a.useOverRelaxation = n != 0;
                }
                case "--no-overrelax" -> a.useOverRelaxation = false;
                case "--minTherm"    -> a.minTherm = Integer.parseInt(need(argv, ++i, s));
                case "--maxTherm"    -> a.maxTherm = Integer.parseInt(need(argv, ++i, s));
                case "--wmax"        -> a.wMax = Double.parseDouble(need(argv, ++i, s));
                case "--qz"          -> a.qz = Double.parseDouble(need(argv, ++i, s));
                case "--initTemp"    -> a.initTemp = Double.parseDouble(need(argv, ++i, s));
                case "--endTemp"     -> a.endTemp = Double.parseDouble(need(argv, ++i, s));
                case "--window"      -> a.window = Double.parseDouble(need(argv, ++i, s));
                case "--domega"      -> a.domega = Double.parseDouble(need(argv, ++i, s));
                case "--nSQW"        -> a.nSQW = Integer.parseInt(need(argv, ++i, s));
                case "--seed-random" -> a.seedRandom = true;
                case "--spins-per-T" -> a.spinsPerT = true;
                case "--targetNeff"  -> a.targetNeff = Double.parseDouble(need(argv, ++i, s));
                case "--gewekeSkip"  -> a.gewekeSkip = Double.parseDouble(need(argv, ++i, s));
                case "--preset"      -> {
                    final String p = need(argv, ++i, s);
                    applyPreset(a, p);
                    a.preset = p;
                }
                case "--decorrelation-log" -> a.decorrelationLog = true;
                case "--memoryMultiplier"  -> a.memoryMultiplier = Double.parseDouble(need(argv, ++i, s));
                case "--field-dir"   -> {
                    String[] xyz = need(argv, ++i, s).split(",");
                    if (xyz.length != 3) throw new IllegalArgumentException("--field-dir attend x,y,z");
                    a.fieldDir = new double[] { Double.parseDouble(xyz[0]),
                                                Double.parseDouble(xyz[1]),
                                                Double.parseDouble(xyz[2]) };
                }
                case "--ac"          -> a.computeAC = true;
                case "--no-ac"       -> a.computeAC = false;
                case "--acReplicas"  -> a.acReplicas = Integer.parseInt(need(argv, ++i, s));
                case "--acH0"        -> a.acH0 = Double.parseDouble(need(argv, ++i, s));
                case "--acPeriods"   -> {
                    final String[] parts = need(argv, ++i, s).split(",");
                    final int[] periods = new int[parts.length];
                    for (int j = 0; j < parts.length; j++) periods[j] = Integer.parseInt(parts[j].trim());
                    a.acPeriods = periods;
                }
                case "--acTransientPeriods" -> a.acTransientPeriods = Integer.parseInt(need(argv, ++i, s));
                case "--acMeasurePeriods"   -> a.acMeasurePeriods = Integer.parseInt(need(argv, ++i, s));
                case "--acTempMax"   -> a.acTempMax = Double.parseDouble(need(argv, ++i, s));
                case "--acTempMin"   -> a.acTempMin = Double.parseDouble(need(argv, ++i, s));
                case "--acTempStep"  -> a.acTempStep = Double.parseDouble(need(argv, ++i, s));
                case "--acFieldDir"  -> {
                    String[] xyz = need(argv, ++i, s).split(",");
                    if (xyz.length != 3) throw new IllegalArgumentException("--acFieldDir attend x,y,z");
                    a.acFieldDir = new double[] { Double.parseDouble(xyz[0]),
                                                  Double.parseDouble(xyz[1]),
                                                  Double.parseDouble(xyz[2]) };
                }
                default -> {
                    if (s.startsWith("--")) throw new IllegalArgumentException("option inconnue " + s);
                    pos.add(s);
                }
            }
        }
        if (pos.size() < 2) throw new IllegalArgumentException("inputDir et outputDir sont obligatoires");
        a.inputDir = pos.get(0);
        a.outputDir = pos.get(1);
        if (pos.size() > 2) a.totalThreads = Integer.parseInt(pos.get(2));
        if (pos.size() > 3) a.threadsPerReplica = Integer.parseInt(pos.get(3));
        if (pos.size() > 4) a.La = Integer.parseInt(pos.get(4));
        if (pos.size() > 5) a.Lb = Integer.parseInt(pos.get(5));
        if (pos.size() > 6) a.hStart = Double.parseDouble(pos.get(6));
        if (pos.size() > 7) a.hEnd = Double.parseDouble(pos.get(7));
        if (pos.size() > 8) a.hStep = Double.parseDouble(pos.get(8));
        if (pos.size() > 9) a.seed = Long.parseLong(pos.get(9));
        if (a.totalThreads < 1) throw new IllegalArgumentException("totalThreads >= 1");
        if (a.threadsPerReplica < 1) throw new IllegalArgumentException("threadsPerReplica >= 1");
        if (a.La < 1 || a.Lb < 1) throw new IllegalArgumentException("La, Lb >= 1");
        if (a.wMax < 0) throw new IllegalArgumentException("--wmax negatif interdit (0 = mode statique)");
        if (!Double.isFinite(a.qz)) {
            throw new IllegalArgumentException("--qz doit etre fini (NaN/Infini rejetes) : " + a.qz);
        }
        // Directions de champ validees INCONDITIONNELLEMENT (comme --qz ci-dessus) : normalize()
        // n'a pas de garde, un vecteur nul ou non fini donnerait des NaN silencieux (replica(),
        // recuit) ; ici l'erreur remonte des le parsing en EXIT_USAGE.
        if (!finiteNonZero(a.fieldDir)) {
            throw new IllegalArgumentException(
                    "--field-dir doit etre un vecteur non nul a composantes finies");
        }
        if (a.acFieldDir != null && !finiteNonZero(a.acFieldDir)) {
            throw new IllegalArgumentException(
                    "--acFieldDir doit etre un vecteur non nul a composantes finies");
        }
        if (!Double.isNaN(a.initTemp) && a.initTemp <= 0) {
            throw new IllegalArgumentException("--initTemp doit etre > 0");
        }
        if (!Double.isNaN(a.endTemp) && a.endTemp <= 0) {
            throw new IllegalArgumentException("--endTemp doit etre > 0");
        }
        if (!Double.isNaN(a.window) && a.window <= 0) {
            throw new IllegalArgumentException("--window doit etre > 0");
        }
        // --domega : deux ecritures d'un meme reglage (W = 2 pi / domega) — la concomitance
        // serait ambigue (laquelle gagne ?), elle est rejetee des le parsing (EXIT_USAGE).
        if (!Double.isNaN(a.domega) && !Double.isNaN(a.window)) {
            throw new IllegalArgumentException(
                    "--domega et --window sont deux ecritures du meme reglage (W = 2*pi/domega)"
                    + " : n'en fournir qu'un");
        }
        if (!Double.isNaN(a.domega) && (!Double.isFinite(a.domega) || a.domega <= 0)) {
            throw new IllegalArgumentException(
                    "--domega doit etre fini et > 0 (resolution visee en energie, hbarre = 1) : "
                            + a.domega);
        }
        // Garde-fou de bande : la bande [0, wmax] doit garder au moins 2 points (pas ~ domega).
        if (!Double.isNaN(a.domega) && !Double.isNaN(a.wMax) && a.wMax > 0
                && a.domega >= a.wMax) {
            throw new IllegalArgumentException("la resolution domega doit etre < wmax"
                    + " (au moins 2 points dans la bande) : domega = " + a.domega
                    + ", wmax = " + a.wMax);
        }
        // --nSQW : nombre d'instantanes du bloc S(Q,w). >= 2 : une moyenne spectrale d'un
        // seul instantane n'a pas de sens ; <= 1 000 000 : garde-fou d'overflow (sqwBlock
        // calcule sweeps = (nSQW + 1) * measurementRateSQW, en int).
        if (a.nSQW != null && (a.nSQW < 2 || a.nSQW > 1_000_000)) {
            throw new IllegalArgumentException(
                    "--nSQW doit etre compris entre 2 et 1 000 000 (instantanes du bloc S(Q,w)"
                    + " ; >= 2 : une moyenne d'un seul instantane n'a pas de sens) : " + a.nSQW);
        }
        if (!Double.isNaN(a.initTemp) && !Double.isNaN(a.endTemp) && a.initTemp <= a.endTemp) {
            throw new IllegalArgumentException("--initTemp (" + a.initTemp
                    + ") doit etre > --endTemp (" + a.endTemp + ")");
        }
        if (!Double.isNaN(a.wMax) && !a.computeSQW) {
            // --wmax ne sert qu'a sqwBlock ; avec --no-sqw il serait silencieusement ignore.
            System.err.println("AVERTISSEMENT : --wmax sans objet avec --no-sqw (S(Q,w) desactive) ; argument ignore.");
        }
        if (a.overrelaxN < 0) throw new IllegalArgumentException("--overrelax >= 0");
        // (2) runUnderField : hStart/hEnd/hStep et leurs validations specifiques ne s'appliquent
        // que si la boucle en champ est active (--field) ; sinon ils sont ignores (nH force a 1,
        // h = 0.0, voir fieldCount).
        if (a.runUnderField) {
            if (!(a.hStep > 0.0)) throw new IllegalArgumentException("hStep > 0");
            if (fieldCount(a) > MAX_FIELDS) {
                throw new IllegalArgumentException("(hEnd - hStart) / hStep = " + fieldCount(a)
                        + " champs, maximum " + MAX_FIELDS + " ; verifier hStep");
            }
        }
        // (3) runAnnealing : desactiver le recuit sans mesure alternative ne produirait aucune
        // observable.
        if (!a.runAnnealing && !a.computeAC) {
            throw new IllegalArgumentException(
                    "--no-annealing necessite --ac (susceptibilite AC)");
        }
        // (11) --ac : bornes du protocole de susceptibilite AC, validees inconditionnellement
        // (meme esprit que --overrelax >= 0 ci-dessus) : une valeur invalide est une erreur
        // d'usage, que --ac soit actif ou non sur cette ligne de commande.
        if (a.acReplicas < 1) throw new IllegalArgumentException("--acReplicas >= 1");
        if (!(a.acH0 > 0.0)) throw new IllegalArgumentException("--acH0 doit etre > 0");
        if (a.acPeriods == null || a.acPeriods.length == 0) {
            throw new IllegalArgumentException("--acPeriods ne doit pas etre vide");
        }
        final Set<Integer> seenAcPeriods = new LinkedHashSet<>();
        for (final int p : a.acPeriods) {
            if (p < 2) {
                throw new IllegalArgumentException("--acPeriods : chaque periode doit etre >= 2, recu " + p);
            }
            if (!seenAcPeriods.add(p)) {
                throw new IllegalArgumentException("--acPeriods : periode dupliquee " + p);
            }
        }
        if (a.acTransientPeriods < 0) throw new IllegalArgumentException("--acTransientPeriods >= 0");
        if (a.acMeasurePeriods < 1) throw new IllegalArgumentException("--acMeasurePeriods >= 1");
        if (!(a.acTempMin > 0.0)) throw new IllegalArgumentException("--acTempMin doit etre > 0");
        if (!(a.acTempMax > a.acTempMin)) {
            throw new IllegalArgumentException("--acTempMax (" + a.acTempMax
                    + ") doit etre > --acTempMin (" + a.acTempMin + ")");
        }
        if (!(a.acTempStep > 0.0)) throw new IllegalArgumentException("--acTempStep doit etre > 0");
        // (9)/(10) : aucun arret anticipe sur stagnation sans detecteur actif. Le client a exclu
        // toute activation implicite : ni --stagnation n'est arme automatiquement ici, ni les
        // arrets ne sont desactives silencieusement -- l'usage est simplement rejete.
        if ((a.stopThermalizationOnStagnation || a.stopProductionOnSaturation) && !a.stagnationDetection) {
            throw new IllegalArgumentException(
                    "--stop-therm-on-stagnation / --stop-prod-on-saturation necessitent un "
                    + "detecteur de stagnation actif : ajouter --stagnation (ou --preset screening)");
        }
        // --domega : ecriture alternative de --window, convertie APRES les validations et
        // AVANT toute construction de ParallelOptions : le reste de la chaine (Delta =
        // pi/w_max, ndt = round(W/Delta) dans sqwBlock) ne voit que la fenetre W = 2 pi/domega.
        if (!Double.isNaN(a.domega)) {
            a.window = fenetreDepuisDeltaOmega(a.domega);
        }
        return a;
    }

    /**
     * Fenetre d'integration W equivalente a une resolution visee {@code domega} en energie :
     * {@code W = 2 pi / domega} (unites hbarre = 1 ; domega et W dans les memes unites que
     * les couplages J de l'input). Package-private et statique pour etre verifiee directement
     * par DomegaTest.
     */
    static double fenetreDepuisDeltaOmega(double domega) {
        return Math.TAU / domega;
    }

    /**
     * Vecteur a composantes finies et de norme^2 &gt; 0 : pre-requis a {@code MathOps.normalize}
     * (sans garde), qui sinon produit des NaN silencieux.
     */
    private static boolean finiteNonZero(double[] v) {
        return Double.isFinite(v[0]) && Double.isFinite(v[1]) && Double.isFinite(v[2])
                && v[0] * v[0] + v[1] * v[1] + v[2] * v[2] > 0.0;
    }

    /**
     * Nombre de valeurs de champ visitees, {@code h_k = hStart + k hStep} pour
     * {@code k = 0 .. nH-1}.
     *
     * <p>Un index entier, et non {@code for (double h = hStart; h < hEnd; h += hStep)} : cette
     * boucle-la accumule l'erreur d'arrondi de {@code hStep} a chaque tour, si bien que le
     * dernier pas tombe soit juste avant, soit juste apres {@code hEnd} selon la valeur
     * binaire exacte de {@code hStep} — le nombre de repliques n'etait donc pas previsible.
     * Avec les valeurs par defaut (0 -&gt; 1 par 0,0125), {@code (hEnd - hStart) / hStep} vaut
     * 79,999... en binaire ; le {@code -1e-9} avant le {@code ceil} absorbe cet epsilon et
     * donne exactement 80 champs, {@code h_79 = 0,9875}.</p>
     *
     * <p>Si {@code a.runUnderField} est faux (defaut, voir {@code --field}), la boucle en champ
     * est desactivee : une seule replique est lancee, a {@code h = 0}, et {@code hStart}/
     * {@code hEnd}/{@code hStep} ne sont pas consultes.</p>
     */
    private static int fieldCount(Args a) {
        if (!a.runUnderField) return 1;
        final double n = (a.hEnd - a.hStart) / a.hStep - 1e-9;
        if (!(n > 0.0)) return 0;
        final double c = Math.ceil(n);
        return c > MAX_FIELDS ? MAX_FIELDS + 1 : (int) c;
    }

    private static String need(String[] argv, int i, String flag) {
        if (i >= argv.length) throw new IllegalArgumentException(flag + " attend une valeur");
        return argv[i];
    }

    // ------------------------------------------------------------------ execution

    /** Une replique qui a echoue : de quoi la relancer sans relire tout le journal. */
    private record Failure(String file, double field, Throwable error) {}

    /** Une replique de recuit a lancer. */
    private record Task(File file, double field) {}

    /** Une replique AC a lancer : fichier, champ statique, indice de replique. */
    private record ACTask(File file, double field, int replicaIdx) {}

    /** Cle de regroupement des repliques AC : un CSV par (fichier, champ statique). */
    private record ACKey(String fileName, double field) {}

    private static int execute(Args a) {
        final long startTime = System.currentTimeMillis();

        final File directory = new File(a.inputDir);
        if (!directory.isDirectory()) {
            System.out.println("Le dossier specifie n'existe pas ou n'est pas un dossier : " + a.inputDir);
            return EXIT_USAGE;
        }
        final File[] files = directory.listFiles();
        if (files == null) {
            System.out.println("Impossible de lister " + a.inputDir);
            return EXIT_USAGE;
        }
        Arrays.sort(files); // ordre de soumission reproductible

        // Le repertoire de sortie est verifie AVANT de lancer quoi que ce soit : sans cela, une
        // faute de frappe ou un montage en lecture seule ne se manifeste qu'a la fin de la
        // premiere replique, apres des heures de calcul jete.
        final File outDir = new File(a.outputDir);
        if (!outDir.exists() && !outDir.mkdirs() && !outDir.isDirectory()) {
            System.out.println("Impossible de creer le repertoire de sortie : " + a.outputDir);
            return EXIT_USAGE;
        }
        if (!outDir.isDirectory()) {
            System.out.println("Le repertoire de sortie n'est pas un dossier : " + a.outputDir);
            return EXIT_USAGE;
        }
        if (!outDir.canWrite()) {
            System.out.println("Le repertoire de sortie n'est pas accessible en ecriture : " + a.outputDir);
            return EXIT_USAGE;
        }

        final int cores = Runtime.getRuntime().availableProcessors();
        final int concurrentReplicas = Math.max(1, a.totalThreads / a.threadsPerReplica);
        if (a.totalThreads > cores) {
            System.out.printf(Locale.US,
                    "AVERTISSEMENT : totalThreads = %d > %d coeurs disponibles ; les threads de "
                    + "CheckerboardMetropolis attendent activement, la sur-souscription degrade "
                    + "fortement le debit.%n", a.totalThreads, cores);
        }
        if (a.threadsPerReplica * concurrentReplicas > a.totalThreads) {
            System.out.printf(Locale.US,
                    "AVERTISSEMENT : threadsPerReplica (%d) x concurrentReplicas (%d) = %d > "
                    + "totalThreads (%d) ; reduire threadsPerReplica.%n",
                    a.threadsPerReplica, concurrentReplicas,
                    a.threadsPerReplica * concurrentReplicas, a.totalThreads);
        }

        final int nH = fieldCount(a);
        final List<Task> tasks = new ArrayList<>();
        final List<ACTask> acTasks = new ArrayList<>();
        // Ordre d'insertion preserve (LinkedHashSet) : l'ordre d'ecriture des CSV AC, et de tout
        // message qui les enumere, reste celui de soumission -- reproductible d'une execution a
        // l'autre malgre le parallelisme des repliques elles-memes.
        final Set<ACKey> acKeys = new LinkedHashSet<>();
        int nFiles = 0;
        for (File file : files) {
            // isFile() : erreur d'usage claire avant la construction (Input est fail-fast), sinon
            // un sous-repertoire nomme bcaoExplor_* finirait en echec de replique.
            if (!file.isFile() || !file.getName().startsWith("bcaoExplor")) continue;
            nFiles++;
            for (int k = 0; k < nH; k++) {
                final double h = a.runUnderField ? a.hStart + k * a.hStep : 0.0;
                // (3) runAnnealing : pas de tache de recuit du tout si --no-annealing (au lieu du
                // no-op de replica() historique) -- --no-annealing exige --ac (voir parse()), donc
                // ce cas ne laisse jamais tasks et acTasks vides simultanement.
                if (a.runAnnealing) tasks.add(new Task(file, h));
                if (a.computeAC) {
                    acKeys.add(new ACKey(file.getName(), h));
                    for (int r = 0; r < a.acReplicas; r++) acTasks.add(new ACTask(file, h, r));
                }
            }
        }

        printOptions(a, cores, concurrentReplicas, nFiles, nH, tasks.size());

        if (tasks.isEmpty() && acTasks.isEmpty()) {
            System.out.println("Aucune replique a lancer (fichiers bcaoExplor* : " + nFiles
                    + ", champs : " + nH + ").");
            return EXIT_USAGE;
        }

        final AtomicBoolean coloringPrinted = new AtomicBoolean(false);
        final ExecutorService executor = Executors.newFixedThreadPool(concurrentReplicas);
        final List<Future<?>> futures = new ArrayList<>(tasks.size());
        for (Task t : tasks) {
            // submit(), pas execute() : execute() confie l'exception au UncaughtExceptionHandler
            // du pool, qui l'imprime au mieux et la perd au pire. Le Future la conserve.
            futures.add(executor.submit(() -> replica(t.file(), t.field(), a, coloringPrinted)));
        }
        // Memes pool et discipline memoire que les taches de recuit ci-dessus : chaque tache AC
        // construit son propre reseau DANS le Callable (voir acReplica), la memoire simultanement
        // engagee reste donc bornee par concurrentReplicas, jamais par acTasks.size().
        final List<Future<List<ACSusceptibility.Sample>>> acFutures = new ArrayList<>(acTasks.size());
        for (ACTask t : acTasks) {
            acFutures.add(executor.submit(
                    () -> acReplica(t.file(), t.field(), t.replicaIdx(), a, coloringPrinted)));
        }
        executor.shutdown();
        try {
            executor.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            executor.shutdownNow();
            System.out.println("Interrompu avant la fin du lot.");
            return EXIT_FAILURE;
        }

        final List<Failure> failures = new ArrayList<>();
        for (int i = 0; i < futures.size(); i++) {
            final Task t = tasks.get(i);
            try {
                futures.get(i).get();
            } catch (ExecutionException e) {
                failures.add(new Failure(t.file().getName(), t.field(),
                        e.getCause() == null ? e : e.getCause()));
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                failures.add(new Failure(t.file().getName(), t.field(), e));
            } catch (RuntimeException e) {
                failures.add(new Failure(t.file().getName(), t.field(), e));
            }
        }

        // ---- Susceptibilite AC : meme mecanique de collecte d'echecs que le recuit ci-dessus
        // (les echecs rejoignent la meme liste "failures"), mais les succes sont regroupes par
        // (fichier, champ) : un CSV n'est ecrit que si TOUTES les repliques AC de ce couple ont
        // reussi.
        final Map<ACKey, List<List<ACSusceptibility.Sample>>> acSuccess = new LinkedHashMap<>();
        final Set<ACKey> acFailedKeys = new LinkedHashSet<>();
        for (int i = 0; i < acFutures.size(); i++) {
            final ACTask t = acTasks.get(i);
            final ACKey key = new ACKey(t.file().getName(), t.field());
            try {
                final List<ACSusceptibility.Sample> samples = acFutures.get(i).get();
                acSuccess.computeIfAbsent(key, k -> new ArrayList<>()).add(samples);
            } catch (ExecutionException e) {
                failures.add(new Failure(t.file().getName(), t.field(),
                        e.getCause() == null ? e : e.getCause()));
                acFailedKeys.add(key);
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                failures.add(new Failure(t.file().getName(), t.field(), e));
                acFailedKeys.add(key);
            } catch (RuntimeException e) {
                failures.add(new Failure(t.file().getName(), t.field(), e));
                acFailedKeys.add(key);
            }
        }

        // Agregation + ecriture CSV sur le THREAD PRINCIPAL (execute() n'est jamais appelee que
        // depuis lui) : un (fichier, champ) dont au moins une replique AC a echoue ne produit pas
        // de CSV, mais n'empeche pas les autres (fichier, champ) d'en produire un.
        int acCsvWritten = 0;
        for (ACKey key : acKeys) {
            if (acFailedKeys.contains(key)) continue;
            final List<List<ACSusceptibility.Sample>> perReplica = acSuccess.get(key);
            final List<ACPoint> points = ACSusceptibility.aggregate(perReplica, a.acH0, key.field());
            final Path csvPath = ACSusceptibilityLog.defaultPath(a.outputDir, key.fileName(), key.field());
            try (ACSusceptibilityLog log = new ACSusceptibilityLog(csvPath)) {
                for (ACPoint p : points) log.write(p);
            }
            acCsvWritten++;
        }

        final long seconds = (System.currentTimeMillis() - startTime) / 1000;
        final int totalTasks = tasks.size() + acTasks.size();
        if (failures.isEmpty()) {
            System.out.println("good job ! " + totalTasks + " repliques terminees"
                    + (acTasks.isEmpty() ? ""
                       : " (" + tasks.size() + " recuit + " + acTasks.size() + " AC)") + ".");
            if (acCsvWritten > 0) {
                System.out.println(acCsvWritten + " fichier(s) ac_susceptibility_*.csv ecrit(s).");
            }
            System.out.println("Tache executee en sec : " + seconds);
            return EXIT_OK;
        }

        System.out.println();
        System.out.printf(Locale.US, "%d replique(s) sur %d ont echoue :%n",
                failures.size(), totalTasks);
        for (Failure f : failures) {
            System.out.printf(Locale.US, "  %-40s h=%.6f : %s%n",
                    f.file(), f.field(), describe(f.error()));
        }
        // La pile de la premiere seulement : mille repliques en echec pour la meme raison
        // rendraient le journal illisible, et c'est la cause qui interesse.
        System.out.println();
        System.out.println("Pile de la premiere erreur :");
        failures.get(0).error().printStackTrace(System.out);
        if (acCsvWritten > 0) {
            System.out.println(acCsvWritten + " fichier(s) ac_susceptibility_*.csv ecrit(s) malgre tout.");
        }
        System.out.println("Tache executee en sec : " + seconds);
        return EXIT_FAILURE;
    }

    private static String describe(Throwable t) {
        final String m = t.getMessage();
        return t.getClass().getSimpleName() + (m == null ? "" : " : " + m);
    }

    /** Recapitulatif des options effectives, avant de lancer quoi que ce soit. */
    private static void printOptions(Args a, int cores, int concurrentReplicas,
                                     int nFiles, int nH, int nTasks) {
        System.out.println("---- AppParallel ----");
        System.out.println("  entree            : " + a.inputDir + " (" + nFiles + " fichier(s) bcaoExplor_*)");
        System.out.println("  sortie            : " + a.outputDir);
        System.out.printf(Locale.US, "  reseau            : %d x %d x 2 = %d sites%n",
                a.La, a.Lb, a.La * a.Lb * HoneycombGeometry.SUBLATTICES);
        // (2) runUnderField : RUPTURE de defaut (voir USAGE) -- sans --field, une seule replique
        // a h=0 est lancee et hStart/hEnd/hStep ne sont ni lus ni valides.
        if (a.runUnderField) {
            System.out.printf(Locale.US,
                    "  champ             : %d valeurs, h = %.6f .. %.6f par %.6f (--field actif)%n",
                    nH, a.hStart, nH > 0 ? a.hStart + (nH - 1) * a.hStep : a.hStart, a.hStep);
        } else {
            System.out.println("  champ             : desactive (1 replique a h=0) -- activer avec --field");
        }
        System.out.printf(Locale.US, "  direction du champ: (%.4f, %.4f, %.4f)%n",
                a.fieldDir[0], a.fieldDir[1], a.fieldDir[2]);
        System.out.printf(Locale.US,
                "  threads           : %d coeurs, totalThreads = %d, threadsPerReplica = %d, "
                + "repliques simultanees = %d%n",
                cores, a.totalThreads, a.threadsPerReplica, concurrentReplicas);
        System.out.println("  repliques         : " + nTasks);
        System.out.println("  graine maitresse  : " + a.seed + (a.seedRandom ? " (aleatoire ; copier cette valeur en argument 10 pour reproduire)" : ""));
        System.out.println("  annealing         : " + (a.runAnnealing
                ? "oui" : "non (mesure alternative : susceptibilite AC, --ac)"));
        final int effectiveOverrelax = a.useOverRelaxation ? a.overrelaxN : 0;
        System.out.println("  sur-relaxation    : " + effectiveOverrelax + " par balayage"
                + (a.useOverRelaxation ? "" : " (desactivee : --no-overrelax ou --overrelax 0)"));
        System.out.println("  relaxation fin T  : " + (a.precessionAtEnd ? "oui" : "non"));
        System.out.println("  refroidissement   : " + (a.adaptiveCooling ? "adaptatif (Huang)" : "geometrique 0.95"));
        System.out.println("  S(Q,w)            : " + (a.computeSQW ? "oui (derniere temperature)" : "non"));
        System.out.println("  w_max (unites Spin): " + (Double.isNaN(a.wMax) ? "historique" : a.wMax));
        if (!Double.isNaN(a.domega)) {
            System.out.println("  resolution visee domega : " + a.domega
                    + " (fenetre W = 2*pi/" + a.domega + " = "
                    + fenetreDepuisDeltaOmega(a.domega) + ")");
        }
        System.out.println("  nSQW              : "
                + (a.nSQW == null ? "40 (historique)" : a.nSQW));
        System.out.println("  qz (r.l.u.)       : " + a.qz + " (composante z du vecteur de diffusion)");
        System.out.println("  T initiale recuit : " + (Double.isNaN(a.initTemp) ? "historique (10)" : a.initTemp));
        System.out.println("  T finale recuit   : " + (Double.isNaN(a.endTemp) ? "historique (0.001)" : a.endTemp));
        System.out.println("  sigma en prod.    : " + (a.freezeSigma ? "gele" : "adaptatif (historique)"));
        if (a.minTherm != null)   System.out.println("  minTherm          : " + a.minTherm);
        if (a.maxTherm != null)   System.out.println("  maxTherm          : " + a.maxTherm);
        if (a.targetNeff != null) System.out.println("  targetNeff        : " + a.targetNeff);
        if (a.gewekeSkip != null) System.out.println("  gewekeSkip        : " + a.gewekeSkip);
        System.out.println("  preset            : " + a.preset);
        if (a.maxProductionSweeps != null)
            System.out.println("  maxProductionSweeps : " + a.maxProductionSweeps);
        if (a.minNeff != null)       System.out.println("  minNeff           : " + a.minNeff);
        if (a.tauMultiplier != null) System.out.println("  tauMultiplier     : " + a.tauMultiplier);
        System.out.println("  stagnation        : " + (a.stagnationDetection ? "active" : "inactive"));
        System.out.println("  arret therm/stagn.: " + (a.stopThermalizationOnStagnation ? "oui" : "non"));
        System.out.println("  arret prod/satur. : " + (a.stopProductionOnSaturation ? "oui" : "non"));
        System.out.println("  decorrelation-log : " + (a.decorrelationLog ? "oui" : "non"));
        if (a.memoryMultiplier != null)
            System.out.println("  memoryMultiplier  : " + a.memoryMultiplier);
        System.out.println("  verbeux           : " + (a.quiet ? "non" : "oui"));
        if (a.computeAC) {
            final double[] acDir = a.acFieldDir != null ? a.acFieldDir : a.fieldDir;
            final int nTAC = acTemperatureCount(a.acTempMax, a.acTempMin, a.acTempStep);
            // Sweeps de MESURE seuls (equilibration exclue, voir ACSusceptibility) : pour une
            // periode P, chaque temperature consomme (transient+measure)*P sweeps de mesure ;
            // sumPeriods = Somme_P (transient+measure)*P, cumulee sur toutes les periodes pour
            // UNE temperature. Le budget total (une replique, toutes periodes, toutes
            // temperatures) est donc nTAC * sumPeriods, puis * acReplicas pour le (fichier, champ).
            long sumPeriods = 0L;
            for (final int p : a.acPeriods) {
                sumPeriods += (long) (a.acTransientPeriods + a.acMeasurePeriods) * p;
            }
            final long measureBudget = (long) nTAC * sumPeriods * a.acReplicas;
            System.out.println("  susceptibilite AC : oui");
            System.out.println("    repliques AC    : " + a.acReplicas);
            System.out.printf(Locale.US, "    H0              : %.6g%n", a.acH0);
            System.out.println("    periodes (sweeps): " + Arrays.toString(a.acPeriods));
            System.out.printf(Locale.US, "    transitoire/mesure : %d / %d periodes%n",
                    a.acTransientPeriods, a.acMeasurePeriods);
            System.out.printf(Locale.US, "    grille T AC     : %.6f .. %.6f par %.6f (%d valeurs)%n",
                    a.acTempMax, a.acTempMin, a.acTempStep, nTAC);
            System.out.printf(Locale.US, "    direction AC    : (%.4f, %.4f, %.4f)%n",
                    acDir[0], acDir[1], acDir[2]);
            System.out.printf(Locale.US,
                    "    budget de mesure (par fichier x champ) : ~%d sweeps (+ equilibrations "
                    + "adaptatives)%n", measureBudget);
        } else {
            System.out.println("  susceptibilite AC : non");
        }
        System.out.println("---------------------");
    }

    /**
     * Nombre de temperatures de la grille AC lineaire descendante {@code tMax, tMax-tStep, ...}
     * (meme formule que {@code ACSusceptibility}, prive la-bas : dupliquee ici pour le seul besoin
     * du recapitulatif de {@link #printOptions}, {@code ACSusceptibility.run} refait le calcul en
     * interne et reste la seule source de verite sur la grille effectivement parcourue).
     */
    private static int acTemperatureCount(double tMax, double tMin, double tStep) {
        return (int) Math.floor((tMax - tMin) / tStep + 1e-9) + 1;
    }

    /** Construit une replique via {@link ReplicaFactory} et la simule. */
    private static void replica(File file, double field, Args a, AtomicBoolean coloringPrinted) {
        final String fileName = file.getName();
        if (!a.runAnnealing) {
            // Inatteignable : parse() exige --ac avec --no-annealing, et runAnnealing = false ne
            // cree aucune tache de recuit (execute()). Garde defensif.
            System.out.println(fileName + " h=" + field + " : --no-annealing, aucun recuit execute.");
            return;
        }
        System.out.println("task " + fileName + " started by " + Thread.currentThread().getName()
                + " field h: " + field);

        final ReplicaFactory.Built built = ReplicaFactory.build(new Input(file.getAbsolutePath()),
                                                               a.La, a.Lb);
        final Lattice lattice = built.lattice();
        final UnitCell uc = built.uc();

        // Champ dans la base de Kitaev, via le tenseur g.
        final double[] H_field = MathOps.scale(field,
                MathOps.matrixXvector(lattice.getGtensor(), MathOps.normalize(a.fieldDir)));
        lattice.setInteractionField(H_field);

        final LatticeColoring coloring = LatticeColoring.of(lattice);
        if (coloringPrinted.compareAndSet(false, true)) {
            System.out.println(coloring.summary());
        }

        final MonteCarlo mc = new MonteCarlo();
        final Spin spin = new Spin(lattice);

        // Graine deterministe par (fichier, champ) : deux repliques distinctes ne partagent
        // jamais la meme trajectoire, et une meme replique est reproductible.
        final long seed = a.seed ^ (fileName.hashCode() * 31L + Double.doubleToLongBits(field));

        final MonteCarlo.ParallelOptions opt = new MonteCarlo.ParallelOptions();
        opt.computeSQW = a.computeSQW;
        opt.overRelaxationPerSweep = a.useOverRelaxation ? a.overrelaxN : 0;
        opt.freezeSigmaInProduction = a.freezeSigma;
        opt.precessionAtEnd = a.precessionAtEnd;
        opt.verbose = !a.quiet;
        opt.wMax = a.wMax;
        opt.initTemp = a.initTemp;
        opt.endTemp = a.endTemp;
        opt.window = a.window;
        opt.qz = a.qz;
        // Meme graine que le sweeper : la configuration initiale devient reproductible elle
        // aussi, donc la replique entiere l'est.
        opt.initialConfigSeed = seed;
        if (!Double.isNaN(a.wMax) && a.wMax > 0) {
            final double fenetre = Double.isNaN(a.window)
                    ? MonteCarlo.SQW_INTEGRATION_WINDOW : a.window;
            final double deltaCible = Math.PI / a.wMax;
            if (deltaCible < spin.getDt()) {
                throw new IllegalArgumentException("--wmax " + a.wMax
                        + " depasse pi/dt (" + (Math.PI / spin.getDt())
                        + ") : l'echantillonnage est deja limite par le pas RK4");
            }
            final double delta = spin.getDt()
                    * Math.max(1, (int) Math.round(deltaCible / spin.getDt()));
            if (delta > fenetre / 2.0) {
                throw new IllegalArgumentException("--wmax " + a.wMax + " trop faible : moins de 2 "
                        + "points spectraux avec la fenetre " + fenetre
                        + " ; minimum ~" + (2.0 * Math.PI / fenetre));
            }
        }
        if (a.adaptiveCooling) opt.cooling = new AdaptiveCoolingSchedule();
        if (a.minTherm != null) opt.thermalization.minSweeps = a.minTherm;
        if (a.maxTherm != null) opt.thermalization.maxSweeps = a.maxTherm;
        if (a.targetNeff != null) opt.thermalization.targetNeffProduction = a.targetNeff;
        if (a.gewekeSkip != null) opt.thermalization.gewekeSkipFraction = a.gewekeSkip;
        if (a.nSQW != null) opt.nSQW = a.nSQW;
        if (a.maxProductionSweeps != null) opt.thermalization.maxProductionSweeps = a.maxProductionSweeps;
        if (a.minNeff != null) opt.thermalization.minNeff = a.minNeff;
        if (a.tauMultiplier != null) opt.thermalization.tauMultiplier = a.tauMultiplier;
        if (a.memoryMultiplier != null) opt.memoryMultiplier = a.memoryMultiplier;

        // ---- Diagnostic de decorrelation et arret par stagnation (voir --preset ci-dessus et
        // --decorrelation-log ci-dessous). opt.stagnation reste null (comportement historique)
        // sauf si --preset screening ou --decorrelation-log l'exigent.
        if (a.stagnationDetection) {
            opt.stagnation = new StagnationDetector.Config();
            opt.stopThermalizationOnStagnation = a.stopThermalizationOnStagnation;
            opt.stopProductionOnSaturation = a.stopProductionOnSaturation;
        }
        if (a.decorrelationLog && opt.stagnation == null) {
            // Mode diagnostic seul (voir USAGE) : le detecteur tourne et alimente le journal,
            // mais aucun des deux commutateurs d'arret anticipe n'est active ci-dessus, donc la
            // trajectoire de la chaine n'est pas modifiee par cette seule option.
            opt.stagnation = new StagnationDetector.Config();
        }

        // Un repertoire de sortie PAR replique : les fichiers produits (EquilibrationLog,
        // DynamicStructureFactor) portent des noms DETERMINISTES (plus d'horodatage ni de
        // compteur statique dans les noms de fichiers, donc plus de course sur le compteur),
        // c'est le repertoire qui porte l'identite du lancement : date-heure, fichier d'entree,
        // champ. Deux repliques ne partagent jamais de repertoire, donc jamais de nom.
        final String uniqueID = new SimpleDateFormat("ddMMyyyy_HHmmss").format(new Date(System.currentTimeMillis()));
        final String outPath = a.outputDir + File.separator + uniqueID + "_" + fileName
                + "_h" + String.format(Locale.US, "%.6f", field);

        // S(Q,w) est ecrit sous le repertoire de la replique (et non sur le chemin historique
        // code en dur dans DynamicStructureFactor) : voir ParallelOptions.sqwOutputDir.
        opt.sqwOutputDir = outPath;
        if (a.spinsPerT) {
            // Instantane des spins a chaque temperature (apres la production, avant la
            // relaxation de fin), nomme avec T pour etre triable par temperature.
            opt.spinSnapshotListener = (T, lat) ->
                    EquilibrationLog.spinsSnapshotAvro(lat, outPath, String.format(Locale.US, "%.5f", T), field);
        }

        // nChunksPerColor fixe (independant de threadsPerReplica) : la trajectoire ne depend que
        // des graines, donc deux runs a nombre de threads different sont comparables bit a bit
        // (tant que threadsPerReplica <= CHUNKS_PER_COLOR ; au-dela on prend threadsPerReplica).
        final int chunks = Math.max(CHUNKS_PER_COLOR, a.threadsPerReplica);
        // StagnationLog est une ressource optionnelle : null quand --decorrelation-log n'est
        // pas demande, ce que le try-with-resources gere nativement (une ressource null n'est
        // simplement pas fermee).
        try (EquilibrationLog log =
                     new EquilibrationLog(EquilibrationLog.defaultPath(outPath));
             StagnationLog stagLog = a.decorrelationLog
                     ? new StagnationLog(StagnationLog.defaultPath(outPath))
                     : null;
             CheckerboardMetropolis sweeper =
                     new CheckerboardMetropolis(lattice, coloring, a.threadsPerReplica, chunks, seed)) {
            opt.log = log;
            opt.stagnationLog = stagLog;
            mc.runParallel(lattice, spin, uc, fileName, H_field, field, sweeper, opt);
        }

        EquilibrationLog.spinOutputAvro(lattice, outPath, field);
        // EquilibrationLog avale les IOException des avro de spins (printStackTrace) : on
        // verifie que les fichiers existent (spins en avro ; les observables sont FUSIONNEES
        // dans le csv d'equilibration), sinon la replique est comptee en echec (code de
        // sortie 1).
        final File outDir = new File(outPath);
        final String[] produced = outDir.list();
        boolean spins = false, obs = false;
        if (produced != null) {
            for (String name : produced) {
                if (name.startsWith("spins_finaux_") && name.endsWith(".avro")) spins = true;
                if (name.equals("equilibration.csv")) obs = true;
            }
        }
        if (!spins || !obs) {
            throw new IllegalStateException("sortie absente dans " + outPath
                    + " (spins_finaux avro=" + spins + ", equilibration.csv=" + obs + ")");
        }

        System.out.println(Thread.currentThread().getName() + " done " + fileName + " h=" + field);
    }

    // ------------------------------------------------------------------ susceptibilite AC

    /**
     * Construit une replique et mesure sa susceptibilite AC ({@link ACSusceptibility#run}), sur
     * la meme recette de construction que {@link #replica} ({@link ReplicaFactory} + {@link Input},
     * champ statique par {@link Lattice#setInteractionField}, coloration, {@code chunks =
     * max(CHUNKS_PER_COLOR, threadsPerReplica)}) -- mais sans passer par {@link MonteCarlo} ni
     * par les ecrivains Avro de {@code EquilibrationLog} : {@link ACSusceptibility} gere
     * elle-meme son equilibration et sa mesure, et ne produit aucun fichier Avro ni journal
     * d'equilibration.
     *
     * @param replicaIdx indice de replique parmi {@code a.acReplicas}, 0-based ; entre dans la
     *                   derivation de {@code seedAC} (voir {@link #splitmix64Finalize}) pour que
     *                   deux repliques AC du meme (fichier, champ) ne partagent jamais leur
     *                   trajectoire.
     * @return les echantillons AC de cette replique, un par couple (periode, temperature), dans
     *         l'ordre produit par {@link ACSusceptibility#run}.
     */
    private static List<ACSusceptibility.Sample> acReplica(File file, double field, int replicaIdx,
                                                            Args a, AtomicBoolean coloringPrinted) {
        final String fileName = file.getName();

        final ReplicaFactory.Built built = ReplicaFactory.build(new Input(file.getAbsolutePath()),
                                                               a.La, a.Lb);
        final Lattice lattice = built.lattice();

        // Champ statique dans la base de Kitaev, via le tenseur g -- meme calcul que replica().
        final double[] hStatic = MathOps.scale(field,
                MathOps.matrixXvector(lattice.getGtensor(), MathOps.normalize(a.fieldDir)));
        lattice.setInteractionField(hStatic);

        final LatticeColoring coloring = LatticeColoring.of(lattice);
        if (coloringPrinted.compareAndSet(false, true)) {
            System.out.println(coloring.summary());
        }

        // Graine : voir splitmix64Finalize. "base" est EXACTEMENT la derivation de replica() (donc
        // la graine que prendrait le recuit pour ce (fichier, champ)) ; l'offset 0xAC... et le
        // multiple de replicaIdx par la constante d'or garantissent que seedAC n'est jamais egale
        // a la graine de recuit, et que deux repliques AC du meme (fichier, champ) sont deux a
        // deux distinctes -- tout en restant une fonction pure des entrees (reproductible).
        final long base = a.seed ^ (fileName.hashCode() * 31L + Double.doubleToLongBits(field));
        final long seedAC = splitmix64Finalize(base + 0xAC00000000000000L
                + replicaIdx * 0x9E3779B97F4A7C15L);

        // Configuration initiale : meme motif que --bench (randomize, rejet de Marsaglia), SANS
        // passer par MonteCarlo (qui utilise Spin.SRG7 -- le meme algorithme, mais couple a une
        // replique de recuit que la susceptibilite AC ne construit pas).
        randomize(lattice, seedAC);

        final int chunks = Math.max(CHUNKS_PER_COLOR, a.threadsPerReplica);

        final ACSusceptibility.Config cfg = new ACSusceptibility.Config();
        cfg.h0 = a.acH0;
        cfg.periodsSweeps = a.acPeriods.clone();
        cfg.transientPeriods = a.acTransientPeriods;
        cfg.measurePeriods = a.acMeasurePeriods;
        cfg.tMax = a.acTempMax;
        cfg.tMin = a.acTempMin;
        cfg.tStep = a.acTempStep;
        // Sur-relaxation EQUILIBRATION SEULEMENT : ACSusceptibility.run coupe deja la sur-relaxation
        // pendant la mesure elle-meme (voir sa javadoc de classe).
        cfg.overRelaxationPerSweep = a.useOverRelaxation ? a.overrelaxN : 0;
        // Meme petit bloc que replica() (minTherm/maxTherm/targetNeff/gewekeSkip) ; pas de
        // maxProductionSweeps/minNeff/tauMultiplier ici, ACSusceptibility n'a pas de phase de
        // production a etendre.
        if (a.minTherm != null) cfg.thermalization.minSweeps = a.minTherm;
        if (a.maxTherm != null) cfg.thermalization.maxSweeps = a.maxTherm;
        if (a.targetNeff != null) cfg.thermalization.targetNeffProduction = a.targetNeff;
        if (a.gewekeSkip != null) cfg.thermalization.gewekeSkipFraction = a.gewekeSkip;
        cfg.verbose = !a.quiet;

        final double[] acDrive = MathOps.matrixXvector(lattice.getGtensor(),
                MathOps.normalize(a.acFieldDir != null ? a.acFieldDir : a.fieldDir));

        try (CheckerboardMetropolis sweeper =
                     new CheckerboardMetropolis(lattice, coloring, a.threadsPerReplica, chunks, seedAC)) {
            return ACSusceptibility.run(lattice, sweeper, hStatic, acDrive, cfg);
        }
    }

    /**
     * Finaliseur (fonction d'avalanche) de SplitMix64 : trois passes xorshift/multiplication
     * (constantes de David Stafford, variante "Mix13", celles de {@code java.util.SplittableRandom})
     * qui dispersent un {@code long} d'entree en une sortie a haute avalanche. Contrairement au
     * generateur SplitMix64 complet, cette methode n'incremente PAS l'etat par la constante d'or
     * avant de melanger : l'appelant ({@link #acReplica}) l'a deja fait explicitement dans le
     * calcul de {@code seedAC} (l'offset {@code 0xAC...} et le multiple de {@code replicaIdx}),
     * donc appliquer ce finaliseur pur et sans effet de bord suffit a obtenir une graine reproductible,
     * distincte de la graine de recuit et deux a deux distincte entre repliques AC.
     */
    private static long splitmix64Finalize(long z) {
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    // ------------------------------------------------------------------ --bench

    private static final String BENCH_USAGE = String.join("\n",
        "Usage :",
        "  java -cp \"bin:lib/*\" AppParallel --bench <La> <Lb> [J1 J2 J3]",
        "  java -cp \"bin:lib/*\" AppParallel --bench <fichierInput> <La> <Lb>");

    /**
     * Dispatche entre les deux formes de {@code --bench} sur la nature du premier argument :
     * un entier est {@code La} (modele honeycomb synthetique), tout le reste est un chemin de
     * fichier d'entree ({@link #benchFromInput}). {@code Integer.parseInt} est essaye en premier
     * et non {@code new File(argv[0]).isFile()} : un nom de fichier purement numerique est bien
     * plus improbable qu'un {@code La} negatif ou mal forme, et cet ordre garde l'erreur d'un
     * entier invalide (« Argument invalide ») plutot que « fichier introuvable ».
     */
    private static int bench(String[] argv) {
        if (argv.length < 2) {
            System.out.println(BENCH_USAGE);
            return EXIT_USAGE;
        }
        final int La;
        try {
            La = Integer.parseInt(argv[0]);
        } catch (NumberFormatException e) {
            return benchFromInput(argv);
        }
        final int Lb;
        final double J1, J2, J3;
        try {
            Lb = Integer.parseInt(argv[1]);
            J1 = argv.length > 2 ? Double.parseDouble(argv[2]) : 1.d;
            J2 = argv.length > 3 ? Double.parseDouble(argv[3]) : 0.2d;
            J3 = argv.length > 4 ? Double.parseDouble(argv[4]) : 0.1d;
        } catch (NumberFormatException e) {
            System.out.println("Argument invalide : " + e.getMessage());
            System.out.println(BENCH_USAGE);
            return EXIT_USAGE;
        }

        final Lattice lattice = honeycomb(La, Lb, J1, J2, J3);
        randomize(lattice, 20260827L);
        final LatticeColoring coloring = LatticeColoring.of(lattice);
        final String modelLine = String.format(Locale.US,
                "honeycomb %dx%dx2 (N = %d), J1=%.4g J2=%.4g J3=%.4g", La, Lb, lattice.size, J1, J2, J3);
        return benchReport(lattice, coloring, modelLine, La, Lb);
    }

    /**
     * Forme {@code --bench <fichierInput> <La> <Lb>} : mesure le debit sur le modele reellement lu
     * dans {@code fichierInput}, construit exactement comme le mode production ({@link ReplicaFactory}).
     */
    private static int benchFromInput(String[] argv) {
        if (argv.length < 3) {
            System.out.println(BENCH_USAGE);
            return EXIT_USAGE;
        }
        final File file = new File(argv[0]);
        // isFile() : erreur d'usage claire avant la construction (Input est fail-fast).
        if (!file.isFile()) {
            System.out.println("Le fichier d'entree n'existe pas ou n'est pas un fichier : " + argv[0]);
            return EXIT_USAGE;
        }
        final int La, Lb;
        try {
            La = Integer.parseInt(argv[1]);
            Lb = Integer.parseInt(argv[2]);
        } catch (NumberFormatException e) {
            System.out.println("Argument invalide : " + e.getMessage());
            System.out.println(BENCH_USAGE);
            return EXIT_USAGE;
        }

        final Lattice lattice = ReplicaFactory.build(new Input(file.getAbsolutePath()), La, Lb).lattice();
        randomize(lattice, 20260827L);
        final LatticeColoring coloring = LatticeColoring.of(lattice);
        final int neighborsPerSite = lattice.getNeighborSites(0).length;
        final String modelLine = String.format(Locale.US, "%s %dx%dx2 (N = %d), %d voisins/site",
                file.getName(), La, Lb, lattice.size, neighborsPerSite);
        return benchReport(lattice, coloring, modelLine, La, Lb);
    }

    /**
     * Tableau de debit (sweeps/s pour {@code nThreads} = 1, 2, 4, ... jusqu'au nombre de coeurs) et
     * recommandation de {@code threadsPerReplica}, commun aux deux formes de {@code --bench}.
     *
     * <p>{@code La}, {@code Lb} ne servent ici qu'a l'exemple de ligne de commande final et au
     * conseil sur la coloration ; le reseau et sa coloration sont deja construits par l'appelant.</p>
     */
    private static int benchReport(Lattice lattice, LatticeColoring coloring, String modelLine,
                                   int La, int Lb) {
        final int cores = Runtime.getRuntime().availableProcessors();

        System.out.println("CPU disponibles : " + cores);
        System.out.println("Modele : " + modelLine);
        System.out.println(coloring.summary());
        System.out.printf(Locale.US, "%d mises a jour par phase de couleur%n",
                lattice.size / coloring.numColors());
        // Conseil purement informatif (jamais d'exception) : sur un reseau nid d'abeille a deux
        // sous-reseaux, La et Lb multiples de 3 permettent souvent une coloration structurelle a
        // 6 couleurs (cf. README § 3) ; au-dela, DSATUR ou un m plus grand degradent l'equilibrage.
        if (coloring.numColors() > 6 && (La % 3 != 0 || Lb % 3 != 0)) {
            System.out.printf(Locale.US,
                    "Conseil : La et Lb multiples de 3 ramenent souvent la coloration a 6 couleurs "
                    + "sur ce type de reseau (%d ici) ; essayer --bench avec La, Lb multiples de 3.%n",
                    coloring.numColors());
        }
        System.out.println();

        System.out.printf(Locale.US, "%-9s %-14s %-12s %-16s%n",
                "threads", "sweeps/s", "vs 1 thread", "maj/thread/phase");
        System.out.println("-".repeat(56));

        double base = Double.NaN, best = -1.d;
        int bestThreads = 1;
        for (int nt = 1; nt <= cores; nt *= 2) {
            final double rate = benchRate(lattice, coloring, nt);
            if (nt == 1) base = rate;
            if (rate > best) { best = rate; bestThreads = nt; }
            System.out.printf(Locale.US, "%-9d %-14.1f %-12s %-16d%n", nt, rate,
                    String.format(Locale.US, "%.2fx", rate / base),
                    lattice.size / coloring.numColors() / nt);
        }

        System.out.println();
        System.out.printf(Locale.US,
                "Recommandation : threadsPerReplica = %d (%.1f sweeps/s, %.2fx)%n",
                bestThreads, best, best / base);
        System.out.printf(Locale.US,
                "  -> avec totalThreads = T, lancer concurrentReplicas = T / %d repliques simultanees%n",
                bestThreads);
        System.out.printf(Locale.US, "  -> par exemple : java -cp \"bin:lib/*\" AppParallel <in> <out> %d %d %d %d%n",
                cores, bestThreads, La, Lb);
        System.out.println("Rappel : il faut ~10^3 mises a jour par thread et par phase de couleur "
                + "pour amortir la barriere ; en dessous, preferer plus de repliques.");
        return EXIT_OK;
    }

    private static double benchRate(Lattice lattice, LatticeColoring coloring, int nThreads) {
        final double T = 1.d, SIGMA = 1.d, MIN_SECONDS = 1.0d;
        try (CheckerboardMetropolis cb = new CheckerboardMetropolis(lattice, coloring, nThreads, nThreads, 1L)) {
            for (int s = 0; s < 50; s++) cb.sweep(T, SIGMA);
            long sweeps = 0L;
            final long t0 = System.nanoTime();
            long elapsed;
            do {
                for (int s = 0; s < 25; s++) cb.sweep(T, SIGMA);
                sweeps += 25;
                elapsed = System.nanoTime() - t0;
            } while (elapsed < MIN_SECONDS * 1e9);
            return sweeps / (elapsed / 1e9);
        }
    }

    // ------------------------------------------------------------------ reseau de test local
    // La geometrie (vecteurs primitifs, sites de base, couches de voisins) vient de
    // HoneycombGeometry, la meme source que celle utilisee pour verifier ce chemin depuis test/ :
    // --bench reste utilisable avec bin/ seul, sans bin_test/ sur le class path, sans dupliquer la
    // boucle de recherche de voisins.

    /** Honeycomb La x Lb x 2, conventions geometriques de {@link HoneycombGeometry}. */
    private static Lattice honeycomb(int La, int Lb, double J1, double J2, double J3) {
        final UnitCell uc = new UnitCell(HoneycombGeometry.a1(), HoneycombGeometry.a2());
        final int b0 = uc.addBasisSite(HoneycombGeometry.basisFractional(0));
        final int b1 = uc.addBasisSite(HoneycombGeometry.basisFractional(1));
        uc.setInteractionGtensor(b0, HoneycombGeometry.identity());
        uc.setInteractionGtensor(b1, HoneycombGeometry.identity());

        final Lattice lattice = new Lattice(
                new int[] { La, Lb, HoneycombGeometry.SUBLATTICES });
        lattice.setKitaevBasis(HoneycombGeometry.identity());
        lattice.setCOB(HoneycombGeometry.identity());
        lattice.setGtensor(HoneycombGeometry.identity());
        lattice.setInteractionField(new double[3]);

        // UnitCell.addInteraction n'ajoute la liste du sous-reseau 1 que si interactions.size()
        // <= 1 : on pre-cree les deux listes pour que le sous-reseau 1 fonctionne meme si le 0
        // est vide.
        while (uc.interactions.size() < HoneycombGeometry.SUBLATTICES) {
            uc.interactions.add(new ArrayList<double[]>());
        }

        final double[] J = { J1, J2, J3 };
        for (int s = 0; s < HoneycombGeometry.SUBLATTICES; s++) {
            for (int shell = 1; shell <= 3; shell++) {
                if (J[shell - 1] == 0.d) continue;
                for (int[] v : HoneycombGeometry.neighbours(s, shell)) {
                    uc.addInteraction(s, v[2], HoneycombGeometry.diag(J[shell - 1]),
                            new byte[] { (byte) v[0], (byte) v[1], (byte) v[2] });
                }
            }
        }
        lattice.createTemplates(uc);
        lattice.setInteractionSites(uc);
        return lattice;
    }

    private static void randomize(Lattice lattice, long seed) {
        final RandomGenerator r = RandomGeneratorFactory.of("Xoshiro256PlusPlus").create(seed);
        final double[] S = lattice.getSpin1Dlattice();
        for (int i = 0; i < lattice.size; i++) {
            double x0, x1, rr;
            do {
                x0 = r.nextDouble() * 2 - 1;
                x1 = r.nextDouble() * 2 - 1;
                rr = x0 * x0 + x1 * x1;
            } while (rr >= 1);
            final double c = 2 * Math.sqrt(1 - rr);
            S[3 * i] = c * x0;
            S[3 * i + 1] = c * x1;
            S[3 * i + 2] = 1 - 2 * rr;
        }
    }
}
