import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

/**
 * Calcule, a partir des journaux CSV d'un run reel, le triplet
 * {@code thermalizationSweeps / thermalizationRate / measurementRate} a utiliser pour du
 * <b>criblage rapide</b> de modeles (voir {@code MonteCarlo}, {@code AdaptiveThermalization},
 * {@code StagnationDetector}) : un resultat valide, pas forcement excellent, obtenu vite,
 * plutot que la thermalisation fixe historique
 * ({@code thermalizationSweeps = 40 000}, {@code thermalizationRate = 80 000},
 * {@code measurementRate = 100}) qui traite toutes les temperatures comme si elles etaient
 * aussi lentes que la transition.
 *
 * <h2>Le probleme : une distribution de tau tres etalee</h2>
 * <p>Sur un recuit de 30x30 spins pres d'une transition de phase, le temps de correlation
 * integre tau (en sweeps) varie de moins d'un sweep loin de la transition a plusieurs
 * centaines de sweeps au point critique (ralentissement critique, tau ~ xi^z). Une seule
 * constante de thermalisation/mesure ne peut donc jamais etre <i>a la fois</i> rapide loin de
 * la transition et suffisante pres d'elle. Cette classe ne cherche pas a resoudre ce dilemme :
 * elle le rend explicite, en choisissant un quantile de la distribution de tau comme objectif
 * ({@link Config#quantile}) et en <b>comptant</b>, plutot qu'en cachant, les temperatures plus
 * lentes que ce quantile qui resteront sous-echantillonnees ({@link Recommendation#nUndersampled()},
 * {@link Recommendation#nUnderThermalized()}). C'est un arbitrage de criblage assume, pas une
 * pretention a la decorrelation universelle.</p>
 *
 * <h2>La formule pivot : de g a tau</h2>
 * <p>Une serie de production echantillonnee tous les {@code Delta} sweeps a une inefficacite
 * statistique {@code g_Delta = 1 + 2 tau_sweeps / Delta} (convention de Chodera, voir
 * {@code EquilibrationDetector}) : c'est simplement la definition de g appliquee a un temps de
 * correlation exprime en <i>echantillons</i> ({@code tau_sweeps / Delta}) plutot qu'en sweeps.
 * On en tire :</p>
 * <pre>
 *   tau_sweeps = (g_Delta - 1) / 2 * Delta
 * </pre>
 * <p>C'est cette formule qui permet de reconstruire un temps de correlation en <b>sweeps</b>,
 * comparable d'une temperature a l'autre, alors que chaque temperature a pu etre echantillonnee
 * a une cadence {@code Delta = measureEvery} differente. Elle est utilisee partout ou
 * {@code tauProductionSweeps} n'est pas fourni directement par le journal (voir
 * {@link #readEquilibrationCsv(Path)}).</p>
 *
 * <h2>Deux journaux, une seule notion de temperature</h2>
 * <p>Le journal d'equilibration ({@code equilibration_*.csv}, ecrit par {@code EquilibrationLog})
 * dit combien de sweeps la thermalisation a consommes et quelle etait sa propre inefficacite
 * statistique. Le journal de decorrelation ({@code decorrelation_*.csv}, ecrit par
 * {@code StagnationLog}) dit, en plus, si la configuration de spins retenue est decorrelee de
 * celle heritee de la temperature precedente (recouvrement q, voir {@code SpinOverlap},
 * {@code StagnationPoint}). {@link #merge(List, List)} fusionne les deux par temperature ; en
 * l'absence du second journal, cette classe fonctionne quand meme (avec
 * {@code tauOverlapSweeps = NaN}, {@code temperaturesDecorrelated = false}) : elle degrade
 * proprement vers le seul critere de production.</p>
 */
public final class SweepCalibrator {

    private SweepCalibrator() { }

    /**
     * Mise en garde valable pour TOUT triplet evalue par cette classe a partir d'un journal
     * donne, quel que soit le scenario ({@link #recommend}, {@link #evaluateFixed}, ou la ligne
     * "historique" du tableau d'arbitrage de {@link #format}) : un texte constant, parce que
     * cette mise en garde ne depend d'aucun parametre, elle est vraie dans tous les cas.
     *
     * <p>{@code tauProductionSweeps} est mesure (ou reconstruit) a partir de la dynamique Monte
     * Carlo qui a produit le journal. Si cette dynamique inclut de la sur-relaxation
     * microcanonique ({@code --overrelax}), tau y est fortement raccourci par rapport a un
     * Metropolis local pur, dont chaque pas ne fait que proposer un flip local. Mesure directe
     * sur ce meme modele a la transition (T = 0.7309767973), sur des series de production de
     * 400 000 sweeps : tau ~= 650-770 sweeps AVEC sur-relaxation (2 sweeps microcanoniques par
     * sweep Metropolis) contre ~900 a 9 500 sweeps SANS, selon le bloc de la serie ; sur
     * l'ensemble de cette serie sans sur-relaxation, le detecteur de Chodera ne trouve meme
     * aucune zone stationnaire (t0 = 384 000 sur 400 000 : la serie n'a pas fini de relaxer).</p>
     *
     * <p>Transposer un triplet d'une dynamique a une autre n'est valide que si tau l'est aussi.
     * Evaluer le triplet historique (ecrit pour un Metropolis local pur) avec les tau de CE
     * journal (mesures avec sur-relaxation) flatte donc ce triplet d'un facteur de l'ordre de
     * 20 a la transition : la ligne "historique" du tableau d'arbitrage sous-estime dans ces
     * proportions le besoin reel du code d'origine, elle ne le disqualifie certainement pas
     * moins que ce que montrent ces chiffres.</p>
     */
    public static final String DYNAMICS_WARNING =
            "ATTENTION dynamique : les tau de ce journal ont ete mesures avec la dynamique qui a "
            + "produit le journal (ici, sur-relaxation microcanonique --overrelax), pas avec un "
            + "Metropolis local pur. Mesure directe sur ce modele a T=0.731 : tau ~= 650-770 "
            + "sweeps AVEC sur-relaxation contre ~900-9500 SANS (Chodera ne trouve meme pas de "
            + "zone stationnaire sur 400000 sweeps de Metropolis pur a cette temperature). "
            + "Transposer un triplet d'une dynamique a une autre n'est valide que si tau l'est "
            + "aussi : la ligne 'historique' ci-dessus, evaluee avec CES tau, flatte le "
            + "Metropolis local pur d'un facteur de l'ordre de 20 a la transition.";

    // ------------------------------------------------------------------------------- Row

    /**
     * Une temperature, telle que relue des journaux (eventuellement fusionnee par {@link #merge}).
     *
     * @param T                          temperature.
     * @param thermSweeps                sweeps de thermalisation reellement consommes.
     * @param t0                         debut de la zone equilibree (Chodera), en sweeps.
     * @param measureEvery               cadence de mesure de production effectivement utilisee
     *                                   par le run (pas la cadence recommandee par cette classe).
     * @param gProduction                g de la serie de production (methode Gamma de Wolff),
     *                                   {@code NaN} si absent du journal.
     * @param neffAchieved               N_eff reellement obtenu sur la production.
     * @param prodSweeps                 sweeps de production reellement effectues.
     * @param tauProductionSweeps        tau_int de la serie de production, en sweeps : fourni
     *                                   directement par le journal de decorrelation, ou
     *                                   reconstruit depuis {@code gProduction} sinon (voir la
     *                                   formule pivot en tete de classe). {@code NaN} si ni l'un
     *                                   ni l'autre n'est disponible.
     * @param tauOverlapSweeps           temps de decroissance du recouvrement de configuration q
     *                                   vers son plateau, en sweeps ({@code NaN} si non mesure :
     *                                   pas de journal de decorrelation, ou premiere temperature
     *                                   du recuit).
     * @param qInf                       plateau du recouvrement q (parametre d'ordre au sens
     *                                   Edwards-Anderson, {@code NaN} si non mesure).
     * @param temperaturesDecorrelated   vrai si le journal de decorrelation atteste que la
     *                                   configuration retenue n'est plus, pour l'essentiel,
     *                                   celle heritee de la temperature precedente. Neutre
     *                                   (false) si non mesure.
     * @param measurementsDecorrelated   vrai si le journal atteste que deux mesures successives
     *                                   de production sont decorrelees. Neutre (false) si non
     *                                   mesure.
     */
    public record Row(double T, long thermSweeps, long t0, int measureEvery, double gProduction,
                       double neffAchieved, long prodSweeps, double tauProductionSweeps,
                       double tauOverlapSweeps, double qInf,
                       boolean temperaturesDecorrelated, boolean measurementsDecorrelated) { }

    // ---------------------------------------------------------------------------- Config

    /** Parametres de la recommandation. Champs publics + setters fluides, comme
     *  {@code AdaptiveThermalization.Config} : un objet de configuration pur, pas de logique. */
    public static final class Config {
        /**
         * Quantile de la distribution de tau_prod retenu comme objectif de criblage. Les
         * temperatures plus lentes que ce quantile sont deliberement sacrifiees (comptees dans
         * {@link Recommendation#nUndersampled()} / {@link Recommendation#nUnderThermalized()})
         * plutot que de dimensionner le triplet entier sur le pire cas (souvent 10 a 100 fois
         * plus lent que la mediane pres d'une transition).
         */
        public double quantile = 0.90;
        /** N_eff vise en production, au quantile retenu. */
        public double targetNeff = 100.0;
        /**
         * K : la thermalisation doit couvrir au moins {@code K * tau_prod} au-dela de t0, pour
         * que la queue de production commence loin du transitoire (meme marge que
         * {@code AdaptiveThermalization.Config#tauMultiplier}).
         */
        public double tauMultiplier = 20.0;
        /**
         * K_q : la thermalisation doit aussi couvrir au moins {@code K_q * tau_overlap}, pour
         * que la memoire de la configuration heritee de la temperature precedente ait decru vers
         * son plateau (voir {@code StagnationPoint}). K_q = 3 laisse une memoire residuelle
         * e^-3 ~ 5 % (Ogielski, PRB 32, 7384 (1985)).
         */
        public double memoryMultiplier = 3.0;
        /** Multiple de tau_prod exige entre deux mesures pour les declarer decorrelees. */
        public double measurementsPerTau = 2.0;

        public Config() { }

        public Config quantile(double v) { this.quantile = v; return this; }
        public Config targetNeff(double v) { this.targetNeff = v; return this; }
        public Config tauMultiplier(double v) { this.tauMultiplier = v; return this; }
        public Config memoryMultiplier(double v) { this.memoryMultiplier = v; return this; }
        public Config measurementsPerTau(double v) { this.measurementsPerTau = v; return this; }

        public Config copy() {
            Config c = new Config();
            c.quantile = quantile;
            c.targetNeff = targetNeff;
            c.tauMultiplier = tauMultiplier;
            c.memoryMultiplier = memoryMultiplier;
            c.measurementsPerTau = measurementsPerTau;
            return c;
        }
    }

    // ----------------------------------------------------------------------- Recommendation

    /**
     * Le triplet demande, ainsi que ce qu'il coute et ce qu'il sacrifie.
     *
     * <p><b>{@link #DYNAMICS_WARNING}</b> s'applique a toute instance de ce record, qu'elle
     * vienne de {@link #recommend} ou de {@link #evaluateFixed} : tous les {@code tau} utilises
     * ici viennent de la dynamique du journal fourni, pas necessairement de celle du triplet
     * evalue.</p>
     *
     * @param thermalizationSweeps    equivalent de la constante historique du meme nom.
     * @param thermalizationRate      = {@code thermalizationSweeps + productionSweeps} :
     *                                equivalent de la constante historique {@code thermalizationRate}
     *                                (l'index de sweep auquel s'arreter, thermalisation comprise).
     * @param measurementRate         equivalent de la constante historique {@code measurementRate}.
     * @param productionSweeps        = {@code thermalizationRate - thermalizationSweeps}.
     * @param tauMedianSweeps         mediane de tau_prod sur les temperatures fournies.
     * @param tauQuantileSweeps       tau_prod au quantile {@link Config#quantile}. Pour
     *                                {@link #recommend}, c'est l'objectif de criblage qui a
     *                                servi a dimensionner le triplet. Pour {@link #evaluateFixed}
     *                                (triplet impose), cette valeur est purement informative :
     *                                elle n'a pas servi a deriver le triplet, qui vient d'ailleurs.
     * @param tauMaxSweeps            tau_prod maximal observe (typiquement a la transition).
     * @param slowestT                temperature ou tau_prod est maximal.
     * @param nTemperatures           nombre de temperatures fournies.
     * @param nUndersampled           nombre de temperatures (calcule sur TOUTES les temperatures
     *                                fournies, pas seulement les retenues) ou
     *                                {@code measurementRate < 2 * tau_prod} : deux mesures
     *                                scalaires successives de production y restent correlees
     *                                avec ce triplet. <b>Ce n'est pas un defaut de justesse</b>
     *                                des lors que g est mesure sur la production elle-meme
     *                                (mode adaptatif, {@code gProduction}) : N_eff et les barres
     *                                d'erreur en tiennent deja compte, on paie seulement des
     *                                mesures redondantes. Voir {@link #verdict()} pour la
     *                                distinction complete avec {@code nUnderThermalized}.
     * @param nUnderThermalized       nombre de temperatures (calcule sur TOUTES les temperatures
     *                                fournies, pas seulement les retenues) ou
     *                                {@code thermalizationSweeps} est inferieur au minimum requis
     *                                pour cette temperature ({@code t0 + tauMultiplier * tau_prod},
     *                                et {@code memoryMultiplier * tauOverlapSweeps} si mesure) :
     *                                la configuration de spins retenue y reste partiellement
     *                                heritee de la temperature precedente avec ce triplet.
     *                                <b>Ceci, contrairement a {@code nUndersampled}, est un vrai
     *                                probleme</b> : la configuration produite n'est pas celle de
     *                                cette temperature.
     * @param sweepsPerReplica        cout total du recuit (sweeps, toutes temperatures) avec ce
     *                                triplet : {@code nTemperatures * thermalizationRate}.
     * @param legacySweepsPerReplica  cout avec le triplet historique fixe (40 000 / 80 000 / 100) :
     *                                {@code nTemperatures * 80 000}.
     * @param adaptiveSweepsPerReplica cout reellement observe dans le journal fourni : somme sur
     *                                les temperatures de {@code thermSweeps + prodSweeps}.
     * @param verdict                 phrase(s) honnetes en francais expliquant ce que ce triplet
     *                                garantit et ce qu'il ne garantit pas, et rappelant la
     *                                distinction physique centrale : une mesure de production
     *                                "sous-echantillonnee" n'est pas une mesure fausse des lors
     *                                que g est mesure sur la production elle-meme (le mode
     *                                adaptatif le fait) -- on paie une redondance, pas un biais.
     *                                Ce qui EST faux, c'est un triplet fixe qui ne mesure jamais
     *                                g et suppose des mesures independantes : il sous-estime
     *                                alors l'erreur d'un facteur sqrt(g) (a la transition de ce
     *                                fichier de reference, g ~ 8.8, soit un facteur ~3 sur
     *                                l'erreur annoncee par le code historique). La vraie exigence
     *                                de decorrelation porte sur la configuration de spins
     *                                conservee et sur deux temperatures successives
     *                                ({@code nUnderThermalized}), pas sur deux mesures scalaires
     *                                de production ({@code nUndersampled}) des lors que g est
     *                                mesure.
     */
    public record Recommendation(
            int thermalizationSweeps,
            int thermalizationRate,
            int measurementRate,
            int productionSweeps,
            double tauMedianSweeps, double tauQuantileSweeps, double tauMaxSweeps,
            double slowestT,
            int nTemperatures,
            int nUndersampled,
            int nUnderThermalized,
            long sweepsPerReplica,
            long legacySweepsPerReplica,
            long adaptiveSweepsPerReplica,
            String verdict) { }

    // -------------------------------------------------------------------------- lecture CSV

    /** Table CSV minimale : index des colonnes par nom, et les lignes de donnees (par champs). */
    private record Table(Map<String, Integer> index, List<String[]> rows) { }

    /**
     * Coupe un fichier CSV en table nom-de-colonne -&gt; index, indifferente a l'ordre des
     * colonnes. Une ligne dont le nombre de champs ne correspond pas a l'en-tete est une erreur
     * fatale (fichier corrompu) : on refuse de deviner, mais on dit exactement quelle ligne pose
     * probleme pour que l'utilisateur puisse aller la regarder.
     *
     * @throws UncheckedIOException si le fichier ne peut pas etre lu.
     * @throws IllegalArgumentException si le fichier est vide, ou si une ligne n'a pas le meme
     *         nombre de champs que l'en-tete (message : numero de ligne 1-based, en-tete comprise).
     */
    private static Table parseCsv(Path csv) {
        List<String> lines;
        try {
            lines = Files.readAllLines(csv, StandardCharsets.UTF_8);
        } catch (IOException e) {
            throw new UncheckedIOException("lecture impossible de " + csv, e);
        }
        if (lines.isEmpty()) throw new IllegalArgumentException("fichier vide : " + csv);
        String[] header = lines.get(0).split(",", -1);
        Map<String, Integer> index = new LinkedHashMap<>();
        for (int i = 0; i < header.length; i++) index.put(header[i], i);
        List<String[]> rows = new ArrayList<>(Math.max(0, lines.size() - 1));
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i);
            if (line.isEmpty()) continue; // derniere ligne vide (fichier termine par un saut de ligne) : toleree
            String[] fields = line.split(",", -1);
            if (fields.length != header.length) {
                throw new IllegalArgumentException("ligne " + (i + 1) + " de " + csv + " : "
                        + fields.length + " champs, " + header.length + " attendus (en-tete)");
            }
            rows.add(fields);
        }
        return new Table(index, rows);
    }

    private static double parseNum(String s) {
        return switch (s) {
            case "NaN" -> Double.NaN;
            case "Inf" -> Double.POSITIVE_INFINITY;
            case "-Inf" -> Double.NEGATIVE_INFINITY;
            default -> Double.parseDouble(s);
        };
    }

    /** Colonne absente -&gt; NaN (jamais d'exception : c'est le contrat documente de l'API). */
    private static double getD(Table t, String[] row, String col) {
        Integer idx = t.index().get(col);
        return idx == null ? Double.NaN : parseNum(row[idx]);
    }

    /** Colonne absente -&gt; {@code neutral} (jamais d'exception). */
    private static long getL(Table t, String[] row, String col, long neutral) {
        Integer idx = t.index().get(col);
        return idx == null ? neutral : Long.parseLong(row[idx]);
    }

    /** Colonne absente -&gt; {@code neutral} (jamais d'exception). */
    private static int getI(Table t, String[] row, String col, int neutral) {
        Integer idx = t.index().get(col);
        return idx == null ? neutral : Integer.parseInt(row[idx]);
    }

    /** Colonne absente -&gt; {@code neutral} (jamais d'exception). {@code "1"} vaut vrai. */
    private static boolean getB(Table t, String[] row, String col, boolean neutral) {
        Integer idx = t.index().get(col);
        return idx == null ? neutral : "1".equals(row[idx]);
    }

    /**
     * Construit un {@link Row} par NOM de colonne : c'est le meme code, qu'il lise le journal
     * d'equilibration ou celui de decorrelation, et il ne suppose jamais un ordre ou une
     * presence de colonne particuliers (voir {@link #readEquilibrationCsv} /
     * {@link #readDecorrelationCsv}).
     */
    private static Row rowFrom(Table t, String[] row) {
        double T = getD(t, row, "T");
        long thermSweeps = getL(t, row, "thermSweeps", 0L);
        long t0 = getL(t, row, "t0", 0L);
        // measureEvery neutre = 1 : une cadence de mesure absente ne doit ni provoquer une
        // exception, ni introduire une division par zero dans la reconstruction de tau ci-dessous.
        int measureEvery = getI(t, row, "measureEvery", 1);
        double gProduction = getD(t, row, "gProduction");
        double neffAchieved = getD(t, row, "neffAchieved");
        long prodSweeps = getL(t, row, "prodSweeps", 0L);

        double tauProductionSweeps = getD(t, row, "tauProductionSweeps");
        if (Double.isNaN(tauProductionSweeps) && !Double.isNaN(gProduction) && measureEvery > 0) {
            // formule pivot, voir le javadoc de tete de classe : tau_sweeps = (g-1)/2 * Delta.
            tauProductionSweeps = (gProduction - 1.0) / 2.0 * measureEvery;
        }

        double tauOverlapSweeps = getD(t, row, "tauOverlapSweeps");
        double qInf = getD(t, row, "qInf");
        boolean temperaturesDecorrelated = getB(t, row, "temperaturesDecorrelated", false);
        boolean measurementsDecorrelated = getB(t, row, "measurementsDecorrelated", false);

        return new Row(T, thermSweeps, t0, measureEvery, gProduction, neffAchieved, prodSweeps,
                tauProductionSweeps, tauOverlapSweeps, qInf, temperaturesDecorrelated, measurementsDecorrelated);
    }

    /**
     * Lit un journal d'equilibration ({@code EquilibrationLog}, colonnes {@code mcIdx,T,
     * thermSweeps,t0,tau,g,neff,zE,zM,prodSweeps,converged,acceptance,sigma,sigmaTherm,ePerSite,
     * eErr,mNorm,seconds,measureEvery,nMeasurements,neffAchieved,gProduction}).
     *
     * <p>Lecture ROBUSTE par nom de colonne : l'ordre des colonnes n'est jamais suppose.
     * {@code tauProductionSweeps} n'existe pas dans ce journal, il est systematiquement
     * reconstruit depuis {@code gProduction} (formule pivot). {@code tauOverlapSweeps},
     * {@code qInf}, {@code temperaturesDecorrelated}, {@code measurementsDecorrelated} n'y
     * existent pas non plus : NaN / false, a completer via {@link #merge} si un journal de
     * decorrelation est disponible.</p>
     */
    public static List<Row> readEquilibrationCsv(Path csv) {
        Table t = parseCsv(csv);
        List<Row> rows = new ArrayList<>(t.rows().size());
        for (String[] fields : t.rows()) rows.add(rowFrom(t, fields));
        return rows;
    }

    /**
     * Lit un journal de decorrelation ({@code StagnationLog}, colonnes {@code mcIdx,T,
     * thermSweeps,prodSweeps,qInf,tauOverlapSweeps,memoryRatio,temperaturesDecorrelated,driftZ,
     * driftPerSweep,stopReason,converged,measureEvery,tauProductionSweeps,measurementsPerTau,
     * measurementsDecorrelated,neffAchieved,neffSaturated,seconds}).
     *
     * <p>Meme mecanique que {@link #readEquilibrationCsv} : lecture par nom de colonne,
     * jamais par position. Ce journal fournit {@code tauProductionSweeps} directement (mesure
     * sur la serie de production elle-meme) : cette valeur est prise telle quelle, la
     * reconstruction par {@code gProduction} ne s'y applique que si elle etait absente (ce
     * journal n'a d'ailleurs pas de colonne {@code gProduction} : {@code Row.gProduction()} y
     * vaut donc NaN, sauf apres {@link #merge} avec le journal d'equilibration).</p>
     */
    public static List<Row> readDecorrelationCsv(Path csv) {
        Table t = parseCsv(csv);
        List<Row> rows = new ArrayList<>(t.rows().size());
        for (String[] fields : t.rows()) rows.add(rowFrom(t, fields));
        return rows;
    }

    /** Deux temperatures sont "la meme" a 1e-9 relatif pres (marge des arrondis d'ecriture CSV,
     *  {@code EquilibrationLog}/{@code StagnationLog} ecrivent en {@code %.10g}). */
    private static boolean sameT(double a, double b) {
        return Math.abs(a - b) <= 1e-9 * Math.max(1.0, Math.abs(a));
    }

    /**
     * Fusionne un journal d'equilibration et un journal de decorrelation par temperature.
     *
     * <p>Pour chaque temperature appariee, les champs propres a la decorrelation
     * ({@code tauOverlapSweeps}, {@code qInf}, {@code temperaturesDecorrelated},
     * {@code measurementsDecorrelated}) viennent du journal de decorrelation ; les autres champs
     * viennent du journal d'equilibration. {@code tauProductionSweeps} prend la valeur du journal
     * de decorrelation quand elle y est finie (mesuree directement sur la production), sinon
     * celle reconstruite depuis le journal d'equilibration.</p>
     *
     * <p><b>Temperatures non appariees</b> : elles sont conservees telles quelles, avec les
     * champs de l'autre journal a leur valeur neutre (NaN / false) plutot qu'exclues. Choix
     * documente : une temperature qui n'a qu'un des deux journaux (run partiel, journal de
     * decorrelation non demande) reste exploitable pour le seul critere que son journal couvre,
     * plutot que de disparaitre silencieusement de la recommandation.</p>
     */
    public static List<Row> merge(List<Row> equilibration, List<Row> decorrelation) {
        List<Row> out = new ArrayList<>(equilibration.size());
        boolean[] used = new boolean[decorrelation.size()];
        for (Row eq : equilibration) {
            int match = -1;
            for (int j = 0; j < decorrelation.size(); j++) {
                if (!used[j] && sameT(eq.T(), decorrelation.get(j).T())) { match = j; break; }
            }
            if (match < 0) {
                out.add(eq);
                continue;
            }
            Row dc = decorrelation.get(match);
            used[match] = true;
            double tauProd = Double.isFinite(dc.tauProductionSweeps())
                    ? dc.tauProductionSweeps() : eq.tauProductionSweeps();
            out.add(new Row(eq.T(), eq.thermSweeps(), eq.t0(), eq.measureEvery(), eq.gProduction(),
                    eq.neffAchieved(), eq.prodSweeps(), tauProd,
                    dc.tauOverlapSweeps(), dc.qInf(),
                    dc.temperaturesDecorrelated(), dc.measurementsDecorrelated()));
        }
        for (int j = 0; j < decorrelation.size(); j++) {
            if (!used[j]) out.add(decorrelation.get(j));
        }
        return out;
    }

    // ----------------------------------------------------------------------------- quantile

    /**
     * Quantile par interpolation lineaire entre rangs (convention dite "type 7", celle de
     * {@code numpy.percentile} et de R par defaut) : sur le tableau trie de taille n, la
     * position fractionnaire est {@code p = q * (n - 1)} et la valeur est interpolee lineairement
     * entre les deux entiers qui l'encadrent. {@code q = 0} renvoie le minimum, {@code q = 1} le
     * maximum, un tableau d'un seul element renvoie cet element quel que soit q. Choisie plutot
     * qu'une convention "rang le plus proche" pour que les chiffres de cette classe se comparent
     * directement a un depouillement fait sous numpy/pandas.
     */
    public static double quantile(double[] values, double q) {
        if (values == null || values.length == 0) throw new IllegalArgumentException("values vide");
        if (!(q >= 0.0 && q <= 1.0)) throw new IllegalArgumentException("q doit etre dans [0,1] : " + q);
        double[] sorted = values.clone();
        Arrays.sort(sorted);
        int n = sorted.length;
        if (n == 1) return sorted[0];
        double pos = q * (n - 1);
        int lo = (int) Math.floor(pos);
        int hi = Math.min(lo + 1, n - 1);
        double frac = pos - lo;
        return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
    }

    // ---------------------------------------------------------------------------- recommend

    /**
     * Calcule le triplet de criblage a partir des temperatures fournies.
     *
     * <h2>Regles de calcul</h2>
     * <ul>
     *   <li>{@code measurementRate = ceil(measurementsPerTau * tauQuantile)}, borne
     *       inferieurement a 1 : c'est la cadence de mesure qui decorrele deux mesures
     *       successives <i>au quantile retenu</i> (donc pour au moins {@code quantile * 100 %}
     *       des temperatures, a l'arrondi pres de {@code ceil}).</li>
     *   <li><b>Temperatures retenues</b> pour dimensionner la thermalisation : celles dont
     *       {@code tauProductionSweeps <= tauQuantile}. Les temperatures plus lentes sont
     *       deliberement sacrifiees : les compter dans le maximum ferait dependre le triplet
     *       entier du pire cas (souvent 100x plus lent que la mediane pres d'une transition),
     *       ruinant l'objectif de criblage rapide. Elles sont comptees dans
     *       {@code nUndersampled} / {@code nUnderThermalized} : sacrifiees, pas cachees.</li>
     *   <li>{@code thermalizationSweeps = ceil(max sur les temperatures retenues de
     *       max(t0 + tauMultiplier*tau, memoryMultiplier*tauOverlap))} ; le second terme du
     *       max n'est applique que si {@code tauOverlapSweeps} est mesure (journal de
     *       decorrelation fusionne), sinon on retombe sur {@code t0 + tauMultiplier*tau} seul.</li>
     *   <li>{@code productionSweeps = ceil(targetNeff * g * measurementRate)} avec
     *       {@code g = 1 + 2*tauQuantile/measurementRate} : c'est directement l'inversion de
     *       {@code N_eff = n/g} au quantile retenu.</li>
     *   <li>{@code thermalizationRate = thermalizationSweeps + productionSweeps}.</li>
     *   <li>{@code nUndersampled} : nombre de temperatures (parmi TOUTES, retenues ou non) ou
     *       {@code measurementRate < 2 * tauProductionSweeps} : deux mesures y restent
     *       correlees avec ce triplet.</li>
     *   <li>{@code nUnderThermalized} : nombre de temperatures (parmi toutes) dont le minimum
     *       requis (meme formule que pour {@code thermalizationSweeps} mais evalue temperature
     *       par temperature) depasse le {@code thermalizationSweeps} retenu.</li>
     *   <li>{@code sweepsPerReplica = nTemperatures * thermalizationRate} ;
     *       {@code legacySweepsPerReplica = nTemperatures * 80 000} ;
     *       {@code adaptiveSweepsPerReplica = somme(thermSweeps + prodSweeps)} observee dans le
     *       journal.</li>
     * </ul>
     *
     * @throws IllegalArgumentException si {@code rows} est vide ou nul, ou si aucune ligne
     *         n'expose de {@code tauProductionSweeps} exploitable (fini).
     */
    public static Recommendation recommend(List<Row> rows, Config cfg) {
        if (rows == null || rows.isEmpty()) throw new IllegalArgumentException("rows vide");
        Config c = cfg == null ? new Config() : cfg;
        int n = rows.size();

        double[] tauFinite = finiteTau(rows);
        double tauMedian = quantile(tauFinite, 0.5);
        double tauQuantile = quantile(tauFinite, c.quantile);
        double tauMax = Arrays.stream(tauFinite).max().orElseThrow();
        Slowest slowest = findSlowest(rows);

        int measurementRate = Math.max(1, (int) Math.ceil(c.measurementsPerTau * tauQuantile));
        double gAtQuantile = 1.0 + 2.0 * tauQuantile / measurementRate;
        int productionSweeps = Math.max(1, (int) Math.ceil(c.targetNeff * gAtQuantile * measurementRate));

        double thermBound = 0.0;
        for (Row r : rows) {
            double tau = r.tauProductionSweeps();
            if (!Double.isFinite(tau) || tau > tauQuantile) continue; // sacrifiee : hors du dimensionnement
            thermBound = Math.max(thermBound, requiredThermSweeps(r, tau, c));
        }
        int thermalizationSweeps = (int) Math.ceil(Math.max(thermBound, 1.0));
        int thermalizationRate = thermalizationSweeps + productionSweeps;

        // Calcule sur TOUTES les temperatures (retenues ou sacrifiees), pas seulement le
        // sous-ensemble retenu pour dimensionner le triplet ci-dessus : sinon le compte serait
        // trompeur, il ne verrait jamais les temperatures qu'on a justement sacrifiees pour
        // aller vite (voir le javadoc de Recommendation). Memes formules que
        // {@link #evaluateFixed}, pour que les deux soient comparables ligne a ligne dans le
        // tableau d'arbitrage de {@link #format}.
        Counts counts = countIssues(rows, measurementRate, thermalizationSweeps, c);

        long sweepsPerReplica = (long) n * thermalizationRate;
        long legacySweepsPerReplica = (long) n * 80_000L;
        long adaptiveSweepsPerReplica = adaptiveCost(rows);

        String verdict = verdict(n, counts.nUndersampled(), counts.nUnderThermalized(),
                slowest.T(), tauMax, slowest.g(), String.format(Locale.US, "quantile %.2f", c.quantile));

        return new Recommendation(thermalizationSweeps, thermalizationRate, measurementRate,
                productionSweeps, tauMedian, tauQuantile, tauMax, slowest.T(), n,
                counts.nUndersampled(), counts.nUnderThermalized(), sweepsPerReplica,
                legacySweepsPerReplica, adaptiveSweepsPerReplica, verdict);
    }

    /**
     * Evalue un triplet {@code thermalizationSweeps / thermalizationRate / measurementRate}
     * IMPOSE, au lieu de le deduire comme {@link #recommend} : utile pour comparer un triplet
     * externe (typiquement le triplet historique fixe, 40 000 / 80 000 / 100) aux memes
     * criteres exacts que ceux utilises pour un triplet recommande. {@code nUndersampled} et
     * {@code nUnderThermalized} sont calcules par {@link #countIssues}, la meme methode que
     * {@link #recommend} utilise en interne : seule la provenance du triplet differe, pas le
     * critere de comptage.
     *
     * <p>{@code tauQuantileSweeps} du resultat vaut tout de meme {@code quantile(tau, cfg.quantile)}
     * (a titre informatif, pour l'affichage), mais ce triplet n'en est PAS derive : contrairement
     * a {@link #recommend}, aucune notion de "temperatures retenues" n'intervient ici, le
     * triplet impose s'applique tel quel a toutes les temperatures.</p>
     *
     * <p><b>Mise en garde indispensable</b> ({@link #DYNAMICS_WARNING}) : les tau de
     * {@code rows} valent pour la dynamique Monte Carlo qui les a mesures. Evaluer un triplet
     * concu pour une AUTRE dynamique (typiquement : le triplet historique, ecrit pour un
     * Metropolis local pur, evalue ici avec des tau mesures sous sur-relaxation) sous-estime
     * fortement le besoin reel de ce triplet. Voir {@link #DYNAMICS_WARNING} pour les chiffres.</p>
     *
     * @throws IllegalArgumentException si {@code rows} est vide/nul, si aucune ligne n'expose de
     *         tau exploitable, si {@code thermalizationSweeps} ou {@code measurementRate} est
     *         &lt;= 0, ou si {@code thermalizationRate <= thermalizationSweeps}.
     */
    public static Recommendation evaluateFixed(List<Row> rows, Config cfg,
            int thermalizationSweeps, int thermalizationRate, int measurementRate) {
        if (rows == null || rows.isEmpty()) throw new IllegalArgumentException("rows vide");
        if (thermalizationSweeps <= 0)
            throw new IllegalArgumentException("thermalizationSweeps doit etre > 0 : " + thermalizationSweeps);
        if (measurementRate <= 0)
            throw new IllegalArgumentException("measurementRate doit etre > 0 : " + measurementRate);
        if (thermalizationRate <= thermalizationSweeps)
            throw new IllegalArgumentException("thermalizationRate (" + thermalizationRate
                    + ") doit etre > thermalizationSweeps (" + thermalizationSweeps + ")");
        Config c = cfg == null ? new Config() : cfg;
        int n = rows.size();

        double[] tauFinite = finiteTau(rows);
        double tauMedian = quantile(tauFinite, 0.5);
        double tauQuantile = quantile(tauFinite, c.quantile); // informatif seulement, voir javadoc
        double tauMax = Arrays.stream(tauFinite).max().orElseThrow();
        Slowest slowest = findSlowest(rows);

        int productionSweeps = thermalizationRate - thermalizationSweeps;

        Counts counts = countIssues(rows, measurementRate, thermalizationSweeps, c);

        long sweepsPerReplica = (long) n * thermalizationRate;
        long legacySweepsPerReplica = (long) n * 80_000L;
        long adaptiveSweepsPerReplica = adaptiveCost(rows);

        String verdict = verdict(n, counts.nUndersampled(), counts.nUnderThermalized(),
                slowest.T(), tauMax, slowest.g(), "triplet impose");

        return new Recommendation(thermalizationSweeps, thermalizationRate, measurementRate,
                productionSweeps, tauMedian, tauQuantile, tauMax, slowest.T(), n,
                counts.nUndersampled(), counts.nUnderThermalized(), sweepsPerReplica,
                legacySweepsPerReplica, adaptiveSweepsPerReplica, verdict);
    }

    /** tauProductionSweeps de {@code rows}, filtre des valeurs non finies. */
    private static double[] finiteTau(List<Row> rows) {
        double[] tauFinite = rows.stream().mapToDouble(Row::tauProductionSweeps)
                .filter(Double::isFinite).toArray();
        if (tauFinite.length == 0)
            throw new IllegalArgumentException("aucune valeur de tauProductionSweeps exploitable dans rows");
        return tauFinite;
    }

    /** Temperature ou tau_prod est maximal, et g de production de cette meme ligne. */
    private record Slowest(double T, double g) { }

    private static Slowest findSlowest(List<Row> rows) {
        double slowestT = Double.NaN;
        double slowestG = Double.NaN;
        double best = Double.NEGATIVE_INFINITY;
        for (Row r : rows) {
            if (Double.isFinite(r.tauProductionSweeps()) && r.tauProductionSweeps() > best) {
                best = r.tauProductionSweeps();
                slowestT = r.T();
                slowestG = r.gProduction();
            }
        }
        return new Slowest(slowestT, slowestG);
    }

    /** max(t0 + K*tau, K_q*tauOverlap) ; le second terme est omis si tauOverlap n'est pas mesure. */
    private static double requiredThermSweeps(Row r, double tau, Config c) {
        double bound = r.t0() + c.tauMultiplier * tau;
        if (Double.isFinite(r.tauOverlapSweeps())) {
            bound = Math.max(bound, c.memoryMultiplier * r.tauOverlapSweeps());
        }
        return bound;
    }

    /** Nombre de temperatures sous-echantillonnees et sous-thermalisees pour un triplet donne
     *  (mesure ou impose, peu importe : critere identique). Toujours sur TOUTES les lignes. */
    private record Counts(int nUndersampled, int nUnderThermalized) { }

    private static Counts countIssues(List<Row> rows, int measurementRate, int thermalizationSweeps, Config c) {
        int nUndersampled = 0;
        int nUnderThermalized = 0;
        for (Row r : rows) {
            double tau = r.tauProductionSweeps();
            if (Double.isFinite(tau) && measurementRate < 2.0 * tau) nUndersampled++;
            double required = requiredThermSweeps(r, Double.isFinite(tau) ? tau : 0.0, c);
            if (thermalizationSweeps < required) nUnderThermalized++;
        }
        return new Counts(nUndersampled, nUnderThermalized);
    }

    /** Cout reellement observe dans le journal : somme(thermSweeps + prodSweeps). */
    private static long adaptiveCost(List<Row> rows) {
        long total = 0L;
        for (Row r : rows) total += r.thermSweeps() + r.prodSweeps();
        return total;
    }

    /**
     * Calcule une {@link Recommendation} par quantile de {@code quantiles}, tout le reste de
     * {@code cfg} etant tenu fixe. C'est ce qui rend testable et affichable (voir
     * {@link #format}) l'arbitrage central de cette classe : deplacer le curseur du quantile de
     * criblage ne change pas seulement {@code measurementRate}, il change tout le triplet et son
     * cout, dans les deux sens (plus rapide en dessous de 0.90, bien plus cher a 1.0).
     *
     * @return un tableau de meme longueur que {@code quantiles}, dans le meme ordre (pas de tri :
     *         l'appelant choisit l'ordre d'affichage).
     * @throws IllegalArgumentException si {@code rows} ou {@code quantiles} est vide ou nul
     *         (propage aussi les exceptions de {@link #recommend}, par exemple un quantile hors
     *         [0,1]).
     */
    public static Recommendation[] scenarios(List<Row> rows, Config cfg, double[] quantiles) {
        if (rows == null || rows.isEmpty()) throw new IllegalArgumentException("rows vide");
        if (quantiles == null || quantiles.length == 0) throw new IllegalArgumentException("quantiles vide");
        Config base = cfg == null ? new Config() : cfg;
        Recommendation[] out = new Recommendation[quantiles.length];
        for (int i = 0; i < quantiles.length; i++) {
            out[i] = recommend(rows, base.copy().quantile(quantiles[i]));
        }
        return out;
    }

    /**
     * Construit le verdict. Point physique central (voir revue) : "sous-echantillonnee"
     * n'est PAS "fausse". Si g est mesure sur la production elle-meme -- ce que fait
     * {@code gProduction} en mode adaptatif -- alors N_eff = n/g et les barres d'erreur qui en
     * decoulent tiennent deja compte de l'auto-correlation residuelle : des mesures correlees
     * ne biaisent pas la moyenne, elles coutent seulement des mesures redondantes (N_eff plus
     * petit que n). Ce qui EST faux, c'est un code qui ne mesure jamais g et suppose des
     * mesures independantes (l'erreur standard historique, {@code sigma/sqrt(n)}) : il
     * sous-estime alors l'erreur d'un facteur sqrt(g), silencieusement. La vraie exigence de
     * decorrelation ne porte donc pas sur {@code nUndersampled} (redondance de mesures
     * scalaires, sans consequence des lors que g est mesure) mais sur {@code nUnderThermalized}
     * (la configuration de spins produite reste partiellement celle heritee de la temperature
     * precedente -- ceci n'a pas de correction a posteriori possible).
     */
    private static String verdict(int n, int nUndersampled, int nUnderTherm,
                                   double slowestT, double tauMax, double slowestG, String scenarioLabel) {
        String biasNote = Double.isFinite(slowestG) && slowestG > 0
                ? String.format(Locale.US,
                        " (a T=%.6f, g=%.1f, un code qui ne mesure pas g sous-estimerait son "
                        + "erreur d'un facteur sqrt(g)=%.1f la)",
                        slowestT, slowestG, Math.sqrt(slowestG))
                : "";
        String principle =
                "Sous-echantillonnee n'est pas synonyme de fausse : des lors que g est mesure sur "
                + "la production elle-meme (mode adaptatif, gProduction), N_eff et les barres "
                + "d'erreur en tiennent deja compte -- on paie une redondance de mesures, pas un "
                + "biais. Ce qui EST faux, c'est un triplet fixe qui ne mesure jamais g et suppose "
                + "des mesures independantes : il sous-estime alors l'erreur d'un facteur sqrt(g)"
                + biasNote + ". La vraie exigence de decorrelation porte sur la configuration de "
                + "spins conservee et sur deux temperatures successives, pas sur deux mesures "
                + "scalaires de production.";

        String status = nUnderTherm == 0
                ? String.format(Locale.US,
                        " Avec ce triplet (%s), aucune des %d temperatures ne reste "
                        + "sous-thermalisee (configuration bien decorrelee de la precedente) ; "
                        + "%d/%d ont des mesures de production redondantes, sans consequence sur "
                        + "l'exactitude. La plus lente est T=%.6f (tau_prod=%.1f sweeps).",
                        scenarioLabel, n, nUndersampled, n, slowestT, tauMax)
                : String.format(Locale.US,
                        " Avec ce triplet (%s), %d/%d temperatures restent "
                        + "sous-thermalisees -- configuration encore partiellement heritee de la "
                        + "precedente, ceci EST un probleme reel -- et %d/%d ont des mesures de "
                        + "production redondantes, sans consequence si g est mesure. La plus lente "
                        + "est T=%.6f (tau_prod=%.1f sweeps).",
                        scenarioLabel, nUnderTherm, n, nUndersampled, n, slowestT, tauMax);

        return principle + status;
    }

    // ------------------------------------------------------------------------------- format

    /**
     * Rapport texte lisible en terminal (aucune ligne ne depasse 100 caracteres) :
     * <ol>
     *   <li>tableau des 10 temperatures les plus lentes par tau_prod decroissant ;</li>
     *   <li>distribution de tau_prod (mediane / quantile retenu / max) ;</li>
     *   <li>bloc pret a coller dans le code historique ;</li>
     *   <li>verdict et comparaison de cout (historique / recommande / adaptatif observe) ;</li>
     *   <li>tableau d'arbitrage : une ligne "historique" (le triplet fixe reel
     *       40 000/80 000/100, evalue par {@link #evaluateFixed} plutot que devine), le triplet
     *       a quantile 0.50, 0.90, 1.00 (voir {@link #scenarios}), et le scenario "adaptatif
     *       observe" (les valeurs reellement utilisees par le controleur, temperature par
     *       temperature -- pas un triplet unique). C'est ce tableau qui repond directement a la
     *       question "quels nombres mettre" en rendant visible, sur la meme colonne de cout, que
     *       couvrir la transition avec un triplet fixe (quantile 1.00) coute davantage que le
     *       triplet historique fixe, alors que le controleur adaptatif, qui n'a besoin d'aucun
     *       triplet, y arrive moins cher que les deux ;</li>
     *   <li>{@link #DYNAMICS_WARNING} : la ligne "historique" ci-dessus est evaluee avec les tau
     *       de CE journal, qui ne sont pas ceux de la dynamique du code historique -- sans cet
     *       avertissement, cette ligne flatterait silencieusement le triplet historique.</li>
     * </ol>
     *
     * <p>Le tableau d'arbitrage est toujours calcule avec une {@link Config} par defaut (seul
     * {@code quantile} varie pour les scenarios a quantile) : il montre l'effet du curseur de
     * criblage isolement, quelle que soit la configuration utilisee pour produire {@code r}
     * lui-meme.</p>
     */
    public static String format(Recommendation r, List<Row> rows) {
        StringBuilder sb = new StringBuilder();
        sb.append("=== SweepCalibrator : triplet de criblage recommande ===\n");

        List<Row> sorted = new ArrayList<>(rows);
        sorted.sort((a, b) -> Double.compare(
                Double.isFinite(b.tauProductionSweeps()) ? b.tauProductionSweeps() : Double.NEGATIVE_INFINITY,
                Double.isFinite(a.tauProductionSweeps()) ? a.tauProductionSweeps() : Double.NEGATIVE_INFINITY));

        sb.append("\nTemperatures les plus lentes (tau_prod decroissant) :\n");
        sb.append(String.format(Locale.US, "%-11s %9s %6s %9s %9s %7s %9s%n",
                "T", "tau_prod", "mEvery", "2*tau", "thermSw", "t0", "neffAch"));
        int shown = Math.min(10, sorted.size());
        for (int i = 0; i < shown; i++) {
            Row row = sorted.get(i);
            sb.append(String.format(Locale.US, "%-11.6f %9.2f %6d %9.2f %9d %7d %9.1f%n",
                    row.T(), row.tauProductionSweeps(), row.measureEvery(),
                    2.0 * row.tauProductionSweeps(), row.thermSweeps(), row.t0(), row.neffAchieved()));
        }

        sb.append(String.format(Locale.US,
                "%ntau_prod (sweeps) : mediane=%.2f  quantile-retenu=%.2f  max=%.2f (T=%.6f)%n",
                r.tauMedianSweeps(), r.tauQuantileSweeps(), r.tauMaxSweeps(), r.slowestT()));

        sb.append("\nBloc pret a coller (remplace les 3 constantes historiques) :\n");
        sb.append("    this.thermalizationSweeps = ").append(r.thermalizationSweeps()).append(";\n");
        sb.append("    this.thermalizationRate = ").append(r.thermalizationRate()).append(";\n");
        sb.append("    this.measurementRate = ").append(r.measurementRate()).append(";\n");

        sb.append('\n');
        appendWrapped(sb, r.verdict(), 100);

        sb.append(String.format(Locale.US,
                "%nCout par replique (sweeps) : historique=%d  recommande=%d  adaptatif-observe=%d%n",
                r.legacySweepsPerReplica(), r.sweepsPerReplica(), r.adaptiveSweepsPerReplica()));
        double accRecommande = (double) r.legacySweepsPerReplica() / Math.max(1, r.sweepsPerReplica());
        double accAdaptatif = (double) r.legacySweepsPerReplica() / Math.max(1, r.adaptiveSweepsPerReplica());
        sb.append(String.format(Locale.US,
                "Acceleration vs historique : recommande=%.1fx  adaptatif-observe=%.1fx%n",
                accRecommande, accAdaptatif));

        appendScenarioTable(sb, rows);

        return sb.toString();
    }

    /**
     * Tableau d'arbitrage : une ligne par quantile de criblage (0.50 / 0.90 / 1.00, Config par
     * defaut sinon) plus une ligne "adaptatif-obs." qui n'est PAS un triplet -- le controleur
     * adaptatif choisit sa propre cadence de mesure et sa propre longueur de thermalisation a
     * chaque temperature (voir {@code AdaptiveThermalization}, {@code MonteCarlo}). Cette ligne
     * affiche donc {@code mediane..max} pour {@code mRate} et {@code thermSw}, "-" pour
     * {@code thermRate} (aucune valeur unique n'existe), et calcule {@code nUndersampled} /
     * {@code nUnderThermalized} <b>ligne a ligne avec les valeurs reellement utilisees</b>
     * (measureEvery et thermSweeps propres a chaque temperature) plutot qu'avec un triplet
     * recommande : c'est la seule maniere de comparer honnetement l'adaptatif a un triplet fixe.
     */
    private static void appendScenarioTable(StringBuilder sb, List<Row> rows) {
        sb.append("\nArbitrage quantile / cout (reponse a \"quels nombres mettre\") :\n");
        sb.append(String.format(Locale.US, "%-16s %9s %11s %9s %13s %9s %9s%n",
                "scenario", "mRate", "thermSw", "thermRate", "cout(sweeps)", "nUnder", "nUnderTh"));

        Config baseCfg = new Config();

        // Point de depart : le triplet historique fixe reel (pas derive, impose), evalue avec
        // les MEMES criteres que les scenarios ci-dessous -- voir DYNAMICS_WARNING plus bas, ce
        // point de comparaison est generereux envers le code historique, pas severe.
        Recommendation historic = evaluateFixed(rows, baseCfg, 40_000, 80_000, 100);
        sb.append(String.format(Locale.US, "%-16s %9d %11d %9d %13d %9s %9s%n",
                "historique", historic.measurementRate(), historic.thermalizationSweeps(),
                historic.thermalizationRate(), historic.sweepsPerReplica(),
                historic.nUndersampled() + "/" + historic.nTemperatures(),
                historic.nUnderThermalized() + "/" + historic.nTemperatures()));

        double[] quantiles = { 0.50, 0.90, 1.00 };
        Recommendation[] scen = scenarios(rows, baseCfg, quantiles);
        for (int i = 0; i < scen.length; i++) {
            Recommendation s = scen[i];
            sb.append(String.format(Locale.US, "%-16s %9d %11d %9d %13d %9s %9s%n",
                    String.format(Locale.US, "quantile %.2f", quantiles[i]),
                    s.measurementRate(), s.thermalizationSweeps(), s.thermalizationRate(),
                    s.sweepsPerReplica(), s.nUndersampled() + "/" + s.nTemperatures(),
                    s.nUnderThermalized() + "/" + s.nTemperatures()));
        }

        // "adaptatif-obs." n'est PAS un triplet : chaque temperature a sa propre cadence de
        // mesure et sa propre longueur de thermalisation (voir MonteCarlo /
        // AdaptiveThermalization). nUndersampled / nUnderThermalized sont donc evalues ligne a
        // ligne avec les valeurs REELLEMENT utilisees a chaque temperature, pas avec un triplet
        // recommande : c'est la seule comparaison honnete avec les lignes a triplet fixe.
        double[] mEvery = rows.stream().mapToDouble(rr -> (double) rr.measureEvery()).toArray();
        double[] thermS = rows.stream().mapToDouble(rr -> (double) rr.thermSweeps()).toArray();
        int adaptiveUndersampled = 0;
        int adaptiveUnderTherm = 0;
        for (Row row : rows) {
            double tau = row.tauProductionSweeps();
            if (Double.isFinite(tau) && row.measureEvery() < 2.0 * tau) adaptiveUndersampled++;
            double required = requiredThermSweeps(row, Double.isFinite(tau) ? tau : 0.0, baseCfg);
            if (row.thermSweeps() < required) adaptiveUnderTherm++;
        }
        sb.append(String.format(Locale.US, "%-16s %9s %11s %9s %13d %9s %9s%n",
                "adaptatif-obs.", medianMaxLabel(mEvery), medianMaxLabel(thermS), "-",
                adaptiveCost(rows), adaptiveUndersampled + "/" + rows.size(),
                adaptiveUnderTherm + "/" + rows.size()));

        sb.append('\n');
        appendWrapped(sb, DYNAMICS_WARNING, 100);
    }

    /** {@code "mediane..max"} (arrondis a l'entier : ce sont des sweeps / cadences entieres). */
    private static String medianMaxLabel(double[] values) {
        if (values.length == 0) return "-";
        long median = Math.round(quantile(values, 0.5));
        long max = Math.round(Arrays.stream(values).max().orElse(Double.NaN));
        return median + ".." + max;
    }

    /**
     * Decoupe {@code text} en lignes d'au plus {@code width} caracteres, sur des frontieres de
     * mots (un mot lui-meme plus long que {@code width} est ecrit tel quel, sans coupure au
     * milieu). Utilise pour que {@link #format} respecte la limite de largeur du rapport meme
     * quand {@link Recommendation#verdict()} est une phrase longue.
     */
    private static void appendWrapped(StringBuilder sb, String text, int width) {
        StringBuilder line = new StringBuilder();
        for (String word : text.split(" ")) {
            if (line.length() > 0 && line.length() + 1 + word.length() > width) {
                sb.append(line).append('\n');
                line.setLength(0);
            }
            if (line.length() > 0) line.append(' ');
            line.append(word);
        }
        if (line.length() > 0) sb.append(line).append('\n');
    }

    // --------------------------------------------------------------------------------- main

    private static final String USAGE =
            "Usage: SweepCalibrator <equilibration.csv> [decorrelation.csv]\n"
            + "         [--target-neff N] [--quantile q] [--tau-multiplier K]\n"
            + "         [--memory-multiplier K] [--measurements-per-tau M]\n"
            + "Lit un ou deux journaux et affiche le triplet de sweeps recommande pour un\n"
            + "criblage rapide et valide (voir le javadoc de SweepCalibrator).";

    public static void main(String[] args) {
        if (args.length < 1) {
            System.err.println(USAGE);
            System.exit(2);
            return;
        }
        Path eqPath = Paths.get(args[0]);
        Path decorPath = null;
        Config cfg = new Config();
        int i = 1;
        if (i < args.length && !args[i].startsWith("--")) {
            decorPath = Paths.get(args[i]);
            i++;
        }
        try {
            for (; i < args.length; i++) {
                String a = args[i];
                switch (a) {
                    case "--target-neff" -> cfg.targetNeff(Double.parseDouble(needValue(args, ++i, a)));
                    case "--quantile" -> cfg.quantile(Double.parseDouble(needValue(args, ++i, a)));
                    case "--tau-multiplier" -> cfg.tauMultiplier(Double.parseDouble(needValue(args, ++i, a)));
                    case "--memory-multiplier" -> cfg.memoryMultiplier(Double.parseDouble(needValue(args, ++i, a)));
                    case "--measurements-per-tau" -> cfg.measurementsPerTau(Double.parseDouble(needValue(args, ++i, a)));
                    default -> throw new IllegalArgumentException("option inconnue : " + a);
                }
            }
        } catch (RuntimeException e) {
            System.err.println("argument invalide : " + e.getMessage());
            System.err.println(USAGE);
            System.exit(2);
            return;
        }

        List<Row> rows;
        try {
            List<Row> eq = readEquilibrationCsv(eqPath);
            rows = decorPath == null ? eq : merge(eq, readDecorrelationCsv(decorPath));
        } catch (RuntimeException e) {
            System.err.println("erreur de lecture : " + e.getMessage());
            System.exit(2);
            return;
        }

        Recommendation rec = recommend(rows, cfg);
        System.out.println(format(rec, rows));
    }

    private static String needValue(String[] args, int idx, String opt) {
        if (idx >= args.length) throw new IllegalArgumentException("valeur manquante pour " + opt);
        return args[idx];
    }
}
