/**
 * Tuple complet du diagnostic de decorrelation et de stagnation pour <b>une</b> temperature,
 * produit par le detecteur de stagnation (voir {@code StagnationDetector}, {@code SpinOverlap}).
 *
 * <p>Il existe pour la meme raison que {@link TemperaturePoint} : sans lui, une vingtaine de
 * nombres circuleraient separement entre le detecteur, {@link StagnationLog#write(StagnationPoint)}
 * et l'appelant, avec le risque permanent d'intervertir deux arguments {@code double} sans que le
 * compilateur ne le voie. Le record nomme chaque champ une fois pour toutes et {@link #extras()}
 * fige l'ordre des colonnes en un seul endroit.</p>
 *
 * <h2>Deux decorrelations distinctes</h2>
 * <p>Le probleme a resoudre est celui du criblage rapide de modeles : on ne cherche pas une
 * configuration parfaitement equilibree, on cherche une configuration <i>valide</i>, obtenue
 * vite. Il faut pour cela s'assurer de deux choses independantes :</p>
 * <ol>
 *   <li><b>decorrelation entre deux temperatures successives</b> ({@code temperaturesDecorrelated}) :
 *       la configuration de spins retenue pour cette temperature ne doit plus etre, pour
 *       l'essentiel, celle heritee de la temperature precedente. On la mesure via le
 *       recouvrement de configuration q(t) = (1/N) sum_i S_i(t).S_i(ref), ou ref est l'etat
 *       herite au debut de la temperature ;</li>
 *   <li><b>decorrelation entre deux mesures successives de production</b>
 *       ({@code measurementsDecorrelated}) : le probleme classique de l'auto-correlation
 *       temporelle d'une serie MCMC, mesure par le temps integre tau_prod de la serie de
 *       production elle-meme.</li>
 * </ol>
 *
 * <h2>Le plateau de q n'est pas zero</h2>
 * <p>Dans une phase ordonnee (ou proche de la transition, ou les domaines sont larges), q(t) ne
 * decroit pas vers 0 : il decroit vers un plateau q_inf &gt; 0, parametre d'ordre au sens
 * d'Edwards-Anderson (cf. Ogielski, PRB 32, 7384 (1985)). Ce plateau existe parce qu'une fraction
 * des spins reste alignee avec l'etat de reference simplement parce que les deux etats
 * appartiennent au meme puits (meme domaine, meme brisure de symetrie locale), pas parce que la
 * chaine n'a pas bouge. Le critere de decorrelation entre temperatures n'est donc pas
 * <i>q -&gt; 0</i> mais <i>q a fini de decroitre vers son plateau</i> : c'est {@code tauOverlapSweeps}
 * qui mesure le temps de cette decroissance, et {@code memoryRatio} qui dit si on lui a laisse
 * assez de sweeps de thermalisation pour l'atteindre. Si en plus {@code qInf} est proche de 1 (chaine
 * gelee : le plateau est atteint mais il vaut presque 1), aucun sweep supplementaire ne
 * decorrelera la chaine de l'etat herite, et {@code stopReason} vaut {@code "STAGNANT_FROZEN"}.</p>
 *
 * @param mcIdx                     index 1-based de la temperature.
 * @param T                         temperature.
 * @param thermSweeps               sweeps de thermalisation reellement effectues.
 * @param prodSweeps                sweeps de production reellement effectues.
 * @param qInf                      plateau du recouvrement de configuration avec l'etat herite de
 *                                  la temperature precedente ; parametre d'ordre au sens
 *                                  Edwards-Anderson, ne tend pas vers 0 en phase ordonnee.
 * @param tauOverlapSweeps          temps de decroissance de ce recouvrement vers son plateau, en
 *                                  sweeps ({@code NaN} si non mesure, par exemple premiere
 *                                  temperature du recuit ou plateau non atteint).
 * @param memoryRatio               thermSweeps / tauOverlapSweeps : combien de temps de memoire de
 *                                  la temperature precedente ont ete effectivement jetes par la
 *                                  thermalisation ({@code NaN} si tauOverlapSweeps l'est).
 * @param temperaturesDecorrelated  vrai si memoryRatio est au moins egal au seuil demande
 *                                  <b>et</b> que le plateau de q a ete atteint.
 * @param driftZ                    z (Geweke-like) de la derive residuelle de l'energie totale en
 *                                  fin de thermalisation : grand |driftZ| = la moyenne bouge encore
 *                                  plus que ne le permet son bruit ({@code NaN} si non mesure).
 * @param driftPerSweep             pente residuelle de l'energie totale, par sweep, a la fin de la
 *                                  thermalisation ({@code NaN} si non mesuree).
 * @param stopReason                raison d'arret de la thermalisation/production, p.ex.
 *                                  {@code "EQUILIBRATED"}, {@code "STAGNANT_PLATEAU"},
 *                                  {@code "STAGNANT_FROZEN"}, {@code "MAX_SWEEPS"}.
 * @param converged                 verdict du controleur d'equilibration
 *                                  ({@code AdaptiveThermalization.converged()}) pour cette
 *                                  temperature : distingue, dans le CSV, un arret parce que
 *                                  l'equilibration a ete detectee d'un arret par stagnation alors
 *                                  que le controleur n'avait pas conclu.
 * @param measureEvery              cadence de mesure retenue en production, en sweeps.
 * @param tauProductionSweeps       tau_int mesure sur la serie de production elle-meme, en sweeps
 *                                  ({@code NaN} si non mesure).
 * @param measurementsPerTau        measureEvery / (2 tauProductionSweeps) : >= 1 signifie que deux
 *                                  mesures successives sont decorrelees ({@code NaN} si
 *                                  tauProductionSweeps l'est).
 * @param measurementsDecorrelated  vrai si measurementsPerTau &gt;= 1.
 * @param neffAchieved              nombre effectif de mesures independantes reellement obtenu sur
 *                                  la production.
 * @param neffSaturated             vrai si les extensions successives de la production ne
 *                                  gagnaient plus de N_eff (rendements decroissants constates).
 * @param seconds                   temps mural passe sur cette temperature.
 * @param hx                        composante x du champ magnetique (base de Kitaev) de la
 *                                  replique ; toujours disponible (valide par runParallel).
 * @param hy                        composante y du champ magnetique (base de Kitaev).
 * @param hz                        composante z du champ magnetique (base de Kitaev).
 */
public record StagnationPoint(
        int mcIdx,
        double T,
        long thermSweeps,
        long prodSweeps,
        double qInf,
        double tauOverlapSweeps,
        double memoryRatio,
        boolean temperaturesDecorrelated,
        double driftZ,
        double driftPerSweep,
        String stopReason,
        boolean converged,
        int measureEvery,
        double tauProductionSweeps,
        double measurementsPerTau,
        boolean measurementsDecorrelated,
        double neffAchieved,
        boolean neffSaturated,
        double seconds,
        double hx,
        double hy,
        double hz) {

    /** Nombre de colonnes renvoyees par {@link #extras()}. */
    public static final int N_EXTRAS = 21;

    /**
     * Validation compacte : ne verifie que ce qui rendrait le point physiquement absurde ou le
     * CSV illisible, pas la coherence fine entre champs (par exemple entre {@code memoryRatio} et
     * {@code temperaturesDecorrelated}), qui releve du detecteur qui construit le point, pas du
     * point lui-meme.
     *
     * <p>{@code NaN} est autorise sur les champs explicitement documentes comme pouvant ne pas
     * etre mesures : {@code qInf}, {@code tauOverlapSweeps}, {@code memoryRatio},
     * {@code driftZ}, {@code driftPerSweep}, {@code tauProductionSweeps},
     * {@code measurementsPerTau}.</p>
     *
     * @throws IllegalArgumentException si T &lt;= 0, si thermSweeps &lt; 0, si prodSweeps &lt; 0,
     *         si measureEvery &lt; 1, ou si stopReason est {@code null} ou vide.
     */
    public StagnationPoint {
        if (T <= 0) throw new IllegalArgumentException("T doit etre > 0 : " + T);
        if (thermSweeps < 0) throw new IllegalArgumentException("thermSweeps doit etre >= 0 : " + thermSweeps);
        if (prodSweeps < 0) throw new IllegalArgumentException("prodSweeps doit etre >= 0 : " + prodSweeps);
        if (measureEvery < 1) throw new IllegalArgumentException("measureEvery doit etre >= 1 : " + measureEvery);
        if (stopReason == null || stopReason.isEmpty())
            throw new IllegalArgumentException("stopReason ne doit pas etre null ou vide : " + stopReason);
    }

    /**
     * Les colonnes numeriques de diagnostic, dans l'ordre documente en tete de classe et repris
     * par {@link StagnationLog#COLUMNS} (les booleens sont ecrits en 0/1, {@code stopReason} en
     * est exclu car il est textuel).
     *
     * <p>Un nouveau tableau a chaque appel : l'appelant le copie et peut le modifier sans
     * consequence.</p>
     */
    public double[] extras() {
        return new double[] {
            mcIdx,                                     // mcIdx
            T,                                          // T
            thermSweeps,                                // thermSweeps
            prodSweeps,                                 // prodSweeps
            qInf,                                       // qInf
            tauOverlapSweeps,                           // tauOverlapSweeps
            memoryRatio,                                // memoryRatio
            temperaturesDecorrelated ? 1.d : 0.d,       // temperaturesDecorrelated
            driftZ,                                     // driftZ
            driftPerSweep,                              // driftPerSweep
            converged ? 1.d : 0.d,                      // converged
            measureEvery,                               // measureEvery
            tauProductionSweeps,                        // tauProductionSweeps
            measurementsPerTau,                         // measurementsPerTau
            measurementsDecorrelated ? 1.d : 0.d,       // measurementsDecorrelated
            neffAchieved,                                // neffAchieved
            neffSaturated ? 1.d : 0.d,                  // neffSaturated
            seconds,                                     // seconds
            hx,                                          // hx
            hy,                                          // hy
            hz                                           // hz
        };
    }
}
