/**
 * Tuple complet du diagnostic d'<b>une</b> temperature du recuit, produit par
 * {@link MonteCarlo#runParallel}.
 *
 * <p>Il existe pour une seule raison : avant, les memes quinze a vingt nombres circulaient
 * separement entre {@code runParallel}, {@link EquilibrationLog#write(TemperaturePoint)} et
 * {@code observablesData}. Une methode a dix-huit parametres {@code double} consecutifs est une
 * invitation permanente a l'interversion silencieuse de deux arguments, et rien dans le
 * compilateur ne l'aurait signalee. Le record nomme chaque champ une fois pour toutes et
 * {@link #extras()} fige l'ordre des colonnes de diagnostic en un seul endroit.</p>
 *
 * <h2>Ordre des colonnes</h2>
 * <p>{@link #extras()} renvoie exactement les colonnes {@code MonteCarlo.COL_THERM} ..
 * {@code MonteCarlo.COL_G_PROD} de {@code observablesData}, dans cet ordre :</p>
 * <ol start="16">
 *   <li>{@code thermSweeps} — sweeps de thermalisation reellement consommes ;</li>
 *   <li>{@code t0Sweep} — sweep a partir duquel la serie est jugee equilibree (Chodera) ;</li>
 *   <li>{@code tauSweeps} — tau_int de la serie de <b>thermalisation</b>, en sweeps ;</li>
 *   <li>{@code g} — inefficacite statistique 1 + 2 tau de la serie de thermalisation ;</li>
 *   <li>{@code neff} — N_eff de la queue de la serie de thermalisation ;</li>
 *   <li>{@code prodSweeps} — sweeps de production <b>totaux</b> (bloc planifie + extensions) ;</li>
 *   <li>{@code converged} — 1 si l'equilibration a ete detectee, 0 si coupee par {@code maxSweeps} ;</li>
 *   <li>{@code seconds} — temps mural passe sur cette temperature ;</li>
 *   <li>{@code measureEvery} — cadence de mesure de production effectivement utilisee ;</li>
 *   <li>{@code nMeasurements} — nombre de mesures de production accumulees ;</li>
 *   <li>{@code sigmaEndOfThermalization} — sigma a la fin de la thermalisation ;</li>
 *   <li>{@code neffAchieved} — N_eff <b>reellement obtenu</b> sur la serie de production ;</li>
 *   <li>{@code gProduction} — g de la serie de production (Wolff), celui qui sert au N_eff ci-dessus.</li>
 * </ol>
 *
 * <h2>Deux tau, deux g, deux N_eff</h2>
 * <p>{@code tauSweeps}, {@code g} et {@code neff} viennent du controleur d'equilibration, donc
 * de la serie <b>de thermalisation</b> : ce sont les nombres qui ont servi a <i>dimensionner</i>
 * la production. {@code gProduction} et {@code neffAchieved} sont mesures <b>a posteriori sur la
 * production elle-meme</b> : ce sont ceux qui disent ce qui a reellement ete livre. Les comparer
 * est le seul moyen de savoir si la prevision etait bonne ; on ne les confond pas, et c'est
 * pourquoi les deux jeux sont stockes.</p>
 *
 * @param mcIdx                    index 1-based de la temperature, comme {@code observablesData[i][0]}.
 * @param T                        temperature.
 * @param thermSweeps              sweeps consommes par la thermalisation adaptative.
 * @param t0Sweep                  debut de la zone equilibree (Chodera), en sweeps.
 * @param tauSweeps                tau_int de la serie de thermalisation, en sweeps.
 * @param g                        g = 1 + 2 tau de la serie de thermalisation.
 * @param neff                     N_eff de la queue de la serie de thermalisation.
 * @param zE                       z de Geweke du canal energie.
 * @param zM                       z de Geweke du canal |M|.
 * @param prodSweeps               sweeps de production totaux (planifie + extensions).
 * @param converged                false si la thermalisation a ete coupee par {@code maxSweeps}.
 * @param acceptance               taux d'acceptation moyen sur toute la production.
 * @param sigma                    sigma a la fin de la production.
 * @param sigmaEndOfThermalization sigma a la fin de la thermalisation ; egal a {@code sigma} bit
 *                                 a bit lorsque {@code freezeSigmaInProduction} est vrai.
 * @param ePerSite                 &lt;E&gt;/N de production.
 * @param eErr                     erreur standard sur &lt;E&gt;/N, corrigee de l'auto-correlation.
 * @param mNorm                    norme de l'aimantation totale moyenne.
 * @param seconds                  temps mural passe sur cette temperature.
 * @param measureEvery             cadence de mesure de production, en sweeps.
 * @param nMeasurements            nombre de mesures de production accumulees.
 * @param neffAchieved             N_eff reellement obtenu = nMeasurements / gProduction.
 * @param gProduction              g de la serie de production (methode Gamma de Wolff).
 * @param hx                       composante x du champ applique (base Kitaev).
 * @param hy                       composante y du champ applique.
 * @param hz                       composante z du champ applique.
 * @param stdE                     3 sigma de la distribution d'energie (largeur, pas une erreur).
 * @param varE                     variance de population de l'energie totale — matiere premiere de
 *                                 C_v/N = varE / (N T^2).
 * @param mx                       &lt;M_x&gt; de production.
 * @param stdMx                    3 sigma de la distribution de M_x.
 * @param my                       &lt;M_y&gt; de production.
 * @param stdMy                    3 sigma de la distribution de M_y.
 * @param mz                       &lt;M_z&gt; de production.
 * @param stdMz                    3 sigma de la distribution de M_z.
 */
public record TemperaturePoint(
        int mcIdx,
        double T,
        long thermSweeps,
        long t0Sweep,
        double tauSweeps,
        double g,
        double neff,
        double zE,
        double zM,
        int prodSweeps,
        boolean converged,
        double acceptance,
        double sigma,
        double sigmaEndOfThermalization,
        double ePerSite,
        double eErr,
        double mNorm,
        double seconds,
        int measureEvery,
        int nMeasurements,
        double neffAchieved,
        double gProduction,
        double hx,
        double hy,
        double hz,
        double stdE,
        double varE,
        double mx,
        double stdMx,
        double my,
        double stdMy,
        double mz,
        double stdMz) {

    /** Nombre de colonnes renvoyees par {@link #extras()}. */
    public static final int N_EXTRAS = 13;

    /**
     * Les colonnes de diagnostic de {@code observablesData}, dans l'ordre documente en tete de
     * classe (indices {@code MonteCarlo.COL_THERM} .. {@code MonteCarlo.COL_G_PROD}).
     *
     * <p>Un nouveau tableau a chaque appel : l'appelant le copie dans {@code combinedData} et
     * peut le modifier sans consequence.</p>
     */
    public double[] extras() {
        return new double[] {
            thermSweeps,               // COL_THERM
            t0Sweep,                   // COL_T0
            tauSweeps,                 // COL_TAU
            g,                         // COL_G
            neff,                      // COL_NEFF
            prodSweeps,                // COL_PROD
            converged ? 1.d : 0.d,     // COL_CONVERGED
            seconds,                   // COL_SECONDS
            measureEvery,              // COL_MEASURE_EVERY
            nMeasurements,             // COL_NMEAS
            sigmaEndOfThermalization,  // COL_SIGMA_THERM
            neffAchieved,              // COL_NEFF_ACHIEVED
            gProduction                // COL_G_PROD
        };
    }
}
