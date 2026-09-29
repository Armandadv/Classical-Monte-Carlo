/**
 * Un point de la courbe de susceptibilite AC : chi'(T, omega) et chi''(T, omega), moyennes sur
 * {@code nReplicas} repliques independantes, produites par {@link ACSusceptibility#aggregate}.
 *
 * <h2>Protocole physique</h2>
 * <p>On mesure la reponse lineaire d'un modele de spins honeycomb a un champ oscillant
 * {@code H(t) = hStatic + h0 . sin(omega t) . acDrive}, {@code t} etant compte en <b>sweeps</b>
 * du moteur Metropolis damier ({@link CheckerboardMetropolis}) : c'est cette dynamique
 * stochastique a bilan detaille (relaxation thermiquement activee), et non une dynamique ODE
 * conservative, qui joue le role d'horloge cinetique. Voir {@link ACSusceptibility} pour le
 * protocole complet (equilibration puis mesure a sigma gele).</p>
 *
 * <h2>Convention de normalisation</h2>
 * <p>La demodulation brute du signal de simulation donne, sur une fenetre de mesure de
 * {@code nMeas} sweeps (un nombre entier de periodes) :</p>
 * <pre>
 *   S'  =  (1/nMeas) . somme_t [ M_par(t) . sin(omega t) ]
 *   S'' = -(1/nMeas) . somme_t [ M_par(t) . cos(omega t) ]
 * </pre>
 * <p>ou {@code M_par(t)} est l'aimantation totale projetee sur la direction du champ oscillant.
 * Au premier ordre en {@code h0}, {@code M_par(t) ~ h0 . [chi' sin(omega t) - chi'' cos(omega t)]
 * . N} (N = nombre de sites), de sorte que {@code S'} et {@code S''} valent toutes deux
 * {@code +chi' . h0 . N / 2} et {@code +chi'' . h0 . N / 2} a un facteur pres : le signe moins
 * introduit dans la definition de {@code S''} ci-dessus compense exactement celui du terme en
 * {@code cos(omega t)} du modele, de sorte que {@code S''} (et donc {@link #chiSecond}) porte
 * directement le signe de {@code chi''} et non son oppose. Une demodulation "naive"
 * {@code somme_t [M_par(t) . cos(omega t)]}, sans ce signe moins, donnerait -chi''.h0.N/2, soit
 * l'oppose de la convention retenue ici.</p>
 * <p><b>Convention de signe :</b> {@link #chiSecond} &gt;= 0 pour une reponse retardee/dissipative,
 * comme sur une mesure experimentale de chi''(T) (les pics de dissipation sont positifs). Voir
 * {@link ACSusceptibility#run} pour l'endroit exact ou ce signe est applique (au moment de la
 * normalisation, sur l'accumulateur brut de la demodulation en cosinus).</p>
 * <p>{@link #chiPrime} et {@link #chiSecond} sont ces quantites <b>divisees par h0</b> (pour
 * obtenir une susceptibilite independante de l'amplitude d'excitation, cf. la section linearite
 * du plan de test) <b>et par le nombre de sites</b> (pour obtenir une grandeur intensive,
 * comparable entre tailles de reseau). Ce facteur constant ne deplace <b>aucun pic</b> de
 * chi'(T) ou chi''(T) : seule l'echelle verticale de la courbe en depend. Voir
 * {@link ACSusceptibility#run} pour la formule exacte de normalisation appliquee.</p>
 *
 * <h2>Decalage champ / mesure</h2>
 * <p>Le champ instantane est ecrit <b>avant</b> le sweep qui l'applique, et l'aimantation est
 * lue <b>apres</b> ce meme sweep : il existe donc un decalage systematique d'un demi-sweep entre
 * la phase du champ et celle de la reponse mesuree. Cela melange une fraction
 * {@code ~sin(pi / periodSweeps)} de chi' dans chi'' (et reciproquement, au meme ordre) ; pour
 * {@code periodSweeps >= 64} ce melange est inferieur a ~5 % et decroit comme
 * {@code 1/periodSweeps}, donc negligeable devant les barres d'erreur statistiques usuelles.</p>
 *
 * @param T                temperature du point.
 * @param omega             pulsation angulaire {@code 2 pi / periodSweeps}, en radians/sweep.
 * @param periodSweeps      periode d'excitation, en sweeps (nombre entier).
 * @param h0                amplitude du champ oscillant utilisee pour produire ce point.
 * @param hStatic           champ statique (metadonnee descriptive du run ; simple scalaire de
 *                          reference transmis tel quel par l'appelant a {@link
 *                          ACSusceptibility#aggregate}, comme le champ de {@code
 *                          EquilibrationLog.defaultPath}).
 * @param nReplicas         nombre de repliques independantes moyennees.
 * @param chiPrime          partie reactive (en phase) de la susceptibilite, normalisee (voir
 *                          ci-dessus), moyennee sur les repliques.
 * @param chiPrimeErr        erreur standard sur {@link #chiPrime} (ecart-type echantillon divise
 *                          par racine(nReplicas) ; 0 si {@code nReplicas == 1}).
 * @param chiSecond         partie dissipative (en quadrature) de la susceptibilite, normalisee,
 *                          moyennee sur les repliques ; convention chiSecond &gt;= 0 pour une
 *                          reponse retardee/dissipative (voir "Convention de normalisation"
 *                          ci-dessus).
 * @param chiSecondErr       erreur standard sur {@link #chiSecond}, meme convention que {@link
 *                          #chiPrimeErr}.
 * @param meanThermSweeps   nombre moyen de sweeps de thermalisation consommes avant la mesure.
 * @param meanAcceptance    taux d'acceptation Metropolis moyen pendant la fenetre de mesure.
 * @param meanSeconds       temps mural moyen (equilibration + mesure) par replique, en secondes.
 */
public record ACPoint(
        double T,
        double omega,
        int periodSweeps,
        double h0,
        double hStatic,
        int nReplicas,
        double chiPrime,
        double chiPrimeErr,
        double chiSecond,
        double chiSecondErr,
        double meanThermSweeps,
        double meanAcceptance,
        double meanSeconds) {
}
