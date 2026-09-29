import java.util.ArrayList;
import java.util.List;

/**
 * Mesure la susceptibilite AC, chi'(omega) et chi''(omega), d'un modele de spins honeycomb par
 * demodulation sous champ oscillant.
 *
 * <h2>Principe physique</h2>
 * <p>On applique un champ {@code H(t) = hStatic + h0 . sin(omega t) . acDrive} et on demodule la
 * projection de l'aimantation sur la direction du champ oscillant. Le temps {@code t} est compte
 * en <b>sweeps</b> du moteur Metropolis damier ({@link CheckerboardMetropolis}) : c'est cette
 * dynamique stochastique <b>a bilan detaille</b> (chaque sweep est une tentative de relaxation
 * thermiquement activee) qui joue le role d'horloge cinetique, et non une dynamique ODE
 * conservative ({@code Spin.solveODE1D}) — celle-ci integre les equations de precession/relaxation
 * deterministes du systeme et n'a pas de notion propre de bruit thermique module par un champ.
 * La pulsation {@code omega = 2 pi / periodSweeps} est donc rigoureusement definie des lors que
 * {@code periodSweeps} est un nombre <b>entier</b> de sweeps, et une fenetre de mesure de longueur
 * un nombre entier de periodes assure l'orthogonalite exacte (sur cette fenetre discrete) des
 * projections sur {@code sin} et {@code cos}.</p>
 * <p><b>Convention de signe :</b> {@link ACPoint#chiSecond} (partie dissipative, en quadrature)
 * est retourne positif pour une reponse retardee/dissipative, comme sur une mesure experimentale
 * de chi''(T) — voir la section "Convention de normalisation" de {@link ACPoint} pour la
 * derivation complete.</p>
 *
 * <h2>Protocole, par temperature</h2>
 * <ol>
 *   <li><b>Equilibration</b> sous le seul champ statique {@code hStatic} : sur-relaxation active
 *       ({@code cfg.overRelaxationPerSweep} passages entre deux sweeps Metropolis), sigma
 *       adaptatif ({@link CheckerboardMetropolis#nextSigma}), duree pilotee par
 *       {@link AdaptiveThermalization} — exactement le motif de
 *       {@code MonteCarlo.runParallel} (boucle de thermalisation, lignes ~955-1020).</li>
 *   <li><b>Mesure</b> sous le champ complet {@code H(t)} : sur-relaxation <b>coupee</b>, sigma
 *       <b>gele</b> a sa valeur de fin d'equilibration. Deux raisons imperatives :
 *       <ul>
 *         <li>l'adaptation de sigma ({@code nextSigma}) viole le bilan detaille — elle change la
 *             loi de transition en cours de mesure, ce qui n'a pas de sens pour une dynamique
 *             cinetique qu'on cherche a interpreter physiquement ;</li>
 *         <li>le taux de rejet Metropolis est module par le champ oscillant lui-meme : laisser
 *             sigma le suivre injecterait artificiellement une reponse a la frequence de mesure,
 *             independante de toute physique du systeme.</li>
 *       </ul>
 *       La sur-relaxation est coupee pour la meme famille de raisons : c'est un mouvement
 *       deterministe et conservatif (reflexion du spin autour de son champ local), donc etranger
 *       a la dynamique cinetique stochastique que l'on utilise comme horloge — l'entrelacer avec
 *       les sweeps Metropolis changerait le temps physique associe a chaque sweep de mesure.</li>
 * </ol>
 *
 * <h2>Balayage en (periode, temperature)</h2>
 * <p>{@link #run} boucle sur les periodes {@code cfg.periodsSweeps} <b>a l'exterieur</b>, et sur
 * une grille de temperatures decroissante <b>a l'interieur</b> : pour chaque periode (chaque
 * omega), on refait une echelle de refroidissement complete depuis {@code cfg.tMax}. Ceci est
 * deliberement redondant en temps de calcul mais indispensable a la correction physique : passer
 * d'une temperature a la suivante <b>sans re-thermaliser depuis zero</b> (field-cooling) fait
 * reposer l'equilibration de chaque temperature sur l'historique de la precedente — un raccourci
 * standard pour balayer un recuit — mais cet historique est celui d'une <i>autre echelle de
 * temps</i> (une autre periode d'excitation) si l'on ne repart pas de {@code tMax} a chaque
 * changement de {@code periodSweeps}. Chaque nouvelle periode efface donc completement l'etat
 * herite en repartant du sommet de la grille de temperature.</p>
 *
 * <h2>Decalage champ / mesure</h2>
 * <p>Le champ instantane {@code H(t)} est ecrit <b>avant</b> le sweep Metropolis qui l'applique
 * (via {@code Lattice.setInteractionField}), et l'aimantation est lue <b>apres</b> ce meme sweep.
 * Il y a donc un decalage systematique d'un demi-sweep entre la phase du champ et celle de la
 * reponse mesuree, qui melange une fraction {@code ~sin(pi / periodSweeps)} de chi' dans chi''
 * (et reciproquement). Ce melange decroit comme {@code 1/periodSweeps} et devient negligeable
 * (moins de ~5 %) des que {@code periodSweeps >= 64} — voir {@link ACPoint} pour le detail du
 * calcul.</p>
 *
 * <h2>Independances</h2>
 * <p>Cette classe ne depend d'aucune classe d'orchestration de recuit
 * ({@code MonteCarlo}, {@code EquilibrationLog}, {@code AppParallel}) : elle ne manipule que
 * {@link Lattice}, {@link CheckerboardMetropolis} et {@link AdaptiveThermalization}.</p>
 */
public final class ACSusceptibility {

    private ACSusceptibility() {}

    /**
     * Parametres du protocole. POJO mutable ; {@link #run} en travaille sur une copie defensive
     * ({@link #copy()}), les modifications apres l'appel n'ont donc aucun effet retroactif.
     */
    public static final class Config {
        /** Amplitude du champ oscillant. */
        public double h0 = 0.02;
        /** Periodes d'excitation explorees, en sweeps (chacune est une echelle de refroidissement complete). */
        public int[] periodsSweeps = { 64, 128, 256 };
        /** Nombre de periodes de regime transitoire ecartees avant d'accumuler la demodulation. */
        public int transientPeriods = 5;
        /** Nombre de periodes de mesure accumulees dans la demodulation. */
        public int measurePeriods = 30;
        /** Temperature de depart de chaque echelle de refroidissement (la plus haute). */
        public double tMax = 2.0;
        /** Temperature d'arret (la grille descend jusqu'a >= tMin - 1e-9). */
        public double tMin = 0.05;
        /** Pas de la grille lineaire de temperature. */
        public double tStep = 0.05;
        /** Sur-relaxations entre deux sweeps Metropolis, EQUILIBRATION SEULEMENT (coupee en mesure). */
        public int overRelaxationPerSweep = 2;
        /** Configuration du controleur d'equilibration adaptative. */
        public AdaptiveThermalization.Config thermalization = new AdaptiveThermalization.Config();
        /** Impression d'une ligne de suivi par (periode, temperature) sur la sortie standard. */
        public boolean verbose = false;

        /** Copie profonde (y compris la configuration de thermalisation imbriquee). */
        public Config copy() {
            Config c = new Config();
            c.h0 = h0;
            c.periodsSweeps = (periodsSweeps == null) ? null : periodsSweeps.clone();
            c.transientPeriods = transientPeriods;
            c.measurePeriods = measurePeriods;
            c.tMax = tMax;
            c.tMin = tMin;
            c.tStep = tStep;
            c.overRelaxationPerSweep = overRelaxationPerSweep;
            c.thermalization = (thermalization == null) ? new AdaptiveThermalization.Config()
                                                         : thermalization.copy();
            c.verbose = verbose;
            return c;
        }
    }

    /**
     * Un point de mesure brut, pour une replique et un couple (temperature, periode) donnes.
     * {@link #chiPrime} et {@link #chiSecond} portent deja la normalisation documentee dans
     * {@link ACPoint} (division par {@code h0} et par le nombre de sites) ; {@link #aggregate}
     * n'a plus qu'a moyenner sur les repliques.
     *
     * @param T              temperature de ce point.
     * @param periodSweeps   periode d'excitation, en sweeps.
     * @param chiPrime       partie reactive normalisee.
     * @param chiSecond      partie dissipative normalisee ; convention chiSecond &gt;= 0 pour une
     *                       reponse retardee/dissipative (voir {@link ACPoint}).
     * @param thermSweeps    sweeps de thermalisation consommes avant la mesure.
     * @param meanAcceptance taux d'acceptation Metropolis moyen sur la fenetre de mesure.
     * @param seconds        temps mural (equilibration + mesure) pour ce point.
     */
    public record Sample(double T, int periodSweeps, double chiPrime, double chiSecond,
                         long thermSweeps, double meanAcceptance, double seconds) { }

    /**
     * Execute le protocole complet (toutes les periodes, toute la grille de temperature) sur une
     * replique, et renvoie un {@link Sample} par couple (periode, temperature).
     *
     * <p>A la sortie — normale ou exceptionnelle, y compris si l'exception survient au milieu de
     * la boucle de mesure d'un point — {@code lattice.getInteractionField()} est restaure a
     * {@code hStatic} : la phase de mesure de chaque (periode, temperature) est entouree d'un
     * {@code try/finally} qui restaure inconditionnellement le champ statique, ce qui rend cet
     * invariant structurel plutot que dependant du chemin d'execution normal.</p>
     *
     * @param lattice  reseau a faire evoluer.
     * @param sweeper  moteur Metropolis damier construit sur {@code lattice} (verifie).
     * @param hStatic  champ statique superpose au champ oscillant, vecteur de longueur 3 dans la
     *                 base de Kitaev (celle de {@code lattice.getSpin1Dlattice()}). Copie
     *                 defensivement : l'appelant peut reutiliser ou muter son tableau ensuite.
     * @param acDrive  direction du champ oscillant, vecteur de longueur 3 non nul dans la meme
     *                 base ; seule sa direction compte, {@link #run} le normalise en interne.
     * @param cfg      parametres du protocole ; copie defensivement, voir {@link Config#copy()}.
     * @return la liste des points de mesure, dans l'ordre periode exterieure / temperature
     *         decroissante interieure.
     * @throws IllegalArgumentException si un argument est invalide (voir {@link #validate}).
     * @throws IllegalStateException si la thermalisation d'un point ne prend pas de decision avant
     *         {@code cfg.thermalization.maxSweeps + cfg.thermalization.sampleEvery} sweeps (meme
     *         garde-fou que {@code MonteCarlo.runParallel}).
     */
    public static List<Sample> run(Lattice lattice, CheckerboardMetropolis sweeper,
                                   double[] hStatic, double[] acDrive, Config cfg) {
        validate(lattice, sweeper, hStatic, acDrive, cfg);
        final Config c = cfg.copy();

        final double[] hStaticLocal = hStatic.clone();
        final double[] acDriveLocal = acDrive.clone();
        final double[] u = normalize(acDriveLocal);

        // Controleur d'equilibration : cree une fois, remis a zero (reset()) avant chaque
        // episode d'equilibration -- comme MonteCarlo.runParallel le fait pour toute une replique
        // (lignes ~919-1020), channels force a 2 (energie, |M|) quel que soit ce que l'appelant a
        // configure ailleurs.
        final AdaptiveThermalization.Config thermCfg = c.thermalization.copy();
        thermCfg.channels = 2;
        final AdaptiveThermalization ctrl = new AdaptiveThermalization(thermCfg);
        final long thermCap = (long) thermCfg.maxSweeps + thermCfg.sampleEvery;

        final int nT = numberOfTemperatures(c);
        final List<Sample> out = new ArrayList<>(c.periodsSweeps.length * nT);

        for (final int period : c.periodsSweeps) {
            final double omega = 2.0 * Math.PI / period;
            final long nTrans = (long) c.transientPeriods * period;
            final long nMeas = (long) c.measurePeriods * period;
            final long nTotal = nTrans + nMeas;

            for (int k = 0; k < nT; k++) {
                final double T = c.tMax - k * c.tStep;
                final long tStartNanos = System.nanoTime();

                // ---- (1) equilibration sous le seul champ statique.
                lattice.setInteractionField(hStaticLocal.clone());
                ctrl.reset();
                double sig = 60.d; // valeur initiale imposee (voir le plan), comme runParallel.
                long thermSweeps = 0L;
                while (!ctrl.isEquilibrated()) {
                    if (thermSweeps >= thermCap) {
                        throw new IllegalStateException("ACSusceptibility : thermalisation sans "
                                + "decision apres " + thermSweeps + " sweeps (periode=" + period
                                + ", T=" + T + ", maxSweeps=" + thermCfg.maxSweeps + ")");
                    }
                    final CheckerboardMetropolis.SweepStats st = sweeper.sweep(T, sig);
                    for (int r = 0; r < c.overRelaxationPerSweep; r++) sweeper.overRelaxationSweep();
                    thermSweeps++;
                    sig = CheckerboardMetropolis.nextSigma(sig, st.rejectionRate());
                    if (ctrl.wantsSample(thermSweeps)) {
                        ctrl.observe(thermSweeps, sweeper.energy(), sweeper.magnetizationNorm());
                    }
                }
                final double sigFrozen = sig;

                // ---- (2) mesure : sur-relaxation coupee, sigma gele, champ H(t) complet. La
                // boucle est entouree d'un try/finally (voir la javadoc de #run) : le champ
                // statique est restaure inconditionnellement en sortie, meme si une exception
                // survient en cours de mesure -- invariant structurel, pas coincidentiel.
                double chiP = 0.d, chiS = 0.d, accAcceptance = 0.d;
                try {
                    final double[] hInst = new double[3];
                    for (long t = 1; t <= nTotal; t++) {
                        final double phase = omega * t;
                        final double drive = Math.sin(phase);
                        hInst[0] = hStaticLocal[0] + c.h0 * drive * acDriveLocal[0];
                        hInst[1] = hStaticLocal[1] + c.h0 * drive * acDriveLocal[1];
                        hInst[2] = hStaticLocal[2] + c.h0 * drive * acDriveLocal[2];
                        lattice.setInteractionField(hInst);

                        final CheckerboardMetropolis.SweepStats st = sweeper.sweep(T, sigFrozen);
                        // Pas de nextSigma, pas d'overRelaxationSweep() : voir la javadoc de classe.

                        if (t > nTrans) {
                            final double mPar = projectedMagnetization(lattice, u);
                            chiP += mPar * drive;
                            chiS += mPar * Math.cos(phase);
                            accAcceptance += st.acceptanceRate();
                        }
                    }
                } finally {
                    // Restauration inconditionnelle du champ statique avant la temperature
                    // suivante (ou avant de rendre la main / propager l'exception).
                    lattice.setInteractionField(hStaticLocal.clone());
                }

                // Normalisation : voir ACPoint. Division par h0 (susceptibilite independante de
                // l'amplitude d'excitation) et par lattice.size (grandeur intensive) ; le facteur
                // ne deplace aucun pic.
                //
                // Signe de chiSecond : la demodulation brute somme_t [mPar(t) . cos(phase)]
                // accumulee dans chiS vaut -chi'' . h0 . lattice.size / 2 (voir la derivation
                // dans ACPoint) -- l'oppose de la convention retenue ici, chiSecond = +chi'' / 2,
                // avec chi'' >= 0 pour une reponse retardee/dissipative comme sur une mesure
                // experimentale. D'ou le signe moins applique ci-dessous, au moment de la
                // normalisation.
                final double norm = (double) nMeas * c.h0 * lattice.size;
                final double chiPrime = chiP / norm;
                final double chiSecond = -chiS / norm;
                final double meanAcceptance = accAcceptance / (double) nMeas;
                final double seconds = (System.nanoTime() - tStartNanos) * 1e-9;

                out.add(new Sample(T, period, chiPrime, chiSecond, thermSweeps, meanAcceptance, seconds));

                if (c.verbose) {
                    System.out.printf(java.util.Locale.US,
                            "  P=%-5d T=%-8.4f therm=%-6d chi'=%-12.6g chi''=%-12.6g acc=%-6.3f %6.2fs%n",
                            period, T, thermSweeps, chiPrime, chiSecond, meanAcceptance, seconds);
                }
            }
        }
        return out;
    }

    /**
     * Regroupe des listes de {@link Sample} produites par plusieurs repliques independantes
     * (memes {@code Config}, donc meme grille de (periode, temperature) dans le meme ordre) en
     * une liste de {@link ACPoint} : moyenne de {@code chiPrime}/{@code chiSecond} sur les
     * repliques, avec erreur standard (ecart-type <b>echantillon</b>, correction de Bessel en
     * {@code n-1}, divise par {@code sqrt(nReplicas)} ; 0 si {@code nReplicas == 1}), et moyennes
     * simples de {@code thermSweeps}, {@code meanAcceptance}, {@code seconds}.
     *
     * @param perReplica une liste par replique, chacune produite par {@link #run} avec la meme
     *                   {@code Config} (donc la meme grille, dans le meme ordre).
     * @param h0         amplitude du champ oscillant, reportee telle quelle dans chaque
     *                   {@link ACPoint} (metadonnee).
     * @param hStatic    valeur de reference du champ statique, reportee telle quelle (metadonnee,
     *                   voir {@link ACPoint#hStatic}).
     * @return un {@link ACPoint} par couple (temperature, periode), dans l'ordre de la premiere
     *         replique.
     * @throws IllegalArgumentException si {@code perReplica} est nul, vide, contient une liste
     *         nulle, ou si les grilles differents d'une replique a l'autre (tailles differentes,
     *         ou temperature/periode qui ne correspondent pas au meme indice).
     */
    public static List<ACPoint> aggregate(List<List<Sample>> perReplica, double h0, double hStatic) {
        if (perReplica == null || perReplica.isEmpty()) {
            throw new IllegalArgumentException("perReplica doit contenir au moins une replique");
        }
        final int nReplicas = perReplica.size();
        final List<Sample> first = perReplica.get(0);
        if (first == null) throw new IllegalArgumentException("perReplica[0] == null");
        final int nPoints = first.size();
        for (int r = 1; r < nReplicas; r++) {
            final List<Sample> rep = perReplica.get(r);
            if (rep == null || rep.size() != nPoints) {
                throw new IllegalArgumentException("grilles incoherentes entre repliques : "
                        + "replique 0 a " + nPoints + " points, replique " + r + " en a "
                        + (rep == null ? "null" : Integer.toString(rep.size())));
            }
        }

        final List<ACPoint> out = new ArrayList<>(nPoints);
        final double[] chiPBuf = new double[nReplicas];
        final double[] chiSBuf = new double[nReplicas];
        for (int i = 0; i < nPoints; i++) {
            final Sample ref = first.get(i);
            double sumChiP = 0.d, sumChiS = 0.d, sumTherm = 0.d, sumAcc = 0.d, sumSec = 0.d;
            for (int r = 0; r < nReplicas; r++) {
                final Sample s = perReplica.get(r).get(i);
                if (s.periodSweeps() != ref.periodSweeps() || Double.compare(s.T(), ref.T()) != 0) {
                    throw new IllegalArgumentException("grilles incoherentes a l'indice " + i
                            + " : replique 0 (T=" + ref.T() + ", periode=" + ref.periodSweeps()
                            + ") vs replique " + r + " (T=" + s.T() + ", periode="
                            + s.periodSweeps() + ")");
                }
                chiPBuf[r] = s.chiPrime();
                chiSBuf[r] = s.chiSecond();
                sumChiP += s.chiPrime();
                sumChiS += s.chiSecond();
                sumTherm += s.thermSweeps();
                sumAcc += s.meanAcceptance();
                sumSec += s.seconds();
            }
            final double meanChiP = sumChiP / nReplicas;
            final double meanChiS = sumChiS / nReplicas;
            final double errChiP = sampleStandardError(chiPBuf, meanChiP);
            final double errChiS = sampleStandardError(chiSBuf, meanChiS);
            final double omega = 2.0 * Math.PI / ref.periodSweeps();

            out.add(new ACPoint(ref.T(), omega, ref.periodSweeps(), h0, hStatic, nReplicas,
                    meanChiP, errChiP, meanChiS, errChiS,
                    sumTherm / nReplicas, sumAcc / nReplicas, sumSec / nReplicas));
        }
        return out;
    }

    /** Ecart-type echantillon (correction de Bessel, n-1) divise par sqrt(n) ; 0 si n &lt;= 1. */
    private static double sampleStandardError(double[] values, double mean) {
        final int n = values.length;
        if (n <= 1) return 0.d;
        double ss = 0.d;
        for (final double v : values) {
            final double d = v - mean;
            ss += d * d;
        }
        final double variance = ss / (n - 1);
        return Math.sqrt(variance / n);
    }

    /**
     * Aimantation totale projetee sur la direction {@code u} (deja normalisee), lue directement
     * sur {@code lattice.getSpin1Dlattice()} — la base de Kitaev, celle dans laquelle
     * {@code hStatic}/{@code acDrive} sont exprimes — et <b>non</b> via
     * {@code CheckerboardMetropolis.magnetization()}, qui renvoie l'aimantation dans une autre
     * base (a*bc, {@code Lattice.getChangeOfBasis()}).
     */
    private static double projectedMagnetization(Lattice lattice, double[] u) {
        final double[] s = lattice.getSpin1Dlattice();
        final int n = lattice.size;
        final double ux = u[0], uy = u[1], uz = u[2];
        double sum = 0.d;
        for (int i = 0; i < n; i++) {
            final int o = 3 * i;
            sum += s[o] * ux + s[o + 1] * uy + s[o + 2] * uz;
        }
        return sum;
    }

    /** Normalise un vecteur 3D ; refuse un vecteur nul ou non fini (direction indefinie). */
    private static double[] normalize(double[] v) {
        final double norm = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        if (!(norm > 0.d) || !Double.isFinite(norm)) {
            throw new IllegalArgumentException(
                    "acDrive doit etre un vecteur non nul et fini : " + java.util.Arrays.toString(v));
        }
        return new double[] { v[0] / norm, v[1] / norm, v[2] / norm };
    }

    /** Nombre de temperatures de la grille lineaire descendante {@code tMax, tMax-tStep, ...}. */
    private static int numberOfTemperatures(Config c) {
        return (int) Math.floor((c.tMax - c.tMin) / c.tStep + 1e-9) + 1;
    }

    /** Validation des bornes du protocole, avant toute allocation ni tout sweep. */
    private static void validate(Lattice lattice, CheckerboardMetropolis sweeper,
                                 double[] hStatic, double[] acDrive, Config cfg) {
        if (lattice == null) throw new IllegalArgumentException("lattice == null");
        if (sweeper == null) throw new IllegalArgumentException("sweeper == null");
        if (sweeper.lattice() != lattice) {
            throw new IllegalArgumentException(
                    "le sweeper ne porte pas sur ce reseau (CheckerboardMetropolis.lattice() != lattice)");
        }
        if (hStatic == null || hStatic.length != 3) {
            throw new IllegalArgumentException("hStatic doit etre un vecteur de longueur 3");
        }
        if (acDrive == null || acDrive.length != 3) {
            throw new IllegalArgumentException("acDrive doit etre un vecteur de longueur 3");
        }
        if (cfg == null) throw new IllegalArgumentException("cfg == null");
        if (!(cfg.h0 > 0.d)) throw new IllegalArgumentException("h0 doit etre > 0 : " + cfg.h0);
        if (cfg.periodsSweeps == null || cfg.periodsSweeps.length == 0) {
            throw new IllegalArgumentException("periodsSweeps ne doit pas etre vide");
        }
        final java.util.Set<Integer> seenPeriods = new java.util.HashSet<>();
        for (final int p : cfg.periodsSweeps) {
            if (p < 2) throw new IllegalArgumentException("periodSweeps doit etre >= 2 : " + p);
            if (!seenPeriods.add(p)) {
                throw new IllegalArgumentException("periodsSweeps contient une periode dupliquee : " + p);
            }
        }
        if (!(cfg.tMin > 0.d)) throw new IllegalArgumentException("tMin doit etre > 0 : " + cfg.tMin);
        if (!(cfg.tMax > cfg.tMin)) {
            throw new IllegalArgumentException(
                    "tMax (" + cfg.tMax + ") doit etre > tMin (" + cfg.tMin + ")");
        }
        if (!(cfg.tStep > 0.d)) throw new IllegalArgumentException("tStep doit etre > 0 : " + cfg.tStep);
        if (cfg.measurePeriods < 1) {
            throw new IllegalArgumentException("measurePeriods doit etre >= 1 : " + cfg.measurePeriods);
        }
        if (cfg.transientPeriods < 0) {
            throw new IllegalArgumentException("transientPeriods doit etre >= 0 : " + cfg.transientPeriods);
        }
    }
}
