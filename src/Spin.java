
import java.util.Date;

import org.apache.commons.math3.util.FastMath;

import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.sampling.distribution.GaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ZigguratSampler;
import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.rng.UniformRandomProvider;

import java.util.random.RandomGenerator;

/**
 * Spin d'un site du reseau : energie locale et totale, variations d'energie mono-spin
 * (getdE/getdE2 : echange + Zeeman seulement) et tirages Monte Carlo (SRG7, Gaussian_move).
 * La dynamique Landau-Lifshitz passe par {@link SpinODE}, integre par le RK4 a pas fixe de
 * {@link Rk4PasFixe} ({@link #Precession1D()}, {@link #solveODE1D(double, double, int, int)}).
 * Le dE du chemin parallele {@link CheckerboardMetropolis} inclut en plus le terme quartique ;
 * ceux d'ici non (voir getdE).
 */
public class Spin {

    UniformRandomProvider rng = RandomSource.XO_SHI_RO_256_PP.create();
    ContinuousSampler GaussSampler = GaussianSampler.of(ZigguratSampler.NormalizedGaussian.of(rng),0.d,1.d);


    final double dt = 0.01; // pas d'integrétion RK et correlation
    /** Pas d'integration RK4, expose pour la derivation du sampling S(Q,w) dans MonteCarlo. */
    public double getDt() { return dt; }
    /**
     * Pas etalonne pour {@link #Precession1D()} par la sonde de derive d'energie (NaN tant
     * que non calibre). Mis en cache : la dynamique de precession ne depend pas de T, le
     * meme pas convient a tous les appels.
     */
    private double precessionDt = Double.NaN;
    private final Lattice lattice;
    private final double[][][] interactionMatrices;


    /** Met en cache, par site, les references des gabarits d'interaction du reseau. */
    public Spin(Lattice lattice) {
        this.lattice = lattice;
        this.interactionMatrices = new double[lattice.size][][];
        for (int k = 0; k < lattice.size; k++) {
            this.interactionMatrices[k] = lattice.getInteractionMatrices(k);
        }
    }


    /** Vecteur unitaire uniforme sur la sphere (rejet : (x0, x1) uniformes dans le disque unite). */
    public double[] SRG7(RandomGenerator random) {
        double x0, x1, r;
        do {
            x0 = random.nextDouble() * 2 - 1;
            x1 = random.nextDouble() * 2 - 1;
            r = x0 * x0 + x1 * x1;
        } while (r >= 1);
        double c = 2 * FastMath.sqrt(1 - r);
        return new double[] { c * x0, c * x1, 1 - 2 * r };
    }

    /** Trois composantes gaussiennes centrees reduites independantes (non normalisees). */
    public double[] Gamma (){
        final double[] randomState = new double[3];
        randomState[0] = GaussSampler.sample();
        randomState[1] = GaussSampler.sample();
        randomState[2] = GaussSampler.sample();

        return randomState;

    }

    // J D Alzate-Cardona etal 2019 J.Phys.:Condens.Matter 31 095802
    /**
     * Mouvement gaussien : nouvel etat du spin du site = normalisation de
     * (ancien + sig * gaussien), donc |S| = 1 impose. Le tirage passe par l'echantillonneur
     * interne ({@link #Gamma()}) ; le parametre random n'est pas utilise.
     */
    public double[] Gaussian_move(int site, double sig, RandomGenerator random) {

        final double[] oldState = lattice.getSpin1D(site);
        final double[] randomState = Gamma();

        final double[] sum = MathOps.vectorizationadd(oldState, MathOps.scale(sig, randomState));
        final double norm = MathOps.norm(sum);
        sum[0] /= norm;
        sum[1] /= norm;
        sum[2] /= norm;
        return sum;
    }

    /**
     * Produit J.Sj pour les gabarits compacts partages avec les noyaux Checkerboard :
     * 2 valeurs (2-PARAM), 3 (diagonale), 5 (Kitaev), sinon pleine 3x3 (triangle superieur
     * lu ; indices 3, 6, 7 ignores par symetrie, comme dans SpinODE).
     */
    public static double[] exchangeEnergySym(double[] M, double[] Sj) {
        final double m00, m01, m02, m11, m12, m22;
        double[] MijSi;
        final double s0 = Sj[0];
        final double s1 = Sj[1];
        final double s2 = Sj[2];

        /*
         * 2-PARAM a.I + b.A : gabarit compact {a, b}
         */
        if (M.length == 2) {
            // 2-PARAM a.I + b.A
            MijSi = new double[] { M[0] * s0 + M[1] * (s1 + s2),
                    M[1] * s0 + M[0] * s1 + M[1] * s2,
                    M[1] * (s0 + s1) + M[0] * s2 };
        }

        /*
         * Matrice diagonale
         */
        else if (M.length == 3) {
            m00 = M[0];
            m11 = M[1];
            m22 = M[2];
            MijSi = new double[] {m00 * s0, m11 * s1 , m22 * s2};
        }

        /*
         * Kitaev : diag(d0, d1, d2) + exactement une paire hors-diagonale, gabarit
         * compact {d0, d1, d2, k, v} (jumeaux exactement egaux m10==m01 etc.) ;
         * termes nuls omis, memes ordres de somme que la forme pleine.
         */
        else if (M.length == 5) {
            final int k = (int) M[3];
            final double v = M[4];
            if (k == 0) {
                MijSi = new double[] { M[0] * s0 + v * s1, v * s0 + M[1] * s1, M[2] * s2 };
            } else if (k == 1) {
                MijSi = new double[] { M[0] * s0 + v * s2, M[1] * s1, v * s0 + M[2] * s2 };
            } else {
                MijSi = new double[] { M[0] * s0, M[1] * s1 + v * s2, v * s1 + M[2] * s2 };
            }
        }

        else {
            m00 = M[0];
            m01 = M[1];
            m02 = M[2];
            m11 = M[4];
            m12 = M[5];
            m22 = M[8];

            MijSi =  new double[] {
                m00 * s0 + m01 * s1 + m02 * s2,
                m01 * s0 + m11 * s1 + m12 * s2,
                m02 * s0 + m12 * s1 + m22 * s2
        };

        }

        return MijSi;
    }

    /**
     * Variation d'energie (echange + Zeeman) du remplacement du spin du site par newState.
     * N'inclut PAS le terme quartique, contrairement au dE de {@link CheckerboardMetropolis} :
     * a ne pas utiliser pour des simulations avec QUARTIC.
     */
    public double getdE(final int site, final double[] newState) {
        double dE = 0.d;
        final double[] oldState = lattice.getSpin1D(site);
        final double[] dSpin = new double[] { newState[0] - oldState[0], newState[1] - oldState[1],
                newState[2] - oldState[2] };
        final int[] interactionSites = lattice.getNeighborSites(site);
        final double[][] interactionMatrices = lattice.getInteractionMatrices(site);

        for (int i = 0; i < interactionSites.length; i++) {
            dE -= MathOps.dot(dSpin,
                    exchangeEnergySym(interactionMatrices[i], lattice.getSpin1D(interactionSites[i]))) ;
        }

        dE -= MathOps.dot(dSpin, lattice.getInteractionField());

        return dE;
    }

    /**
     * Comme {@link #getdE} (sans terme quartique), mais factorise les voisins consecutifs
     * partageant la meme instance de matrice : J.(somme des spins) une fois par groupe au
     * lieu de J.Sj par voisin — mathematiquement equivalent, moins de produits.
     */
    public double getdE2(final int site, final double[] newState) {
        double dE = 0.d;

        final double[] oldState = lattice.getSpin1D(site);

        final double[] dSpin = new double[] { newState[0] - oldState[0], newState[1] - oldState[1],
                newState[2] - oldState[2] };

        final int[] interactionSites = lattice.getNeighborSites(site);
        final double[][] interactionMatrices = lattice.getInteractionMatrices(site);

        int i = 0;

        while (i < interactionSites.length) {
            // On récupère la matrice d'interaction courante
            double[] currentMatrix = interactionMatrices[i];

            double[] spinSum = new double[] { 0.d, 0.d, 0.d };

            // Tant que la prochaine matrice est identique alors on factorise les spins
            do {
                double[] neighborSpin = lattice.getSpin1D(interactionSites[i]);
                spinSum[0] += neighborSpin[0];
                spinSum[1] += neighborSpin[1];
                spinSum[2] += neighborSpin[2];

                i++;
            } while (i < interactionSites.length && currentMatrix == interactionMatrices[i]);
            double[] factorizedEnergy = exchangeEnergySym(currentMatrix, spinSum);
            dE -= MathOps.dot(dSpin, factorizedEnergy);
        }
        dE -= MathOps.dot(dSpin, lattice.getInteractionField());
        return dE;
    }

    /**
     * Energie du site : echange sur les voisins d'interaction declares (chaque liaison
     * compte une fois, garde site > voisin) + Zeeman. Le terme quartique n'y figure pas :
     * getEnergy l'ajoute a part via {@link #getQuarticEnergy()}.
     */
    public double getEnergyLocal(int site) {
        double energy = 0.d;
        double[] s0 = lattice.getSpin1D(site);
        int[] interactionSites = lattice.getNeighborSites(site);
        double[][] interactionMatrices = lattice.getInteractionMatrices(site);
        double[] JkjxSj = new double[3];

        for (int i = 0; i < interactionSites.length; i++) {
            if (site > interactionSites[i]) {
                JkjxSj = MathOps.vectorizationadd(JkjxSj,
                        exchangeEnergySym(interactionMatrices[i],
                                lattice.getSpin1D(interactionSites[i])));
            }
        }
        energy -= MathOps.dot(JkjxSj, s0);
        energy -= MathOps.dot(s0, lattice.getInteractionField());
        return energy;
    }

    /** Energie totale : somme des energies locales, plus le terme quartique s'il est actif. */
    public double getEnergy() {
        double energy = 0.d;
        for (int site = 0; site < lattice.size; site++) {
            energy += getEnergyLocal(site);
        }
        if (lattice.hasQuartic()) {
            energy += getQuarticEnergy();
        }
        return energy;
    }

    /**
     * Partie quartique premiers voisins : {@code E4 = -b somme_liaisons (S_i^g S_j^g)^2},
     * chaque liaison comptee une fois (garde site > voisin, meme convention que
     * {@link #getEnergyLocal}). Renvoie 0 si le terme n'est pas actif.
     */
    public double getQuarticEnergy() {
        if (!lattice.hasQuartic()) {
            return 0.d;
        }
        final double b = lattice.getQuarticB();
        double e4 = 0.d;
        for (int site = 0; site < lattice.size; site++) {
            final double[] s0 = lattice.getSpin1D(site);
            final int[] neighbors = lattice.getNeighborSites(site);
            for (int i = 0; i < neighbors.length; i++) {
                if (site > neighbors[i]) {
                    final int g = lattice.getQuarticGamma(site, i);
                    final double[] s1 = lattice.getSpin1D(neighbors[i]);
                    final double p = s0[g] * s1[g];
                    e4 += p * p;
                }
            }
        }
        return -b * e4;
    }


    /**
     * Intègre la dynamique de précession (Rk4PasFixe, 4 évaluations de dérivée par pas,
     * sans framework) et sauvegarde ndt échantillons normalisés de S(t) pour le canal SQW.
     *
     * <p>Deux sondes de diagnostic sont imprimées à chaque appel (générales, aucune
     * hypothèse physique) : dérive de norme max (pas de temps trop grand ?) et
     * décorrélation &lt;S(0),S(tf)&gt;/(3N) (fenêtre assez longue ?).</p>
     *
     * <p>Nombre de pas : l'ancien code commons-math avançait à pas fixe {@code dt} en
     * plafonnant le temps final (dernier pas clampé) ; le nouveau calcule
     * {@code nSteps = (int)((tf - t0) / dt)} (troncature). Pour tous les appelants actuels
     * {@code (tf - t0)} est un multiple exact de {@code dt}, donc la troncature ne change
     * rien (quotient déjà entier, validé par décompilation de commons-math 3.6.1) ;
     * l'échantillonnage {@code index 0-based % measurementRate == 0} est lui identique
     * (le StepHandler sauvegardait aussi quand {@code saveCounter % rate == 0},
     * {@code saveCounter} valant l'index 0-based du pas).</p>
     */
    public double[] solveODE1D(final double t0, final double tf,
                               final int measurementRate, final int ndt) {
        final int nSteps = Math.max(1, (int) ((tf - t0) / dt));
        final double[] dSdt = new double[ndt * lattice.size * 3];
        final double[] spins = lattice.getSpin1Dlattice();
        final SpinODE ode = new SpinODE(lattice);
        final Rk4PasFixe rk4 = new Rk4PasFixe((t, s, dS) -> ode.computeDerivatives(t, s, dS));
        final int[] compteur = new int[1];
        rk4.integre(t0, tf, nSteps, spins, (t, s, nPasFait) -> {
            if (nPasFait % measurementRate == 0 && compteur[0] < ndt) {
                for (int i = 0; i < lattice.size; i++) {
                    final double norm = MathOps.norm(s, i * 3);
                    final int index = compteur[0] * lattice.size * 3 + i * 3;
                    dSdt[index] = s[i * 3] / norm;
                    dSdt[index + 1] = s[i * 3 + 1] / norm;
                    dSdt[index + 2] = s[i * 3 + 2] / norm;
                }
                compteur[0]++;
            }
            return false;
        });
        return dSdt;
    }


    /**
     * Précession d'essai sous la dynamique pure (conservative), intégrée par le RK4 maison
     * à pas fixe. Les deux paramètres d'intégration ne sont pas figés, ils sont pilotés
     * par les sondes :
     * <ul>
     *   <li>sonde 1 — dérive d'énergie relative sur une fenêtre de base : au-delà de
     *       {@code deriveMax}, dt est divisé par deux et l'essai recommence (la précession
     *       exacte conserve E). Étalonnage une seule fois et mis en cache : la dynamique
     *       ne dépend pas de la température ;</li>
     *   <li>sonde 2 — décorrélation ⟨S(0), S(t)⟩/(3N) : l'intégration se poursuit par
     *       fenêtres de {@code tf0} (sans repartir de l'état initial) tant que la
     *       corrélation résiduelle dépasse {@code corrMax}, dans la limite de
     *       {@code fenetresMax} fenêtres. À basse T la corrélation peut plafonner sans
     *       jamais décorréler (dynamique quasi périodique) : le plafond borne alors le
     *       coût.</li>
     * </ul>
     */
    public void Precession1D() {
        final double deriveMax = 1e-6d;   // sonde 1 : dérive d'énergie relative tolérée
        final double corrMax = 0.05d;     // sonde 2 : corrélation résiduelle tolérée
        final double dtMin = 1e-4d;       // plancher d'étalonnage : au-delà, dt est gardé tel quel
        final double tf0 = 5.0d;         // fenêtre de base (valeur historique)
        final int fenetresMax = 16;        // plafond sonde 2 : tf <= tf0 * fenetresMax

        final SpinODE ode = new SpinODE(lattice);
        final Rk4PasFixe rk4 = new Rk4PasFixe((t, s, dS) -> ode.computeDerivatives(t, s, dS));

        final double[] initial = lattice.getSpin1Dlattice().clone();
        final double[] s = lattice.getSpin1Dlattice();
        final double eAvant = getEnergy();

        // Sonde 1 : étalonnage de dt sur une fenêtre de base (une fois pour la vie du Spin).
        if (Double.isNaN(precessionDt)) {
            double essai = 0.05d;
            while (true) {
                System.arraycopy(initial, 0, s, 0, s.length);
                rk4.integre(0.0d, tf0, Math.max(1, (int) Math.round(tf0 / essai)), s, null);
                final double derive =
                        Math.abs(getEnergy() - eAvant) / Math.max(Math.abs(eAvant), 1.0);
                if (derive <= deriveMax || essai <= dtMin) { precessionDt = essai; break; }
                essai /= 2.0d;
            }
        }

        // Intégration par fenêtres de tf0 ; sonde 2 (décorrélation) à chaque frontière.
        System.arraycopy(initial, 0, s, 0, s.length);
        double corr = 1.0d;
        int fenetres = 0;
        while (true) {
            rk4.integre(fenetres * tf0, (fenetres + 1) * tf0,
                    Math.max(1, (int) Math.round(tf0 / precessionDt)), s, null);
            fenetres++;
            for (int i = 0; i < lattice.size; i++) {
                final int b = i * 3;
                final double norm = MathOps.norm(s, b);
                s[b] /= norm;
                s[b + 1] /= norm;
                s[b + 2] /= norm;
            }
            corr = 0.0d;
            for (int i = 0; i < initial.length; i++) corr += initial[i] * s[i];
            corr /= initial.length;
            if (Math.abs(corr) <= corrMax || fenetres >= fenetresMax) break;
        }

        lattice.setSpin1Dlattice(s);
    }
}
