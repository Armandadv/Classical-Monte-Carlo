/**
 * Schema de refroidissement adaptatif pour un recuit Monte Carlo : le pas en temperature
 * se resserre automatiquement la ou les fluctuations d'energie (donc la chaleur specifique)
 * sont grandes, c'est-a-dire au voisinage de la transition de phase.
 *
 * <h2>Derivation</h2>
 * <p>Huang, Romeo &amp; Sangiovanni-Vincentelli (1986), <i>An efficient general cooling
 * schedule for simulated annealing</i>, IEEE ICCAD-86, p. 381-384, proposent</p>
 * <pre>
 *   T_{k+1} = T_k exp(-lambda T_k / sigma_E(T_k))
 * </pre>
 * <p>ou sigma_E est l'ecart-type de l'energie a l'equilibre a la temperature T_k. L'idee
 * est de garder constant le "recouvrement" entre les distributions d'energie de deux
 * temperatures consecutives : la variation de l'energie moyenne d'un pas doit rester une
 * fraction fixe de la largeur de la distribution.</p>
 *
 * <p>Pour un systeme de N sites, sigma_E est extensive :
 * sigma_E = sqrt(Var(E)) = T sqrt(N c_v) avec c_v = Var(E) / (N T^2) la chaleur specifique
 * <b>par site</b> (en unites k_B = 1). En reportant :</p>
 * <pre>
 *   ratio = T_{k+1}/T_k = exp(-lambda T / (T sqrt(N c_v))) = exp(-lambda / sqrt(N c_v))
 * </pre>
 * <p>Ce ratio depend de la taille du systeme, ce qui est genant : le meme lambda donnerait
 * des schemas tres differents pour L = 8 et L = 32. On utilise donc la forme
 * <b>par site</b> (equivalente a un lambda re-echelonne lambda' = lambda / sqrt(N)) :</p>
 * <pre>
 *   ratio = exp(-lambda / sqrt(c_v)),      c_v = Var(E_total) / (N T^2)
 * </pre>
 * <p>Avec lambda ~ 0.05 : c_v = 1 -&gt; ratio = 0.9512 (equivalent au 0.95 geometrique
 * historiquement utilise dans ce code) ; c_v = 4 (pic de C_v) -&gt; 0.9753, pas deux fois
 * plus fin ; c_v = 0.25 (regime paramagnetique ou gele) -&gt; 0.9048, on descend plus vite.
 * Le ratio est ensuite borne a [rMin, rMax] pour eviter tout blocage ou tout saut.</p>
 *
 * <p>C'est la version discrete de l'idee de <b>vitesse thermodynamique constante</b> de
 * Salamon, Nulton, Harland, Pedersen, Ruppeiner &amp; Liao (1988), <i>Simulated annealing
 * with constant thermodynamic speed</i>, Comput. Phys. Commun. <b>49</b>, 423-428 :
 * dT/dt proportionnel a -T / (tau sqrt(C)), qui prescrit egalement des pas plus petits
 * la ou la chaleur specifique est grande.</p>
 *
 * <p><b>Attention</b> : avec des pas adaptatifs le nombre de temperatures parcourues n'est
 * <i>pas connu a l'avance</i> (contrairement au schema geometrique ou
 * nMcs = floor(log(Tend/Tinit)/log(dT))). Les structures de sortie doivent donc etre des
 * listes a croissance dynamique, et les barres de progression doivent etre exprimees en
 * log(T) plutot qu'en nombre d'iterations.</p>
 */
public final class AdaptiveCoolingSchedule {

    private final double lambda;
    private final double rMin;
    private final double rMax;
    private final double geometricFallback;

    private double lastRatio;
    private double lastCv;

    /** Valeurs par defaut : lambda = 0.05, rMin = 0.85, rMax = 0.99, repli geometrique 0.95. */
    public AdaptiveCoolingSchedule() {
        this(0.05, 0.85, 0.99, 0.95);
    }

    /**
     * @param lambda            parametre de Huang (adimensionne) ; plus il est grand, plus
     *                          les pas sont fins. 0.05 reproduit ~0.95 pour c_v = 1.
     * @param rMin              ratio minimal autorise (pas le plus grossier).
     * @param rMax              ratio maximal autorise (pas le plus fin).
     * @param geometricFallback ratio utilise lorsque c_v n'est pas exploitable
     *                          (NaN, &lt;= 0, T = 0, N = 0).
     * @throws IllegalArgumentException si {@code lambda <= 0}, si la fenetre de ratios ne
     *         verifie pas {@code 0 < rMin <= rMax < 1}, ou si
     *         {@code geometricFallback} n'est pas dans {@code ]0, 1[} (les valeurs non finies
     *         sont rejetees par les memes comparaisons).
     *
     *         <p>Ce n'est pas du pedantisme : tout ratio &gt;= 1 rend la suite
     *         T_{k+1} = T_k . ratio non decroissante, et la boucle
     *         {@code while (T > endTemp)} de {@code MonteCarlo.runParallel} ne peut alors
     *         <b>jamais</b> se terminer &mdash; le recuit tourne indefiniment a temperature
     *         constante ou croissante. Un ratio &lt;= 0 est tout aussi fatal (temperature
     *         nulle ou negative), et lambda &lt;= 0 inverse l'adaptation (pas plus
     *         <i>grossiers</i> la ou C_v est grand). Le seul moment ou ces erreurs sont
     *         detectables a coup sur est la construction.</p>
     */
    public AdaptiveCoolingSchedule(double lambda, double rMin, double rMax, double geometricFallback) {
        if (!(lambda > 0.0)) {
            throw new IllegalArgumentException("lambda doit etre > 0 : " + lambda);
        }
        if (!(rMin > 0.0) || !(rMin <= rMax) || !(rMax < 1.0)) {
            throw new IllegalArgumentException(
                    "il faut 0 < rMin <= rMax < 1 (sinon la boucle de recuit ne termine pas) : rMin="
                            + rMin + ", rMax=" + rMax);
        }
        if (!(geometricFallback > 0.0) || !(geometricFallback < 1.0)) {
            throw new IllegalArgumentException(
                    "geometricFallback doit etre dans ]0, 1[ : " + geometricFallback);
        }
        this.lambda = lambda;
        this.rMin = rMin;
        this.rMax = rMax;
        this.geometricFallback = geometricFallback;
        this.lastRatio = geometricFallback;
        this.lastCv = Double.NaN;
    }

    /**
     * Temperature suivante.
     *
     * @param varEnergyTotal variance (de population) de l'energie <b>totale</b> mesuree a la
     *                       temperature T pendant la production.
     * @param nSites         nombre de sites du reseau.
     * @return T * ratio, avec ratio = clamp(exp(-lambda / sqrt(c_v)), rMin, rMax) et
     *         c_v = varEnergyTotal / (nSites T^2) ; ratio = geometricFallback si c_v n'est
     *         pas un reel strictement positif.
     */
    public double next(double T, double varEnergyTotal, int nSites) {
        double cv = Double.NaN;
        if (nSites > 0 && Double.isFinite(T) && T > 0.0 && Double.isFinite(varEnergyTotal)) {
            cv = varEnergyTotal / (nSites * T * T);
        }
        lastCv = cv;

        double ratio;
        if (!(cv > 0.0) || Double.isNaN(cv) || !Double.isFinite(cv)) {
            ratio = geometricFallback;
        } else {
            ratio = Math.exp(-lambda / Math.sqrt(cv));
            if (!Double.isFinite(ratio)) ratio = geometricFallback;
            if (ratio < rMin) ratio = rMin;
            if (ratio > rMax) ratio = rMax;
        }
        lastRatio = ratio;
        return T * ratio;
    }

    /** Dernier ratio T_{k+1}/T_k effectivement applique (apres bornage). */
    public double lastRatio() { return lastRatio; }

    /** Derniere chaleur specifique par site utilisee (NaN si non exploitable). */
    public double lastCv() { return lastCv; }

    public double lambda() { return lambda; }
    public double rMin() { return rMin; }
    public double rMax() { return rMax; }
    public double geometricFallback() { return geometricFallback; }
}
