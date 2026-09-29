/**
 * Recouvrement de configuration (overlap d'Edwards-Anderson) entre l'etat courant d'un
 * reseau de spins et une reference capturee anterieurement.
 *
 * <pre>
 *   q(t) = (1/N) sum_i S_i(t) . S_i(ref)
 * </pre>
 * <p>avec des spins unitaires, donc q dans [-1, 1] et q(0) = 1 par construction. Cette
 * observable ne mesure pas la relaxation vers l'equilibre thermique (c'est le role de
 * {@code EquilibrationDetector} sur l'energie ou |M|) mais l'oubli de la configuration
 * <b>heritee</b> d'une temperature precedente : dans une phase ordonnee ou frustree, q ne
 * tend pas vers 0 mais decroit vers un plateau q_inf non nul (parametre d'ordre
 * d'Edwards-Anderson). Voir spec_stagnation.md, section "Observable centrale", et :</p>
 * <ul>
 *   <li><b>Ogielski, A. T.</b> (1985). <i>Dynamics of three-dimensional Ising spin
 *       glasses in thermal equilibrium.</i> Phys. Rev. B <b>32</b>, 7384. &mdash; q(t)
 *       comme diagnostic d'equilibration pour des systemes frustres, decroissance vers
 *       un plateau non nul.</li>
 *   <li><b>Marinari, E. &amp; Parisi, G.</b> &mdash; utilisation de l'overlap comme
 *       parametre d'ordre et diagnostic de decorrelation en Monte Carlo de verres de
 *       spin.</li>
 * </ul>
 *
 * <p>Usage attendu (voir {@code StagnationDetector}) : capturer la configuration heritee
 * en debut de palier de temperature, puis suivre {@code overlap(lattice)} au fil des
 * sweeps de thermalisation. {@link #analyse} extrait ensuite de la serie q(t) un temps de
 * decroissance tau_q et le plateau q_inf : thermSweeps &gt;= K . tau_q (K ~ 3-5) garantit
 * une memoire residuelle &lt; 2 % ; q_inf &gt;= qFrozen (~0.9) signale une chaine gelee
 * (aucun sweep supplementaire ne la decorrelera de l'etat herite).</p>
 *
 * <p>Cette classe est volontairement independante de {@code EquilibrationDetector} pour
 * la partie "etat + reseau" (elle en reutilise seulement les utilitaires statistiques
 * dans {@link #analyse}) : elle ne connait que des vecteurs de spins bruts, jamais
 * l'energie ou l'aimantation.</p>
 */
public final class SpinOverlap {

    private final int nSites;
    private double[] reference;   // longueur 3*nSites ; null si aucune reference capturee

    /**
     * @param nSites nombre de sites du reseau (&gt;= 1).
     * @throws IllegalArgumentException si {@code nSites < 1}.
     */
    public SpinOverlap(int nSites) {
        if (nSites < 1) {
            throw new IllegalArgumentException("nSites doit etre >= 1, recu " + nSites);
        }
        this.nSites = nSites;
    }

    /** Nombre de sites configure a la construction. */
    public int nSites() { return nSites; }

    /**
     * Copie defensive de la configuration courante du reseau comme reference.
     *
     * <p>{@code lattice.getSpin1Dlattice()} expose le tableau <b>vivant</b> des spins (pas
     * de copie, pas d'allocation par site) : capturer sans copier laisserait la reference
     * suivre silencieusement les modifications ulterieures du reseau, ce qui viderait la
     * notion meme de "configuration heritee". D'ou la copie ici.</p>
     *
     * @throws IllegalArgumentException si la taille du reseau ne correspond pas a
     *         {@code nSites}.
     */
    public void capture(Lattice lattice) {
        double[] s = lattice.getSpin1Dlattice();
        checkLength(s.length, "lattice");
        ensureBuffer();
        System.arraycopy(s, 0, reference, 0, 3 * nSites);
    }

    /**
     * Copie defensive de la reference depuis un tableau brut {@code [Sx0,Sy0,Sz0,Sx1,...]}
     * de longueur {@code 3*nSites}.
     *
     * @throws IllegalArgumentException si {@code spins1D} est null ou de longueur incorrecte.
     */
    public void capture(double[] spins1D) {
        if (spins1D == null) throw new IllegalArgumentException("spins1D null");
        checkLength(spins1D.length, "spins1D");
        ensureBuffer();
        System.arraycopy(spins1D, 0, reference, 0, 3 * nSites);
    }

    /** Vrai si une reference a ete capturee (et non oubliee par {@link #clear()}). */
    public boolean hasReference() { return reference != null; }

    /** Oublie la reference courante ({@link #hasReference()} redevient false). */
    public void clear() { reference = null; }

    /**
     * q = (1/N) sum_i S_i . S_i^ref entre la configuration courante du reseau et la
     * reference capturee. Boucle unique O(3N), aucune allocation (pas d'appel a
     * {@code getSpin1D(site)}, qui alloue un {@code double[3]} par site).
     *
     * @throws IllegalStateException si aucune reference n'a ete capturee.
     * @throws IllegalArgumentException si la taille du reseau ne correspond pas a
     *         {@code nSites}.
     */
    public double overlap(Lattice lattice) {
        requireReference();
        double[] s = lattice.getSpin1Dlattice();
        checkLength(s.length, "lattice");
        return dot(s, reference, nSites);
    }

    /**
     * q entre {@code spins1D} et la reference capturee.
     *
     * @throws IllegalStateException si aucune reference n'a ete capturee.
     * @throws IllegalArgumentException si {@code spins1D} est null ou de longueur incorrecte.
     */
    public double overlap(double[] spins1D) {
        requireReference();
        if (spins1D == null) throw new IllegalArgumentException("spins1D null");
        checkLength(spins1D.length, "spins1D");
        return dot(spins1D, reference, nSites);
    }

    /** Copie defensive de la reference courante (pour les tests) ; longueur {@code 3*nSites}. */
    public double[] reference() {
        requireReference();
        return reference.clone();
    }

    /**
     * Recouvrement statique entre deux configurations brutes {@code a} et {@code b}, sans
     * aucun etat. Meme semantique que {@link #overlap(double[])} mais reutilisable sans
     * instance (p.ex. pour comparer deux configurations independantes dans les tests).
     *
     * @throws IllegalArgumentException si {@code nSites < 1}, si {@code a} ou {@code b}
     *         est null, ou si leur longueur n'est pas exactement {@code 3*nSites}.
     */
    public static double overlap(double[] a, double[] b, int nSites) {
        if (nSites < 1) throw new IllegalArgumentException("nSites doit etre >= 1, recu " + nSites);
        if (a == null || b == null) throw new IllegalArgumentException("a et b ne doivent pas etre null");
        int need = 3 * nSites;
        if (a.length != need) {
            throw new IllegalArgumentException("a.length=" + a.length + " != 3*nSites=" + need);
        }
        if (b.length != need) {
            throw new IllegalArgumentException("b.length=" + b.length + " != 3*nSites=" + need);
        }
        return dot(a, b, nSites);
    }

    private static double dot(double[] a, double[] b, int nSites) {
        double sum = 0.0;
        int n3 = 3 * nSites;
        for (int i = 0; i < n3; i++) {
            sum += a[i] * b[i];
        }
        return sum / nSites;
    }

    private void ensureBuffer() {
        if (reference == null) reference = new double[3 * nSites];
    }

    private void checkLength(int length, String what) {
        int need = 3 * nSites;
        if (length != need) {
            throw new IllegalArgumentException(
                    what + " de longueur " + length + " incoherente avec nSites=" + nSites
                            + " (attendu " + need + ")");
        }
    }

    private void requireReference() {
        if (reference == null) {
            throw new IllegalStateException("aucune reference capturee : appeler capture(...) d'abord");
        }
    }

    // ------------------------------------------------------------------
    // Analyse de la decroissance q(t) -> q_inf
    // ------------------------------------------------------------------

    /**
     * Resultat de l'analyse d'une serie de recouvrement q(t) sur une fenetre {@code [from, to)}.
     *
     * @param qInf      plateau estime : moyenne de q sur la derniere fraction
     *                  {@code plateauFraction} de la fenetre.
     * @param tauSamples temps de decroissance vers le plateau, en <b>echantillons</b>
     *                  (interpole lineairement, pas quantifie au pas d'echantillonnage) ;
     *                  0 si la serie est deja au plateau des {@code from} ;
     *                  {@code NaN} si le seuil de decroissance n'est jamais atteint.
     * @param plateaued vrai si la premiere et la seconde moitie de la fenetre de plateau
     *                  ont des moyennes compatibles (voir {@link #analyse}) ; toujours
     *                  false si {@code tauSamples} est NaN (la decroissance n'a pas fini).
     * @param nUsed     nombre d'echantillons de la fenetre {@code [from, to)} effectivement
     *                  analyses ({@code to - from}, ou 0 si la fenetre est vide).
     */
    public record Decay(double qInf, double tauSamples, boolean plateaued, int nUsed) { }

    /** {@link #analyse(double[], int, int, double)} avec {@code plateauFraction = 0.5}. */
    public static Decay analyse(double[] q, int from, int to) {
        return analyse(q, from, to, 0.5);
    }

    /**
     * Estime le plateau q_inf et le temps de decroissance tau_q d'une serie q(t) sur
     * {@code [from, to)}.
     *
     * <h3>q_inf</h3>
     * <p>Moyenne de q sur la derniere fraction {@code plateauFraction} de la fenetre
     * (par defaut la seconde moitie).</p>
     *
     * <h3>tau_q</h3>
     * <p>Premier indice t (en echantillons, relatif a {@code from}) tel que
     * {@code q[t] - qInf <= (q[from] - qInf) / e}, avec interpolation lineaire entre les
     * deux points encadrants pour ne pas quantifier tau au pas d'echantillonnage
     * (important quand tau &lt; quelques pas d'echantillonnage, cf. la meme remarque pour
     * {@code AdaptiveThermalization.tauSweepsFromSamples}). Si {@code q[from] - qInf <= 0}
     * la serie est deja au plateau : {@code tauSamples = 0}. Si le seuil n'est jamais
     * atteint dans la fenetre, {@code tauSamples = NaN} et {@code plateaued = false}
     * (une decroissance qui n'a pas atteint 1/e de son ecart initial au plateau n'est de
     * toute evidence pas terminee).</p>
     *
     * <h3>plateaued</h3>
     * <p>La fenetre de plateau (les derniers {@code plateauFraction . (to - from)}
     * echantillons) est coupee en deux moities egales ; {@code plateaued} est vrai si
     * leurs moyennes different de moins d'une erreur standard combinee
     * {@code sqrt(se1^2 + se2^2)}, {@code se} etant {@link EquilibrationDetector#standardError}
     * de chaque moitie. C'est le meme principe que le test de Geweke
     * ({@link EquilibrationDetector#gewekeZ}) applique specifiquement a la fenetre de
     * plateau plutot qu'a toute la queue : on ne demande pas que q soit devenu constant,
     * seulement que sa valeur recente ne bouge plus au-dela du bruit statistique.</p>
     *
     * @param plateauFraction fraction de la fenetre utilisee pour q_inf et le test de
     *                        plateau, dans (0, 1] ; valeur non finie ou hors bornes -&gt; 0.5.
     * @return {@code nUsed = 0} et toutes les autres valeurs a leur defaut si la fenetre
     *         est vide.
     */
    public static Decay analyse(double[] q, int from, int to, double plateauFraction) {
        int n = to - from;
        if (n <= 0) return new Decay(0.0, Double.NaN, false, 0);

        double frac = (Double.isFinite(plateauFraction) && plateauFraction > 0.0 && plateauFraction <= 1.0)
                ? plateauFraction : 0.5;
        int plateauLen = Math.max(1, (int) Math.round(frac * n));
        if (plateauLen > n) plateauLen = n;
        int plateauFrom = to - plateauLen;

        double qInf = EquilibrationDetector.mean(q, plateauFrom, to);

        double tauSamples = tauFromThreshold(q, from, to, qInf);

        boolean plateaued;
        if (Double.isNaN(tauSamples)) {
            // Le seuil de 1/e n'a pas ete atteint : la decroissance n'est pas finie, le
            // plateau ne peut donc pas etre considere comme atteint quel que soit ce que
            // dit la comparaison de moyennes ci-dessous (qui pourrait accidentellement
            // "voir" un palier transitoire au milieu d'une decroissance lente).
            plateaued = false;
        } else {
            plateaued = plateauStable(q, plateauFrom, to, plateauLen);
        }

        return new Decay(qInf, tauSamples, plateaued, n);
    }

    private static double tauFromThreshold(double[] q, int from, int to, double qInf) {
        double target = q[from] - qInf;
        // Tolerance relative (et non target <= 0.0 strict) : qInf est la moyenne flottante
        // d'une sous-fenetre de q, calculee par sommation puis division. Pour une serie
        // logiquement constante (p.ex. 45 copies exactement egales de 0.95), cette moyenne
        // ne retombe generalement PAS bit-a-bit sur la valeur commune (0.95 n'est pas exact
        // en binaire ; la somme de n copies approchees derive de quelques ulps), si bien que
        // target peut valoir un epsilon strictement positif (~1e-16) au lieu de 0 exactement
        // -- et ce MEME epsilon affecte alors tous les points de la fenetre, qui ne peuvent
        // donc jamais franchir le seuil : sans cette tolerance, la boucle ci-dessous
        // renverrait a tort NaN / plateaued=false sur une serie pourtant parfaitement
        // stable. 1e-12 absorbe tres largement le bruit d'arrondi (n termes, erreur relative
        // ~n.eps_machine, eps_machine ~ 2e-16) tout en restant negligeable devant tout signal
        // physique de decroissance.
        double eps = 1e-12 * Math.max(1.0, Math.abs(qInf));
        if (target <= eps) return 0.0;   // deja au plateau (ou en dessous, aux arrondis pres)

        double thresh = qInf + target / Math.E;
        for (int t = from; t < to; t++) {
            if (q[t] <= thresh) {
                if (t == from) return 0.0;
                double qPrev = q[t - 1];
                double qCurr = q[t];
                double denom = qCurr - qPrev;
                double f = (denom != 0.0) ? (thresh - qPrev) / denom : 0.0;
                if (f < 0.0) f = 0.0;
                if (f > 1.0) f = 1.0;
                return (t - 1 - from) + f;
            }
        }
        return Double.NaN;   // jamais atteint dans la fenetre
    }

    private static boolean plateauStable(double[] q, int plateauFrom, int to, int plateauLen) {
        int half = plateauLen / 2;
        if (half < 1) return true;   // fenetre trop courte pour etre coupee : rien a contredire

        // Court-circuit numerique : si la fenetre est constante a une tolerance flottante
        // pres (typiquement une serie MC deja parfaitement gelee), on ne passe pas par la
        // comparaison de moyennes ci-dessous. Elle serait sinon numeriquement instable dans
        // ce regime : deux valeurs bit-a-bit identiques repetees n1 puis n2 fois donnent des
        // moyennes flottantes legerement differentes (n1 != n2 => arrondis de sommation
        // differents), et l'erreur standard correspondante (variance residuelle de l'ordre
        // de l'epsilon machine au carre, inflatee par une auto-correlation g mal definie sur
        // un signal degenere) n'a plus le bon ordre de grandeur pour arbitrer entre les deux.
        double lo = q[plateauFrom], hi = q[plateauFrom];
        for (int i = plateauFrom + 1; i < to; i++) {
            double v = q[i];
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
        double scale = Math.max(1.0, Math.max(Math.abs(lo), Math.abs(hi)));
        if ((hi - lo) <= 1e-9 * scale) return true;

        int mid = plateauFrom + half;
        double meanFirst = EquilibrationDetector.mean(q, plateauFrom, mid);
        double meanSecond = EquilibrationDetector.mean(q, mid, to);
        double seFirst = EquilibrationDetector.standardError(q, plateauFrom, mid);
        double seSecond = EquilibrationDetector.standardError(q, mid, to);
        double se = Math.sqrt(seFirst * seFirst + seSecond * seSecond);
        // se peut etre NaN (valeur non finie dans q) : la comparaison est alors false,
        // ce qui refuse prudemment de declarer un plateau sur des donnees invalides.
        return Math.abs(meanFirst - meanSecond) <= se;
    }
}
