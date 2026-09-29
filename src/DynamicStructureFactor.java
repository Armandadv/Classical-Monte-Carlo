import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Locale;

import org.apache.avro.Schema;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;
import org.apache.commons.math3.util.FastMath;
import org.jtransforms.fft.DoubleFFT_1D;
import org.jtransforms.fft.DoubleFFT_2D;

/**
 * Structure dynamique de spin S(Q,ω). Pour chaque instantané S(t) reçu : FFT spatiale 2D par
 * sous-réseau, FFT temporelle avec fenêtre gaussienne par sous-réseau — une seule fois par
 * (qx, qy), les séries S0/S1 ne dépendant pas de (h, k) (linéarité de la FFT) —, combinaison
 * spectrale inter-sous-réseaux par le facteur de phase e^{-i2πQ·(r_B−r_A)}, projection sur le
 * tenseur g puis ⊥ au vecteur de diffusion (garde Q ≈ 0 : direction neutron indéfinie, S gardé
 * tel quel), accumulation en moyenne courante. L'Avro {@code structure_dipolaire_h<champ>.avro}
 * n'est écrit qu'au (nSQW+1)-ième appel de {@link #processSQW} — d'où la reconstruction de
 * l'instance à chaque bloc S(Q,w) (cf. {@code MonteCarlo.sqwBlock}) : le compteur interne ne se
 * remet pas à zéro. Les propriétés du schéma (N, Delta, Qzz, Window, DeltaOmega) tracent les paramètres de mesure.
 * Équivalent mathématiquement pour tout (Na, Nb) ; arrondis différents de l'implémentation
 * historique (~1e-13 rel., cf. BenchDSF).
 *
 * <p>Pairing des blocs miroirs : pour qzz rigoureusement nul, la relation exacte
 * I(2−h, 2−k, (Na−qx)%Na, (Nb−qy)%Nb, ω) = I(h, k, qx, qy, (time−ω)%time) (spins réels →
 * hermiticité spatiale : le facteur de phase et la projection ⊥Q se conjuguent cohéremment)
 * permet de ne calculer que les blocs (h, k) lexicographiquement maximaux et de remplir les
 * autres par les fréquences négatives de leur miroir (passe de remplissage en fin de
 * {@link #processSQW}). Deux restrictions validées :
 * (R1) les colonnes/lignes de Nyquist paires (qx = Na/2 si Na pair ; qy = Nb/2 si Nb pair)
 * sont exclues du pairing — g(Na/2) = +1/2 auto-miroir, ambiguïté ±1/2 mod G : la base
 * bi-atomique n'est pas périodique sous Q → Q+G — et restent calculées dans les 9 blocs ;
 * (R2) le pairing entier est désactivé si qzz ≠ 0 (composante z non mirrorée dans Qkitaev :
 * Q̂(−Q) ≠ −Q̂(Q) sinon ; option --qz d'AppParallel). Le pairing requiert qzz == 0
 * (composante z du vecteur de diffusion) ; deltaZ (offset z des sous-réseaux) est quelconque —
 * à qzz = 0 le terme de phase z s'annule exactement (0.0 × deltaZ = +0.0).</p>
 *
 * <p>FFT spatiale empaquetée (÷2 FFT) : les deux sous-réseaux d'une même composante partagent
 * un seul {@code complexForward} (partie réelle = sous-réseau 0, partie imaginaire = sous-réseau 1,
 * toutes deux réelles — {@code FFT(a + i·b) = FFT(a) + i·FFT(b)}) ; la séparation se fait au
 * store par paire de bins miroirs, le bin miroir écrit conjugué bit à bit. Miroir des FFT
 * temporelles (÷2 FFT) : les séries du point miroir (qxm, qym) étant les conjuguées exactes de
 * celles du point (hermiticité spatiale du packing), leur FFT temporelle est reconstruite sans
 * FFT par conjugué-renversé — seul le point (qx, qy) ≥lex son miroir est calculé. Équivalences
 * mathématiques vérifiées (arrondis différents, ~1e-13 rel., cf. BenchDSF).</p>
 *
 * <p>Zero-padding : interpole (les bins existants préservés), ne réduit pas le coût FFT —
 * outil d'analyse, pas d'optimisation.</p>
 *
 * <p>FFT spatiale réelle exclue (sonde RealFFTProbe) : realForwardFull/complexForward =
 * 1.01 (12×12), 0.87 (30×30), 0.47 (6×4), 0.79 (3×3) — pas de gain sur les tailles
 * non-puissances-de-2 ; le packed realForward est réservé aux puissances de 2. Décision
 * documentée pour ne pas y revenir.</p>
 */
public class DynamicStructureFactor {
    private double[][] S_Qw;
    private int totalQ;
    private final int time;
    private final int Na, Nb, H, K, numSpins;
    private double qzz;
    private int Nmesures;
    private final double delta; // pas d'échantillonnage temporel = dt*measurementRate (ħ/meV)
    private Lattice lattice;
    /**
     * Répertoire racine des fichiers S(Q,w). Défaut : répertoire courant. AppParallel le
     * redirige vers le répertoire de la réplique via {@link #setOutputDirPath(String)}.
     */
    private String outputDirPath = ".";

    /** Redirige la sortie S(Q,w). */
    public void setOutputDirPath(String dir) {
        if (dir == null || dir.isEmpty()) throw new IllegalArgumentException("outputDirPath vide");
        this.outputDirPath = dir.endsWith(File.separator) ? dir : dir + File.separator;
    }

    /**
     * Fixe la composante z du vecteur de diffusion, en r.l.u. (defaut 0 : honeycomb 2D, cf.
     * constructeur) ; enregistree dans l'Avro comme propriete "Qzz" (traceabilite de l'option
     * --qz d'AppParallel). A appeler AVANT la premiere accumulation : changer qzz en cours de
     * run basculerait le pairing (cf. javadoc de classe) et melangerait deux conventions dans
     * la moyenne courante.
     */
    public void setQzz(double qzz) {
        this.qzz = qzz;
    }

    public String getOutputDirPath() { return outputDirPath; }

    /** Lecture de l'accumulation pour BenchDSF uniquement. */
    double[][] getS_Qw() { return S_Qw; }

    /**
     * Bilans par étage de {@link #processSQW} pour BenchDSF (nanosecondes, remis à zéro à
     * chaque appel) : {@code spatialNanos} = boucle FFT2D packée + stores séparés (de la
     * préparation des buffers à la fin des stores fftSpatial) ; {@code temporalNanos} =
     * étape temporelle des points CALCULÉS seulement (les FFT1D — le miroir (qxm, qym) est
     * reconstruit sans FFT, cf. miroir des FFT temporelles) ; {@code projectionNanos} =
     * projection G hoistée + boucle (h, k) du point (combine + ⊥Q + moyenne) + reconstruction
     * conjugué-renversé du miroir + boucle (h, k) du miroir ; {@code fillNanos} = passe de
     * remplissage des blocs miroirs. Champs package-private lus directement par BenchDSF
     * (plus simple que des getters) ; ~2 appels {@code System.nanoTime} par étage et par
     * (qx, qy) — overhead négligeable.
     */
    long spatialNanos, temporalNanos, projectionNanos, fillNanos;

    public DynamicStructureFactor(int H, int K, int time, double delta, Lattice lattice) {
        this.lattice = lattice;
        this.Na = lattice.lengthA;
        this.Nb = lattice.lengthB;
        this.time = time;
        this.delta = delta;
        this.numSpins = lattice.size;
        this.H = H;
        this.K = K;
        // Reseau honeycomb 2D : pas de composante l du vecteur de diffusion.
        // (Ancienne valeur BaCoAsO : 4.67.)
        this.qzz = 0.d;
        this.totalQ = Na * H * Nb * K;
        this.S_Qw = new double[totalQ][time];
        for (int Q = 0; Q < totalQ; Q++) {
            Arrays.fill(S_Qw[Q], 0.0);
        }
        this.Nmesures = 0;
    }

    private double getQValue(int value, final int length) {
        return (value <= length / 2) ? value / (double) (length) : -(length - value) / (double) (length);
    }

    /**
     * Accumule un instantané S(t) dans la moyenne courante de S(Q,ω). L'écriture Avro n'a lieu
     * qu'au (nSQW+1)-ième appel (contrainte de {@code MonteCarlo.sqwBlock}).
     *
     * @param filein paramètre historique : n'entre plus dans le nom du fichier.
     */
    public void processSQW(UnitCell uc, double[] S_t_flat, final String filein, final int nSQW, final double B) {

        // Bilans par étage (BenchDSF) : reset à chaque appel.
        spatialNanos = 0L;
        temporalNanos = 0L;
        projectionNanos = 0L;
        fillNanos = 0L;

        final int time_flat = (int) S_t_flat.length / (numSpins * 3);
        if (time_flat != time) {
            throw new IllegalStateException("instantanés reçus (" + time_flat
                    + ") != fenêtre du constructeur (" + time + ")");
        }
        final double[][] Gtensor = lattice.getGtensor();
        final double factorX = 1.d / Math.sqrt(lattice.lengthA * lattice.lengthB);
        final double factorT = 1.d / Math.sqrt(time);

        /*
         * Base directe de la maille, en cartésien et en unité a : {a/a, b/a, c/a} — le
         * troisieme vecteur vaut c/a (c/a = 1 revient à un réseau strictement 2D, pas
         * d'anisotropie hors plan). Valeur BaCoAsO a T = 4.2 K (p87 Cyrille) :
         * c/a = 23.25/4.997 ≈ 4.6528, utilisée pour l'interprétation de --qz (r.l.u.
         * de c*). Convention identique au legacy (src/legacy/Spin_legacy_getCorrelation.txt,
         * lignes 25-28).
         */
        final double[][] basis = new double[][] { { 1.d, 0.d, 0.d }, { 0.d, 1.d, 0.d }, { 0.d, 0.d, 23.25/4.997 } };
        final double[][] kitaev_basis_reverse = MathOps.inverse(MathOps.mult2D(basis, lattice.getChangeOfBasis()));

        // Initialisation des objets FFT JTransform.
        // Convention DoubleFFT_2D(rows, cols) : rows = j (a2, Nb), cols = i (a1, Na) : élément
        // (r, c) à 2*(r*Na + c) ; le stockage transposé ci-dessous restaure [qx ↔ i/a1][qy ↔ j/a2].
        DoubleFFT_2D fft2D = new DoubleFFT_2D(Nb, Na);
        DoubleFFT_1D fft1D = new DoubleFFT_1D(time);

        // Pré-allocate buffer pour les FFT spatiales
        double[] fftBuffer = new double[Na * Nb * 2]; // Réel et imaginaire
        // Pré-allocate buffer pour les FFT temporelles
        double[] fftTimeBuffer = new double[time * 2]; // Réel et imaginaire

        // Allocation de la mémoire pour les résultats des FFT spatiales : quatre tableaux
        // primitifs (Re/Im par sous-réseau) au lieu d'objets Complex (empreinte mémoire
        // réduite) ; remplis par le packing ÷2 FFT ci-dessous — équivalence mathématique,
        // arrondis différents de l'implémentation historique (cf. BenchDSF ; sonde
        // PackingProbe : écart max ~6e-15 sur (6,4), (3,3), (12,12)).
        // +1 Uniquement si Na et Nb paires prend en compte que FFT(z=n/2) x(f = +1/(2
        // Delta),-1/(2 Delta))
        final double[][][][] fftSpatial0Re = new double[time][3][Na][Nb];
        final double[][][][] fftSpatial0Im = new double[time][3][Na][Nb];
        final double[][][][] fftSpatial1Re = new double[time][3][Na][Nb];
        final double[][][][] fftSpatial1Im = new double[time][3][Na][Nb];

        // Transformation de Fourier spatiale : PACKING ÷2 FFT. Les DEUX sous-réseaux d'une
        // même composante sont empaquetés dans un seul buffer complexe — partie réelle =
        // S(i, j, t) sous-réseau 0, partie imaginaire = S(i, j, t) sous-réseau 1, toutes deux
        // RÉELLES : FFT(a + i·b) = FFT(a) + i·FFT(b) pour a, b réels — un complexForward par
        // (t, comp) au lieu de deux. Lecture directe de S_t_flat : l'ancien
        // S_t_reshaped[t][site][xyz] valait exactement S_t_flat[(t * numSpins + site) * 3 + xyz]
        // (identité prouvée).
        // Étage spatial (bilan BenchDSF) : de la préparation des buffers à la fin des stores.
        final long tSpatial = System.nanoTime();
        for (int t = 0; t < time; t++) {
            for (int xyz = 0; xyz < 3; xyz++) {
                // Préparation du buffer packé : cellule (i, j) → position (row=j, col=i) ;
                // fftBuffer[2*(j*Na+i)] = sous-réseau 0 (réel), fftBuffer[2*(j*Na+i)+1] =
                // sous-réseau 1 (imaginaire).
                for (int j = 0; j < Nb; j++) {
                    for (int i = 0; i < Na; i++) {
                        final int site0 = lattice.siteIndex(i, j, 0);
                        final int site1 = lattice.siteIndex(i, j, 1);
                        fftBuffer[2 * (j * Na + i)] = S_t_flat[(t * numSpins + site0) * 3 + xyz];
                        fftBuffer[2 * (j * Na + i) + 1] = S_t_flat[(t * numSpins + site1) * 3 + xyz];
                    }
                }

                // FFT 2D in-place du buffer packé (les deux sous-réseaux d'un coup)
                fft2D.complexForward(fftBuffer);

                // Séparation au store, par PAIRE de bins miroirs (k, km) lue UNE fois (règle
                // lex (f2, f1) ≥ (f2m, f1m) — couvre chaque paire et les auto-miroirs une
                // fois, les DEUX bins écrits à chaque fois). Z = buffer packé aux deux bins :
                // F0(k) = (Z(k) + conj(Z(km)))/2 (sous-réseau 0) et F1(k) = (Z(k) −
                // conj(Z(km)))/(2i) (sous-réseau 1) — formules validées empiriquement par la
                // sonde PackingProbe. Lecture (f1=y, f2=x) : écrite en [f2][f1] — transposée
                // pour [qx↔i][qy↔j]. factorX s'applique aux 4 sorties.
                for (int f1 = 0; f1 < Nb; f1++) {
                    final int f1m = (Nb - f1) % Nb;
                    for (int f2 = 0; f2 < Na; f2++) {
                        final int f2m = (Na - f2) % Na;
                        if (f2 < f2m || (f2 == f2m && f1 < f1m)) continue; // paire déjà traitée
                        final int kBin = 2 * (f1 * Na + f2);
                        final int kmBin = 2 * (f1m * Na + f2m);
                        final double zreK = fftBuffer[kBin];
                        final double zimK = fftBuffer[kBin + 1];
                        final double zreM = fftBuffer[kmBin];
                        final double zimM = fftBuffer[kmBin + 1];
                        final double f0re = (zreK + zreM) * 0.5;
                        final double f0im = (zimK - zimM) * 0.5;
                        final double f1re = (zimK + zimM) * 0.5;
                        final double f1im = (zreM - zreK) * 0.5;
                        fftSpatial0Re[t][xyz][f2][f1] = f0re * factorX;
                        fftSpatial0Im[t][xyz][f2][f1] = f0im * factorX;
                        fftSpatial1Re[t][xyz][f2][f1] = f1re * factorX;
                        fftSpatial1Im[t][xyz][f2][f1] = f1im * factorX;
                        if (f2m != f2 || f1m != f1) {
                            // Miroir distinct : conjugué EXACT bit à bit (A(−k) = conj(A(k)),
                            // par construction de la séparation).
                            fftSpatial0Re[t][xyz][f2m][f1m] = f0re * factorX;
                            fftSpatial0Im[t][xyz][f2m][f1m] = -(f0im * factorX);
                            fftSpatial1Re[t][xyz][f2m][f1m] = f1re * factorX;
                            fftSpatial1Im[t][xyz][f2m][f1m] = -(f1im * factorX);
                        }
                    }
                }
            }
        }
        spatialNanos += System.nanoTime() - tSpatial;

        // Allocation de la mémoire pour les données combinées après application du
        // facteur de phase
        final int NaNb = Na * Nb;

        // delta = position(B) − position(A) en coordonnées fractionnaires (le code historique
        // prenait B direct, correct seulement si A = origine ; z = 0.1733… − 0.16 = 0.0133…
        // sur la fixture BaCoAsO). À qzz = 0 (défaut), qzz·deltaZ = +0.0 exact : chemin
        // bit-identique au 2D historique.
        final double deltaX = uc.basis.get(1)[0] - uc.basis.get(0)[0];
        final double deltaY = uc.basis.get(1)[1] - uc.basis.get(0)[1];
        final double deltaZ = uc.basis.get(1)[2] - uc.basis.get(0)[2];

        // Tableaux spectraux hoistés hors des boucles (qx, qy) — alloués UNE fois par
        // processSQW, intégralement réécrits à chaque usage (aucune valeur n'est lue avant
        // d'avoir été réécrite dans l'itération courante) :
        // - x0/x1 : FFT temporelle de chaque sous-réseau (les séries S0/S1 ne dépendent pas
        //   de (h, k) → une FFT temporelle par (qx, qy) au lieu de H·K, cf. linéarité) ;
        // - g0/g1 : projection G de ces FFT (G est la même matrice pour tous les Q → hoistée
        //   elle aussi, cf. linéarité) ;
        // - g0m/g1m : projection G du point MIROIR (qxm, qym), reconstruite par
        //   conjugué-renversé depuis g0/g1 (pas de FFT ni de mat-vec pour le miroir, cf.
        //   commentaire d'ouverture de la boucle (qx, qy)) — alloués aussi UNE fois ;
        // - reS/imS/reS_perp/imS_perp : vecteurs de travail du bin ω (ex-alloués par Q/ω).
        final double[][] x0Re = new double[3][time];
        final double[][] x0Im = new double[3][time];
        final double[][] x1Re = new double[3][time];
        final double[][] x1Im = new double[3][time];
        final double[][] g0Re = new double[3][time];
        final double[][] g0Im = new double[3][time];
        final double[][] g1Re = new double[3][time];
        final double[][] g1Im = new double[3][time];
        final double[][] g0mRe = new double[3][time];
        final double[][] g0mIm = new double[3][time];
        final double[][] g1mRe = new double[3][time];
        final double[][] g1mIm = new double[3][time];
        final double[] reS = new double[3];
        final double[] imS = new double[3];
        final double[] reS_perp = new double[3];
        final double[] imS_perp = new double[3];

        // ---- Pairing des blocs miroirs : précalcul -----------------------------------------
        // Le miroir exige la composante z du vecteur de diffusion rigoureusement nulle :
        // Q̂(−Q) ≠ −Q̂(Q) sinon (composante z non mirrorée dans Qkitaev) — un --qz non nul
        // désactive le pairing entier (R2).
        final boolean pairingActif = (qzz == 0.0);
        // (R1) Colonnes/lignes de Nyquist paires exclues du pairing : g(Na/2) = +1/2 est
        // auto-miroir en indice mais ambigu ±1/2 mod G — la base bi-atomique n'est pas
        // périodique sous Q → Q+G. Ces points restent calculés dans les 9 blocs.
        final boolean[] nyqX = new boolean[Na];
        for (int qx = 0; qx < Na; qx++) {
            nyqX[qx] = (Na % 2 == 0 && qx == Na / 2);
        }
        final boolean[] nyqY = new boolean[Nb];
        for (int qy = 0; qy < Nb; qy++) {
            nyqY[qy] = (Nb % 2 == 0 && qy == Nb / 2);
        }
        // Blocs (h, k) remplis = miroir (2−h, 2−k) dans la grille [0,H)×[0,K) ET (h, k) lex
        // strictement plus petit que son miroir (les blocs maximaux sont calculés, les
        // minimaux remplis par leur miroir). À H=K=3 : remplis {(0,0),(0,1),(0,2),(1,0)}.
        // Garantie du prédicat : hc = 2−h ∈ [0,H) et kc = 2−k ∈ [0,K) quand il rend vrai,
        // et le miroir d'un bloc rempli est toujours calculé (jamais rempli lui-même).
        final boolean[][] blocRempli = new boolean[H][K];
        for (int h = 0; h < H; h++) {
            for (int k = 0; k < K; k++) {
                blocRempli[h][k] = blocRempli(h, k, H, K);
            }
        }

        // Transformation de Fourier temporelle et calcul de S(Q, ω). Boucles (qx, qy)
        // externes, boucles (h, k) internes : l'étape temporelle ne dépend pas de (h, k).
        //
        // MIROIR DES FFT TEMPORELLES (÷2 FFT) : pour des spins réels, les séries spatiales du
        // point miroir (qxm, qym) = ((Na−qx)%Na, (Nb−qy)%Nb) sont les CONJUGUÉES EXACTES de
        // celles du point (qx, qy) — hermiticité spatiale de la FFT packée : les bins miroirs
        // sont écrits conjugués bit à bit (étage spatial) ; la fenêtre gaussienne, réelle
        // diagonale, commute avec la conjugaison (bit à bit, sonde MirrorProbe). La FFT
        // temporelle du miroir se reconstruit donc SANS FFT par conjugué-renversé :
        // X_m(ω) = conj(X((time−ω)%time)) — sonde MirrorProbe : exact (bit-exact à time=8,
        // ~8e-15 à time=100). G est réel → g0m/g1m = conjugués-renversés de g0/g1 (le produit
        // matrice-vecteur réel commute avec la conjugaison, bit à bit) : ni x0m/x1m ni mat-vec
        // pour le miroir. Règle de calcul : le point (qx, qy) est CALCULÉ ssi (qx, qy) ≥lex
        // (qxm, qym) — règle lex correcte pour toute parité (PAS de qx > Na/2) ; l'involution
        // (miroir du miroir = le point) garantit que chaque point est calculé exactement une
        // fois, directement ou par son miroir. Si le point est strictement ≥lex son miroir
        // distinct : l'étape temporelle du point UNE fois, puis SA boucle (h, k) avec ses
        // g0/g1, puis la reconstruction et la boucle (h, k) du miroir avec ses qxxm/qyym, son
        // facteur de phase et son Qkitaev propres (via getQValue, pas de formule manuelle).
        // L'ordre des Q accumulés change (les Q du miroir intercalés) : sans effet, S_Qw[Q]
        // est indépendant par Q. Le pairing des blocs (h, k) lit S_Qw fini et n'interagit pas
        // (nyqX[qxm] == nyqX[qx] : même motif de skip des deux côtés).
        for (int qx = 0; qx < Na; qx++) {
            double qxx = getQValue(qx, Na) - 1.d; // shift en X pour centré le 0 à l'origine
            final int qxm = (Na - qx) % Na;
            final double qxxm = getQValue(qxm, Na) - 1.d;
            for (int qy = 0; qy < Nb; qy++) {
                final int qym = (Nb - qy) % Nb;
                final boolean miroirDistinct = (qxm != qx) || (qym != qy);
                if (miroirDistinct && (qx < qxm || (qx == qxm && qy < qym))) {
                    continue; // lex-minimal : ce point sera calculé par son miroir
                }
                double qyy = getQValue(qy, Nb) - 1.d;
                final double qyym = getQValue(qym, Nb) - 1.d;

                // ---- Étape temporelle : UNE fois par point calculé (le miroir est reconstruit)
                final long tTemporal = System.nanoTime();
                for (int comp = 0; comp < 3; comp++) {
                    // Sous-réseau 0 : série temporelle → fenêtre → FFT → x0Re/x0Im
                    for (int t = 0; t < time; t++) {
                        fftTimeBuffer[2 * t] = fftSpatial0Re[t][comp][qx][qy];
                        fftTimeBuffer[2 * t + 1] = fftSpatial0Im[t][comp][qx][qy];
                    }
                    // Fenetrage avant la FFT temporelle (sans objet pour time == 1 :
                    // en mode statique la fenetre ecraserait l'unique echantillon) ;
                    // fenêtre réelle diagonale, passe à travers la combinaison spectrale.
                    if (time > 1) {
                        MathOps.applyGaussianWindow(fftTimeBuffer, time);
                    }
                    // FFT 1D in-place sur fftTimeBuffer
                    fft1D.complexForward(fftTimeBuffer);
                    for (int omega = 0; omega < time; omega++) {
                        x0Re[comp][omega] = fftTimeBuffer[2 * omega] * factorT;
                        x0Im[comp][omega] = fftTimeBuffer[2 * omega + 1] * factorT;
                    }
                    // Sous-réseau 1 : idem depuis fftSpatial1Re/Im → x1Re/x1Im
                    for (int t = 0; t < time; t++) {
                        fftTimeBuffer[2 * t] = fftSpatial1Re[t][comp][qx][qy];
                        fftTimeBuffer[2 * t + 1] = fftSpatial1Im[t][comp][qx][qy];
                    }
                    if (time > 1) {
                        MathOps.applyGaussianWindow(fftTimeBuffer, time);
                    }
                    fft1D.complexForward(fftTimeBuffer);
                    for (int omega = 0; omega < time; omega++) {
                        x1Re[comp][omega] = fftTimeBuffer[2 * omega] * factorT;
                        x1Im[comp][omega] = fftTimeBuffer[2 * omega + 1] * factorT;
                    }
                }
                temporalNanos += System.nanoTime() - tTemporal;

                // ---- Projection G hoistée : G est la même matrice pour tous les Q, donc
                // G·(X0 + φ·X1) = G·X0 + φ·(G·X1) — même linéarité que la FFT temporelle.
                // Les produits matrice-vecteur sortent de la boucle (h, k) : 4 par (qx, qy, ω)
                // au lieu de 2 par (Q, ω), et zéro allocation dans les boucles internes.
                // Chronométrée dans l'étage projection (bilan BenchDSF) : c'est une projection,
                // et chez old le mat-vec G est dans la boucle ω (comparaison cohérente).
                final long tProjection = System.nanoTime();
                for (int comp = 0; comp < 3; comp++) {
                    for (int omega = 0; omega < time; omega++) {
                        g0Re[comp][omega] = Gtensor[comp][0] * x0Re[0][omega]
                                + Gtensor[comp][1] * x0Re[1][omega] + Gtensor[comp][2] * x0Re[2][omega];
                        g0Im[comp][omega] = Gtensor[comp][0] * x0Im[0][omega]
                                + Gtensor[comp][1] * x0Im[1][omega] + Gtensor[comp][2] * x0Im[2][omega];
                        g1Re[comp][omega] = Gtensor[comp][0] * x1Re[0][omega]
                                + Gtensor[comp][1] * x1Re[1][omega] + Gtensor[comp][2] * x1Re[2][omega];
                        g1Im[comp][omega] = Gtensor[comp][0] * x1Im[0][omega]
                                + Gtensor[comp][1] * x1Im[1][omega] + Gtensor[comp][2] * x1Im[2][omega];
                    }
                }

                // ---- Combinaison spectrale des projections G + projection par Q : boucle
                // (h, k) du point courant (extraite dans accumuleBlocs, réutilisée telle
                // quelle pour le miroir).
                accumuleBlocs(qx, qy, qxx, qyy, g0Re, g0Im, g1Re, g1Im, kitaev_basis_reverse,
                        deltaX, deltaY, deltaZ, pairingActif, blocRempli, nyqX, nyqY,
                        reS, imS, reS_perp, imS_perp);

                if (miroirDistinct) {
                    // Reconstruction du miroir par conjugué-renversé (PAS de FFT) : cf.
                    // commentaire d'ouverture ; chronométrée dans l'étage projection.
                    for (int comp = 0; comp < 3; comp++) {
                        for (int omega = 0; omega < time; omega++) {
                            final int om = (time - omega) % time;
                            g0mRe[comp][omega] = g0Re[comp][om];
                            g0mIm[comp][omega] = -g0Im[comp][om];
                            g1mRe[comp][omega] = g1Re[comp][om];
                            g1mIm[comp][omega] = -g1Im[comp][om];
                        }
                    }
                    accumuleBlocs(qxm, qym, qxxm, qyym, g0mRe, g0mIm, g1mRe, g1mIm,
                            kitaev_basis_reverse, deltaX, deltaY, deltaZ, pairingActif,
                            blocRempli, nyqX, nyqY, reS, imS, reS_perp, imS_perp);
                }
                projectionNanos += System.nanoTime() - tProjection;
            }
        }

        // ---- Passe de remplissage des blocs miroirs ----------------------------------------
        // I(2−h, 2−k, (Na−qx)%Na, (Nb−qy)%Nb, ω) = I(h, k, qx, qy, (time−ω)%time) — exact
        // (spins réels → hermiticité spatiale : le facteur de phase et la projection ⊥Q se
        // conjuguent cohéremment) ; points de Nyquist pairs exclus (R1). Le miroir
        // (hc, kc) = (2−h, 2−k) d'un bloc rempli est dans la grille (garanti par le
        // précalcul) et toujours calculé (jamais rempli lui-même : lex plus grand). La copie
        // écrase la moyenne courante du bloc rempli par celle de la source — identique, car
        // la relation vaut par instantané et la moyenne est linéaire. Points de Nyquist des
        // blocs remplis : calculés normalement ci-dessus, non écrasés ici.
        // Bracket du timer À L'INTÉRIEUR du if : fillNanos == 0 ssi pairing désactivé (R2),
        // sémantique exploitée par le check de garde qzz.
        if (pairingActif) {
            final long tFill = System.nanoTime();
            for (int h = 0; h < H; h++) {
                for (int k = 0; k < K; k++) {
                    if (!blocRempli[h][k]) continue;
                    final int hc = 2 - h, kc = 2 - k;
                    for (int qx = 0; qx < Na; qx++) {
                        if (nyqX[qx]) continue;
                        final int qxm = (Na - qx) % Na;
                        for (int qy = 0; qy < Nb; qy++) {
                            if (nyqY[qy]) continue;
                            final int qym = (Nb - qy) % Nb;
                            final int Qcible = (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);
                            final int Qsource = (hc * NaNb * K) + (kc * NaNb) + (qxm * Nb + qym);
                            for (int omega = 0; omega < time; omega++) {
                                S_Qw[Qcible][omega] = S_Qw[Qsource][(time - omega) % time];
                            }
                        }
                    }
                }
            }
            fillNanos += System.nanoTime() - tFill;
        }

        Nmesures += 1;

        if (Nmesures > nSQW) {
            // Matrice changement de base a,b vers x,y une BON
            double[][] P_xy = { { 1.,0.}, { 1./Math.sqrt(3.), 2./Math.sqrt(3.) } };
            try {

               // Résolution réelle par bin (hbarre = 1), calculée AVANT le schéma : celui-ci
               // la référence (propriété DeltaOmega), avec la durée intégrée effective
               // Window = time*delta (= W après l'arrondi de ndt dans sqwBlock).
               final double omegaPerBin = Math.TAU / (time * delta); // meV par bin (hbarre = 1)
                Schema schema = new Schema.Parser().parse("{" +
                        "  \"type\": \"record\"," +
                        "  \"name\": \"SQW\"," +
                        "  \"N\": " + time + "," +
                        "  \"Delta\": " + delta + "," +
                        "  \"Qzz\": " + qzz + "," +
                        "  \"Window\": " + (time * delta) + "," +
                        "  \"DeltaOmega\": " + omegaPerBin + "," +
                        "  \"fields\": [" +
                        "    {\"name\": \"k_x\", \"type\": \"float\"}," +
                        "    {\"name\": \"k_y\", \"type\": \"float\"}," +
                        "    {\"name\": \"w\", \"type\": \"int\"}," +
                        "    {\"name\": \"w_meV\", \"type\": \"float\"}," +
                        "    {\"name\": \"S_QW\", \"type\": \"float\"}" +
                        "  ]" +
                        "}");

                // filein n'entre plus dans le nom : le repertoire porte l'identite (fusion nommage 26/09/2026).
                File filePath = new File(outputDirPath);
                File file = new File(filePath,
                        "structure_dipolaire_h" + String.format(Locale.US, "%.6f", B) + ".avro");

                if (!filePath.exists()) {
                    filePath.mkdirs();  // Créer le répertoire s'il n'existe pas
                }

               DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(schema);
               DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter);
               dataFileWriter.setCodec(CodecFactory.snappyCodec());
               dataFileWriter.create(schema, file);
               for (int h = 0; h < H; h++) {
                   double hklx = h;
                   for (int k = 0; k < K; k++) {
                       double hkly = k;
                       for (int qx = 0; qx < Na; qx++) {
                           double qxx = getQValue(qx, Na) -1.d;
                           for (int qy = 0; qy < Nb; qy++) {
                               int Q = (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);
                               double qyy = getQValue(qy, Nb) -1.d;
                               double[] Q_xy = MathOps.matrixXvector(P_xy, new double[] { qxx + hklx, qyy + hkly });

                               for (int w = 0; w <= time / 2; w++) {   // on jette les energies négative
                                   GenericRecord record = new GenericData.Record(schema);
                                   record.put("k_x", (float) Q_xy[0]);
                                   record.put("k_y", (float) Q_xy[1]);
                                   record.put("w", (int) w);
                                   record.put("w_meV", (float) (omegaPerBin * w));
                                   record.put("S_QW", (float) S_Qw[Q][w]);
                                   dataFileWriter.append(record);
                               }
                           }
                       }
                   }
               }
               dataFileWriter.close();
           } catch (IOException e) {
               e.printStackTrace();
           }
       }
    }

    /**
     * Boucle (h, k) d'un point (qx, qy) : combinaison spectrale des projections G hoistées
     * (g0/g1) avec le facteur de phase e^{-i2πQ·δ}, projection ⊥ Qkitaev (garde Q ≈ 0 :
     * direction neutron indéfinie, S gardé tel quel) et moyenne courante dans S_Qw. Extraite
     * de processSQW pour être réutilisée TELLE QUELLE par le point (qx, qy) et par son miroir
     * temporel (qxm, qym) — chacun avec ses propres g0/g1 (calculés ou reconstruits
     * conjugués-renversés) et ses propres qxx/qyy (via getQValue), d'où facteur de phase et
     * Qkitaev corrects pour les deux. Les blocs de pairing (blocRempli, hors points de
     * Nyquist pairs R1) sont sautés : remplis par la passe de remplissage en fin de
     * processSQW — même motif de skip pour le point et son miroir (nyqX[qxm] == nyqX[qx]).
     * Zéro allocation en dehors des deux vecteurs Qkitaev par bloc (comme avant l'extraction).
     */
    private void accumuleBlocs(final int qx, final int qy, final double qxx, final double qyy,
            final double[][] g0Re, final double[][] g0Im, final double[][] g1Re, final double[][] g1Im,
            final double[][] kitaev_basis_reverse, final double deltaX, final double deltaY,
            final double deltaZ, final boolean pairingActif, final boolean[][] blocRempli,
            final boolean[] nyqX, final boolean[] nyqY,
            final double[] reS, final double[] imS, final double[] reS_perp, final double[] imS_perp) {
        final int NaNb = Na * Nb;
        for (int h = 0; h < H; h++) {
            double hklx = h;
            for (int k = 0; k < K; k++) {
                double hkly = k;
                // Pairing des blocs miroirs : le bin (Q, ω) entier est skippé — il
                // sera rempli par la passe de remplissage (fréquences négatives de
                // son miroir). Blocs non remplis et points de Nyquist pairs (R1) :
                // calculés normalement.
                if (pairingActif && blocRempli[h][k] && !nyqX[qx] && !nyqY[qy]) continue;
                // Indice Q
                int Q = (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);

                // Calcul du facteur de phase : e^{-i(2PI*Q.r_alpha)} exprimé en r.l.u.
                double Qr = (qxx + hklx) * deltaX + (qyy + hkly) * deltaY + qzz * deltaZ;
                // e^{iθ} avec θ = -TAU*Qr, calculé directement — identique au bit à
                // l'ancien new Complex(0, θ).exp() (FastMath.exp(0) = 1.0, ×1.0 exact).
                final double theta = -Math.TAU * Qr;
                final double pfRe = FastMath.cos(theta);
                final double pfIm = FastMath.sin(theta);

                // Pré-calcul QKitaev (dans la base Kitaev) et sa norme pour ce Q
                double[] Qkitaev = MathOps.matrixXvector(kitaev_basis_reverse,
                        new double[] { qxx + hklx, qyy + hkly, qzz});
                // (ancienne conversion ×23.25/4.997, déplacée INVERSEMENT dans basis ligne ~165
                // via l'inverse de la base — NE PAS réactiver ici : double comptage (c/a)²)
                double QKitaev_norm_sq = Qkitaev[0] * Qkitaev[0] + Qkitaev[1] * Qkitaev[1]
                        + Qkitaev[2] * Qkitaev[2];

                // Combinaison spectrale des projections G (G·X0 + φ·G·X1 — linéarité,
                // équivalence mathématique, arrondis différents de l'implémentation
                // historique, cf. BenchDSF), puis ⊥Q — seule la projection perpendiculaire
                // dépend du point Q (la direction Q̂ change).
                for (int omega = 0; omega < time; omega++) {

                    // Sx(omega) = reS[0] + i imS[0], etc. — vecteur combiné déjà projeté G
                    reS[0] = g0Re[0][omega] + (g1Re[0][omega] * pfRe - g1Im[0][omega] * pfIm);
                    imS[0] = g0Im[0][omega] + (g1Im[0][omega] * pfRe + g1Re[0][omega] * pfIm);
                    reS[1] = g0Re[1][omega] + (g1Re[1][omega] * pfRe - g1Im[1][omega] * pfIm);
                    imS[1] = g0Im[1][omega] + (g1Im[1][omega] * pfRe + g1Re[1][omega] * pfIm);
                    reS[2] = g0Re[2][omega] + (g1Re[2][omega] * pfRe - g1Im[2][omega] * pfIm);
                    imS[2] = g0Im[2][omega] + (g1Im[2][omega] * pfRe + g1Re[2][omega] * pfIm);

                    /*
                     * *****************************Neutron part ****************************************
                     */

                    // Calcul du produit scalaire pour la partie réelle
                    double dotRe = reS[0] * Qkitaev[0] + reS[1] * Qkitaev[1] + reS[2] * Qkitaev[2];
                    // Garde-fou qzz = 0 (reseau 2D) : au bin Gamma, Q = 0 exactement et
                    // la projection perpendiculaire diviserait par |Q|^2 = 0 (0/0 = NaN).
                    // Convention Q = 0 : direction neutron indefinie, on garde S tel quel.
                    if (QKitaev_norm_sq > 1e-12) {
                        reS_perp[0] = reS[0] - (dotRe / QKitaev_norm_sq) * Qkitaev[0];
                        reS_perp[1] = reS[1] - (dotRe / QKitaev_norm_sq) * Qkitaev[1];
                        reS_perp[2] = reS[2] - (dotRe / QKitaev_norm_sq) * Qkitaev[2];
                    } else {
                        reS_perp[0] = reS[0];
                        reS_perp[1] = reS[1];
                        reS_perp[2] = reS[2];
                    }

                    // Calcul du produit scalaire pour la partie imaginaire
                    double dotIm = imS[0] * Qkitaev[0] + imS[1] * Qkitaev[1] + imS[2] * Qkitaev[2];
                    if (QKitaev_norm_sq > 1e-12) {
                        imS_perp[0] = imS[0] - (dotIm / QKitaev_norm_sq) * Qkitaev[0];
                        imS_perp[1] = imS[1] - (dotIm / QKitaev_norm_sq) * Qkitaev[1];
                        imS_perp[2] = imS[2] - (dotIm / QKitaev_norm_sq) * Qkitaev[2];
                    } else {
                        imS_perp[0] = imS[0];
                        imS_perp[1] = imS[1];
                        imS_perp[2] = imS[2];
                    }

                    /*
                     * projection des ODS perp au vecteur de diffusion
                     */

                    double perp_mag_sq = 0.0;
                    for (int i = 0; i < 3; i++) {
                        perp_mag_sq += (reS_perp[i] * reS_perp[i] + imS_perp[i] * imS_perp[i]);
                    }

                    double sum = perp_mag_sq / numSpins;

                    S_Qw[Q][omega] = (Nmesures * S_Qw[Q][omega] + sum) / (Nmesures + 1);

                }
            }
        }
    }

    /**
     * Bloc (h, k) « rempli » par pairing miroir : vrai ssi le miroir (2−h, 2−k) est dans la
     * grille [0,H)×[0,K) ET (h, k) est lexicographiquement strictement plus petit que son
     * miroir (les blocs maximaux sont calculés, les minimaux remplis par leur miroir) ;
     * faux sinon — miroir hors grille (bloc calculé : rien à copier) ou auto-miroir
     * (h, k) = (2−h, 2−k) (bloc calculé). Garantit que hc = 2−h ∈ [0,H) et kc = 2−k ∈ [0,K)
     * quand le résultat est vrai.
     */
    private static boolean blocRempli(final int h, final int k, final int H, final int K) {
        final int hc = 2 - h;
        final int kc = 2 - k;
        if (hc < 0 || hc >= H || kc < 0 || kc >= K) {
            return false; // miroir hors grille : bloc calculé
        }
        return (h < hc) || (h == hc && k < kc);
    }

}