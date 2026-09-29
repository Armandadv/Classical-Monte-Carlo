import java.io.File;
import java.io.IOException;

import org.apache.avro.Schema;
import org.apache.avro.file.CodecFactory;
import org.apache.avro.file.DataFileWriter;
import org.apache.avro.generic.GenericData;
import org.apache.avro.generic.GenericDatumWriter;
import org.apache.avro.generic.GenericRecord;
import org.apache.avro.io.DatumWriter;
import org.jtransforms.fft.DoubleFFT_2D;

/**
 * Facteur de structure de DIMERES (canal de liaison) :
 *
 * <pre>   S_dimer(Q) = (1/N) * &lt; | somme_b (tau_b - tauMean) * e^{-i Q . R_b} |^2 &gt;</pre>
 *
 * avec tau_b = (S_i^g * S_j^g)^2 la variable de dimere de la liaison b = (i, j)_g
 * (g = etiquette Kitaev de la liaison), placee au MILIEU de la liaison R_b, et N
 * le nombre de sites. tauMean = moyenne de tau_b sur toutes les liaisons
 * (3*Na*Nb) de LA configuration : le fond uniforme de densite de dimeres est
 * retire AVANT la FFT, le correlateur est donc CONNECTE. Reference pour la
 * convention (correlateur connecte du champ de dimeres, trois orientations de
 * liaisons, normalisation 1/N) : Eq. (2) de Z. Yan et al., Phys. Rev. Lett.
 * (2022), arXiv:2202.11100 ; cf. aussi DimerSFRecompute pour le recalcul
 * standalone a partir d'un instantane spins_T*.avro. C'est la variable dans laquelle le pinch point du liquide
 * de dimeres (point RK de l'ensemble KBQ b pur) est le plus net : contrairement
 * au canal couleur {@link ColorStructureFactor} (densite de couleur sur les
 * SITES), il n'y a pas de termes intra-sous-reseau qui diluent l'anisotropie.
 *
 * <h2>Mise en oeuvre</h2>
 * Chaque site A (sous-reseau 0) porte exactement une liaison de chaque couleur :
 * tau_g est rangee sur la grille de la cellule du site A, une grille par couleur.
 * La FFT de chaque grille donne tau_g(q) aux positions des sites A ; le milieu de
 * liaison ajoute un facteur de phase par couleur :
 *
 * <pre>   phase_g = exp( -i TAU * ( a_idx * (fx + dx_g)/2 + b_idx * (fy + dy_g)/2 ) )</pre>
 *
 * ou (fx, fy) = position fractionnaire du site B dans la maille ({@code uc.basis.get(1)}),
 * (dx_g, dy_g) = decalage de cellule de la liaison g, et (a_idx, b_idx) = (qxx + h,
 * qyy + k) en r.l.u. — meme convention -TAU que le DSF. Les trois couleurs sont
 * combinees COHEREMMENT (une seule somme complexe), puis l'intensite est moyenee
 * sur les nSQW+1 configurations (moyenne glissante) et ecrite en Avro a la
 * derniere temperature (cf. MonteCarlo.sqwBlock).
 *
 * <h2>Conventions partagees avec ColorStructureFactor / DSF</h2>
 * Ordre vrai des sites : site = 2*(a + Na*b) + sr (cf. Lattice.siteIndex, a
 * rapide suivant A1, b lent suivant A2) ; buffer FFT en Nb lignes (b) x Na
 * colonnes (a), DoubleFFT_2D(Nb, Na), bin (qx, qy) relu a (qy*Na + qx) — correct
 * pour Na et Nb quelconques. Indice Q lineaire identique au DSF. Pas de
 * normalisation sqrt(Na*Nb) (echelle propre, non comparable au DSF). Pas de
 * composante z (reseau honeycomb 2D, qzz = deltaZ = 0 ; le delta 3D généralisé vit dans
 * DynamicStructureFactor).
 *
 * <h2>Ancre de validation physique</h2>
 * Sur la variete des appariements parfaits (T -&gt; 0, b pur), somme_b tau_b = N/2
 * exactement POUR CHAQUE configuration, donc tauMean = (N/2)/(3N/2) = 1/3.
 * La soustraction annule la composante uniforme TOTALE du champ : S_dimer = 0
 * exactement aux points ou les trois phases de milieu de liaison valent 1, en
 * pratique (h,k,qx,qy) = (1,1,0,0) c'est-a-dire aIdx = bIdx = 0 (au lieu de
 * (N/2)^2 / N = N/4 = 288 pour N = 1152 sans soustraction). Aux AUTRES points
 * Gamma replis (aIdx, bIdx entiers mais phases -1 selon la couleur), il reste
 * un signal d'imbalance entre couleurs — pas le fond de densite. Hors points
 * Gamma replies (aIdx ou bIdx non entier), la somme de phase annule le terme
 * uniforme et le spectre est identique a l'ancienne convention. Les
 * extinctions systematiques en M (zeros du facteur de forme) sont attendues,
 * comme dans le canal couleur.
 */
public class DimerStructureFactor {

    // Liaisons premiers voisins du site A (sous-reseau 0) : decalage de cellule
    // (dx, dy) et couleur Kitaev g (0 = x, 1 = y, 2 = z). Convention de l'ordre
    // des liaisons du fichier d'entree (formes A/B/C -> y/x/z), identique a
    // QUARTIC_GAMMA_BY_FORM du generateur Python : (-1,-1)->A->y, (-1,0)->B->x,
    // (0,0)->C->z.
    private static final int[][] BONDS = {
            { 0, 0, 2 },      // ( 0, 0) -> C -> z
            { -1, 0, 0 },     // (-1, 0) -> B -> x
            { -1, -1, 1 },    // (-1,-1) -> A -> y
    };

    private final int H, K, Na, Nb, NaNb, numSpins, totalQ;
    private final double[] S_dimer;      // intensite moyennee par Q
    private int Nmesures = 0;
    private String outputDirPath = null; // null = repertoire courant (writeAvro)
    private double qzz = 0.d;            // composante z du vecteur de diffusion (r.l.u.) ; traceabilite seule (honeycomb 2D)

    /**
     * Alloue les accumulateurs pour une grille de H x K zones : dimensions Na, Nb et
     * nombre de sites lus dans le lattice ; S_dimer(Q) part de zero.
     */
    public DimerStructureFactor(int H, int K, Lattice lattice) {
        this.H = H;
        this.K = K;
        this.Na = lattice.lengthA;
        this.Nb = lattice.lengthB;
        this.NaNb = Na * Nb;
        this.numSpins = lattice.size;
        this.totalQ = NaNb * H * K;
        this.S_dimer = new double[totalQ];
    }

    /**
     * Index lineaire Q d'un point (h, k, qx, qy), identique au DSF et utilise a
     * l'ecriture Avro : source unique, plus aucune formule dupliquee.
     */
    private int qIndex(int h, int k, int qx, int qy) {
        return (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);
    }

    /** Redirige la sortie Avro (rejette null et vide, comme le DSF). */
    public void setOutputDirPath(String dir) {
        if (dir == null || dir.isEmpty()) throw new IllegalArgumentException("outputDirPath vide");
        this.outputDirPath = dir.endsWith(File.separator) ? dir : dir + File.separator;
    }

    /**
     * Fixe la composante z du vecteur de diffusion (r.l.u.). Traceabilite seule : inertie
     * physique (deltaZ = 0) hors DynamicStructureFactor — le delta 3D généralisé (z des
     * sous-réseaux lu depuis l'input) vit dans DynamicStructureFactor ; la propriete Qzz des
     * Avro reflète l'option --qz d'AppParallel.
     */
    public void setQzz(double qzz) {
        this.qzz = qzz;
    }

    /**
     * Coordonnee reelle d'un bin FFT (r.l.u.) : value/length si value &lt;= length/2,
     * sinon -(length - value)/length, soit [-0.5, 0.5].
     */
    private double getQValue(int value, final int length) {
        return (value <= length / 2) ? value / (double) (length) : -(length - value) / (double) (length);
    }

    /**
     * Accumule une configuration : tau_g = (S_A^g * S_B^g)^2 sur les 3 grilles de
     * liaisons (cellule du site A), soustraction du fond tauMean (correlateur
     * connecte), FFT par couleur, combinaison coherente avec les phases de
     * milieu de liaison, moyenne glissante.
     *
     * @param uc        maille (position fractionnaire du site B : basis.get(1))
     * @param S_t_flat  spins aplaties ; seuls les N premiers triplets (premier
     *                  instantane) sont lus.
     */
    public void processConfig(UnitCell uc, double[] S_t_flat) {

        final double fx = uc.basis.get(1)[0];    // position fractionnaire du site B
        final double fy = uc.basis.get(1)[1];

        // ---- (1) grilles tau_g sur les cellules des sites A, soustraction du fond ----
        // Ordre vrai des sites : site = 2*(a + Na*b) + sr ; le site B voisin de A(a, b)
        // par la liaison (dx, dy) est en cellule (a + dx mod Na, b + dy mod Nb), sr = 1.
        final double[][] tau = new double[3][NaNb];
        double sommeTau = 0.0;
        for (int bond = 0; bond < 3; bond++) {
            final int dx = BONDS[bond][0];
            final int dy = BONDS[bond][1];
            final int g = BONDS[bond][2];

            for (int b = 0; b < Nb; b++) {
                for (int a = 0; a < Na; a++) {
                    final int siteA = 2 * (a + Na * b);
                    // cellule du site B, avec conditions periodiques
                    final int a2 = Math.floorMod(a + dx, Na);
                    final int b2 = Math.floorMod(b + dy, Nb);
                    final int siteB = 2 * (a2 + Na * b2) + 1;
                    final double t = (S_t_flat[siteA * 3 + g] * S_t_flat[siteB * 3 + g]);
                    final double t2 = t * t;           // tau = (S_A^g * S_B^g)^2
                    tau[bond][b * Na + a] = t2;
                    sommeTau += t2;
                }
            }
        }

        // Fond de dimeres retire AVANT la FFT : tauMean = moyenne de tau sur les
        // 3*Na*Nb liaisons de LA configuration -> correlateur connecte ; la
        // composante uniforme (le pied en Gamma, N/4 a T->0 sur l'appariement
        // parfait) ne masque plus le pinch point.
        final double tauMean = sommeTau / (3.0 * NaNb);

        final DoubleFFT_2D fft2D = new DoubleFFT_2D(Nb, Na);
        final double[] fftBuffer = new double[2 * NaNb];
        final double[][] reT = new double[3][NaNb];
        final double[][] imT = new double[3][NaNb];

        for (int bond = 0; bond < 3; bond++) {
            for (int lin = 0; lin < NaNb; lin++) {
                fftBuffer[2 * lin] = tau[bond][lin] - tauMean;
                fftBuffer[2 * lin + 1] = 0.0;
            }
            fft2D.complexForward(fftBuffer);
            // Bin (qx, qy) relu a (qy * Na + qx) : cf. convention de remplissage
            for (int qy = 0; qy < Nb; qy++) {
                for (int qx = 0; qx < Na; qx++) {
                    final int lin = qy * Na + qx;
                    reT[bond][lin] = fftBuffer[2 * lin];
                    imT[bond][lin] = fftBuffer[2 * lin + 1];
                }
            }
        }

        // ---- (2) Depliage multi-zones : boucles et index identiques au DSF ----
        for (int h = 0; h < H; h++) {
            double hklx = h;
            for (int k = 0; k < K; k++) {
                double hkly = k;
                for (int qx = 0; qx < Na; qx++) {
                    double qxx = getQValue(qx, Na) - 1.d;
                    for (int qy = 0; qy < Nb; qy++) {
                        double qyy = getQValue(qy, Nb) - 1.d;
                        final int Q = qIndex(h, k, qx, qy);
                        final int lin = qy * Na + qx;

                        // indices r.l.u. du Q courant (a_idx, b_idx), comme dans la phase DSF
                        final double aIdx = qxx + hklx;
                        final double bIdx = qyy + hkly;

                        // ---- (3) somme coherente des 3 couleurs de liaison ----
                        double sumRe = 0.0, sumIm = 0.0;
                        for (int bond = 0; bond < 3; bond++) {
                            final int dx = BONDS[bond][0];
                            final int dy = BONDS[bond][1];
                            // milieu de liaison : ((fx + dx)/2, (fy + dy)/2) en fractionnaire
                            final double phaseAngle = -Math.TAU
                                    * (aIdx * (fx + dx) / 2.0 + bIdx * (fy + dy) / 2.0);
                            final double phRe = Math.cos(phaseAngle);
                            final double phIm = Math.sin(phaseAngle);
                            final double tre = reT[bond][lin];
                            final double tim = imT[bond][lin];
                            sumRe += tre * phRe - tim * phIm;
                            sumIm += tre * phIm + tim * phRe;
                        }
                        double intensite = (sumRe * sumRe + sumIm * sumIm) / numSpins;

                        // Moyenne glissante sur les configurations
                        S_dimer[Q] = (Nmesures * S_dimer[Q] + intensite) / (Nmesures + 1);
                    }
                }
            }
        }

        Nmesures += 1;
    }

    /**
     * Ecrit le fichier Avro fusionne (une ligne par Q) sous le nom COMPLET fourni, compose
     * par l'appelant (ex. {@code structure_dimere_h0.000000.avro}).
     * Appele une seule fois, apres la boucle des nSQW+1 configurations.
     */
    public void writeAvro(String fileName) {
        try {
            Schema schema = new Schema.Parser().parse("{" +
                    "  \"type\": \"record\"," +
                    "  \"name\": \"DimerStructure\"," +
                    "  \"Qzz\": " + qzz + "," +
                    "  \"fields\": [" +
                    "    {\"name\": \"h\", \"type\": \"int\"}," +
                    "    {\"name\": \"k\", \"type\": \"int\"}," +
                    "    {\"name\": \"qx\", \"type\": \"int\"}," +
                    "    {\"name\": \"qy\", \"type\": \"int\"}," +
                    "    {\"name\": \"k_x\", \"type\": \"float\"}," +
                    "    {\"name\": \"k_y\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_dimer\", \"type\": \"float\"}" +
                    "  ]" +
                    "}");

            String dir = (outputDirPath == null) ? "." : outputDirPath;
            File filePath = new File(dir);
            if (!filePath.exists()) {
                filePath.mkdirs();
            }
            File file = new File(filePath, fileName);

            DatumWriter<GenericRecord> datumWriter = new GenericDatumWriter<>(schema);
            DataFileWriter<GenericRecord> dataFileWriter = new DataFileWriter<>(datumWriter);
            dataFileWriter.setCodec(CodecFactory.snappyCodec());
            dataFileWriter.create(schema, file);

            // Mapping k_x/k_y : identique au DSF (P_xy sur (qxx + h, qyy + k))
            double[][] P_xy = { { 1., 0. }, { 1. / Math.sqrt(3.), 2. / Math.sqrt(3.) } };
            for (int h = 0; h < H; h++) {
                double hklx = h;
                for (int k = 0; k < K; k++) {
                    double hkly = k;
                    for (int qx = 0; qx < Na; qx++) {
                        double qxx = getQValue(qx, Na) - 1.d;
                        for (int qy = 0; qy < Nb; qy++) {
                            double qyy = getQValue(qy, Nb) - 1.d;
                            final int Q = qIndex(h, k, qx, qy);
                            double[] Q_xy = MathOps.matrixXvector(P_xy,
                                    new double[] { qxx + hklx, qyy + hkly });

                            GenericRecord record = new GenericData.Record(schema);
                            record.put("h", (int) h);
                            record.put("k", (int) k);
                            record.put("qx", (int) qx);
                            record.put("qy", (int) qy);
                            record.put("k_x", (float) Q_xy[0]);
                            record.put("k_y", (float) Q_xy[1]);
                            record.put("S_dimer", (float) S_dimer[Q]);
                            dataFileWriter.append(record);
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
