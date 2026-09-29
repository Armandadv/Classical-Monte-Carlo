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
 * Facteur de structure de COULEUR (canal quadrupolaire), invariant sous les
 * retournements Z2 du modele KBQ :
 *
 *   S_col(Q) = (1/numSpins) * Sum_a < | Sum_sites n_i^a * e^{-i Q . r_i} |^2 >
 *
 * avec n_i^a = (S_i^a)^2 : les CARRES des composantes de spin, et non les spins
 * eux-memes. C'est le canal quadrupolaire : aucun signe de spin n'y survit,
 * contrairement au DSF magnetique de DynamicStructureFactor.
 *
 * Depliage multi-zones identique au DSF : memes boucles (h, k, qx, qy), meme indice
 * lineaire Q = (h*Na*Nb*K) + (k*Na*Nb) + (qx*Nb + qy), et MEME phase inter-sous-reseaux
 * copiee verbatim de DSF.processSQW (convention -TAU incluse ; reseau honeycomb 2D :
 * qzz = 0 et deltaZ = 0, les sous-reseaux sont coplanaires ; le delta 3D généralisé vit
 * dans DynamicStructureFactor). Statique : une configuration
 * par mesure ; l'intensite
 * est moyenee sur les nSQW+1 configurations (moyenne glissante dans processConfig) puis
 * ecrite en Avro par writeAvro (appele une seule fois, a la derniere temperature,
 * cf. MonteCarlo.sqwBlock).
 *
 * CONVENTION FFT (a ne PAS confondre avec le DSF) : l'ordre vrai des sites est
 * site = 2*(a + Na*b) + sr (cf. Lattice.siteIndex : a = index rapide, direction
 * A1 ; b = index lent, direction A2). Le buffer complexe packed est donc organise en
 * Nb lignes (direction b) x Na colonnes (direction a), rempli en row-major (b*Na + a)
 * et transforme par DoubleFFT_2D(Nb, Na) ; le bin (qx, qy) est relu a l'index
 * (qy*Na + qx). Convention directe, correcte pour Na et Nb quelconques — desormais
 * la meme convention que le DSF (aligne sur siteIndex, correct pour tout Na/Nb).
 *
 * Contrairement au DSF, PAS de normalisation sqrt(Na*Nb) (leur factorX) : les
 * echelles Color et SQW ne sont pas comparables.
 */
public class ColorStructureFactor {

    private final int H, K, Na, Nb, numSpins, totalQ;
    private final double[] S_color;      // intensite moyennee par Q
    private int Nmesures = 0;
    private String outputDirPath = null; // null = repertoire courant
    private double qzz = 0.d;            // composante z du vecteur de diffusion (r.l.u.) ; honeycomb 2D (ancienne valeur BaCoAsO : 4.67)

    /**
     * Alloue les accumulateurs pour une grille de H x K zones : dimensions Na, Nb et
     * nombre de sites lus dans le lattice ; S_color(Q) part de zero.
     */
    public ColorStructureFactor(int H, int K, Lattice lattice) {
        this.H = H;
        this.K = K;
        this.Na = lattice.lengthA;
        this.Nb = lattice.lengthB;
        this.numSpins = lattice.size;
        this.totalQ = Na * H * Nb * K;
        this.S_color = new double[totalQ];
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
     * Accumule une configuration dans la moyenne glissante de S_col(Q).
     *
     * @param uc       maille unite ; sert uniquement a la phase inter-sous-reseaux
     *                 (positions relatives des 2 sous-reseaux, uc.basis.get(1)).
     * @param S_t_flat spins aplaties [instantane][site][3] ; seuls les N premiers triplets
     *                 (premier instantane) sont lus, ELEVES AU CARRE :
     *                 n[site][a] = (S_t_flat[site*3 + a])^2. Ce sont des CARRES, source
     *                 unique de verite du canal quadrupolaire (invariant Z2).
     */
    public void processConfig(UnitCell uc, double[] S_t_flat) {

        final int NaNb = Na * Nb;

        // ---- (1) FFT 2D des n^a = (S^a)^2, par composante et par sous-reseau sr ----
        // Ordre vrai des sites (cf. Lattice.siteIndex) : site = 2*(a + Na*b) + sr,
        // a = index rapide (direction A1, Na valeurs), b = index lent (direction A2, Nb
        // valeurs). Buffer organise en Nb lignes (b) x Na colonnes (a), row-major
        // (b*Na + a), transforme par DoubleFFT_2D(Nb, Na) ; cf. javadoc de classe.
        final DoubleFFT_2D fft2D = new DoubleFFT_2D(Nb, Na);
        final double[] fftBuffer = new double[2 * NaNb];
        final double[][] reS0 = new double[3][NaNb];
        final double[][] imS0 = new double[3][NaNb];
        final double[][] reS1 = new double[3][NaNb];
        final double[][] imS1 = new double[3][NaNb];

        for (int comp = 0; comp < 3; comp++) {
            for (int sr = 0; sr < 2; sr++) {
                // lignes = direction b (Nb), colonnes = direction a (Na) :
                for (int b = 0; b < Nb; b++) {
                    for (int a = 0; a < Na; a++) {
                        final int site = 2 * (a + Na * b) + sr;
                        final double s = S_t_flat[site * 3 + comp];
                        fftBuffer[2 * (b * Na + a)] = s * s;   // n^comp = (S^comp)^2
                        fftBuffer[2 * (b * Na + a) + 1] = 0.0; // partie imaginaire
                    }
                }
                fft2D.complexForward(fftBuffer);
                // Bin (qx, qy) : ligne qy (frequence selon b), colonne qx (frequence selon a),
                // relu a l'index (qy * Na + qx).
                for (int qy = 0; qy < Nb; qy++) {
                    for (int qx = 0; qx < Na; qx++) {
                        final int lin = qy * Na + qx;
                        if (sr == 0) {
                            reS0[comp][lin] = fftBuffer[2 * lin];
                            imS0[comp][lin] = fftBuffer[2 * lin + 1];
                        } else {
                            reS1[comp][lin] = fftBuffer[2 * lin];
                            imS1[comp][lin] = fftBuffer[2 * lin + 1];
                        }
                    }
                }
            }
        }

        // ---- (2) Depliage multi-zones : boucles et phase copiees du DSF ----
        // Phase inter-sous-reseaux copiee verbatim de DSF.processSQW (convention -TAU incluse) :
        final double deltaX = uc.basis.get(1)[0];
        final double deltaY = uc.basis.get(1)[1];
        final double deltaZ = 0.d; // canal statique insensible au z : offset inter-sous-reseaux nul codé ici ; le delta 3D généralisé vit dans DynamicStructureFactor

        for (int h = 0; h < H; h++) {
            double hklx = h;
            for (int k = 0; k < K; k++) {
                double hkly = k;
                for (int qx = 0; qx < Na; qx++) {
                    double qxx = getQValue(qx, Na) - 1.d;
                    for (int qy = 0; qy < Nb; qy++) {
                        double qyy = getQValue(qy, Nb) - 1.d;
                        // Indice Q, identique au DSF
                        int Q = (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);
                        // Bin FFT stocke a l'index (qy * Na + qx) : cf. convention de remplissage
                        final int lin = qy * Na + qx;

                        double Qr = (qxx + hklx) * deltaX + (qyy + hkly) * deltaY + qzz * deltaZ;
                        double phaseAngle = -Math.TAU * Qr;
                        double phaseRe = Math.cos(phaseAngle);
                        double phaseIm = Math.sin(phaseAngle);

                        // ---- (3) Combinaison des sous-reseaux puis intensite ----
                        double intensite = 0.0;
                        for (int a = 0; a < 3; a++) {
                            final double s0re = reS0[a][lin];
                            final double s0im = imS0[a][lin];
                            final double s1re = reS1[a][lin];
                            final double s1im = imS1[a][lin];
                            // S1 multiplie par le facteur de phase e^{-i TAU Qr}
                            final double s1phRe = s1re * phaseRe - s1im * phaseIm;
                            final double s1phIm = s1re * phaseIm + s1im * phaseRe;
                            final double reComb = s0re + s1phRe;
                            final double imComb = s0im + s1phIm;
                            intensite += reComb * reComb + imComb * imComb;
                        }
                        intensite /= numSpins;

                        // Moyenne glissante sur les configurations
                        S_color[Q] = (Nmesures * S_color[Q] + intensite) / (Nmesures + 1);
                    }
                }
            }
        }

        Nmesures += 1;
    }

    /**
     * Ecrit le facteur de structure de couleur moyenne en Avro (codec snappy, comme le DSF).
     * Un enregistrement par (h, k, qx, qy) ; k_x/k_y = composantes (h,k)+(qx,qy) tournees
     * dans la base cartesienne x,y via P_xy (copiee du DSF).
     *
     * @param fileName NOM COMPLET du fichier (sans repertoire), compose par l'appelant
     *                 (ex. {@code structure_couleur_h0.000000.avro}) ; ecrit dans le repertoire
     *                 fixe par {@link #setOutputDirPath} (defaut : repertoire courant).
     */
    public void writeAvro(String fileName) {
        final int NaNb = Na * Nb;
        // Matrice changement de base a,b vers x,y (copiee du DSF)
        double[][] P_xy = { { 1., 0. }, { 1. / Math.sqrt(3.), 2. / Math.sqrt(3.) } };
        try {
            Schema schema = new Schema.Parser().parse("{" +
                    "  \"type\": \"record\"," +
                    "  \"name\": \"ColorStructure\"," +
                    "  \"Qzz\": " + qzz + "," +
                    "  \"fields\": [" +
                    "    {\"name\": \"h\", \"type\": \"int\"}," +
                    "    {\"name\": \"k\", \"type\": \"int\"}," +
                    "    {\"name\": \"qx\", \"type\": \"int\"}," +
                    "    {\"name\": \"qy\", \"type\": \"int\"}," +
                    "    {\"name\": \"k_x\", \"type\": \"float\"}," +
                    "    {\"name\": \"k_y\", \"type\": \"float\"}," +
                    "    {\"name\": \"S_color\", \"type\": \"float\"}" +
                    "  ]" +
                    "}");

            // null (jamais configure) = repertoire courant
            File filePath = new File((outputDirPath == null) ? "." : outputDirPath);
            File file = new File(filePath, fileName);

            if (!filePath.exists()) {
                filePath.mkdirs();  // creer le repertoire s'il n'existe pas
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
                        double qxx = getQValue(qx, Na) - 1.d;
                        for (int qy = 0; qy < Nb; qy++) {
                            int Q = (h * NaNb * K) + (k * NaNb) + (qx * Nb + qy);
                            double qyy = getQValue(qy, Nb) - 1.d;
                            double[] Q_xy = MathOps.matrixXvector(P_xy, new double[] { qxx + hklx, qyy + hkly });

                            GenericRecord record = new GenericData.Record(schema);
                            record.put("h", h);
                            record.put("k", k);
                            record.put("qx", qx);
                            record.put("qy", qy);
                            record.put("k_x", (float) Q_xy[0]);
                            record.put("k_y", (float) Q_xy[1]);
                            record.put("S_color", (float) S_color[Q]);
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
