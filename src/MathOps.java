import java.util.Arrays;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.FastMath;

/** Outils mathematiques statiques : vecteurs 3D, matrices 3x3, produits, normes. */
public final class MathOps {

    private MathOps() {
        // Classe utilitaire : aucune instance, appeler les methodes via MathOps.scale(...).
    }

    /** Produit de deux matrices carrees de meme taille. */
    public static double[][] mult2D(final double[][] A, final double[][] B) {
        int size = A.length;
        double[][] R = new double[size][size];
        for (int i = 0; i < size; i++) {
            for (int k = 0; k < size; k++) {
                double Aik = A[i][k];
                for (int j = 0; j < size; j++) {
                    R[i][j] += Aik * B[k][j];
                }
            }
        }
        return R;
    }

    /** Inverse d'une matrice carree (via commons-math). */
    public static double[][] inverse(double[][] A){
        RealMatrix matrix = MatrixUtils.createRealMatrix(A);
        RealMatrix inverse = MatrixUtils.inverse(matrix);
        return inverse.getData();
    }

    /** Transposee d'une matrice. */
    // https://stackoverflow.com/questions/15449711/transpose-double-matrix-with-a-java-function
    public static double[][] transpose(final double[][] A) {
        final double[][] T = new double[A[0].length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                T[j][i] = A[i][j];
            }
        }
        return T;
    }

    /** Difference de deux vecteurs. */
    public static double[] subtract(final double[] array1, final double[] array2) {
        double[] result = new double[array1.length];
        for (int i = 0; i < array1.length; i++) {
            result[i] = array1[i] - array2[i];
        }
        return result;
    }

    /** Norme euclidienne d'un vecteur. */
    public static double norm(final double[] array1) {
        double sumOfSquares = 0.;
        for (double component : array1)
            sumOfSquares = Math.fma(component, component, sumOfSquares);
        return FastMath.sqrt(sumOfSquares);
    }

    /** Norme euclidienne des 3 composantes situees a partir de startIndex. */
    public static double norm(final double[] array, final int startIndex) {
        double sumOfSquares = 0.;
        for (int i = startIndex; i < startIndex + 3; i++) {
            sumOfSquares = Math.fma(array[i], array[i], sumOfSquares);
        }
        return FastMath.sqrt(sumOfSquares);
    }

    /** Vecteur normalise (copie : l'entree n'est pas modifiee). */
    public static double[] normalize(final double[] array1) {
        double[] result = Arrays.copyOf(array1, array1.length);
        double norme = norm(result);
        for (int i = 0; i < result.length; i++) {
            result[i] /= norme;
        }
        return result;
    }

    /** Produit scalaire de deux vecteurs. */
    public static double dot(final double[] array1, final double[] array2) {
        double dotProd = 0.;
        for (int i = 0; i < array1.length; i++)
            dotProd = Math.fma(array1[i], array2[i], dotProd);
        return dotProd;
    }

    /** Produit nombre x vecteur. */
    public static double[] scale(final double number, final double[] array1) {
        double[] result = new double[array1.length];
        for (int i = 0; i < result.length; i++) {
            result[i] = number * array1[i];
        }
        return result;
    }

    /** Somme de deux vecteurs. */
    public static double[] add(final double[] array1, final double[] array2) {
        double[] result = new double[array1.length];
        for (int i = 0; i < array1.length; i++) {
            result[i] = array1[i] + array2[i];
        }
        return result;
    }

    /** Somme de deux vecteurs 3D, codee a plat. */
    public static double[] vectorizationadd(final double[] array1, final double[] array2) {
        return new double[] { array1[0] + array2[0], array1[1] + array2[1], array1[2] + array2[2] };
    }

    /** Produit matrice x vecteur. */
    public static double[] matrixXvector(final double[][] matrix, final double[] vector) {
        final int rows = matrix.length;
        double[] result = new double[rows];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < vector.length; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }
        return result;
    }

    /** Carre. */
    public static double sq(double a) {
        return a * a;
    }

    /**
     * Applique un fenêtrage gaussien in-place sur les données FFT temporelles.
     *
     * Fenêtre AUTOMATIQUE : sigma = 0.375 * Ntime, règlé pour que la largeur à
     * mi-hauteur du pic élastique (signal constant fenêtré) tombe sur ~1 bin de
     * la grille en énergie (FWHM = 2.355/sigma_T = 1.07 * 2*pi/Ntime), quelle
     * que soit la résolution --wmax/--domega choisie : le pic de Bragg ne bave
     * plus en énergie.
     *
     * @param fftTimeBuffer Tableau de données FFT temporelles intercalées [Re, Im,
     *                      Re, Im, ...].
     * @param Ntime         Nombre de pas de temps.
     */
    static void applyGaussianWindow(double[] fftTimeBuffer, int Ntime) {
        final double sigma = 0.75 * Ntime / 2.0d;
        final double mean = Ntime / 2.0;
        for (int t = 0; t < Ntime; t++) {
            double window = FastMath.exp(-0.5 * MathOps.sq((t - mean) / sigma));
            fftTimeBuffer[2 * t] *= window; // Partie réelle
            fftTimeBuffer[2 * t + 1] *= window; // Partie imaginaire
        }
    }

}
