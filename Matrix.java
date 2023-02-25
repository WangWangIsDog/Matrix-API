/* *****************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix file1.txt file2.txt
 *  Dependencies: In.java StdOut.java
 *
 *  % java Matrix matrix1.txt matrix2.txt
 *
 **************************************************************************** */

import edu.princeton.cs.algs4.In;
import edu.princeton.cs.algs4.StdOut;

public class Matrix {
    public static void main(String[] args) {

        // unit test for vector x vector
        // String[] data1 = new In(args[0]).readAllStrings();
        // String[] data2 = new In(args[1]).readAllStrings();

        // int N = Integer.parseInt(data1[0]);
        // int M = Integer.parseInt(data2[0]);
        //
        // if (N != M) throw new IllegalArgumentException("Two vectors must have equal sizes.");
        //
        // double[] firstVector = new double[N];
        // for (int i = 0; i < N; i++)
        //     firstVector[i] = Double.parseDouble(data1[i + 1]);
        //
        // double[] secondVector = new double[N];
        // for (int i = 0; i < N; i++)
        //     secondVector[i] = Double.parseDouble(data2[i + 1]);
        //
        // StdOut.println(dot(firstVector, secondVector));

        // unit test for matrix x matrix
        // String[] data1 = new In(args[0]).readAllStrings();
        // String[] data2 = new In(args[1]).readAllStrings();
        //
        // int N1 = Integer.parseInt(data1[0]);
        // int M1 = Integer.parseInt(data1[1]);
        // int N2 = Integer.parseInt(data2[0]);
        // int M2 = Integer.parseInt(data2[1]);
        //
        // if (M1 != N2) throw new IllegalArgumentException("Two matrices don`t match");
        //
        // double[][] firstMatrix = new double[N1][M1];
        // for (int i = 0; i < N1; i++)
        //     for (int j = 0; j < M1; j++)
        //         firstMatrix[i][j] = Double.parseDouble(data1[i * M1 + j + 2]);
        //
        // double[][] secondMatrix = new double[N2][M2];
        // for (int i = 0; i < N2; i++)
        //     for (int j = 0; j < M2; j++)
        //         secondMatrix[i][j] = Double.parseDouble(data2[i * M2 + j + 2]);
        //
        // double[][] result = mult(firstMatrix, secondMatrix);
        //
        // for (int i = 0; i < result.length; i++)
        //     for (int j = 0; j < result[0].length; j++) {
        //         StdOut.print(result[i][j] + " ");
        //         if (j == result[0].length - 1)
        //             StdOut.println();
        //     }

        // unit test for transpose
        // String[] data = new In(args[0]).readAllStrings();
        //
        // int N = Integer.parseInt(data[0]);
        // int M = Integer.parseInt(data[1]);
        // double[][] matrix = new double[N][M];
        // for (int i = 0; i < N; i++)
        //     for (int j = 0; j < M; j++)
        //         matrix[i][j] = Double.parseDouble(data[i * M + j + 2]);
        //
        // double[][] result = transpose(matrix);
        // for (int i = 0; i < M; i++)
        //     for (int j = 0; j < N; j++) {
        //         StdOut.print(result[i][j] + " ");
        //         if (j == N - 1)
        //             StdOut.println();
        //     }

        // unit test for matrix x vector
        // String[] data1 = new In(args[0]).readAllStrings();
        // String[] data2 = new In(args[1]).readAllStrings();
        //
        // int N1 = Integer.parseInt(data1[0]);
        // int M1 = Integer.parseInt(data1[1]);
        // int N2 = Integer.parseInt(data2[0]);
        //
        // if (M1 != N2) throw new IllegalArgumentException("com.lixing.algorithm.matrix.Matrix and vector are not matched");
        //
        // double[][] matrix = new double[N1][M1];
        // for (int i = 0; i < N1; i++)
        //     for (int j = 0; j < M1; j++)
        //         matrix[i][j] = Double.parseDouble(data1[i * M1 + j + 2]);
        //
        // double[] vector = new double[N2];
        // for (int i = 0; i < N2; i++)
        //     vector[i] = Double.parseDouble(data2[i + 1]);
        //
        // double[] result = mult(matrix, vector);
        // for (int i = 0; i < result.length; i++)
        //     StdOut.println(result[i]);

        // unit test for vector x matrix
        String[] data1 = new In(args[0]).readAllStrings();
        String[] data2 = new In(args[1]).readAllStrings();

        int M1 = Integer.parseInt(data1[0]);
        int N2 = Integer.parseInt(data2[0]);
        int M2 = Integer.parseInt(data2[1]);

        if (M1 != N2) throw new IllegalArgumentException("Vector and matrix don`t match");

        double[] vector = new double[M1];
        for (int i = 0; i < M1; i++)
            vector[i] = Double.parseDouble(data1[i + 1]);

        double[][] matrix = new double[N2][M2];
        for (int i = 0; i < N2; i++)
            for (int j = 0; j < M2; j++)
                matrix[i][j] = Double.parseDouble(data2[i * M2 + j + 2]);

        double[] result = mult(vector, matrix);
        for (int i = 0; i < M2; i++) {
            StdOut.print(result[i] + " ");
            if (i == M2 - 1)
                StdOut.println();
        }
    }

    // dot product of two vectors
    public static double dot(double[] a, double[] b) {

        double result = 0.0;
        for (int i = 0; i < a.length; i++) {
            result += a[i] * b[i];
        }
        return result;
    }

    // multiplication between two matrices
    public static double[][] mult(double[][] a, double[][] b) {

        int N = a.length;
        int M = b[0].length;
        int K = a[0].length;
        double[][] result = new double[N][M];

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                for (int k = 0; k < K; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return result;
    }

    // transpose of matrix
    public static double[][] transpose(double[][] a) {

        int N = a.length;
        int M = a[0].length;

        double[][] result = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                result[i][j] = a[j][i];

        return result;
    }

    // multiplication between matrix and vector
    public static double[] mult(double[][] a, double[] b) {

        int N = a.length;
        int M = a[0].length;

        double[] result = new double[N];
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
                result[i] += a[i][j] * b[j];

        return result;
    }

    // multiplication between vector and matrix
    public static double[] mult(double[] a, double[][] b) {

        int M1 = a.length;
        int M2 = b[0].length;

        double[] result = new double[M1];
        for (int i = 0; i < M2; i++)
            for (int j = 0; j < M1; j++)
                result[i] += a[j] * b[j][i];

        return result;
    }

}
