/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package predictionmodel;

import Jama.Matrix;

/**
 *
 * @author Ana
 */
public class Particle {

    double[] mu;
    double[][] A;
    double[][] sigma;
    double L;
    double[][] X;

    public Particle(double[][] sigma, double L, double[][] X) {
        this.sigma = sigma;
        this.L = L;
        this.X = X;
    }

    public Particle(double[] miu, double[][] sigma, double L, double[][] X, double[][] A) {
        this.mu = miu;
        this.sigma = sigma;
        this.L = L;
        this.X = X;
        this.A = A;
    }

    public Particle() {
    }

    public double[][] getA() {
        return A;
    }

    public void setA(double[][] A) {
        this.A = A;
    }

    public double[][] getX() {
        return X;
    }

    public void setX(double[][] X) {
        this.X = X;
    }

    public double[] getMu() {
        return mu;
    }

    public void setMu(double[] mu) {
        this.mu = mu;
    }

    public double[][] getSigma() {
        return sigma;
    }

    public void setSigma(double[][] sigma) {
        this.sigma = sigma;
    }

    public double getL() {
        return L;
    }

    public void setL(double L) {
        this.L = L;
    }

    /*sum of Xi * Xi_transpose*/
    public double[][] multiplySecond(double[][] first, double[][] second) {
        int d = first.length;
        double[][] matrix = new double[d][d];

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {

                matrix[i][j] = first[i][0] * second[0][j];
            }
        }
        return matrix;
    }

    /*sum of 2 matrices*/
    public double[][] add(double[][] A, double[][] B) {
        int m = A.length;
        int n = A[0].length;
        double[][] C = new double[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[i][j] = A[i][j] + B[i][j];
            }
        }
        return C;
    }

    /*the transpose of a matrix*/
    public double[][] transpose(double[][] A) {
        int m = A.length;
        int n = A[0].length;
        double[][] C = new double[n][m];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                C[j][i] = A[i][j];
            }
        }
        return C;
    }

    /*the calcul of the low triangular matrix, using Cholesky decomposition*/
    public double[][] calculTriangularMatrixL(double[][] sigmaScale) {
        double[][] matrixL;
        Matrix M = new Matrix(sigmaScale);
        matrixL = M.chol().getL().getArrayCopy();
        return matrixL;
    }

    /*generation of mu vector using Box Muller method*/
    public double[] generateMu(int size) {
        BoxMuller BM = new BoxMuller();
        double[] vector = BM.boxMuller(size);
        return vector;
    }

    public void readVector(double[] Z) {
        for (int i = 0; i < Z.length; i++) {
            System.out.print(" " + Z[i] +"\n");
        }
    }

    /*perturbation of the encoded vector mu*/
    public double[] perturbMu(double[] mu, int aux) {
        double[] perturbedMu = new double[mu.length];
        double val = 0.0;
        for (int i = 0; i < mu.length; i++) {
            BoxMuller BM = new BoxMuller();
            val = BM.boxMuller(1)[0];
            perturbedMu[i] = mu[i] + aux * val;
        }
        return perturbedMu;
    }

    /*d = dimension of the space, n = degrees of freedom*/
    public double[][] generateParticleVector(int d, double[][] scaleMatrix, double[] mu, int n) {
        double[][] particleVector = new double[d][n];
        double[][] lowTriangularMatrix = calculTriangularMatrixL(scaleMatrix);

        for (int j = 0; j < n; j++) {
            BoxMuller BM = new BoxMuller();
            double[] vector = BM.boxMuller(d);
            double[] product = multiplyMV(lowTriangularMatrix, vector);

            for (int i = 0; i < d; i++) {
                particleVector[i][j] = mu[i] + product[i];
            }
        }

        return particleVector;
    }

    /*calculate sigma parameter for particle*/
    public double[][] calculateSigma(double[][] XM, int n) {
        int d = XM.length;
        double[][] sigma = new double[d][d];
        for (int i = 0; i < n; i++) {
            //        System.out.println("pasul  " + i);
            double[][] x = new double[d][1];
            double[][] xt = new double[1][d];
            for (int j = 0; j < d; j++) {
                x[j][0] = XM[j][i];
            }
            xt = transpose(x);
            double[][] product = multiplySecond(x, xt);
            sigma = add(sigma, product);
        }
        return sigma;
    }

    /*perturb the encoded vector X used to calculate sigma, with value a that influences the strength of the perturbation*/
    public double[][] perturbX(double[][] X, double[][] Y, double a, int d, int n) {
        double[][] perturbedX = new double[d][n];
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < n; j++) {
                perturbedX[i][j] = Math.sqrt(a) * X[i][j] + Math.sqrt(1 - a) * Y[i][j];
            }
        }
        return perturbedX;
    }

    /*product scalar with matrix*/
    public double[][] multipySM(double a, double[][] A) {
        double[][] B = new double[A.length][A[0].length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                B[i][j] = a * A[i][j];
            }
        }
        return B;
    }

    public double[] multiplyMV(double[][] A, double[] x) {
        int p = x.length;

        //System.out.println("A after calcul low triangular matrix       00"+A[0][0] +"  01" +A[0][1]);
        double[] y = new double[p];
        for (int i = 0; i < p; i++) {
            y[i] = A[0][0] * x[i];
        }

        return y;

    }

    /*multiplication of two matrices*/
    public double[][] multiplyMM(double[][] A, double[][] B) {
        int mA = A.length;
        int nA = A[0].length;
        int mB = B.length;
        int nB = B[0].length;
        if (nA != mB) {
            throw new RuntimeException("Illegal matrix dimensions.");
        }
        double[][] C = new double[mA][nB];
        for (int i = 0; i < mA; i++) {
            for (int j = 0; j < nB; j++) {
                for (int k = 0; k < nA; k++) {
                    C[i][j] += (A[i][k] * B[k][j]);
                }
            }
        }
        return C;
    }

    public double[][] scaleMatrix(int d) {
       /* double[][] scaleSigma = new double[d][d];
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                if (i == j) {
                    scaleSigma[i][i] = 20;
                } else {
                    scaleSigma[i][j] = 0;
                }
            }
        }*/

        double[][] scaleSigma = {{10000, 0}, {0, 10000}};

      /*  System.out.println("Scale Matrix");
        readMatrix(scaleSigma);
        */return scaleSigma;
    }

    public void readMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j] + "   ");
            }
            System.out.println();
        }
    }

    /*the calculation of the likelihood of a Particle*/
    public double calculL(double[][] samples, double[][] sigma, double[] mu) {
        Matrix sigmaM = new Matrix(sigma);
        int p = sigma.length;
        int N = samples[0].length;
        double det = sigmaM.det();
        double[][] inverse = sigmaM.inverse().getArrayCopy();

        double val1 = -(N * p) / 2.0;
        double val2 = -N / 2.0;
        double val3 = 0.0;
        double piDouble = 2 * Math.PI;

        double[] pow = new double[N];

        for (int i = 0; i < N; i++) {
            double[] firstSubstraction = new double[p];

            for (int j = 0; j < p; j++) {
                firstSubstraction[j] = samples[j][i] - mu[j];
            }

            double[] firstMultiplication = multiplyVM(firstSubstraction, inverse);
            double secondMultiplication = dot(firstMultiplication, firstSubstraction);
            pow[i] = secondMultiplication / 2.0;
        }

        for (int i = 0; i < N; i++) {
            val3 += pow[i];
        }

        double Lvalue = val1 * Math.log(piDouble) + val2 * Math.log(det) - val3;
        return Lvalue;
    }

    public double dot(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new RuntimeException("Illegal vector dimensions.");
        }
        double sum = 0.0;
        for (int i = 0; i < x.length; i++) {
            sum += x[i] * y[i];
        }
        return sum;
    }

    public double[] multiplyVM(double[] x, double[][] A) {
        int m = A.length;
        int n = A[0].length;
        if (x.length != m) {
            throw new RuntimeException("Illegal matrix dimensions.");
        }
        double[] y = new double[n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < m; i++) {
                y[j] += (A[i][j] * x[i]);
            }
        }
        return y;
    }

}
