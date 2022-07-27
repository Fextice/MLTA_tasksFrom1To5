package ru.vsu.cs.course1;

class Task {
    public static int N = 5; // кол-во столбцов

    public static void getConf(double[][] mt, double[][] temp, int p, int q, int n) { // p, q - вычеркнутые строка и столбец (номера)
        int i = 0, j = 0;
        for (int row = 0; row < n; row++) { //перебирает строки
            for (int col = 0; col < n; col++) { //перебирает столбцы
                if (row != p && col != q) {
                    temp[i][j++] = mt[row][col]; //заполняет матрицу без вычеркнутых элементов
                    if (j == n - 1) {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    public static double recursionOfAMatrix(double[][] mt, int n) { // нахождения детерминанта рекурсией
        int D = 0; //определитель
        if (n == 1) {
            return mt[0][0]; //элементы строки
        }
        double[][] temp = new double[N][N];
        int sign = 1; // -1 в степени (i+j)
        for (int f = 0; f < n; f++) {
            getConf(mt, temp, 0, f, n); //mt[0][f] - элементы первой строки (они идут по очереди)
            D += sign * mt[0][f] * recursionOfAMatrix(temp, n - 1); //алгебраическое дополнение (матрица после вычёркивания)
            sign = -sign;
        }
        return D;
    }

    static void adj(double[][] mt, double[][] adj) // нахождение сопряженной матрицы
    {
        if (N == 1) {
            adj[0][0] = 1;
            return;
        }
        int sign = 1;
        double[][] temp = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                getConf(mt, temp, i, j, N);
                sign = ((i + j) % 2 == 0) ? 1 : -1;
                adj[j][i] = (sign) * (recursionOfAMatrix(temp, N - 1));
            }
        }
    }

    static double[][] inv(double[][] mt) { // нахождение развернутой матрицы
        double det = recursionOfAMatrix(mt, N);
        double[][] adj = new double[N][N];
        adj(mt, adj);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                mt[i][j] = adj[i][j] / det; // нахождение развернутой по формуле "развернутая матрица(mt) = adj(mt)/det(mt)"
            }
        }
        return mt;
    }

	 public static double[] gaussian(double[][] A, double[] B) {
        for (int p = 0; p < N; p++) {
            int max = p;
            for (int i = p + 1; i < N; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }
            double[] temp = A[p]; A[p] = A[max]; A[max] = temp;
            double t = B[p]; B[p] = B[max]; B[max] = t;
            for (int i = p + 1; i < N; i++) {
                double alpha = A[i][p] / A[p][p];
                B[i] -= alpha * B[p];
                for (int j = p; j < N; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }
        double[] x = new double[N];
        for (int i = N - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < N; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (B[i] - sum) / A[i][i];
        }
        return x;
    }

    public static double[] kramer(double[][] A, double[] B)
    {
        double[][] temp = new double[N][N];
        double[] x = new double[N];
        for(int i=0;i<N;i++)
        {
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    if(k == i)
                        temp[j][k] = B[j];
                    else
                        temp[j][k] = A[j][k];
                }
            }
            x[i]=recursionOfAMatrix(temp,N)/recursionOfAMatrix(A,N);
        }
        return x;
    }

    public static double[][] mul(double[][] FirstMt, double[][] SecondMt) {

        int FirstMtRows = FirstMt.length;
        int FirstMtColumns = FirstMt[0].length;
        int SecondMtRows = SecondMt.length;
        int SecondMtColumns = SecondMt[0].length;

        if (FirstMtColumns != SecondMtRows) {
            throw new IllegalArgumentException("Cтрочки: " + FirstMtColumns + " не равняются колонкам " + SecondMtRows + ".");
        }

        double[][] MulMatrix = new double[FirstMtRows][SecondMtColumns];
        for (int i = 0; i < FirstMtRows; i++) {
            for (int j = 0; j < SecondMtColumns; j++) {
                MulMatrix[i][j] = 0.000;
            }
        }

        for (int i = 0; i < FirstMtRows; i++) { // aRow
            for (int j = 0; j < SecondMtColumns; j++) { // bColumn
                for (int k = 0; k < FirstMtColumns; k++) { // aColumn
                    MulMatrix[i][j] += FirstMt[i][k] * SecondMt[k][j];
                }
            }
        }
        return MulMatrix;
    }
}
