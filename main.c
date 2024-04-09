#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 51
#define M 51
#define K 60001


typedef double* Vector;
typedef double** Matrix;
typedef double*** Matrix3D;

float B(float x, float t)
{
    return 0.5; /* 15 * (x-0.5) / (1 - pow(abs(15 * (x-0.5)), 2))*/
}

float C(float x, float t)
{
    return 0.5;
}

Matrix3D init_zero_matrix3d(int x_size, int y_size, int z_size) {
    Matrix3D matrix = malloc(sizeof(Matrix) * z_size);
    if (matrix == NULL) {
        return NULL;
    }

    for (int k = 0; k < z_size; k++)
    {
        matrix[k] = malloc(sizeof(Vector) * y_size);
        if (matrix[k] == NULL) {
            return NULL;
        }
        for (int j = 0; j < y_size; j++){
            matrix[k][j] = malloc(sizeof(float) * x_size);
            if (matrix[k][j] == NULL) {
                return NULL;
            }
        }
    }

    for (int k = 0; k < z_size; k++)
        for (int j = 0; j < y_size; j++)
            for (int i = 0; i < x_size; i++)
                matrix[k][j][i] = 0.0;

    return matrix;
}

Vector init_zero_vector(int size) {
    Vector vector = malloc(sizeof(float) * size);
    if (vector == NULL) {
        return NULL;
    }
    for (int i = 0; i < size; i++)
        vector[i] = 0;

    return vector;
}

Matrix init_zero_matrix(int x_size, int y_size) {
    Matrix matrix = malloc(sizeof(float*) * y_size);
    if (matrix == NULL) {
        return NULL;
    }

    for (int j = 0; j < y_size; j++){
        matrix[j] = malloc(sizeof(float) * x_size);
        if (matrix[j] == NULL) {
            return NULL;
        }
    }

    for (int j = 0; j < y_size; j++)
        for (int i = 0; i < x_size; i++)
            matrix[j][i] = 0.0;

    return matrix;
}

float Kef(float phi)
{
    return phi* phi * phi;
}

int main()
{
    float pa = 1.0;
    float ps0 = 1.0;
    float rhof = 2.0;
    float rhos = 2600.0;
    double tau0 = 1.0 / 31104000.0;
    double k0 = 1.0 / 10000000000000.0;
    float alpha = 101368;
    int L = 2000;
    int H = 50;
    float g = 9.81;
    double rhof_g = rhof * g * H / alpha;

    float tau =  1.0 / (K - 1);
    float hx = 1.0 / (N - 1);
    float hz = 1.0 / (M - 1);

    double alph = tau0 * L * L / (k0 * alpha);

    Matrix3D phi0 = init_zero_matrix3d(N, M, K);
    Matrix3D phi0_12 = init_zero_matrix3d(N, M, K);

    for (int j = 0; j < M; j++)
    {
        for (int i = 0; i < N; i++)
        {
            phi0[0][j][i] = 0.82;
            phi0_12[0][j][i] = 0.82;
        }
    }

    Vector alphx = init_zero_vector(N);
    Vector betax = init_zero_vector(N);
    Vector alphz = init_zero_vector(M);
    Vector betaz = init_zero_vector(M);

    for (int k = 0; k < K-1; k++)
    {
        for (int j = 1; j < M-1; j++)
        {
            alphx[0] = 0;
            betax[0] = 0.82;
            for (int i = 0; i < N-1; i++)
            {
                float BB = B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k);
                float A = B(i * hx, tau * k) * tau / (4 * hx);
                float G = A;
                float P = 1;
                float F = phi0[k][j][i] - C(i * hx, tau * k) * tau / (4 * hz) * (phi0[k][j+1][i] - phi0[k][j-1][i]) - phi0[k][j][i] * tau / (4*hx) * BB;
                alphx[i+1] = G / (P - alphx[i]*A);
                betax[i+1] = (F + A*betax[i]) / (P - alphx[i]*A);
            }

            phi0_12[k+1][j][N-1] = 0.82; /* betax[N-1]/(1 - alphx[N-1]); */
            for (int i = N-2; i >= 0; i--)
                phi0_12[k + 1][j][i] = alphx[i+1] * phi0_12[k+1][j][i+1] + betax[i+1];
        }

        for (int i = 1; i < N-1; i++)
        {
            alphz[0] = 0;
            betaz[0] = 0.82;
            for (int j = 0; j < M-1; j++)
            {
                float BB = B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k);
                float A = C(i * hx, tau * k) * tau / (4 * hz);
                float G = A;
                float P = 1;
                float F = phi0_12[k][j][i] - B(i * hx, tau * k) * tau / (4 * hx) * (phi0_12[k][j][i+1] - phi0_12[k][j][i-1]) - phi0[k][j][i] * tau / (4*hx) * BB;
                alphz[j+1] = G / (P - alphz[j] * A);
                betaz[j+1] = (F + A * betaz[j]) / (P - alphz[j] * A);
            }

            phi0[k+1][M-1][i] = betaz[M-1]/(1 - alphz[M-1]);
            for (int j = M-2; j >= 0; j--)
                phi0[k + 1][j][i] = alphz[j+1] * phi0[k+1][j+1][i] + betaz[j+1];
        }
    }

    float v[N];
    for (int i = 0; i < N; i++)
    {
        v[i] = 0.5 * exp(-pow((20.0 * i * hx - 10.0), 10.0));
    }

    Matrix3D pf = init_zero_matrix3d(N, M, K);
    Matrix3D ps = init_zero_matrix3d(N, M, K);

    Matrix integrate_H = init_zero_matrix(N, K);
    Matrix integrate_H2 = init_zero_matrix(N, K);
    Matrix3D integrate_Z = init_zero_matrix3d(N, M, K);
    Matrix3D integrate_Z2 = init_zero_matrix3d(N, M, K);

    /* float ps_integration_H[K][N];
       float ps_integration_Z = init_zero_matrix3d(N, M, K); */

    for (int k = 0; k < K; k++)
    {
        for (int i = 0; i < N; i++)
        {
            float sum_h = 0;
            float sum_h2 = 0;
            for (int j = 0; j < M; j++)
            {
                sum_h += hz / Kef(phi0[k][j][i]);
                sum_h2 += j * hz / Kef(phi0[k][j][i]);
            }
            integrate_H[k][i] = sum_h;
            integrate_H2[k][i] = sum_h2;
            for (int j = 0; j < M; j++)
            {
                float sum_z = 0;
                float sum_z2 =0;
                for (int jj = 0; jj < j; jj++)
                {
                    sum_z = hz / Kef(phi0[k][jj][i]);
                    sum_z2 = jj * hz / Kef(phi0[k][jj][i]);
                }
                integrate_Z[k][j][i] = sum_z;
                integrate_Z2[k][j][i] = sum_z2;
            }
        }
    }

    for (int k = 0; k < K-1; k++)
    {
        float vt = exp(-pow((10 * tau * k - 0.5), 10));
        for (int j = 0; j < M; j++)
        {
            for (int i = 0; i < N; i++)
            {
                float BB = (B((i + 1) * hx, tau * k) - B((i - 1) * hx, tau * k)) / (2 * hx);
                pf[k][j][i] = pa + rhof_g  * (1 - j * hz) + alph * BB * (integrate_Z2[k][j][i] - integrate_H2[k][i]) - \
                            alph * (phi0[k][j][i] * v[i] * vt + BB) * (integrate_Z[k][j][i] - integrate_H[k][i]);
            }
        }
    }

    FILE* pf_file = fopen("./pf_file.txt", "w");
    if (pf_file == NULL)
        return -1;

    for (int j = 0; j < M; j++)
        for (int i = 0; i < N; i++)
            fprintf(pf_file, "%f ", pf[100][j][i]);

    fprintf(pf_file, "%c", '\n');

    fclose(pf_file);
    printf("\tSuccessful!\n");

    return 0;
}
