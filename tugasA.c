// Tugas Pemrograman A

// Kelompok 1:
// Andrew Kristofer Jian		NPM. 2206059673
// Hanif Nur Ilham Sanjaya		NPM. 2206059692

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 1000

// Struktur data
typedef struct {
    double year[MAX_POINTS];
    double internet[MAX_POINTS];
    double population[MAX_POINTS];
    int n;
} FullDataset;

// Function prototype
int read_full_csv(const char *filename, FullDataset *ds);
void newton_interpolate(double *x, double *y, int n, double xq, double *result);
void polynomial_regression(double *x, double *y, int n, int degree, double *coeffs);
double eval_polynomial(const double *coeffs, int degree, double x);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s Data_Tugas_Pemrograman_A.csv\n", argv[0]);
        return EXIT_FAILURE;
    }

    FullDataset ds;
    if (read_full_csv(argv[1], &ds) != 0) {
        fprintf(stderr, "Error: gagal membaca file '%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    // Array bantu
    double *years = ds.year;
    double *inet  = ds.internet;
    double *pop   = ds.population;
    int n = ds.n;

    // 1. Estimasi data hilang via Newton (menggunakan 4 titik terdekat)
    double missing_years[] = {2005, 2006, 2015, 2016};
    int nm = sizeof(missing_years)/sizeof(double);
    printf("Estimasi data hilang menggunakan Newton Interpolation:\n");
    for (int i = 0; i < nm; i++) {
        double yq = missing_years[i];
        // Pilih 4 titik terdekat (2 sebelum dan 2 sesudah jika memungkinkan)
        double xs[4], ys_pop[4], ys_inet[4];
        int indices[4];
        // Inisialisasi indices dengan nilai yang tidak valid
        for (int j = 0; j < 4; j++) indices[j] = -1;
        // Cari 4 titik terdekat
        for (int j = 0, count = 0; j < n && count < 4; j++) {
            if (fabs(years[j] - yq) < 5 && years[j] != yq) { // Ambil tahun dalam jarak 5 tahun
                indices[count++] = j;
            }
        }
        // Jika kurang dari 4 titik, ambil titik pertama yang tersedia
        for (int j = 0, k = 0; k < 4 && j < n; j++) {
            int used = 0;
            for (int m = 0; m < 4; m++) {
                if (indices[m] == j) {
                    used = 1;
                    break;
                }
            }
            if (!used && indices[k] == -1) indices[k++] = j;
        }
        // Salin data ke xs, ys_pop, ys_inet
        for (int k = 0; k < 4; k++) {
            if (indices[k] >= 0 && indices[k] < n) {
                xs[k] = years[indices[k]];
                ys_pop[k] = pop[indices[k]];
                ys_inet[k] = inet[indices[k]];
            } else {
                // Fallback: gunakan titik pertama jika indeks tidak valid
                xs[k] = years[0];
                ys_pop[k] = pop[0];
                ys_inet[k] = inet[0];
            }
        }
        double est_pop, est_inet;
        newton_interpolate(xs, ys_pop, 4, yq, &est_pop);
        newton_interpolate(xs, ys_inet, 4, yq, &est_inet);
        printf("  Tahun %.0f: Internet=%.4f%%, Populasi=%.0f\n", yq, est_inet, est_pop);
    }

    // 2. Regresi polinomial derajat 2 (derajat 3: ada masalah overfitting)
    int deg = 2;
    double coeff_pop[4], coeff_inet[4];
    polynomial_regression(years, pop, n, deg, coeff_pop);
    polynomial_regression(years, inet, n, deg, coeff_inet);

    printf("\nPersamaan polinomial (degree=%d):\n", deg);
    printf("  Populasi: y = ");
    for (int i = deg; i >= 0; i--) {
        printf("%.6g x^%d %c ", coeff_pop[i], i, (i>0?'+':' '));
    }
    printf("\n");

    printf("  Internet: y = ");
    for (int i = deg; i >= 0; i--) {
        printf("%.6g x^%d %c ", coeff_inet[i], i, (i>0?'+':' '));
    }
    printf("\n");

    // 3. Prediksi masa depan
    double pop2030  = eval_polynomial(coeff_pop, deg, 2030);
    double inet2035 = eval_polynomial(coeff_inet, deg, 2035);
    printf("\nPrediksi:\n");
    printf("  Populasi 2030 = %.0f\n", pop2030);
    printf("  Internet 2035 = %.4f%%\n", inet2035);

    // 4. Tulis data untuk plotting (tahun,internet,populasi)
    FILE *fp = fopen("plot_data.csv","w");
    if (fp) {
        fprintf(fp, "year,internet,pop\n");
        for (int i = 0; i < n; i++) {
            fprintf(fp, "%.0f,%.4f,%.0f\n", years[i], inet[i], pop[i]);
        }
        fclose(fp);
        printf("\nFile 'plot_data.csv' siap untuk plotting di Python Colab.\n");
    }

    return EXIT_SUCCESS;
}

int read_full_csv(const char *filename, FullDataset *ds) {
    FILE *f = fopen(filename, "r");
    if (!f) return -1;
    char line[256];
    ds->n = 0;
    // Lewati header
    if (!fgets(line, sizeof(line), f)) { fclose(f); return -1; }
    while (fgets(line, sizeof(line), f) && ds->n < MAX_POINTS) {
        double yr, iuser, pop;
        if (sscanf(line, "%lf,%lf,%lf", &yr, &iuser, &pop) == 3) {
            ds->year[ds->n]       = yr;
            ds->internet[ds->n]   = iuser;
            ds->population[ds->n] = pop;
            ds->n++;
        }
    }
    fclose(f);
    return 0;
}

void newton_interpolate(double *x, double *y, int n, double xq, double *result) {
    // Alokasi untuk tabel divided differences
    double **dd = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        dd[i] = (double *)malloc(n * sizeof(double));
    }
    // Inisialisasi kolom pertama dengan y
    for (int i = 0; i < n; i++) {
        dd[i][0] = y[i];
    }
    // Hitung divided differences
    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dd[i][j] = (dd[i+1][j-1] - dd[i][j-1]) / (x[i+j] - x[i]);
        }
    }
    // Evaluasi polinomial Newton
    double res = dd[0][0];
    double term = 1.0;
    for (int i = 1; i < n; i++) {
        term *= (xq - x[i-1]);
        res += dd[0][i] * term;
    }
    *result = res;
    // Bebaskan memori
    for (int i = 0; i < n; i++) {
        free(dd[i]);
    }
    free(dd);
}

void polynomial_regression(double *x, double *y, int n, int degree, double *coeffs) {
    int d = degree;
    int m = d + 1;

    // Alokasi dinamis untuk array X
    double *X = (double *)malloc((2*d + 1) * sizeof(double));
    if (!X) {
        fprintf(stderr, "Error: gagal alokasi memori untuk X\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < 2*d + 1; i++) {
        X[i] = 0;
        for (int j = 0; j < n; j++) X[i] += pow(x[j], i);
    }

    // Alokasi dinamis untuk matriks B
    double **B = (double **)malloc(m * sizeof(double *));
    if (!B) {
        fprintf(stderr, "Error: gagal alokasi memori untuk B\n");
        free(X);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < m; i++) {
        B[i] = (double *)malloc((m + 1) * sizeof(double));
        if (!B[i]) {
            fprintf(stderr, "Error: gagal alokasi memori untuk B[%d]\n", i);
            for (int k = 0; k < i; k++) free(B[k]);
            free(B);
            free(X);
            exit(EXIT_FAILURE);
        }
    }

    // Inisialisasi matriks B
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            B[i][j] = X[i + j];
        }
        B[i][m] = 0;
        for (int j = 0; j < n; j++) B[i][m] += y[j] * pow(x[j], i);
    }

    // Eliminasi Gauss
    for (int i = 0; i < m; i++) {
        double div = B[i][i];
        if (div == 0) {
            fprintf(stderr, "Error: pembagian oleh nol pada eliminasi Gauss\n");
            for (int k = 0; k < m; k++) free(B[k]);
            free(B);
            free(X);
            exit(EXIT_FAILURE);
        }
        for (int j = i; j <= m; j++) B[i][j] /= div;
        for (int k = 0; k < m; k++) {
            if (k != i) {
                double f = B[k][i];
                for (int j = i; j <= m; j++) B[k][j] -= f * B[i][j];
            }
        }
    }

    // Salin koefisien
    for (int i = 0; i < m; i++) coeffs[i] = B[i][m];

    // Bebaskan memori
    for (int i = 0; i < m; i++) free(B[i]);
    free(B);
    free(X);
}

double eval_polynomial(const double *coeffs, int degree, double x) {
    double s = 0.0;
    for (int i = 0; i <= degree; i++) s += coeffs[i] * pow(x, i);
    return s;
}