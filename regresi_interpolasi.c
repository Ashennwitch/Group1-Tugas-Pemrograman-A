// Skeleton Code for Tugas Pemrograman A (Curve Fitting & Interpolation) in C
// Group of 3 students – implement regression, interpolation, CSV I/O, and plotting setup

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// --- Constants and Types --------------------------------------------------
#define MAX_POINTS 1000    // maximum number of data points per series

typedef struct {
    double x[MAX_POINTS];
    double y[MAX_POINTS];
    int n;
} DataSeries;

// --- Function Prototypes --------------------------------------------------

// 1. CSV I/O
int read_csv(const char *filename, DataSeries *ds);
int write_csv(const char *filename, const DataSeries *ds);

// 2. Regression (Least Squares)
void linear_regression(const DataSeries *ds, double *m, double *c);
void polynomial_regression(const DataSeries *ds, int degree, double *coeffs);

// 3. Interpolation
double newton_interpolate(const DataSeries *ds, double x);
double lagrange_interpolate(const DataSeries *ds, double x);

// 4. Plotting helper (generates Gnuplot script)
int generate_gnuplot_script(const char *datafile, const char *output_png);

// 5. Utility
void usage(const char *progname);

// --- Main ------------------------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // 1. Read data
    DataSeries ds;
    if (read_csv(argv[1], &ds) != 0) {
        fprintf(stderr, "Error: cannot read CSV file '%s'\n", argv[1]);
        return EXIT_FAILURE;
    }

    // 2. Handle missing values: detect gaps and choose method
    //    Example: fill missing by regression or interpolation
    //    (Placeholder: assume no missing points for now)

    // 3. Perform linear regression
    double m, c;
    linear_regression(&ds, &m, &c);
    printf("Linear fit: y = %.6f x + %.6f\n", m, c);

    // 4. Perform polynomial regression (e.g., degree 2)
    int degree = 2;
    double coeffs[10] = {0}; // supports up to degree 9
    polynomial_regression(&ds, degree, coeffs);
    printf("Poly fit (deg=%d): ", degree);
    for (int i = 0; i <= degree; i++) {
        printf("a%d=%.6f ", i, coeffs[i]);
    }
    printf("\n");

    // 5. Interpolate at a test point
    double x_query = ds.x[0];
    double y_newton = newton_interpolate(&ds, x_query);
    double y_lagrange = lagrange_interpolate(&ds, x_query);
    printf("Interpolation at x=%.3f: Newton=%.6f, Lagrange=%.6f\n", x_query, y_newton, y_lagrange);

    // 6. Write fitted values for plotting
    DataSeries fit;
    fit.n = ds.n;
    for (int i = 0; i < ds.n; i++) {
        fit.x[i] = ds.x[i];
        fit.y[i] = m * ds.x[i] + c;  // use linear model here
    }
    write_csv("fit_output.csv", &fit);

    // 7. Generate Gnuplot script
    generate_gnuplot_script("fit_output.csv", "fit_plot.png");
    printf("Generated plot script 'plot_fit.gp' and data 'fit_output.csv'\n");

    return EXIT_SUCCESS;
}

// --- Function Stubs --------------------------------------------------------

int read_csv(const char *filename, DataSeries *ds) {
    FILE *fp = fopen(filename, "r");
    if (!fp) return -1;
    ds->n = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp) && ds->n < MAX_POINTS) {
        if (line[0] == '#') continue; // skip comments
        double x, y;
        if (sscanf(line, "%lf , %lf", &x, &y) == 2) {
            ds->x[ds->n] = x;
            ds->y[ds->n] = y;
            ds->n++;
        }
    }
    fclose(fp);
    return 0;
}

int write_csv(const char *filename, const DataSeries *ds) {
    FILE *fp = fopen(filename, "w");
    if (!fp) return -1;
    for (int i = 0; i < ds->n; i++) {
        fprintf(fp, "%.6f,%.6f\n", ds->x[i], ds->y[i]);
    }
    fclose(fp);
    return 0;
}

void linear_regression(const DataSeries *ds, double *m, double *c) {
    double sumx=0, sumy=0, sumxy=0, sumx2=0;
    int n = ds->n;
    for (int i = 0; i < n; i++) {
        sumx  += ds->x[i];
        sumy  += ds->y[i];
        sumxy += ds->x[i] * ds->y[i];
        sumx2 += ds->x[i] * ds->x[i];
    }
    *m = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx*sumx);
    *c = (sumy - (*m)*sumx) / n;
}

void polynomial_regression(const DataSeries *ds, int degree, double *coeffs) {
    // TODO: build and solve normal equations matrix of size (degree+1)x(degree+1)
    // Placeholder: zero coefficients
    for (int i = 0; i <= degree; i++) coeffs[i] = 0.0;
}

double newton_interpolate(const DataSeries *ds, double x) {
    // TODO: compute divided differences and evaluate Newton form
    return 0.0;
}

double lagrange_interpolate(const DataSeries *ds, double x) {
    // TODO: evaluate Lagrange polynomials
    return 0.0;
}

int generate_gnuplot_script(const char *datafile, const char *output_png) {
    FILE *gp = fopen("plot_fit.gp", "w");
    if (!gp) return -1;
    fprintf(gp,
        "set terminal pngcairo size 800,600\n"
        "set output '%s'\n"
        "set title 'Data vs Linear Fit'\n"
        "set xlabel 'X'\n"
        "set ylabel 'Y'\n"
        "plot '%s' using 1:2 with points pt 7 title 'Original', \
"
        "     'fit_output.csv' using 1:2 with lines lw 2 title 'Linear Fit'\n",
        output_png, datafile
    );
    fclose(gp);
    return 0;
}

void usage(const char *progname) {
    printf("Usage: %s data.csv\n", progname);
    printf("  data.csv : input CSV with two columns X,Y\n");
}
