/*
 * main.c
 * Sobel edge detection for PGM/PPM input images.
 *
 * Compile:  gcc -Wall -Wextra -O2 -o sobel main.c -lm
 * Usage:    ./sobel
 * Outputs:  sobel.pgm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    int rows;
    int cols;
    double *data;
} Matrix;

typedef Matrix Image;

static Matrix make_sobel_x_kernel(void);
static Matrix make_sobel_y_kernel(void);

Matrix create_matrix(int rows, int cols) {
    Matrix m;
    m.rows = rows;
    m.cols = cols;
    m.data = (double *)calloc((size_t)(rows * cols), sizeof(double));
    if (!m.data) {
        fprintf(stderr, "create_matrix: out of memory");
        exit(EXIT_FAILURE);
    }
    return m;
}

void free_matrix(Matrix *m) {
    free(m->data);
    m->data = NULL;
    m->rows = m->cols = 0;
}

static inline double mat_get(const Matrix *m, int r, int c) {
    return m->data[r * m->cols + c];
}

static inline void mat_set(const Matrix *m, int r, int c, double v) {
    m->data[r * m->cols + c] = v;
}

/*--------------------------------------------------
 *          INPUT IMAGE I/O .PPM OR .PGM
 *--------------------------------------------------*/
static void pgm_skip(FILE *fp) {
    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '#') {
            while ((c = fgetc(fp)) != EOF && c != '\n') {}
        } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            /* skip */
        } else {
            ungetc(c, fp);
            break;
        }
    }
}

Matrix pgm_load(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        fprintf(stderr, "pgm_load: cannot open '%s'\n", filename);
        exit(EXIT_FAILURE);
    }

    char magic[3];
    if (!fgets(magic, sizeof(magic), fp)) {
        fprintf(stderr, "pgm_load: bad header\n");
        exit(EXIT_FAILURE);
    }
    int type = magic[1] - '0';
    if (type != 2 && type != 3 && type != 5 && type != 6) {
        fprintf(stderr, "pgm_load: unsupported format P%d\n", type);
        exit(EXIT_FAILURE);
    }
    int binary = (type == 5 || type == 6);
    int colour = (type == 3 || type == 6);

    pgm_skip(fp);
    int cols, rows, maxval;
    if (fscanf(fp, "%d", &cols) != 1 ||
        (pgm_skip(fp), fscanf(fp, "%d", &rows) != 1) ||
        (pgm_skip(fp), fscanf(fp, "%d", &maxval) != 1)) {
        fprintf(stderr, "pgm_load: malformed header\n");
        exit(EXIT_FAILURE);
    }
    fgetc(fp);

    Matrix m = create_matrix(rows, cols);
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double grey;
            if (colour) {
                double rv, gv, bv;
                if (binary) {
                    rv = (double)fgetc(fp);
                    gv = (double)fgetc(fp);
                    bv = (double)fgetc(fp);
                } else {
                    int ri = 0, gi = 0, bi = 0;
                    if (fscanf(fp, "%d %d %d", &ri, &gi, &bi) != 3) {
                        ri = gi = bi = 0;
                    }
                    rv = ri; gv = gi; bv = bi;
                }
                grey = 0.299 * rv + 0.587 * gv + 0.114 * bv;
            } else {
                if (binary) {
                    grey = (double)fgetc(fp);
                } else {
                    int v = 0;
                    if (fscanf(fp, "%d", &v) != 1) v = 0;
                    grey = (double)v;
                }
            }
            mat_set(&m, r, c, grey);
        }
    }
    fclose(fp);
    return m;
}

void pgm_save(const char *filename, const Matrix *m) {
    double mn = mat_get(m, 0, 0), mx = mn;
    for (int i = 0; i < m->rows * m->cols; i++) {
        if (m->data[i] < mn) mn = m->data[i];
        if (m->data[i] > mx) mx = m->data[i];
    }

    FILE *fp = fopen(filename, "wb");
    if (!fp) {
        fprintf(stderr, "pgm_save: cannot open '%s'\n", filename);
        exit(EXIT_FAILURE);
    }
    fprintf(fp, "P5\n%d %d\n255\n", m->cols, m->rows);

    double range = (mx > mn) ? (mx - mn) : 1.0;
    for (int r = 0; r < m->rows; r++) {
        for (int c = 0; c < m->cols; c++) {
            int v = (int)round((mat_get(m, r, c) - mn) / range * 255.0);
            if (v < 0) v = 0;
            if (v > 255) v = 255;
            fputc(v, fp);
        }
    }
    fclose(fp);
}

/*--------------------------------------------------
 *              CONVOLVE FUNCTION
 *--------------------------------------------------*/
Matrix convolve(Matrix m1, Matrix m2) {
    Matrix out = create_matrix(m1.rows, m1.cols);
    int ar = (m2.rows - 1) / 2;
    int ac = (m2.cols - 1) / 2;

    for (int r = 0; r < m1.rows; r++) {
        for (int c = 0; c < m1.cols; c++) {
            int top = r - ar, left = c - ac;
            int bottom = top + m2.rows - 1, right = left + m2.cols - 1;

            if (top < 0 || left < 0 || bottom >= m1.rows || right >= m1.cols)
                continue;

            double sum = 0.0;
            for (int fr = 0; fr < m2.rows; fr++)
                for (int fc = 0; fc < m2.cols; fc++)
                    sum += mat_get(&m1, top + fr, left + fc) * mat_get(&m2, fr, fc);

            mat_set(&out, r, c, sum);
        }
    }
    return out;
}

/*--------------------------------------------------
 *                SOBEL FILTER
 *--------------------------------------------------*/
static void scale_matrix_to_255(Matrix *m) {
    double mn = m->data[0], mx = m->data[0];
    for (int i = 1; i < m->rows * m->cols; i++) {
        if (m->data[i] < mn) mn = m->data[i];
        if (m->data[i] > mx) mx = m->data[i];
    }
    if (mx <= mn) {
        for (int i = 0; i < m->rows * m->cols; i++) m->data[i] = 0.0;
        return;
    }
    for (int i = 0; i < m->rows * m->cols; i++)
        m->data[i] = (m->data[i] - mn) / (mx - mn) * 255.0;
}

static Matrix make_sobel_x_kernel(void) {
    Matrix k = create_matrix(3, 3);
    static const double values[9] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
    for (int i = 0; i < 9; i++) k.data[i] = values[i];
    return k;
}

static Matrix make_sobel_y_kernel(void) {
    Matrix k = create_matrix(3, 3);
    static const double values[9] = { 1, 2, 1, 0, 0, 0, -1, -2, -1 };
    for (int i = 0; i < 9; i++) k.data[i] = values[i];
    return k;
}

Image sobel(Image img) {
    Matrix kx = make_sobel_x_kernel();
    Matrix ky = make_sobel_y_kernel();
    Matrix gx = convolve(img, kx);
    Matrix gy = convolve(img, ky);
    Matrix out = create_matrix(img.rows, img.cols);

    for (int r = 0; r < img.rows; r++)
        for (int c = 0; c < img.cols; c++) {
            double dx = mat_get(&gx, r, c);
            double dy = mat_get(&gy, r, c);
            mat_set(&out, r, c, sqrt(dx * dx + dy * dy));
        }

    scale_matrix_to_255(&out);
    free_matrix(&kx);
    free_matrix(&ky);
    free_matrix(&gx);
    free_matrix(&gy);
    return out;
}

/*--------------------------------------------------
 *                  MAIN
 *--------------------------------------------------*/
int main(void) {
    char inputFilename[] = "DW.ppm";
    char sobelFilename[] = "sobel.pgm";

    Image img = pgm_load(inputFilename);
    Image sobel_img = sobel(img);
    pgm_save(sobelFilename, &sobel_img);

    free_matrix(&img);
    free_matrix(&sobel_img);

    printf("Sobel image saved to '%s'\n", sobelFilename);
    return EXIT_SUCCESS;
}