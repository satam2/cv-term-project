/*
 * main.c
 * Implements a Matrix type, grayscale image I/O, convolution, Sobel, and
 * Canny edge detection for PGM/PPM input images.
 *
 * Compile:  gcc -Wall -Wextra -O2 -o main main.c -lm
 * Usage:    ./main
 * Outputs:  canny_no_gaussian.pgm
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// define Matrix struct
typedef struct {
    int rows;
    int cols;
    double *data;
} Matrix;

typedef Matrix Image;

Matrix make_gaussian_kernel_5x5(void);
static Matrix make_sobel_x_kernel(void);
static Matrix make_sobel_y_kernel(void);
static Matrix apply_hysteresis_threshold(Matrix suppressed, double low_threshold, double high_threshold);
Image canny(Image img);

// zero-initialized Matrix
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
/* Skip whitespace and comments in a PGM header */
static void pgm_skip(FILE *fp) {
    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '#') {
            while ((c = fgetc(fp)) != EOF && c != '\n') {
            }
        } else if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            /* skip */
        } else {
            ungetc(c, fp);
            break;
        }
    }
}

/* Load a PGM (P2/P5) or PPM (P3/P6) into a greyscale Matrix.
 * PPM colour pixels are converted to greyscale using the standard
 * luminance formula: Y = 0.299*R + 0.587*G + 0.114*B            */
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
    int type = magic[1] - '0'; /* 2=P2, 3=P3, 5=P5, 6=P6 */
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
    /* consume exactly one whitespace byte after maxval */
    fgetc(fp);

    Matrix m = create_matrix(rows, cols);

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double grey;
            if (colour) {
                /* Read R, G, B and convert to greyscale */
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
                    rv = ri;
                    gv = gi;
                    bv = bi;
                }
                grey = 0.299 * rv + 0.587 * gv + 0.114 * bv;
            } else {
                if (binary) {
                    grey = (double)fgetc(fp);
                } else {
                    int v = 0;
                    if (fscanf(fp, "%d", &v) != 1) {
                        v = 0;
                    }
                    grey = (double)v;
                }
            }
            mat_set(&m, r, c, grey);
        }
    }
    fclose(fp);
    return m;
}

/* Save a Matrix as a P5 (binary) PGM.
 * Values are first scaled so that the full [0,255] range is used. */
void pgm_save(const char *filename, const Matrix *m) {
    /* Find min and max */
    double mn = mat_get(m, 0, 0), mx = mn;
    for (int i = 0; i < m->rows * m->cols; i++) {
        if (m->data[i] < mn) {
            mn = m->data[i];
        }
        if (m->data[i] > mx) {
            mx = m->data[i];
        }
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
            if (v < 0) {
                v = 0;
            }
            if (v > 255) {
                v = 255;
            }
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
    int ar = (m2.rows - 1) / 2; // anchor row
    int ac = (m2.cols - 1) / 2; // anchor col

    for (int r = 0; r < m1.rows; r++) {
        for (int c = 0; c < m1.cols; c++) {
            /* Keep the filter fully inside the image; untouched cells stay 0. */
            int top = r - ar;
            int left = c - ac;
            int bottom = top + m2.rows - 1;
            int right = left + m2.cols - 1;

            if (top < 0 || left < 0 || bottom >= m1.rows || right >= m1.cols) {
                continue;
            }

            double sum = 0.0;
            for (int fr = 0; fr < m2.rows; fr++) {
                for (int fc = 0; fc < m2.cols; fc++) {
                    sum += mat_get(&m1, top + fr, left + fc) * mat_get(&m2, fr, fc);
                }
            }

            mat_set(&out, r, c, sum);
        }
    }

    return out;
}

/*--------------------------------------------------
 *                MEDIAN FILTER
 *--------------------------------------------------*/
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    return (da > db) - (da < db);
}

Matrix median_filter(Matrix m1, Matrix m2) {
    Matrix out = create_matrix(m1.rows, m1.cols);
    int ar = (m2.rows - 1) / 2; // anchor row
    int ac = (m2.cols - 1) / 2; // anchor col
    int n = m2.rows * m2.cols;
    double *buf = (double *)malloc((size_t)n * sizeof(double));

    if (!buf) {
        fprintf(stderr, "median_filter: out of memory\n");
        exit(1);
    }

    for (int r = 0; r < m1.rows; r++) {
        for (int c = 0; c < m1.cols; c++) {
            int top = r - ar;
            int left = c - ac;
            int bottom = top + m2.rows - 1;
            int right = left + m2.cols - 1;

            if (top < 0 || left < 0 || bottom >= m1.rows || right >= m1.cols) {
                continue;
            }

            int k = 0;
            for (int fr = 0; fr < m2.rows; fr++) {
                for (int fc = 0; fc < m2.cols; fc++) {
                    buf[k++] = mat_get(&m1, top + fr, left + fc);
                }
            }

            qsort(buf, (size_t)n, sizeof(double), cmp_double);
            mat_set(&out, r, c, buf[n / 2]); /* lower-median for even n */
        }
    }

    free(buf);
    return out;
}

/*--------------------------------------------------
 *                CANNY EDGE DETECTOR
 *--------------------------------------------------*/
Image canny(Image img) {
    Matrix kx = make_sobel_x_kernel();
    Matrix ky = make_sobel_y_kernel();
    Matrix gx = convolve(img, kx);
    Matrix gy = convolve(img, ky);
    Matrix magnitude = create_matrix(img.rows, img.cols);
    Matrix suppressed = create_matrix(img.rows, img.cols);
    double max_value = 0.0;

    for (int r = 0; r < img.rows; r++) {
        for (int c = 0; c < img.cols; c++) {
            double dx = mat_get(&gx, r, c);
            double dy = mat_get(&gy, r, c);
            double mag = sqrt(dx * dx + dy * dy);
            mat_set(&magnitude, r, c, mag);
        }
    }

    for (int r = 1; r < img.rows - 1; r++) {
        for (int c = 1; c < img.cols - 1; c++) {
            double dx = mat_get(&gx, r, c);
            double dy = mat_get(&gy, r, c);
            double angle = atan2(dy, dx) * 180.0 / M_PI;
            double mag = mat_get(&magnitude, r, c);
            double n1 = 0.0;
            double n2 = 0.0;

            if (angle < 0.0) {
                angle += 180.0;
            }

            if ((angle >= 0.0 && angle < 22.5) || angle >= 157.5) {
                n1 = mat_get(&magnitude, r, c - 1);
                n2 = mat_get(&magnitude, r, c + 1);
            } else if (angle < 67.5) {
                n1 = mat_get(&magnitude, r - 1, c + 1);
                n2 = mat_get(&magnitude, r + 1, c - 1);
            } else if (angle < 112.5) {
                n1 = mat_get(&magnitude, r - 1, c);
                n2 = mat_get(&magnitude, r + 1, c);
            } else {
                n1 = mat_get(&magnitude, r - 1, c - 1);
                n2 = mat_get(&magnitude, r + 1, c + 1);
            }

            /* Keep only local maxima along the gradient direction. */
            if (mag >= n1 && mag >= n2) {
                mat_set(&suppressed, r, c, mag);
                if (mag > max_value) {
                    max_value = mag;
                }
            }
        }
    }

    Matrix out = create_matrix(img.rows, img.cols);
    if (max_value > 0.0) {
        double high_threshold = 0.1 * max_value;
        double low_threshold = 0.05 * high_threshold;
        /* Keep weak responses only when they connect to a strong edge. */
        Matrix hysteresis = apply_hysteresis_threshold(suppressed, low_threshold, high_threshold);
        free_matrix(&out);
        out = hysteresis;
    }

    free_matrix(&kx);
    free_matrix(&ky);
    free_matrix(&gx);
    free_matrix(&gy);
    free_matrix(&magnitude);
    free_matrix(&suppressed);
    return out;
}

static Matrix apply_hysteresis_threshold(Matrix suppressed, double low_threshold, double high_threshold) {
    Matrix out = create_matrix(suppressed.rows, suppressed.cols);
    int total = suppressed.rows * suppressed.cols;
    int *stack = (int *)malloc((size_t)total * sizeof(int));
    int top = 0;

    if (!stack) {
        fprintf(stderr, "apply_hysteresis_threshold: out of memory\n");
        exit(EXIT_FAILURE);
    }

    for (int r = 0; r < suppressed.rows; r++) {
        for (int c = 0; c < suppressed.cols; c++) {
            if (mat_get(&suppressed, r, c) >= high_threshold) {
                mat_set(&out, r, c, 255.0);
                stack[top++] = r * suppressed.cols + c;
            }
        }
    }

    while (top > 0) {
        int index = stack[--top];
        int r = index / suppressed.cols;
        int c = index % suppressed.cols;

        for (int dr = -1; dr <= 1; dr++) {
            for (int dc = -1; dc <= 1; dc++) {
                int nr;
                int nc;

                if (dr == 0 && dc == 0) {
                    continue;
                }

                nr = r + dr;
                nc = c + dc;
                if (nr < 0 || nr >= suppressed.rows || nc < 0 || nc >= suppressed.cols) {
                    continue;
                }

                if (mat_get(&out, nr, nc) == 255.0) {
                    continue;
                }

                if (mat_get(&suppressed, nr, nc) >= low_threshold) {
                    mat_set(&out, nr, nc, 255.0);
                    stack[top++] = nr * suppressed.cols + nc;
                }
            }
        }
    }

    free(stack);
    return out;
}

/*--------------------------------------------------
 *                  HELPERS
 *--------------------------------------------------*/
/* 3x3 normalised averaging (box) filter */
Matrix make_averaging_kernel_3x3(void) {
    Matrix k = create_matrix(3, 3);
    for (int i = 0; i < 9; i++) {
        k.data[i] = 1.0 / 9.0;
    }
    return k;
}

/* 5x5 all-ones matrix used as the "shape" for the median window */
Matrix make_median_kernel_5x5(void) {
    Matrix k = create_matrix(5, 5);
    for (int i = 0; i < 25; i++) {
        k.data[i] = 1.0;
    }
    return k;
}

Matrix make_gaussian_kernel_5x5(void) {
    Matrix k = create_matrix(5, 5);
    static const double values[25] = {
         2,  4,  5,  4,  2,
         4,  9, 12,  9,  4,
         5, 12, 15, 12,  5,
         4,  9, 12,  9,  4,
         2,  4,  5,  4,  2
    };

    for (int i = 0; i < 25; i++) {
        k.data[i] = values[i] / 159.0;
    }
    return k;
}

static Matrix make_sobel_x_kernel(void) {
    Matrix k = create_matrix(3, 3);
    static const double values[9] = {
        -1, 0, 1,
        -2, 0, 2,
        -1, 0, 1
    };

    for (int i = 0; i < 9; i++) {
        k.data[i] = values[i];
    }
    return k;
}

static Matrix make_sobel_y_kernel(void) {
    Matrix k = create_matrix(3, 3);
    static const double values[9] = {
         1,  2,  1,
         0,  0,  0,
        -1, -2, -1
    };

    for (int i = 0; i < 9; i++) {
        k.data[i] = values[i];
    }
    return k;
}

/* Load one image and write its Canny edge map to disk. */
void edgeDetection(char *inputFilename, char *cannyFilename) {
    Image img = pgm_load(inputFilename);

    Matrix gaussian = make_gaussian_kernel_5x5();
    Image blurred = convolve(img, gaussian);

    Image canny_img = canny(blurred);

    pgm_save(cannyFilename, &canny_img);

    free_matrix(&img);
    free_matrix(&canny_img);
    free_matrix(&gaussian);
    free_matrix(&blurred);
}

/*--------------------------------------------------
 *                  MAIN
 *--------------------------------------------------*/
int main(void) {
    char inputFilename[] = "blur_noise_and_contrast.ppm";
    char cannyFilename[] = "canny_with_gaussian.pgm";

    edgeDetection(inputFilename, cannyFilename);
    printf("Canny image saved to '%s'\n", cannyFilename);
    return EXIT_SUCCESS;
}
