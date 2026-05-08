#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "netpbm.h"

#define FEATURES 6
#define FILTER_SIZE 5

typedef struct {
    double f[FEATURES];
    int row;
    int col;
    int label;
} DataPoint;

typedef struct {
    int row;
    int col;
} Point;

Image labelDarkComponents(Image inputImg, int threshold, int minSize) {
    int h = inputImg.height;
    int w = inputImg.width;

    Image output = createImage(h, w);

    int **labels = malloc(h * sizeof(int *));
    for (int r = 0; r < h; r++) {
        labels[r] = calloc(w, sizeof(int));
    }

    Point *queue = malloc(h * w * sizeof(Point));

    int colors[8][3] = {
        {255, 0, 0},
        {0, 255, 0},
        {0, 0, 255},
        {255, 255, 0},
        {255, 0, 255},
        {0, 255, 255},
        {255, 128, 0},
        {128, 0, 255}
    };

    int currentLabel = 1;

    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {

            /*
               Dark pixel condition.
               Lower threshold = only very dark lakes.
               Higher threshold = lakes + shadows.
            */
            if (inputImg.map[r][c].i > threshold || labels[r][c] != 0) {
                continue;
            }

            int front = 0;
            int back = 0;

            queue[back++] = (Point){r, c};
            labels[r][c] = currentLabel;

            int count = 0;

            while (front < back) {
                Point p = queue[front++];
                count++;

                int dr[4] = {-1, 1, 0, 0};
                int dc[4] = {0, 0, -1, 1};

                for (int k = 0; k < 4; k++) {
                    int rr = p.row + dr[k];
                    int cc = p.col + dc[k];

                    if (rr < 0 || rr >= h || cc < 0 || cc >= w) {
                        continue;
                    }

                    if (labels[rr][cc] != 0) {
                        continue;
                    }

                    if (inputImg.map[rr][cc].i <= threshold) {
                        labels[rr][cc] = currentLabel;
                        queue[back++] = (Point){rr, cc};
                    }
                }
            }

            if (count < minSize) {
                for (int i = 0; i < back; i++) {
                    labels[queue[i].row][queue[i].col] = -1;
                }
            } else {
                currentLabel++;
            }
        }
    }

    /*
       Draw result.
       Dark components become colored.
       Everything else becomes original grayscale.
    */
    for (int r = 0; r < h; r++) {
        for (int c = 0; c < w; c++) {
            int label = labels[r][c];

            if (label > 0) {
                int colorIndex = label % 8;

                setPixel(
                    output,
                    r,
                    c,
                    colors[colorIndex][0],
                    colors[colorIndex][1],
                    colors[colorIndex][2],
                    255
                );
            } else {
                int gray = inputImg.map[r][c].i;
                setPixel(output, r, c, gray, gray, gray, gray);
            }
        }
    }

    for (int r = 0; r < h; r++) {
        free(labels[r]);
    }

    free(labels);
    free(queue);

    return output;
}

void makeLawsFilter(const int v1[5], const int v2[5], int filter[5][5]) {
    for (int r = 0; r < 5; r++) {
        for (int c = 0; c < 5; c++) {
            filter[r][c] = v1[r] * v2[c];
        }
    }
}

double applyFilterAt(Image img, int row, int col, int filter[5][5]) {
    double sum = 0.0;
    int offset = FILTER_SIZE / 2;

    for (int fr = 0; fr < FILTER_SIZE; fr++) {
        for (int fc = 0; fc < FILTER_SIZE; fc++) {
            int rr = row + fr - offset;
            int cc = col + fc - offset;

            if (rr >= 0 && rr < img.height && cc >= 0 && cc < img.width) {
                sum += img.map[rr][cc].i * filter[fr][fc];
            }
        }
    }

    return sum;
}

Matrix computeLawsEnergy(Image img, int filter[5][5], int windowSize) {
    Matrix filtered = createMatrix(img.height, img.width);
    Matrix energy = createMatrix(img.height, img.width);

    int radius = windowSize / 2;

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            filtered.map[r][c] = applyFilterAt(img, r, c, filter);
        }
    }

    for (int r = 0; r < img.height; r++) {
        for (int c = 0; c < img.width; c++) {
            double sum = 0.0;
            int count = 0;

            for (int wr = -radius; wr <= radius; wr++) {
                for (int wc = -radius; wc <= radius; wc++) {
                    int rr = r + wr;
                    int cc = c + wc;

                    if (rr >= 0 && rr < img.height && cc >= 0 && cc < img.width) {
                        sum += fabs(filtered.map[rr][cc]);
                        count++;
                    }
                }
            }

            energy.map[r][c] = sum / count;
        }
    }

    deleteMatrix(filtered);
    return energy;
}

double distanceSquared(double *a, double *b) {
    double sum = 0.0;

    for (int i = 0; i < FEATURES; i++) {
        double d = a[i] - b[i];
        sum += d * d;
    }

    return sum;
}

void normalizeFeatures(DataPoint *points, int n) {
    for (int f = 0; f < FEATURES; f++) {
        double minVal = points[0].f[f];
        double maxVal = points[0].f[f];

        for (int i = 1; i < n; i++) {
            if (points[i].f[f] < minVal) minVal = points[i].f[f];
            if (points[i].f[f] > maxVal) maxVal = points[i].f[f];
        }

        double range = maxVal - minVal;
        if (range == 0.0) range = 1.0;

        for (int i = 0; i < n; i++) {
            points[i].f[f] = (points[i].f[f] - minVal) / range;
        }
    }
}

void kmeans(DataPoint *points, int n, int k, int maxIters) {
    double centroids[k][FEATURES];

    srand((unsigned int)time(NULL));

    for (int cluster = 0; cluster < k; cluster++) {
        int index = rand() % n;

        for (int f = 0; f < FEATURES; f++) {
            centroids[cluster][f] = points[index].f[f];
        }
    }

    for (int iter = 0; iter < maxIters; iter++) {
        int changed = 0;

        for (int i = 0; i < n; i++) {
            int bestLabel = 0;
            double bestDist = distanceSquared(points[i].f, centroids[0]);

            for (int cluster = 1; cluster < k; cluster++) {
                double d = distanceSquared(points[i].f, centroids[cluster]);

                if (d < bestDist) {
                    bestDist = d;
                    bestLabel = cluster;
                }
            }

            if (points[i].label != bestLabel) {
                points[i].label = bestLabel;
                changed = 1;
            }
        }

        double sums[k][FEATURES];
        int counts[k];

        for (int cluster = 0; cluster < k; cluster++) {
            counts[cluster] = 0;

            for (int f = 0; f < FEATURES; f++) {
                sums[cluster][f] = 0.0;
            }
        }

        for (int i = 0; i < n; i++) {
            int label = points[i].label;
            counts[label]++;

            for (int f = 0; f < FEATURES; f++) {
                sums[label][f] += points[i].f[f];
            }
        }

        for (int cluster = 0; cluster < k; cluster++) {
            if (counts[cluster] == 0) continue;

            for (int f = 0; f < FEATURES; f++) {
                centroids[cluster][f] = sums[cluster][f] / counts[cluster];
            }
        }

        if (!changed) break;
    }
}

Image segmentTexture(Image inputImg, int segments) {
    int height = inputImg.height;
    int width = inputImg.width;
    int n = height * width;

    int L5[5] = { 1,  4,  6,  4,  1 };
    int E5[5] = {-1, -2,  0,  2,  1 };
    int S5[5] = {-1,  0,  2,  0, -1 };
    int R5[5] = { 1, -4,  6, -4,  1 };

    int E5L5[5][5];
    int L5E5[5][5];
    int E5E5[5][5];
    int S5S5[5][5];

    makeLawsFilter(E5, L5, E5L5);
    makeLawsFilter(L5, E5, L5E5);
    makeLawsFilter(E5, E5, E5E5);
    makeLawsFilter(S5, S5, S5S5);

    Matrix energy1 = computeLawsEnergy(inputImg, E5L5, 15);
    Matrix energy2 = computeLawsEnergy(inputImg, L5E5, 15);
    Matrix energy3 = computeLawsEnergy(inputImg, E5E5, 15);
    Matrix energy4 = computeLawsEnergy(inputImg, S5S5, 15);

    DataPoint *points = malloc(n * sizeof(DataPoint));

    int index = 0;

    for (int r = 0; r < height; r++) {
        for (int c = 0; c < width; c++) {
            points[index].f[0] = energy1.map[r][c];
            points[index].f[1] = energy2.map[r][c];
            points[index].f[2] = energy3.map[r][c];
            points[index].f[3] = energy4.map[r][c];

            points[index].f[4] = inputImg.map[r][c].i / 255.0;
            points[index].f[5] = 0.5 * ((double)c / width);

            points[index].row = r;
            points[index].col = c;
            points[index].label = 0;

            index++;
        }
    }

    normalizeFeatures(points, n);
    kmeans(points, n, segments, 50);

    Image output = createImage(height, width);

    int colors[4][3] = {
    {255, 0, 0},
    {0, 180, 0},
    {0, 0, 255},
    {255, 255, 0}
    };

    for (int i = 0; i < n; i++) {
        int label = points[i].label % segments;
        int r = points[i].row;
        int c = points[i].col;

        int gray = inputImg.map[r][c].i;

    int red   = (gray + colors[label][0]) / 2;
    int green = (gray + colors[label][1]) / 2;
    int blue  = (gray + colors[label][2]) / 2;

    setPixel(output, r, c, red, green, blue, gray);
    }

    deleteMatrix(energy1);
    deleteMatrix(energy2);
    deleteMatrix(energy3);
    deleteMatrix(energy4);
    free(points);

    return output;
}


int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Usage:\n");
        printf("  %s input.ppm output.ppm segments\n", argv[0]);
        return 1;
    }

    Image input = readImage(argv[1]);
    int segments = atoi(argv[3]);

    Image output = segmentTexture(input, segments);

    writeImage(output, argv[2]);

    deleteImage(input);
    deleteImage(output);

    return 0;
}


// Directions to run:
// cd Texture
// gcc main.c netpbm.c -o main -lm   
// .main DW.ppm texture_output.ppm 4
