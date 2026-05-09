/*
 * fracture_detect.cpp
 *
 * Detects linear fractures and X-crossings in satellite/aerial rock imagery.
 * Written in C-style (manual structs, malloc, qsort) using OpenCV 4.x C++ API
 * — the old OpenCV C API (IplImage etc.) was removed in OpenCV 4.
 *
 * Pipeline:
 *   1. Grayscale conversion
 *   2. Gaussian blur           suppress rocky texture noise
 *   3. Canny edge detection    find linear edge features
 *   4. Probabilistic Hough     extract line segments
 *   5. Keep top-N longest      filter to dominant fractures
 *   6. X-crossing detection    pairs with sufficient angular difference
 *                              whose infinite extensions meet inside the image
 *   7. Cluster nearby crossings merge duplicates within a radius
 *
 * Build:
 *   g++ fracture_detect.cpp -o fracture_detect \
 *       $(pkg-config --cflags --libs opencv4) -lm -std=c++11
 * 
 *   clang++ fracture_detect.cpp -o fracture_detect -std=c++11 
 *       $(pkg-config --cflags --libs opencv4)
 *
 * Usage:
 *   ./fracture_detect <input_image> <output_prefix>
 *   e.g. ./fracture_detect DW.ppm out
 *   produces: out_edges.png  out_lines.png  out_crossings.png
 */

#include <opencv2/opencv.hpp>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <algorithm>

/* ── Tunable parameters ─────────────────────────────────────────────────── */
#define BLUR_KERNEL       9      /* Gaussian kernel size (must be odd)        */
#define CANNY_LOW         50.0   /* Canny lower threshold                     */
#define CANNY_HIGH        130.0  /* Canny upper threshold                     */
#define HOUGH_THRESHOLD   110    /* Accumulator votes needed                  */
#define HOUGH_MIN_LEN     150.0  /* Minimum segment length in pixels          */
#define HOUGH_MAX_GAP     25.0   /* Max gap within a segment                  */
#define CROSS_ANGLE_MIN   50.0   /* Min crossing angle (degrees)              */
#define CROSS_ANGLE_MAX   85.0   /* Max crossing angle (degrees)              */
#define CROSS_RADIUS      60.0   /* Cluster radius for merging crossings      */
#define TOP_N_LINES       80     /* Keep only the N longest segments          */
/* ─────────────────────────────────────────────────────────────────────────── */


/* ── Segment ─────────────────────────────────────────────────────────────── */
typedef struct {
    float x1, y1, x2, y2, length;
} Segment;

/* ── Crossing ────────────────────────────────────────────────────────────── */
typedef struct {
    float x, y;
    int   count;
} Cross;

/* ── Helpers ─────────────────────────────────────────────────────────────── */
static float seg_length(float x1, float y1, float x2, float y2) {
    return sqrtf((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
}

/* Angle of segment in [0, 180) degrees */
static float seg_angle(const Segment *s) {
    float a = (float)(atan2f(s->y2-s->y1, s->x2-s->x1) * 180.0f / (float)M_PI);
    if (a <   0.0f) a += 180.0f;
    if (a >= 180.0f) a -= 180.0f;
    return a;
}

/* Infinite-line intersection. Returns 1 and fills ox/oy if not parallel. */
static int line_intersect(const Segment *A, const Segment *B,
                           float *ox, float *oy) {
    float x1=A->x1, y1=A->y1, x2=A->x2, y2=A->y2;
    float x3=B->x1, y3=B->y1, x4=B->x2, y4=B->y2;
    float d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
    if (fabsf(d) < 1e-6f) return 0;
    float t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / d;
    *ox = x1 + t*(x2-x1);
    *oy = y1 + t*(y2-y1);
    return 1;
}

/* qsort comparator: descending length */
static int cmp_desc(const void *a, const void *b) {
    float la = ((const Segment *)a)->length;
    float lb = ((const Segment *)b)->length;
    return (lb > la) - (lb < la);
}

/* ── Main ────────────────────────────────────────────────────────────────── */
int main(int argc, char **argv) {

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_image> <output_prefix>\n", argv[0]);
        return 1;
    }
    const char *input  = argv[1];
    const char *prefix = argv[2];

    /* Output filenames */
    char fname_edges[512], fname_lines[512], fname_cross[512];
    snprintf(fname_edges, sizeof(fname_edges), "%s_edges.png",     prefix);
    snprintf(fname_lines, sizeof(fname_lines), "%s_lines.png",     prefix);
    snprintf(fname_cross, sizeof(fname_cross), "%s_crossings.png", prefix);

    /* ── 1. Load ─────────────────────────────────────────────────────────── */
    cv::Mat src = cv::imread(input, cv::IMREAD_COLOR);
    if (src.empty()) { fprintf(stderr, "Cannot open: %s\n", input); return 1; }
    int W = src.cols, H = src.rows;
    printf("Image: %d x %d\n", W, H);

    /* ── 2. Grayscale ────────────────────────────────────────────────────── */
    cv::Mat gray;
    cv::cvtColor(src, gray, cv::COLOR_BGR2GRAY);

    /* ── 3. Gaussian blur ────────────────────────────────────────────────── */
    cv::Mat blurred;
    cv::GaussianBlur(gray, blurred, cv::Size(BLUR_KERNEL, BLUR_KERNEL), 0);

    /* ── 4. Canny ────────────────────────────────────────────────────────── */
    cv::Mat edges;
    cv::Canny(blurred, edges, CANNY_LOW, CANNY_HIGH);
    cv::imwrite(fname_edges, edges);
    printf("Saved: %s\n", fname_edges);

    /* ── 5. Probabilistic Hough ──────────────────────────────────────────── */
    std::vector<cv::Vec4i> hlines;
    cv::HoughLinesP(edges, hlines,
                    1, CV_PI/180,
                    HOUGH_THRESHOLD,
                    HOUGH_MIN_LEN,
                    HOUGH_MAX_GAP);

    int n_hough = (int)hlines.size();
    printf("Hough segments found: %d\n", n_hough);

    /* ── 6. Convert to plain C array of Segment ──────────────────────────── */
    Segment *segs = (Segment *)malloc(n_hough * sizeof(Segment));
    if (!segs) { fprintf(stderr, "OOM\n"); return 1; }

    for (int i = 0; i < n_hough; i++) {
        segs[i].x1 = (float)hlines[i][0];
        segs[i].y1 = (float)hlines[i][1];
        segs[i].x2 = (float)hlines[i][2];
        segs[i].y2 = (float)hlines[i][3];
        segs[i].length = seg_length(segs[i].x1, segs[i].y1,
                                    segs[i].x2, segs[i].y2);
    }

    /* ── 7. Sort by length, keep top N ──────────────────────────────────── */
    qsort(segs, n_hough, sizeof(Segment), cmp_desc);
    int n_segs = (n_hough < TOP_N_LINES) ? n_hough : TOP_N_LINES;
    printf("Segments kept (top-%d): %d\n", TOP_N_LINES, n_segs);

    /* ── 8. Draw kept lines ───────────────────────────────────────────────── */
    cv::Mat vis_lines = src.clone();
    for (int i = 0; i < n_segs; i++) {
        cv::line(vis_lines,
                 cv::Point((int)segs[i].x1, (int)segs[i].y1),
                 cv::Point((int)segs[i].x2, (int)segs[i].y2),
                 cv::Scalar(0, 220, 255), 2, cv::LINE_AA);
    }
    cv::imwrite(fname_lines, vis_lines);
    printf("Saved: %s\n", fname_lines);

    /* ── 9. Find X-crossings ─────────────────────────────────────────────── */
    int max_raw = n_segs * n_segs;
    float *raw_x = (float *)malloc(max_raw * sizeof(float));
    float *raw_y = (float *)malloc(max_raw * sizeof(float));
    if (!raw_x || !raw_y) { fprintf(stderr, "OOM\n"); return 1; }

    int n_raw = 0;
    for (int i = 0; i < n_segs; i++) {
        for (int j = i+1; j < n_segs; j++) {
            float a1   = seg_angle(&segs[i]);
            float a2   = seg_angle(&segs[j]);
            float diff = fabsf(a1 - a2);
            if (diff > 90.0f) diff = 180.0f - diff;

            if (diff < CROSS_ANGLE_MIN || diff > CROSS_ANGLE_MAX) continue;

            float ix, iy;
            if (!line_intersect(&segs[i], &segs[j], &ix, &iy)) continue;

            /* Must land inside the image */
            if (ix < 0 || ix >= W || iy < 0 || iy >= H) continue;

            raw_x[n_raw] = ix;
            raw_y[n_raw] = iy;
            n_raw++;
        }
    }
    printf("Raw crossings: %d\n", n_raw);

    /* ── 10. Cluster crossings ───────────────────────────────────────────── */
    Cross *crosses = (Cross *)malloc(n_raw * sizeof(Cross));
    int   *used    = (int   *)calloc(n_raw, sizeof(int));
    if (!crosses || !used) { fprintf(stderr, "OOM\n"); return 1; }
    int n_crosses = 0;

    for (int i = 0; i < n_raw; i++) {
        if (used[i]) continue;
        float sx = raw_x[i], sy = raw_y[i];
        int   cnt = 1;
        used[i] = 1;
        for (int j = i+1; j < n_raw; j++) {
            if (used[j]) continue;
            float dx = raw_x[i] - raw_x[j];
            float dy = raw_y[i] - raw_y[j];
            if (sqrtf(dx*dx + dy*dy) < CROSS_RADIUS) {
                sx += raw_x[j]; sy += raw_y[j];
                cnt++; used[j] = 1;
            }
        }
        crosses[n_crosses].x     = sx / cnt;
        crosses[n_crosses].y     = sy / cnt;
        crosses[n_crosses].count = cnt;
        n_crosses++;
    }
    printf("Merged crossings: %d\n", n_crosses);

    /* ── 11. Draw crossings ──────────────────────────────────────────────── */
    cv::Mat out = vis_lines.clone();
    for (int i = 0; i < n_crosses; i++) {
        cv::Point c((int)crosses[i].x, (int)crosses[i].y);
        /* Red circle */
        cv::circle(out, c, 14, cv::Scalar(0,0,255), 2, cv::LINE_AA);
        /* Yellow crosshair */
        cv::line(out, cv::Point(c.x-11,c.y), cv::Point(c.x+11,c.y),
                 cv::Scalar(0,255,255), 2, cv::LINE_AA);
        cv::line(out, cv::Point(c.x,c.y-11), cv::Point(c.x,c.y+11),
                 cv::Scalar(0,255,255), 2, cv::LINE_AA);
    }
    cv::imwrite(fname_cross, out);
    printf("Saved: %s\n", fname_cross);

    /* ── Cleanup ─────────────────────────────────────────────────────────── */
    free(segs); free(raw_x); free(raw_y); free(crosses); free(used);

    return 0;
}