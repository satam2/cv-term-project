# Fracture & X-Crossing Detector

Detects linear fractures and X-crossings in satellite or aerial rock imagery using OpenCV.

## Build

```bash
clang++ fracture_detect.cpp -o fracture_detect -std=c++11 $(pkg-config --cflags --libs opencv4)
```

## Run

```bash
./fracture_detect <input_image> <output_prefix>
```

Example:
```bash
./fracture_detect DW.ppm output
```

Produces three output PNGs:
- `output_edges.png` — raw Canny edge map
- `output_lines.png` — detected fracture lines drawn in cyan
- `output_crossings.png` — lines + red circles at X-junctions

---

## Tunable Parameters

All parameters are defined at the top of `fracture_detect.cpp`. Current value below is what we started with. Later on changed these parameters to see new results.

| Parameter | Default | What it controls | Lower value | Higher value |
|---|---|---|---|---|
| `BLUR_KERNEL` | `9` | Gaussian blur kernel size (must be odd) | More detail, more noise | Smoother, less noise |
| `CANNY_LOW` | `50.0` | Canny weak edge threshold | Detects more edges | Misses faint edges |
| `CANNY_HIGH` | `130.0` | Canny strong edge threshold | More edges found | Only strongest edges |
| `HOUGH_THRESHOLD` | `150` | Votes needed to declare a line | More lines found | Fewer, stronger lines |
| `HOUGH_MIN_LEN` | `150.0` | Minimum segment length in pixels | Catches shorter fractures | Only long fractures |
| `HOUGH_MAX_GAP` | `25.0` | Max gap allowed within one segment | Segments stay broken | Connects dotted lines |
| `CROSS_ANGLE_MIN` | `35.0` | Min angle between two lines to count as X (degrees) | Allows shallow crossings | Only steep crossings |
| `CROSS_ANGLE_MAX` | `90.0` | Max angle between two lines to count as X (degrees) | Excludes near-perpendicular | Includes all angles |
| `CROSS_RADIUS` | `60.0` | Clustering radius — merges nearby crossings (pixels) | More separate X markers | Fewer, broader clusters |
| `TOP_N_LINES` | `80` | Keep only the N longest Hough segments | Fewer dominant lines | More lines including subtle ones |

---
