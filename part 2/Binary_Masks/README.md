# Binary Mask (Stage 3)

## Modifications
- Reads `../CLAHE_Images/clahe.png`.
- Applies **Sato tubeness ridge filter** (`skimage.filters.sato`) to highlight elongated dark line-like structures.
- Thresholds the ridge response by **percentile** rather than Otsu (the response distribution is heavily right-skewed, not bimodal, so Otsu was too lenient).
- Cleans the mask with morphological **closing** to bridge tiny gaps along a fracture, then **`remove_small_objects`** to drop residual vegetation specks.
- Also writes `ridge_response.png` as a debug preview of the raw filter output.

## Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| `sigmas` | `range(1, 4)` | Fractures in `DW.ppm` are roughly 1–4 px wide; scanning this range catches both thin and slightly thicker lineaments. |
| `black_ridges` | `True` | Fractures appear as dark lines on a lighter background. |
| `PCT` (percentile threshold) | `97.0` | Keeps only the strongest ~3% of ridge responses. Lower (`95`) leaked vegetation; higher (`99`) lost real fractures. |
| `closing` disk radius | `1` | Bridges 1-pixel gaps without fattening lines. |
| `remove_small_objects` `max_size` | `80` | Drops short noisy components while keeping real (longer) fracture fragments. `connectivity=2` (8-connectivity). |

## Outcome
`binary_mask.png` is a clean foreground/background image where white pixels trace the fracture lineaments. Most blobby vegetation and uniform-texture regions are correctly suppressed, while linear features survive as connected (often slightly broken) line segments. The mask is dense enough to skeletonize but sparse enough that the skeleton will not be dominated by noise.
