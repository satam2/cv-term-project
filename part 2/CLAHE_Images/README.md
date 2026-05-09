# CLAHE (Stage 2)

## Modifications
- Reads `../Median_Filter_Images/median.png` and applies Contrast Limited Adaptive Histogram Equalization (`cv2.createCLAHE`).
- Run **after** the median filter so noise that would otherwise be boosted has already been suppressed.

## Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| `clipLimit` | `2.0` | Standard moderate clipping. Higher values (e.g. `4.0`) over-amplified background texture and made vegetation compete with fractures; lower values (`1.0`) left fractures too faint. |
| `tileGridSize` | `(8, 8)` | Tiles are small enough to enhance local fracture contrast in shaded/illuminated regions independently, but large enough to avoid blocky tile artifacts at boundaries. |

## Outcome
`clahe.png` shows clearly higher local contrast: dark fracture lineaments are now visibly darker than their immediate neighborhood, regardless of whether they sit in a bright or shadowed part of the scene. This makes the next stage (Sato ridge filter + thresholding) much more reliable, since fractures now have a stronger local signal-to-background ratio.
