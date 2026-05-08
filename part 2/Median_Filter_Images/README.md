# Median Filter (Stage 1)

## Modifications
- Reads `../DW.ppm`, converts BGR → grayscale, applies a single-pass median filter via `cv2.medianBlur`.
- This stage runs **before** CLAHE so contrast enhancement does not amplify salt-and-pepper noise.

## Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| `ksize` | `5` | Large enough to kill isolated speckle and small vegetation texture; small enough to keep linear fracture edges sharp. `ksize=3` left too much speckle; `ksize=7` started to smear thin fractures. |

## Outcome
Output `median.png` is visually nearly identical to the source but with isolated bright/dark pixels and small high-frequency vegetation texture removed. Thin linear features (the fractures we care about) are preserved, which is exactly the input shape CLAHE needs.
