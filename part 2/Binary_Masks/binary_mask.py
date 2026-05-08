"""
Stage 3: Binary mask of fracture features.

Reads ../CLAHE_Images/clahe.png (denoised + contrast-enhanced) and
produces a binary mask where foreground (255) == "this pixel belongs
to a fracture-like linear feature".

Approach (uses scikit-image + OpenCV libraries, as requested):
  1. Sato tubeness / ridge filter (skimage.filters.sato) tuned for dark
     ridges (black_ridges=True) and for ridge widths typical of fractures.
     This boosts elongated line-like structures and suppresses blobby
     vegetation.
  2. Percentile threshold of the ridge response (top ~3% of pixels) to
     pull out only the strongest line-like evidence. Otsu was too lenient
     here because the distribution is heavily right-skewed rather than
     bimodal.
  3. Morphological close to bridge tiny gaps along a fracture.
  4. remove_small_objects to drop residual vegetation specks.

Output: binary_mask.png (uint8, 0 / 255).
"""
from pathlib import Path
import cv2
import numpy as np
from skimage.filters import sato
from skimage.morphology import remove_small_objects, closing, disk

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "CLAHE_Images" / "clahe.png"
DST = HERE / "binary_mask.png"
RIDGE_DST = HERE / "ridge_response.png"

img = cv2.imread(str(SRC), cv2.IMREAD_GRAYSCALE)
if img is None:
    raise FileNotFoundError(f"Could not read {SRC}")

img_f = img.astype(np.float32) / 255.0

# Sato is a Hessian-based ridge filter. sigmas controls the line widths it
# responds to (in pixels). Fractures here are roughly 1-4 px wide, so we
# scan a small range. black_ridges=True targets dark linear features on a
# lighter background, which matches typical lineament appearance.
ridge = sato(img_f, sigmas=range(1, 4), black_ridges=True)
ridge_norm = ridge / (ridge.max() + 1e-12)

cv2.imwrite(str(RIDGE_DST), (ridge_norm * 255).astype(np.uint8))

# Top-percentile threshold: keep only the strongest 3% ridge responses.
# Tunable knob -- raise PCT to be stricter, lower to keep more.
PCT = 97.0
thr = np.percentile(ridge_norm, PCT)
mask = ridge_norm > thr

mask = closing(mask, disk(1))

# Drop short noise components. max_size in skimage 0.26+ removes objects
# of size <= max_size (replaces deprecated min_size with same semantics).
mask = remove_small_objects(mask, max_size=80, connectivity=2)

out = (mask.astype(np.uint8)) * 255
cv2.imwrite(str(DST), out)
print(f"[BinaryMask] wrote {DST}")
print(f"[BinaryMask] threshold @ {PCT:.1f}th percentile = {thr:.4f}")
print(f"[BinaryMask] foreground pixels: {int(mask.sum())} "
      f"({100.0 * mask.mean():.2f}% of image)")
print(f"[BinaryMask] ridge-response preview: {RIDGE_DST}")
