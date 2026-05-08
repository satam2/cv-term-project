"""
Stage 1: Median filter for noise reduction.

Reads ../DW.ppm directly, converts to grayscale, and applies a median
filter to remove salt-and-pepper noise and small high-frequency texture
(e.g. isolated vegetation pixels) while preserving the edges of linear
fracture features. Running median BEFORE CLAHE prevents the contrast
enhancement step from amplifying speckle noise.

Output: median.png (grayscale, in this folder).
"""
from pathlib import Path
import cv2

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "DW.ppm"
DST = HERE / "median.png"

img_bgr = cv2.imread(str(SRC), cv2.IMREAD_COLOR)
if img_bgr is None:
    raise FileNotFoundError(f"Could not read {SRC}")

gray = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2GRAY)

# ksize=5 is a good balance: kills isolated speckle / small vegetation
# texture without smearing the linear fractures we want to preserve.
denoised = cv2.medianBlur(gray, ksize=5)

cv2.imwrite(str(DST), denoised)
print(f"[Median] wrote {DST}  shape={denoised.shape}")
