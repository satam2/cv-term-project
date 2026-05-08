"""
Stage 2: Median filter for noise reduction.

Reads ../CLAHE_Images/clahe.png and applies a median filter to remove
salt-and-pepper noise and small high-frequency texture (e.g. isolated
vegetation pixels) while preserving the edges of linear fracture features.

Output: median.png (grayscale, in this folder).
"""
from pathlib import Path
import cv2

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "CLAHE_Images" / "clahe.png"
DST = HERE / "median.png"

img = cv2.imread(str(SRC), cv2.IMREAD_GRAYSCALE)
if img is None:
    raise FileNotFoundError(f"Could not read {SRC}")

# ksize=5 is a good balance: kills isolated speckle / small vegetation
# texture without smearing the linear fractures we want to preserve.
denoised = cv2.medianBlur(img, ksize=5)

cv2.imwrite(str(DST), denoised)
print(f"[Median] wrote {DST}  shape={denoised.shape}")
