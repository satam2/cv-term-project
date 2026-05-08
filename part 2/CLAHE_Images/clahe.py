"""
Stage 2: CLAHE preprocessing.

Reads ../Median_Filter_Images/median.png (already denoised) and applies
Contrast Limited Adaptive Histogram Equalization to enhance local
contrast so faint linear fracture features become more separable from
background. Running CLAHE AFTER the median filter avoids amplifying
salt-and-pepper noise.

Output: clahe.png (grayscale, in this folder).
"""
from pathlib import Path
import cv2

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "Median_Filter_Images" / "median.png"
DST = HERE / "clahe.png"

gray = cv2.imread(str(SRC), cv2.IMREAD_GRAYSCALE)
if gray is None:
    raise FileNotFoundError(f"Could not read {SRC}")

clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
enhanced = clahe.apply(gray)

cv2.imwrite(str(DST), enhanced)
print(f"[CLAHE] wrote {DST}  shape={enhanced.shape}")
