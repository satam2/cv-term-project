"""
Stage 1: CLAHE preprocessing.

Reads ../DW.ppm, converts to grayscale, applies Contrast Limited Adaptive
Histogram Equalization to enhance local contrast so faint linear fracture
features become more separable from background.

Output: clahe.png (grayscale, in this folder).
"""
from pathlib import Path
import cv2

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "DW.ppm"
DST = HERE / "clahe.png"

img_bgr = cv2.imread(str(SRC), cv2.IMREAD_COLOR)
if img_bgr is None:
    raise FileNotFoundError(f"Could not read {SRC}")

gray = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2GRAY)

clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
enhanced = clahe.apply(gray)

cv2.imwrite(str(DST), enhanced)
print(f"[CLAHE] wrote {DST}  shape={enhanced.shape}")
