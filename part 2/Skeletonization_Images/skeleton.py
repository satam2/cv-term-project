"""
Stage 4: Skeletonization.

Reads ../Binary_Masks/binary_mask.png and reduces every connected fracture
region to a 1-pixel-wide centerline using skimage.morphology.skeletonize
(Zhang-Suen / Lee thinning). The result is the topology of the fracture
network: a graph of line segments suitable for vectorization.

Outputs:
  skeleton.png         -- white-on-black skeleton.
  skeleton_overlay.png -- skeleton drawn in red over the original image,
                          for sanity-checking placement.
"""
from pathlib import Path
import cv2
import numpy as np
from skimage.morphology import skeletonize

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "Binary_Masks" / "binary_mask.png"
ORIG = HERE.parent / "DW.ppm"
DST = HERE / "skeleton.png"
OVERLAY = HERE / "skeleton_overlay.png"

mask = cv2.imread(str(SRC), cv2.IMREAD_GRAYSCALE)
if mask is None:
    raise FileNotFoundError(f"Could not read {SRC}")

bool_mask = mask > 127
skel = skeletonize(bool_mask)

skel_u8 = (skel.astype(np.uint8)) * 255
cv2.imwrite(str(DST), skel_u8)

# Red overlay on the original satellite image to show where fractures land.
orig = cv2.imread(str(ORIG), cv2.IMREAD_COLOR)
overlay = orig.copy()
overlay[skel] = (0, 0, 255)  # BGR red
cv2.imwrite(str(OVERLAY), overlay)

print(f"[Skeleton] wrote {DST}")
print(f"[Skeleton] skeleton pixels: {int(skel.sum())}")
print(f"[Skeleton] overlay: {OVERLAY}")
