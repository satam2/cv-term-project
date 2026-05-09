r"""
Stage 5: Vectorize the skeleton with Freeman 8-direction chain codes.

Reads ../Skeletonization_Images/skeleton.png and:

  1. Splits the skeleton at junction pixels (>= 3 skeleton neighbors) so
     each remaining connected component is a simple polyline with two
     endpoints (or a closed loop).
  2. Walks each polyline from an endpoint, recording the direction taken
     at every step using Freeman 8-connectivity codes.

     Direction lookup (image coords, +y = down):
           3  2  1
           4  *  0
           5  6  7

     Note: y grows downward in image coordinates, so direction 6 is
     "down" and direction 2 is "up".
  3. Filters out tiny segments (<= MIN_LEN pixels) -- these are mostly
     residual noise rather than real fractures.

Outputs:
  chain_codes.txt              -- per-segment: id, length, start (x,y),
                                  end (x,y), and the code string.
  chain_code_visualization.png -- each fracture in a distinct color over
                                  the original image, with a green dot at
                                  the start point of each chain.
  chain_code_skeleton.png      -- same colored segments on a black bg.
"""
from pathlib import Path
import cv2
import numpy as np
from scipy.ndimage import convolve, label

HERE = Path(__file__).resolve().parent
SKEL_SRC = HERE.parent / "Skeletonization_Images" / "skeleton.png"
ORIG = HERE.parent / "DW.ppm"
TXT_DST = HERE / "chain_codes.txt"
VIS_DST = HERE / "chain_code_visualization.png"
SKEL_VIS_DST = HERE / "chain_code_skeleton.png"

MIN_LEN = 10  # discard segments shorter than this many pixels

# (dy, dx) -> Freeman direction code (image coords: +y is down).
DIR_CODE = {
    (0, 1): 0, (-1, 1): 1, (-1, 0): 2, (-1, -1): 3,
    (0, -1): 4, (1, -1): 5, (1, 0): 6, (1, 1): 7,
}
NEIGHBOR_OFFSETS = list(DIR_CODE.keys())

skel = cv2.imread(str(SKEL_SRC), cv2.IMREAD_GRAYSCALE)
if skel is None:
    raise FileNotFoundError(f"Could not read {SKEL_SRC}")
skel_b = (skel > 127).astype(np.uint8)

# Count skeleton neighbors for every skeleton pixel.
kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]], dtype=np.uint8)
nbr_count = convolve(skel_b, kernel, mode="constant", cval=0)
nbr_count = nbr_count * skel_b  # zero out non-skeleton positions

# A junction has >=3 skeleton neighbors. Removing them splits the
# skeleton into simple chains we can walk linearly.
junctions = nbr_count >= 3
chains = skel_b.astype(bool) & ~junctions

# Connected components of the broken-up skeleton.
structure = np.ones((3, 3), dtype=np.uint8)
labels, n_components = label(chains, structure=structure)
print(f"[ChainCode] components after junction split: {n_components}")


def neighbors(y, x, mask):
    h, w = mask.shape
    for dy, dx in NEIGHBOR_OFFSETS:
        ny, nx = y + dy, x + dx
        if 0 <= ny < h and 0 <= nx < w and mask[ny, nx]:
            yield ny, nx, dy, dx


def trace_component(comp_mask):
    """Walk a single connected polyline; return (path, codes)."""
    ys, xs = np.where(comp_mask)
    if len(ys) == 0:
        return [], []

    # Prefer an endpoint (exactly 1 neighbor inside this component) as the
    # starting pixel. If the component is a closed loop, fall back to any
    # pixel.
    start = None
    for y, x in zip(ys, xs):
        deg = sum(1 for _ in neighbors(y, x, comp_mask))
        if deg == 1:
            start = (y, x)
            break
    if start is None:
        start = (int(ys[0]), int(xs[0]))

    path = [start]
    codes = []
    visited = np.zeros_like(comp_mask, dtype=bool)
    visited[start] = True
    cur = start
    while True:
        nxt = None
        for ny, nx, dy, dx in neighbors(cur[0], cur[1], comp_mask):
            if not visited[ny, nx]:
                nxt = (ny, nx, dy, dx)
                break
        if nxt is None:
            break
        ny, nx, dy, dx = nxt
        codes.append(DIR_CODE[(dy, dx)])
        path.append((ny, nx))
        visited[ny, nx] = True
        cur = (ny, nx)
    return path, codes


# Distinct colors for each kept segment (HSV around the wheel).
def color_for(i, total):
    hue = int(180 * i / max(total, 1))
    hsv = np.uint8([[[hue, 220, 255]]])
    bgr = cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR)[0, 0]
    return int(bgr[0]), int(bgr[1]), int(bgr[2])


orig = cv2.imread(str(ORIG), cv2.IMREAD_COLOR)
vis = orig.copy()
skel_vis = np.zeros_like(orig)

segments = []
for comp_id in range(1, n_components + 1):
    comp_mask = labels == comp_id
    if comp_mask.sum() < MIN_LEN:
        continue
    path, codes = trace_component(comp_mask)
    if len(codes) < MIN_LEN - 1:
        continue
    segments.append((path, codes))

print(f"[ChainCode] kept segments (len >= {MIN_LEN}): {len(segments)}")

with open(TXT_DST, "w") as f:
    f.write("# Freeman 8-direction chain codes\n")
    f.write("# Direction lookup (image coords, +y = down):\n")
    f.write("#   3 2 1\n")
    f.write("#   4 . 0\n")
    f.write("#   5 6 7\n\n")
    for i, (path, codes) in enumerate(segments):
        sy, sx = path[0]
        ey, ex = path[-1]
        color = color_for(i, len(segments))
        for (y, x) in path:
            vis[y, x] = color
            skel_vis[y, x] = color
        cv2.circle(vis, (int(sx), int(sy)), 3, (0, 255, 0), -1)
        cv2.circle(skel_vis, (int(sx), int(sy)), 3, (0, 255, 0), -1)
        f.write(f"segment {i:03d}  length={len(codes)}  ")
        f.write(f"start=({sx},{sy})  end=({ex},{ey})\n")
        f.write("  codes: " + "".join(str(c) for c in codes) + "\n\n")

cv2.imwrite(str(VIS_DST), vis)
cv2.imwrite(str(SKEL_VIS_DST), skel_vis)
print(f"[ChainCode] wrote {TXT_DST}")
print(f"[ChainCode] wrote {VIS_DST}")
print(f"[ChainCode] wrote {SKEL_VIS_DST}")
