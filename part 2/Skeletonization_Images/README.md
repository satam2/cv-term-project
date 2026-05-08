# Skeletonization (Stage 4)

## Modifications
- Reads `../Binary_Masks/binary_mask.png`, binarizes it with a `> 127` threshold, and reduces every connected fracture region to a 1-pixel-wide centerline using `skimage.morphology.skeletonize` (Zhang-Suen / Lee thinning).
- Also produces `skeleton_overlay.png`: the skeleton drawn in red on top of the original `DW.ppm` for visual sanity-checking of placement.

## Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| Binarization threshold | `> 127` | Mask is already 0/255 from Stage 3; 127 is the safe midpoint. |
| `skeletonize` algorithm | default (Zhang-Suen / Lee) | Topology-preserving thinning suitable for line networks. No tuning knobs needed. |
| Overlay color | `(0, 0, 255)` BGR (red) | Maximum contrast against the satellite imagery's gray/green palette. |

## Outcome
`skeleton.png` collapses each fracture region into a single-pixel-wide curve, exposing the **topology** of the fracture network (a graph of segments and junctions). The overlay confirms the skeleton lies along the visible lineaments rather than along noise. This 1-px representation is the input the chain-code stage requires — without it, Freeman 8-connectivity tracing would be ambiguous.
