# Chain Code Vectorization (Stage 5)

## Modifications
- Reads `../Skeletonization_Images/skeleton.png`.
- Detects **junction pixels** (skeleton pixels with ≥ 3 skeleton neighbors via 3×3 convolution) and removes them so the remaining components are simple polylines with two endpoints (or closed loops).
- Connected-component labels each polyline (`scipy.ndimage.label`, 8-connectivity).
- For each component, finds an endpoint (degree-1 pixel) as the starting point and walks neighbor-by-neighbor, recording the **Freeman 8-direction code** at every step.
- Filters out short fragments (residual noise rather than real fractures).
- Renders two visualizations: colored segments over the original image (`chain_code_visualization.png`) and over a black background (`chain_code_skeleton.png`), with a green dot marking each chain's start point. Writes per-segment id, length, start/end coordinates, and the full code string to `chain_codes.txt`.

## Parameters
| Parameter | Value | Reason |
|-----------|-------|--------|
| Junction definition | `nbr_count >= 3` | A skeleton pixel with 3+ skeleton neighbors is a branch point; removing them guarantees each remaining component is a simple chain. |
| Connectivity | 8-neighborhood (`3×3` ones structure) | Matches Freeman 8-direction encoding. |
| `MIN_LEN` | `10` px | Discards very short fragments (mostly noise / spur tips). Lowering it inflated the segment count with non-fracture clutter; raising it dropped real short fractures. |
| Direction encoding | Freeman 8 with image coords (+y down): `0`=E, `1`=NE, `2`=N, `3`=NW, `4`=W, `5`=SW, `6`=S, `7`=SE | Standard Freeman convention adapted for image (y-down) coordinates. |
| Segment colors | HSV hue swept around the wheel | Distinct colors per segment for visual separation. |

## Outcome
The pipeline produces a compact **vector representation** of the fracture network: each fracture is a single starting `(x, y)` plus a string of digits 0–7. `chain_codes.txt` holds these codes for downstream analysis (length statistics, orientation histograms, etc.), and the visualizations confirm the chains follow the visible lineaments end-to-end. Splitting at junctions before tracing means each chain is unambiguous (no branching choices during the walk), at the cost of representing branching fractures as multiple segments meeting at a former junction pixel — an acceptable tradeoff for clean, well-defined chain codes.
