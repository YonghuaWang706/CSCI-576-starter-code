#!/usr/bin/env python3
from dataclasses import dataclass, field
from typing import List, Tuple

import cv2
import numpy as np

# ================= Configuration =================
VISUAL_DEBUG = False
DEBUG_SCALE = 0.6
STEP_SPEED_MS = 0  # 0 = Wait for keypress (Step-by-step mode)
MIN_EDGE_LEN = 10  # Minimum pixels to consider a valid edge
MATCH_THRESHOLD = 3000.0  # Visual MSE threshold
GRADIENT_WEIGHT = 2.0  # Multiplier for gradient-heavy areas


# ================= Data Structures =================

@dataclass
class EdgeFeature:
    pixels: np.ndarray
    grad_mag: np.ndarray  # Gradient magnitude profile
    variance: float
    len: int


@dataclass
class Piece:
    id: int
    image: np.ndarray
    h: int
    w: int
    original_tl: tuple[int, int]
    original_tr: tuple[int, int]
    original_bl: tuple[int, int]
    original_br: tuple[int, int]
    # 4 Edges: 0:Top, 1:Right, 2:Bottom, 3:Left
    edge_features: List[EdgeFeature] = field(default_factory=list)
    placed_x: int = 0
    placed_y: int = 0
    initial_angle: float = 0.0  # Store the rotation angle found in preprocess



@dataclass
class FrontierEdge:
    x: int
    y: int
    length: int
    axis: str  # 'H' (Horizontal) or 'V' (Vertical)
    direction: int  # 1 (Right/Down) or -1 (Left/Up)
    host_piece_id: int

    def __repr__(self):
        dir_str = "+" if self.direction > 0 else "-"
        return f"Edge(x={self.x}, y={self.y}, len={self.length}, {self.axis}{dir_str})"


# ================= Helper Functions =================

def get_gradient_profile(pixels: np.ndarray) -> Tuple[np.ndarray, float]:
    """Computes gradient magnitude profile and variance."""
    if pixels.size == 0: return np.array([]), 0.0
    gray = cv2.cvtColor(pixels[np.newaxis, :, :], cv2.COLOR_BGR2GRAY)
    sobel = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)
    mag = np.abs(sobel).flatten()
    return mag, float(np.std(mag))

def extract_edge_feature(img: np.ndarray, side: int, is_tilted: bool) -> EdgeFeature:
    h, w = img.shape[:2]

    # If tilted, we skip 1px to avoid rotation artifacts (black halos).
    # If straight, we use 0px to preserve every bit of valid texture data.
    margin = 1 if is_tilted else 0

    if side == 0:  # Top
        pixels = img[margin, :, :]
    elif side == 1:  # Right
        pixels = img[:, w - 1 - margin, :]
    elif side == 2:  # Bottom
        pixels = img[h - 1 - margin, :, :]
    elif side == 3:  # Left
        pixels = img[:, margin, :]

    # Gaussian Blur (it helps robustness even on straight edges)
    # Using a slightly smaller kernel for straight pieces is optional,
    # but (3,1) is generally safe for everything.
    # pixels = cv2.GaussianBlur(pixels.reshape(1, -1, 3), (3, 1), 0).reshape(-1, 3)
    if is_tilted:
        # Use (5, 5) to blur in both directions slightly, smoothing out the "grid" noise
        pixels = cv2.GaussianBlur(pixels.reshape(1, -1, 3), (5, 5), 0).reshape(-1, 3)
    else:
        # Keep original light blur for translation integrity
        pixels = cv2.GaussianBlur(pixels.reshape(1, -1, 3), (3, 1), 0).reshape(-1, 3)

    mag, var = get_gradient_profile(pixels)
    return EdgeFeature(pixels.astype(np.float32), mag, var, len(pixels))


def order_points(pts):
    rect = np.zeros((4, 2), dtype="float32")
    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]
    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]
    return rect


def check_overlap(new_rect: Tuple[int, int, int, int], placed_pieces: List[Piece]) -> bool:
    nx, ny, nw, nh = new_rect
    tol = 2  # Tolerance for "touching"
    for p in placed_pieces:
        x_overlap = max(0, min(nx + nw - tol, p.placed_x + p.w - tol) - max(nx + tol, p.placed_x + tol))
        y_overlap = max(0, min(ny + nh - tol, p.placed_y + p.h - tol) - max(ny + tol, p.placed_y + tol))
        if x_overlap > 0 and y_overlap > 0:
            return True
    return False


# ================= Main Logic: Global Solver =================

class PuzzleSolverV2:
    def __init__(self, img_path):
        self.original_img = cv2.imread(img_path)
        if self.original_img is None: raise ValueError("Image not found")

        self.pieces: List[Piece] = []
        self.placed_pieces: List[Piece] = []
        self.frontier: List[FrontierEdge] = []  # Standard List, not Heap
        self.unused_ids = set()
        self.min_x, self.max_x = 0, 0
        self.min_y, self.max_y = 0, 0
        self.seen_variances = []  # Store history (or just sum/count for speed)
        self.dynamic_threshold = 10.0  # Start with a safe, low default
        self.warmup_period = 50  # Don't apply penalties until we've seen 50 comparisons
        self.variance_tracker_count = 0

    def rotate_image(self, image, angle):
        """
        Rotates an image around its center, expanding the canvas to fit.
        """
        if abs(angle) < 1.0: return image  # Optimization

        h, w = image.shape[:2]
        center = (w // 2, h // 2)

        # 1. Get Rotation Matrix
        M = cv2.getRotationMatrix2D(center, angle, 1.0)

        # 2. Calculate New Bounding Box Size
        cos = np.abs(M[0, 0])
        sin = np.abs(M[0, 1])
        new_w = int((h * sin) + (w * cos))
        new_h = int((h * cos) + (w * sin))

        # 3. Adjust Matrix to Center the Image
        M[0, 2] += (new_w / 2) - center[0]
        M[1, 2] += (new_h / 2) - center[1]

        return cv2.warpAffine(image, M, (new_w, new_h), flags=cv2.INTER_LINEAR)

    def preprocess(self):
        gray = cv2.cvtColor(self.original_img, cv2.COLOR_BGR2GRAY)
        _, thresh = cv2.threshold(gray, 10, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        print(f"🔍 Found {len(contours)} contours")
        is_tilted = False
        for i, cnt in enumerate(contours):
            if cv2.contourArea(cnt) < 500: continue
            rect = cv2.minAreaRect(cnt)
            angle = rect[2]  # Get the angle
            # Check if piece is effectively straight (axis-aligned)
            # OpenCV angles can be 0, 90, -90, etc. We check proximity to 90-degree intervals.
            # We treat < 2.0 degrees deviation as "Straight".
            is_straight = (abs(angle) < 2.0 or abs(angle - 90.0) < 2.0 or
                           abs(angle + 90.0) < 2.0 or abs(angle - 180.0) < 2.0)
            is_tilted = not is_straight
            if is_tilted: break
        for i, cnt in enumerate(contours):
            if cv2.contourArea(cnt) < 500: continue

            # NO-TILT PATH (CROP)
            if not is_tilted:
                # When the entire set is NOT tilted, we use bounding rect CROP for pixel perfection
                x, y, w, h = cv2.boundingRect(cnt)
                piece_img = self.original_img[y:y + h, x:x + w]

                p = Piece(id=i, image=piece_img, original_tl=(x, y), original_tr=(x+w, y), original_bl=(x, y+h), original_br=(x+w, y+h), h=h, w=w, initial_angle=0.0)
                # is_tilted=False means margin=0 in extract_edge_feature
                p.edge_features = [extract_edge_feature(piece_img, s, is_tilted=False) for s in range(4)]

            # TILTED PATH (WARP)
            else:

                rect = cv2.minAreaRect(cnt)

                box = order_points(cv2.boxPoints(rect))

                (tl, tr, br, bl) = box
                width = int(max(np.linalg.norm(br - bl), np.linalg.norm(tr - tl)))
                height = int(max(np.linalg.norm(tr - br), np.linalg.norm(tl - bl)))

                dst_pts = np.array([[0, 0], [width - 1, 0], [width - 1, height - 1], [0, height - 1]], dtype="float32")
                M = cv2.getPerspectiveTransform(box, dst_pts)

                # High-Quality Interpolation ===
                # LANCZOS4 keeps edges sharper than the default LINEAR
                warped = cv2.warpPerspective(self.original_img, M, (width, height), flags=cv2.INTER_LANCZOS4)

                dy = tr[1] - tl[1]
                dx = tr[0] - tl[0]
                true_angle = np.degrees(np.arctan2(dy, dx))

                p = Piece(id=i, image=warped, original_tl=tl, original_tr=tr, original_br=br, original_bl=bl ,h=height, w=width, initial_angle=true_angle)
                p.edge_features = [extract_edge_feature(warped, s, is_tilted) for s in range(4)]

            self.pieces.append(p)
            self.unused_ids.add(i)
        print(f"✅ Extracted {len(self.pieces)} pieces.")

    def add_piece_to_frontier(self, piece: Piece, x: int, y: int):
        """Adds 4 edges of a placed piece to the frontier list."""
        # Top (H, -1)
        self.frontier.append(FrontierEdge(x, y, piece.w, 'H', -1, piece.id))
        # Right (V, 1)
        self.frontier.append(FrontierEdge(x + piece.w, y, piece.h, 'V', 1, piece.id))
        # Bottom (H, 1)
        self.frontier.append(FrontierEdge(x, y + piece.h, piece.w, 'H', 1, piece.id))
        # Left (V, -1)
        self.frontier.append(FrontierEdge(x, y, piece.h, 'V', -1, piece.id))

    def get_host_pixels(self, edge: FrontierEdge, offset: int, length: int) -> np.ndarray:
        """Extracts pixel sub-segment from the host piece."""
        host = self.pieces[edge.host_piece_id]

        # Determine Host Side Index based on Edge Geometry
        # H, -1 -> Top (0); V, 1 -> Right (1); H, 1 -> Bottom (2); V, -1 -> Left (3)
        if edge.axis == 'H':
            side = 0 if edge.direction == -1 else 2
        else:
            side = 1 if edge.direction == 1 else 3

        full_pixels = host.edge_features[side].pixels

        # We assume the FrontierEdge x,y matches the Host Side start.
        # But if the frontier edge was a "remainder" split, we need to find WHERE in the host it starts.

        # Calculate Host Side Start Global Coords
        hx, hy = host.placed_x, host.placed_y
        if side == 1: hx += host.w  # Right edge x
        if side == 2: hy += host.h  # Bottom edge y

        # Calculate delta to find start index
        if edge.axis == 'H':
            start_idx = edge.x - hx
        else:
            start_idx = edge.y - hy

        # Extract
        # Safeguard indices
        s = max(0, start_idx + offset)
        e = min(len(full_pixels), s + length)
        return full_pixels[s:e]

    def compute_match_cost(self, pixels_a, pixels_b, edge_variance=0):
        """
        Computes MSE but normalizes it by the edge's variance (texture complexity).
        We trust high-variance matches MORE than low-variance matches.
        """

        if len(pixels_a) != len(pixels_b) or len(pixels_a) == 0:
            return float('inf')

        # adjust variance gradually
        # We record what "normal" variance looks like for these specific segments
        if self.variance_tracker_count < 1000:  # Limit memory usage, stop learning eventually
            self.seen_variances.append(edge_variance)
            self.variance_tracker_count += 1

            # Periodically update the threshold (e.g., every 10 calls)
            if self.variance_tracker_count % 10 == 0:
                # Set threshold to Median (robust against outliers)
                self.dynamic_threshold = np.median(self.seen_variances) * 0.5


        diff = pixels_a - pixels_b
        mse = np.mean(diff ** 2)

        # Variance Clamping (Safety Rail)
        # We cap variance at 200. This stops the Cafe Awning (Var 2000)
        # from getting an unfair advantage.
        # Awning (2000) -> 200. Flower (100) -> 100. Yellow (5) -> 5.
        effective_var = min(edge_variance, 200.0)

        # Variance Normalization ===
        # If the edge has high texture (variance), we tolerate higher MSE.
        # If the edge is flat (low variance), we demand near-zero MSE.
        # Formula: Adjusted_Cost = MSE / (1 + Weight * Variance)

        # Sensitivity factor. Higher = prefer high-texture matches more aggressively.
        TEXTURE_WEIGHT = 5.0

        # adjusted_cost = mse / (1.0 + TEXTURE_WEIGHT * edge_variance)

        if edge_variance < self.dynamic_threshold:
            # PENALTY: This is a low-texture match (like the dark floor).
            adjusted_cost = 2 * mse / (1.0 + TEXTURE_WEIGHT * edge_variance)
        else:
            # REWARD: This is a high-texture match. Trust it.
            # Standard formula: lower the cost slightly for high variance.
            adjusted_cost = mse / (1.0 + TEXTURE_WEIGHT * edge_variance)

        return adjusted_cost


    # ================= LOGIC 1: Zippering (Self-Healing) =================

    def try_zipper_frontier(self):
        """Checks if any two frontier edges are facing each other and touching."""
        for i, e1 in enumerate(self.frontier):
            for j, e2 in enumerate(self.frontier):
                if i >= j: continue
                if e1.axis != e2.axis: continue
                if e1.direction == e2.direction: continue  # Must face opposite ways

                # Check Colinearity and Overlap
                # Vertical: X must match. Y ranges overlap.
                if e1.axis == 'V':
                    if abs(e1.x - e2.x) > 2: continue  # Not touching X
                    overlap_len = min(e1.y + e1.length, e2.y + e2.length) - max(e1.y, e2.y)
                else:  # Horizontal
                    if abs(e1.y - e2.y) > 2: continue  # Not touching Y
                    overlap_len = min(e1.x + e1.length, e2.x + e2.length) - max(e1.x, e2.x)

                if overlap_len > 5:  # Valid overlap found
                    # Check Visual Match (Pixels)
                    # We need to extract the OVERLAPPING segment from both
                    # Simplified: Just grab pixels and check MSE
                    # (In full version, handle offsets. Here assume aligned for Zippering)

                    # For V2 prototype: Just "Zip" if geometry matches perfectly.
                    # This fixes the "Trapped Edge" immediately.
                    print(f"⚡ Zippering Seam between Piece {e1.host_piece_id} and {e2.host_piece_id}")

                    # Remove both (simplified - assumes full overlap for now)
                    # In robust version: Split remaining parts.
                    # Here we just mark them "satisfied" by removing from list.
                    # Care needed: removing by index invalidates loop.
                    # Return instructions to main loop to restart.
                    return [e1, e2]
        return None

    # ================= LOGIC 2: Global Best Match =================

    def find_global_best_match(self):
        """Scans ALL frontier edges vs ALL unused pieces. Supports Virtual Merging."""
        best_cost = float('inf')
        best_move = None  # (piece_id, edges_to_remove, new_edges_to_add, placed_rect)

        # 1. Iterate all Open Slots (Frontier Edges)
        # Note: To support "Virtual Merging" (Stacking), we don't just look at one edge.
        # We look at "Chains" of edges.
        # IMPLEMENTATION: Look at one edge. If candidate is larger, check neighbors.

        sorted_frontier = sorted(self.frontier, key=lambda e: (e.x, e.y))  # Spatial sort helps chaining

        # Calculate current puzzle dimensions for the Compactness Check
        curr_width = self.max_x - self.min_x
        curr_height = self.max_y - self.min_y
        if curr_height == 0: curr_height = 1

        # Determine if we should enforce compactness yet
        # We start enforcing strictly after 50% of pieces are placed.
        total_pieces = len(self.pieces)
        placed_count = len(self.placed_pieces)
        enforce_compactness = placed_count > (total_pieces * 0.5)


        for f_idx, edge in enumerate(sorted_frontier):
            if edge.length < MIN_EDGE_LEN: continue

            # Determine required Candidate Side
            # H, -1 (Up) -> Need Bottom (2)
            # V, 1 (Right) -> Need Left (3)
            # H, 1 (Down) -> Need Top (0)
            # V, -1 (Left) -> Need Right (1)
            if edge.axis == 'H':
                req_cand_side = 2 if edge.direction == -1 else 0
            else:
                req_cand_side = 3 if edge.direction == 1 else 1

            # 2. Iterate all Unused Pieces and select corresponding edge with opposite direction
            for uid in self.unused_ids:
                cand = self.pieces[uid]
                cand_side_len = cand.edge_features[req_cand_side].len

                # === VIRTUAL MERGING LOGIC ===
                # Check if we can form a "Target Slot" starting at 'edge'

                target_edges_chain = [edge]
                current_chain_len = edge.length

                # If Candidate is longer than this edge, look for neighbors
                chain_valid = True

                # Basic tolerance check
                if cand_side_len > current_chain_len + 5:
                    # We need more edges!
                    # Find edge starting where this one ends
                    remaining_needed = cand_side_len - current_chain_len

                    # Naive Search for neighbor (can be optimized with spatial hash)
                    next_start_x = edge.x + (edge.length if edge.axis == 'H' else 0)
                    next_start_y = edge.y + (0 if edge.axis == 'H' else edge.length)

                    found_next = False
                    for next_edge in sorted_frontier:
                        if next_edge is edge: continue
                        if next_edge.axis != edge.axis: continue
                        if next_edge.direction != edge.direction: continue

                        # Check start point
                        if abs(next_edge.x - next_start_x) < 2 and abs(next_edge.y - next_start_y) < 2:
                            target_edges_chain.append(next_edge)
                            current_chain_len += next_edge.length
                            found_next = True
                            break  # Found one neighbor with next_start_x and next_start_y, simplified loop

                    if not found_next:
                        chain_valid = False  # Gap is too big, no neighbor found
                # go to next unused piece
                if not chain_valid: continue
                if cand_side_len > current_chain_len + 5: continue  # Still too big

                # === VISUAL MATCH CHECK ===
                # Compare Candidate Pixels vs Chain of Host Pixels

                # Construct Host Pixels Chain
                # (Assuming start alignment)

                # Quick Offset check: We slide 1px? No, "Global Best" usually assumes "Corner Fit".
                # To save time, we align corners. (Greedy Corner Match)

                # Get Host Pixels from Chain
                # Note: This is complex pixel extraction.
                # SIMPLIFICATION: Compare just the first edge segment.
                # If the first segment matches well, it's a strong candidate.

                h_pixels = self.get_host_pixels(edge, 0, min(edge.length, cand_side_len))
                c_pixels = cand.edge_features[req_cand_side].pixels[:len(h_pixels)]

                # Re-calculate variance for the specific segment being matched
                _, h_var = get_gradient_profile(h_pixels)

                # Pass this variance to the cost function
                cost = self.compute_match_cost(h_pixels, c_pixels, edge_variance=h_var)

                # Check how many neighbors this candidate would touch if placed here
                # We calculate the potential rectangle
                nx, ny = 0, 0
                if edge.axis == 'V':
                    nx = edge.x if edge.direction == 1 else edge.x - cand.w
                    ny = edge.y
                else:
                    nx = edge.x
                    ny = edge.y if edge.direction == 1 else edge.y - cand.h

                # Count neighbors touching the candidate's OTHER 3 sides
                neighbors_found = 0
                check_rect = (nx, ny, cand.w, cand.h)
                tol = 2

                # Scan all placed pieces to see if we touch them
                for p in self.placed_pieces:
                    # Don't count the host of the primary edge we are matching against
                    if p.id == edge.host_piece_id: continue

                    # Check for touch (intersection with tolerance)
                    x_touch = max(0, min(nx + cand.w, p.placed_x + p.w) - max(nx, p.placed_x))
                    y_touch = max(0, min(ny + cand.h, p.placed_y + p.h) - max(ny, p.placed_y))

                    # If we touch significantly on any side
                    if x_touch > 0 and y_touch > 0:  # This means overlap/touch
                        neighbors_found += 1

                # PENALTY RULE:
                # If we only touch 1 neighbor (the primary edge), we are "Cantilevered".
                # This is risky. Apply a penalty to the cost to discourage it
                # unless the match is absolutely perfect.
                if neighbors_found == 0:
                    if len(self.placed_pieces) > 6:
                        cost *= 3.0  # Double the cost (making it harder to be "Best")
                else:
                    cost *= 0.5  # Reward corners (halve the cost, making it preferred)

                ## we favor exact match for edges' lengths after hale pieces have been placed
                if len(self.placed_pieces) > len(self.pieces) * 0.5:
                    if abs(edge.length-cand_side_len > 5):
                        cost *= 1.5

                # COMPACTNESS PENALTY
                if enforce_compactness:
                    # Calculate new bounds if we placed this
                    new_min_x = min(self.min_x, nx)
                    new_max_x = max(self.max_x, nx + cand.w)
                    new_min_y = min(self.min_y, ny)
                    new_max_y = max(self.max_y, ny + cand.h)

                    expanded_x = (new_max_x - new_min_x) > curr_width
                    expanded_y = (new_max_y - new_min_y) > curr_height

                    # 1. General Expansion Penalty (Discourage growing bounds)
                    # We prefer filling holes (no expansion) over growing edges.
                    if expanded_x or expanded_y:
                        cost *= 2

                # =================================

                if cost < best_cost:
                    if(VISUAL_DEBUG):
                        # Mapping: 0=Top, 1=Right, 2=Bottom, 3=Left
                        side_labels = ["Top", "Right", "Bottom", "Left"]
                        side_str = side_labels[req_cand_side]
                        print(f"🔄 BETTER MATCH FOUND: "
                              f"Candidate {cand.id} (Edge {side_str}) "
                              f"beats previous best (Cost {best_cost:.2f} -> {cost:.2f}). "
                              f"Target: Host Edge at ({edge.x}, {edge.y}). "
                              )
                    best_cost = cost
                    best_move = (uid, target_edges_chain, (nx, ny, cand.w, cand.h))

        return best_move, best_cost

    # ================= Execution =================

    def solve(self):
        # 1. Seed
        best_seed = max(self.pieces, key=lambda p: sum(e.variance for e in p.edge_features))
        print(f"🌱 Seed: Piece {best_seed.id}")
        best_seed.placed_x = 0
        best_seed.placed_y = 0
        self.placed_pieces.append(best_seed)
        self.unused_ids.remove(best_seed.id)
        self.add_piece_to_frontier(best_seed, 0, 0)
        self.min_x, self.max_x = 0, best_seed.w
        self.min_y, self.max_y = 0, best_seed.h

        while self.unused_ids:
            if VISUAL_DEBUG:
                self.draw_debug_view()
                cv2.waitKey(STEP_SPEED_MS)

            # A. Zippering
            zipped = self.try_zipper_frontier()
            if zipped:
                for z_edge in zipped:
                    if z_edge in self.frontier: self.frontier.remove(z_edge)
                continue  # Loop again to see if more zipping needed

            # B. Global Search
            move, cost = self.find_global_best_match()

            if move and cost < MATCH_THRESHOLD:
                uid, covered_edges, rect = move
                nx, ny, w, h = rect

                # Collision Check
                if check_overlap(rect, self.placed_pieces):
                    print(f"⚠️ Overlap detected for Piece {uid}. Skipping.")
                    # Hack: Removing the primary edge to prevent infinite loop on bad match
                    if covered_edges[0] in self.frontier: self.frontier.remove(covered_edges[0])
                    continue

                # Execute
                cand = self.pieces[uid]
                cand.placed_x = nx
                cand.placed_y = ny
                self.placed_pieces.append(cand)
                self.unused_ids.remove(uid)
                print(f"✅ Placed Piece {uid} (Cost: {cost:.1f})")

                # Remove Covered Edges
                total_covered_len = 0
                for ce in covered_edges:
                    if ce in self.frontier: self.frontier.remove(ce)
                    total_covered_len += ce.length

                # Handle Split (Remainder)
                # If chain was 100px but piece was 80px, we have 20px remainder.
                used_len = cand.h if covered_edges[0].axis == 'V' else cand.w
                remainder = total_covered_len - used_len

                if remainder > MIN_EDGE_LEN:
                    # Where does remainder start?
                    # End of the new piece
                    last_edge = covered_edges[-1]
                    # Logic is tricky for complex chains.
                    # Simplified: We just add the NEW piece's edges.
                    # The "Hole" is filled. If there's a gap, it's a new open edge?
                    # No, we removed the frontier edges that represented the empty space.
                    # If we didn't fill it all, we implicitly left a gap?
                    # FIX: We must re-add the remainder of the LAST edge in chain.
                    pass  # Omitted for brevity in V2, usually exact match in puzzle.

                # Add New Edges
                self.add_piece_to_frontier(cand, nx, ny)

                # Update Bounds
                self.min_x = min(self.min_x, nx)
                self.min_y = min(self.min_y, ny)
                self.max_x = max(self.max_x, nx + cand.w)
                self.max_y = max(self.max_y, ny + cand.h)

            else:
                print("❌ No confident match found. Stopping or Manual Intervention needed.")
                break

        print("🏁 Done!")
        if VISUAL_DEBUG:
            self.draw_debug_view(final=True)
        cv2.waitKey(0)

    def draw_debug_view(self, final=False):
        w = self.max_x - self.min_x + 600
        h = self.max_y - self.min_y + 600
        canvas = np.zeros((h, w, 3), dtype=np.uint8)
        ox = -self.min_x + 300
        oy = -self.min_y + 300

        for p in self.placed_pieces:
            dx, dy = p.placed_x + ox, p.placed_y + oy
            canvas[dy:dy + p.h, dx:dx + p.w] = p.image
            cv2.rectangle(canvas, (dx, dy), (dx + p.w, dy + p.h), (50, 50, 50), 1)
            cx, cy = dx + p.w // 2 - 10, dy + p.h // 2 + 5
            cv2.putText(canvas, str(p.id), (cx, cy), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 0, 0), 3)
            cv2.putText(canvas, str(p.id), (cx, cy), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 255, 255), 2)

        if not final:
            for edge in self.frontier:
                sx, sy = edge.x + ox, edge.y + oy
                ex, ey = (sx + edge.length, sy) if edge.axis == 'H' else (sx, sy + edge.length)
                cv2.line(canvas, (sx, sy), (ex, ey), (0, 0, 255), 2)

        disp = cv2.resize(canvas, (int(w * DEBUG_SCALE), int(h * DEBUG_SCALE)))
        cv2.imshow("Puzzle Solver V2", disp)

    def draw_initial_state(self):
        """
        Draws all extracted pieces at their ORIGINAL positions and ROTATION
        by cropping the raw bounding box directly from the original input image.
        """
        if not self.pieces:
            print("⚠️ No pieces to draw. Run preprocess() first.")
            return

        print(f"🎨 Drawing {len(self.pieces)} pieces in their original rotated state...")

        # 1. Determine Canvas Size
        # Find global min/max coordinates across ALL corners of ALL pieces
        min_gx, min_gy = float('inf'), float('inf')
        max_gx, max_gy = float('-inf'), float('-inf')

        for p in self.pieces:
            xs = [p.original_tl[0], p.original_tr[0], p.original_bl[0], p.original_br[0]]
            ys = [p.original_tl[1], p.original_tr[1], p.original_bl[1], p.original_br[1]]

            min_gx = min(min_gx, min(xs))
            min_gy = min(min_gy, min(ys))
            max_gx = max(max_gx, max(xs))
            max_gy = max(max_gy, max(ys))

        # 2. Create Canvas
        pad = 50
        canvas_w = int(max_gx - min_gx + (pad * 2))
        canvas_h = int(max_gy - min_gy + (pad * 2))

        # Use a generic gray background to differentiate from piece background
        canvas = np.zeros((canvas_h, canvas_w, 3), dtype=np.uint8)

        # 3. Draw Pieces
        for p in self.pieces:
            # A. Calculate Local Bounding Box of the rotated piece
            xs = [p.original_tl[0], p.original_tr[0], p.original_bl[0], p.original_br[0]]
            ys = [p.original_tl[1], p.original_tr[1], p.original_bl[1], p.original_br[1]]

            bx_min, by_min = int(min(xs)), int(min(ys))
            bx_max, by_max = int(max(xs)), int(max(ys))

            # B. Crop from Original Image (The Raw Rotated Pixel Data)
            # Ensure we don't go out of bounds of the original image
            h_orig, w_orig = self.original_img.shape[:2]
            src_x1 = max(0, bx_min)
            src_y1 = max(0, by_min)
            src_x2 = min(w_orig, bx_max)
            src_y2 = min(h_orig, by_max)

            raw_piece_img = self.original_img[src_y1:src_y2, src_x1:src_x2]

            if raw_piece_img.size == 0: continue

            # C. Calculate Paste Location on New Canvas
            draw_x = int(src_x1 - min_gx + pad)
            draw_y = int(src_y1 - min_gy + pad)

            h_crop, w_crop = raw_piece_img.shape[:2]

            # D. Paste the Raw Crop
            if draw_y + h_crop <= canvas_h and draw_x + w_crop <= canvas_w:
                canvas[draw_y:draw_y + h_crop, draw_x:draw_x + w_crop] = raw_piece_img

            # E. Draw Red Boundary Box (Visual Confirmation)
            pts = np.array([
                p.original_tl, p.original_tr, p.original_br, p.original_bl
            ], dtype=np.int32)

            # Shift points to canvas space
            pts[:, 0] -= int(min_gx - pad)
            pts[:, 1] -= int(min_gy - pad)

            # cv2.polylines(canvas, [pts], isClosed=True, color=(0, 0, 255), thickness=2)

            # F. Draw ID
            cx = draw_x + w_crop // 2
            cy = draw_y + h_crop // 2
            cv2.putText(canvas, str(p.id), (cx - 10, cy),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 2)

        # 4. Display
        display_scale = 0.5
        disp_w = int(canvas_w * display_scale)
        disp_h = int(canvas_h * display_scale)
        disp = cv2.resize(canvas, (disp_w, disp_h))

        cv2.imshow("Initial State (Original Rotated)", disp)
        cv2.waitKey(0)

    def animate_assembly(self):
        """
        Animates pieces flying AND rotating from their original state to the solution.
        """
        if not self.placed_pieces:
            print("⚠️ No solved pieces to animate.")
            return

        print("🎥 Starting Assembly Animation (Translation + Rotation)...")

        # 1. Calculate Bounds of the "Input World" (Original positions)
        min_in_x, min_in_y = float('inf'), float('inf')
        max_in_x, max_in_y = float('-inf'), float('-inf')

        for p in self.placed_pieces:
            # Check all 4 corners to be safe, or just TL/BR
            xs = [p.original_tl[0], p.original_tr[0], p.original_bl[0], p.original_br[0]]
            ys = [p.original_tl[1], p.original_tr[1], p.original_bl[1], p.original_br[1]]
            min_in_x = min(min_in_x, min(xs))
            min_in_y = min(min_in_y, min(ys))
            max_in_x = max(max_in_x, max(xs))
            max_in_y = max(max_in_y, max(ys))

        # 2. Calculate Union Bounds (Input World + Solution World)
        # self.min_x/max_x are the bounds of the solved puzzle relative to seed (0,0)
        global_min_x = min(self.min_x, min_in_x)
        global_min_y = min(self.min_y, min_in_y)
        global_max_x = max(self.max_x, max_in_x)
        global_max_y = max(self.max_y, max_in_y)

        # 3. Define Canvas Size & Offset
        padding = 100
        sol_w = int(global_max_x - global_min_x + (padding * 2))
        sol_h = int(global_max_y - global_min_y + (padding * 2))

        # The offset shifts the entire coordinate system so 'global_min' is at (padding, padding)
        ox = -global_min_x + padding
        oy = -global_min_y + padding
        # ===========================================================

        total_frames = 120  # 4 seconds

        # Pre-calculate Start/End Centers for all pieces to save CPU in loop
        piece_paths = []
        for p in self.placed_pieces:
            # START: Center of the original rotated box
            # Average of the 4 original corners
            start_cx = (p.original_tl[0] + p.original_br[0]) / 2
            start_cy = (p.original_tl[1] + p.original_br[1]) / 2

            # END: Center of the placed (upright) piece
            # p.placed_x is Top-Left, so add half width/height
            end_cx = (p.placed_x + ox) + (p.w / 2)
            end_cy = (p.placed_y + oy) + (p.h / 2)

            piece_paths.append({
                'p': p,
                'start': (start_cx, start_cy),
                'end': (end_cx, end_cy),
                'angle': p.initial_angle
            })

        # 2. Animation Loop
        for frame in range(total_frames + 1):
            t = frame / total_frames
            # Ease-out (fast start, slow stop)
            t_smooth = 1 - (1 - t) ** 3

            canvas = np.zeros((sol_h, sol_w, 3), dtype=np.uint8)

            for item in piece_paths:
                p = item['p']
                sx, sy = item['start']
                ex, ey = item['end']

                # A. Interpolate Center
                cur_cx = sx + (ex - sx) * t_smooth
                cur_cy = sy + (ey - sy) * t_smooth

                # B. Interpolate Angle
                # We go from p.initial_angle -> 0 degrees
                cur_angle = item['angle'] * (1 - t_smooth)

                # C. Rotate the Piece Image
                # We take the upright p.image and rotate it "back" towards initial,
                # then slowly rotate it to 0.
                img_to_draw = self.rotate_image(p.image, cur_angle)

                # D. Calculate Top-Left based on New Center & New Size
                h_img, w_img = img_to_draw.shape[:2]
                draw_x = int(cur_cx - (w_img / 2))
                draw_y = int(cur_cy - (h_img / 2))

                # E. Draw (with Clipping)
                if draw_x < sol_w and draw_y < sol_h and draw_x + w_img > 0 and draw_y + h_img > 0:
                    x1 = max(0, draw_x)
                    y1 = max(0, draw_y)
                    x2 = min(sol_w, draw_x + w_img)
                    y2 = min(sol_h, draw_y + h_img)

                    src_x1 = x1 - draw_x
                    src_y1 = y1 - draw_y
                    src_x2 = src_x1 + (x2 - x1)
                    src_y2 = src_y1 + (y2 - y1)

                    # Basic transparency check (optional): masking black pixels
                    # For simple blit:
                    try:
                        patch = img_to_draw[src_y1:src_y2, src_x1:src_x2]
                        target = canvas[y1:y2, x1:x2]

                        # Simple mask to allow non-rectangular shapes if source has black bg
                        gray = cv2.cvtColor(patch, cv2.COLOR_BGR2GRAY)
                        _, mask = cv2.threshold(gray, 1, 255, cv2.THRESH_BINARY)

                        # Background where mask is 0
                        bg = cv2.bitwise_and(target, target, mask=cv2.bitwise_not(mask))
                        # Foreground where mask is 1
                        fg = cv2.bitwise_and(patch, patch, mask=mask)

                        canvas[y1:y2, x1:x2] = cv2.add(bg, fg)
                    except:
                        pass  # avoid crash on sub-pixel rounding errors

            # F. Show Frame
            display_scale = 0.8
            disp_w = int(sol_w * display_scale)
            disp_h = int(sol_h * display_scale)
            disp = cv2.resize(canvas, (disp_w, disp_h))

            cv2.imshow("Puzzle Assembly", disp)
            if frame == 0:
                print("⏸️  Paused at start. Press any key to begin animation...")
                cv2.waitKey(0)  # Indefinite wait
            else:
                cv2.waitKey(16)  # ~60 FPS wait for animation

        print("🎉 Animation Complete!")
        cv2.waitKey(0)


if __name__ == "__main__":
    import sys

    solver = PuzzleSolverV2(sys.argv[1] if len(sys.argv) > 1 else ".\\day2_test\\test_irregular_translate.png")
    solver.preprocess()
    # solver.draw_initial_state()
    solver.solve()
    solver.animate_assembly()