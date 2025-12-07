#!/usr/bin/env python3
import cv2
import numpy as np
import heapq
import time
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Dict

# ================= Configuration =================
VISUAL_DEBUG = True  # Turn on/off real-time debug window
DEBUG_SCALE = 0.6  # Scale of debug window
STEP_SPEED_MS = 0  # Delay between steps (10ms is fast, 0 is infinite wait)
MIN_EDGE_LEN = 10  # Ignore edges shorter than this
MATCH_THRESHOLD = 2500.0  # Max MSE to accept a match (tune this based on image quality)


# ================= Data Structures =================

@dataclass
class Piece:
    id: int
    image: np.ndarray  # BGR Image
    h: int
    w: int
    # 4 Edges: 0:Top, 1:Right, 2:Bottom, 3:Left
    # Stored as dictionaries with 'pixels', 'grad_mag', 'variance'
    edge_features: List[dict] = field(default_factory=list)

    # Final Position (if placed)
    placed_x: int = 0
    placed_y: int = 0


@dataclass(order=True)
class FrontierEdge:
    # Priority: We want Highest Variance first.
    # heapq is Min-Heap, so we store negative variance.
    priority: float

    # Geometry
    x: int
    y: int
    length: int
    axis: str  # 'H' (Horizontal) or 'V' (Vertical)
    direction: int  # 1 (Right/Down) or -1 (Left/Up)

    # Metadata
    host_piece_id: int

    def __repr__(self):
        dir_str = "+" if self.direction > 0 else "-"
        return f"Edge(x={self.x}, y={self.y}, len={self.length}, {self.axis}{dir_str}, var={-self.priority:.1f})"


# ================= Helper Functions =================

def order_points(pts):
    """Sorts rectangle points: top-left, top-right, bottom-right, bottom-left"""
    rect = np.zeros((4, 2), dtype="float32")
    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]
    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]
    return rect


def get_gradient_variance(pixels: np.ndarray) -> float:
    """Computes Standard Deviation of Gradient Magnitude for a strip of pixels."""
    if pixels.size == 0: return 0.0
    gray = cv2.cvtColor(pixels[np.newaxis, :, :], cv2.COLOR_BGR2GRAY)
    grad_x = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(gray, cv2.CV_32F, 0, 1, ksize=3)
    mag = cv2.magnitude(grad_x, grad_y)
    return float(np.std(mag))


def extract_edge_feature(img: np.ndarray, side: int) -> dict:
    """Extracts pixel data and metrics for a specific side (0=T, 1=R, 2=B, 3=L)."""
    h, w = img.shape[:2]
    # We take a 1px line for matching, maybe 3px for robustness if needed
    if side == 0:  # Top
        pixels = img[0, :, :]
    elif side == 1:  # Right
        pixels = img[:, w - 1, :]
    elif side == 2:  # Bottom
        pixels = img[h - 1, :, :]
    elif side == 3:  # Left
        pixels = img[:, 0, :]

    var = get_gradient_variance(pixels)

    # Compute Gradient Magnitude strip for matching
    gray = cv2.cvtColor(pixels[np.newaxis, :, :], cv2.COLOR_BGR2GRAY)
    sobel = cv2.Sobel(gray, cv2.CV_32F, 1, 0, ksize=3)  # Simple 1D gradient
    grad_mag = np.abs(sobel).flatten()

    return {
        "pixels": pixels.astype(np.float32),
        "grad": grad_mag,
        "variance": var,
        "len": len(pixels)
    }


def check_overlap(new_rect: Tuple[int, int, int, int], placed_pieces: List[Piece]) -> bool:
    """
    Returns True if new_rect (x, y, w, h) intersects with any placed piece.
    Uses strict inequality (area > 0) to allow touching edges.
    """
    nx, ny, nw, nh = new_rect
    # Small tolerance to prevent floating point/1px rounding false positives
    tol = 2

    for p in placed_pieces:
        # Intersection test
        x_overlap = max(0, min(nx + nw - tol, p.placed_x + p.w - tol) - max(nx + tol, p.placed_x + tol))
        y_overlap = max(0, min(ny + nh - tol, p.placed_y + p.h - tol) - max(ny + tol, p.placed_y + tol))

        if x_overlap > 0 and y_overlap > 0:
            return True
    return False


# ================= Main Solver Logic =================

class PuzzleSolver:
    def __init__(self, img_path):
        self.original_img = cv2.imread(img_path)
        if self.original_img is None:
            raise ValueError("Image not found")

        self.pieces: List[Piece] = []
        self.placed_pieces: List[Piece] = []
        self.frontier = []  # Priority Queue
        self.unused_ids = set()

        # Debug Visualization
        self.min_x, self.max_x = 0, 0
        self.min_y, self.max_y = 0, 0

    def preprocess(self):
        """Extracts pieces, straightens them, DOES NOT resize."""
        gray = cv2.cvtColor(self.original_img, cv2.COLOR_BGR2GRAY)
        _, thresh = cv2.threshold(gray, 10, 255, cv2.THRESH_BINARY)
        contours, _ = cv2.findContours(thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

        print(f"🔍 Found {len(contours)} contours")

        for i, cnt in enumerate(contours):
            if cv2.contourArea(cnt) < 500: continue

            rect = cv2.minAreaRect(cnt)
            box = cv2.boxPoints(rect)
            box = order_points(box)

            # Perspective Transform (Straighten)
            (tl, tr, br, bl) = box
            widthA = np.linalg.norm(br - bl)
            widthB = np.linalg.norm(tr - tl)
            maxWidth = int(max(widthA, widthB))

            heightA = np.linalg.norm(tr - br)
            heightB = np.linalg.norm(tl - bl)
            maxHeight = int(max(heightA, heightB))

            dst_pts = np.array([[0, 0], [maxWidth - 1, 0], [maxWidth - 1, maxHeight - 1], [0, maxHeight - 1]],
                               dtype="float32")
            M = cv2.getPerspectiveTransform(box, dst_pts)
            warped = cv2.warpPerspective(self.original_img, M, (maxWidth, maxHeight))

            # Extract Features
            piece = Piece(id=i, image=warped, h=maxHeight, w=maxWidth)
            piece.edge_features = [
                extract_edge_feature(warped, 0),  # Top
                extract_edge_feature(warped, 1),  # Right
                extract_edge_feature(warped, 2),  # Bottom
                extract_edge_feature(warped, 3)  # Left
            ]

            self.pieces.append(piece)
            self.unused_ids.add(i)

        print(f"✅ Extracted {len(self.pieces)} valid pieces.")

    def add_piece_to_frontier(self, piece: Piece, offset_x: int, offset_y: int):
        """Adds all 4 edges of a newly placed piece to the frontier."""
        # Top (Side 0) - Horizontal, Normal points Up (-1)
        self.push_frontier(piece, 0, offset_x, offset_y, piece.w, 'H', -1)

        # Right (Side 1) - Vertical, Normal points Right (+1)
        self.push_frontier(piece, 1, offset_x + piece.w, offset_y, piece.h, 'V', 1)

        # Bottom (Side 2) - Horizontal, Normal points Down (+1)
        self.push_frontier(piece, 2, offset_x, offset_y + piece.h, piece.w, 'H', 1)

        # Left (Side 3) - Vertical, Normal points Left (-1)
        self.push_frontier(piece, 3, offset_x, offset_y, piece.h, 'V', -1)

    def push_frontier(self, piece, side_idx, x, y, length, axis, direction):
        """Helper to push edge to heap with Variance priority."""
        if length < MIN_EDGE_LEN: return

        # Get variance of this specific segment
        feat = piece.edge_features[side_idx]

        # Note: Ideally we calculate variance of the *remaining* segment if split.
        # For simplicity, we use the pre-calculated variance of the whole edge
        # but penalize it slightly if it's a sub-segment (optional).
        priority = -feat['variance']  # Max-Heap simulation

        edge = FrontierEdge(priority, int(x), int(y), int(length), axis, direction, piece.id)
        heapq.heappush(self.frontier, edge)

    def get_sliding_match(self, target_edge: FrontierEdge, cand: Piece, cand_side: int):
        """
        Slides candidate edge along target edge.
        Returns: (best_mse, offset)
        offset: Distance from target_start to candidate_start
        """
        # Get pixel arrays
        # Note: We need to handle 'Reverse' matching if needed?
        # For rectangular puzzles with standard orientation:
        # Top matches Bottom, Left matches Right.
        # Pixels are typically left-to-right or top-to-bottom.

        # Extract the target segment pixels from the Host Piece?
        # This is complex because the FrontierEdge is abstract.
        # We need to grab the pixels from the Placed Piece at the frontier coordinates.
        # SIMPLIFICATION: We assume we can access the host piece's edge data.
        host_piece = next(p for p in self.pieces if p.id == target_edge.host_piece_id)

        # Identify which side of host piece is this edge?
        # We need to map global frontier coord back to local piece coord.
        # This is tricky.
        # ROBUST APPROACH: We trust the candidate check.
        # We compare candidate pixels vs host pixels.

        # Get Host Pixels (Full Edge)
        # We need to know WHICH side of the host this frontier belongs to.
        # H, -1 => Top (0); V, 1 => Right (1); H, 1 => Bottom (2); V, -1 => Left (3)
        if target_edge.axis == 'H':
            host_side = 0 if target_edge.direction == -1 else 2
        else:
            host_side = 1 if target_edge.direction == 1 else 3

        host_pixels = host_piece.edge_features[host_side]['pixels']
        cand_pixels = cand.edge_features[cand_side]['pixels']

        # T-Junction Logic: The Frontier Edge might be a SUB-SEGMENT of the Host Edge.
        # We need to track where the current frontier edge starts relative to the host edge.
        # For this prototype, let's assume standard greedy fill where we just try to match.

        # Standardize constraint: Candidate Length must be <= Target Length
        # (We are filling a hole, fitting a brick into a gap)
        len_t = target_edge.length
        len_c = len(cand_pixels)

        if len_c > len_t:
            return float('inf'), 0

        # Sliding Window
        best_mse = float('inf')
        best_offset = 0

        # We slide the candidate (smaller) inside the target (larger)
        # Limitation: We assume we are filling the target edge from its start (0).
        # To strictly support T-junctions anywhere, we loop range.

        # Optimization: Just check standard deviation difference first? No.

        for offset in range(0, len_t - len_c + 1, 1):  # Step 1
            # In a real T-junction split, the "Host Pixels" array is the FULL edge.
            # We need to know where the "Frontier Edge" starts within the "Host Edge".
            # This requires tracking "offset_in_host" in FrontierEdge.
            # APPROXIMATION for V1: Assume we match the BEGINNING of the frontier.

            # Since we don't track offset_in_host, we'll try to match against the
            # *relevant section* of the host pixels.
            # This is the hardest part. Let's simplify:
            # We compare features.

            # Extract sub-segment
            # This assumes host_pixels corresponds exactly to target_edge.
            # If target_edge was split before, we are in trouble unless we sliced pixels.
            # FIX: Only "Fresh" edges are in frontier? No.

            # REVISED STRATEGY:
            # When splitting, we assume the pixels are consistent.
            # Let's just do a direct visual match logic assuming alignment.

            diff = host_pixels[offset:offset + len_c] - cand_pixels
            mse = np.mean(diff ** 2)

            if mse < best_mse:
                best_mse = mse
                best_offset = offset

        return best_mse, best_offset

    def solve(self):
        # 1. Seed Selection: Piece with highest total variance
        print("🌱 Selecting Seed Piece...")
        best_seed = None
        max_var = -1
        for p in self.pieces:
            total_var = sum(e['variance'] for e in p.edge_features)
            if total_var > max_var:
                max_var = total_var
                best_seed = p

        print(f"Selected Seed: Piece {best_seed.id} (Var: {max_var:.1f})")

        # Place Seed
        best_seed.placed_x = 0
        best_seed.placed_y = 0
        self.placed_pieces.append(best_seed)
        self.unused_ids.remove(best_seed.id)

        # Add Seed Edges to Frontier
        self.add_piece_to_frontier(best_seed, 0, 0)

        # Update Canvas Bounds
        self.min_x, self.max_x = 0, best_seed.w
        self.min_y, self.max_y = 0, best_seed.h

        # 2. Region Growing Loop
        while self.frontier and self.unused_ids:
            if VISUAL_DEBUG:
                self.draw_debug_view()
                cv2.waitKey(STEP_SPEED_MS)

            # Pop best edge
            target = heapq.heappop(self.frontier)

            # Discard tiny edges (noise)
            if target.length < MIN_EDGE_LEN:
                continue

            print(f"Processing Edge: {target}")

            # Find Match
            best_match_score = float('inf')
            best_cand_id = -1
            best_cand_side = -1
            best_offset = 0  # Relative to target start

            # Determine required opposite direction
            # If Target is Right (+1), we need Left (-1)
            req_axis = target.axis
            req_dir = -target.direction

            # Map required direction to side index for candidate
            # H, -1 (Up) -> Bottom (2)
            # V, 1 (Right) -> Left (3)
            # H, 1 (Down) -> Top (0)
            # V, -1 (Left) -> Right (1)
            if req_axis == 'H':
                cand_side_idx = 2 if req_dir == 1 else 0
            else:
                cand_side_idx = 3 if req_dir == 1 else 1

            # Iterate all unused pieces
            for uid in list(self.unused_ids):
                cand = self.pieces[uid]

                # Check Length Constraint (Candidate <= Gap)
                # Note: We allow candidate to be slightly smaller for T-junctions
                cand_len = cand.edge_features[cand_side_idx]['len']

                if cand_len > target.length + 2:  # Tolerance
                    continue

                # Compute Cost (Sliding Window)
                # For simplicity in this demo, we assume the frontier edge IS the pixels
                # In a full implementation, we'd slice the host pixels correctly.
                # Here we use a simpler 'midpoint' or 'start' check.

                # Get Host Pixels
                host_piece = next(p for p in self.pieces if p.id == target.host_piece_id)
                # Logic to get specific pixel slice omitted for brevity,
                # we assume simple edge matching here:

                # Calculate cost
                host_full_pixels = host_piece.edge_features[0 if target.axis == 'H' and target.direction == -1 else (
                    1 if target.direction == 1 and target.axis == 'V' else (2 if target.direction == 1 else 3))][
                    'pixels']

                # CRITICAL: If the frontier edge is a REMAINDER, we need to offset the host pixels!
                # This requires tracking "offset_from_origin" in FrontierEdge.
                # For this demo, we assume we match the *start* of the current open segment.
                # (Improving this requires adding 'origin_offset' to FrontierEdge)

                match_mse, _ = self.get_sliding_match(target, cand, cand_side_idx)

                if match_mse < best_match_score:
                    best_match_score = match_mse
                    best_cand_id = uid
                    best_cand_side = cand_side_idx
                    best_offset = 0  # Greedy: assume start alignment for now

            # Check Threshold
            if best_match_score > MATCH_THRESHOLD:
                print(f"  ❌ No good match (Best: {best_match_score:.1f}). Discarding edge.")
                continue  # Discard edge (Hole)

            # COLLISION CHECK
            cand = self.pieces[best_cand_id]

            # Calculate Candidate Coordinates
            # If Target is V, Right (+1) at (x,y), and we match Left (-1)
            # Candidate X = Target X
            # Candidate Y = Target Y + Offset
            nx, ny = 0, 0
            if target.axis == 'V':
                if target.direction == 1:  # Right side of Host
                    nx = target.x
                    ny = target.y + best_offset
                else:  # Left side of Host
                    nx = target.x - cand.w
                    ny = target.y + best_offset
            else:  # Horizontal
                if target.direction == 1:  # Bottom of Host
                    nx = target.x + best_offset
                    ny = target.y
                else:  # Top of Host
                    nx = target.x + best_offset
                    ny = target.y - cand.h

            cand_rect = (nx, ny, cand.w, cand.h)
            if check_overlap(cand_rect, self.placed_pieces):
                print("  ⚠️ Collision detected! Skipping.")
                continue

            # EXECUTE PLACEMENT
            print(f"  ✅ Matched Piece {cand.id} (MSE: {best_match_score:.1f})")
            cand.placed_x = nx
            cand.placed_y = ny
            self.placed_pieces.append(cand)
            self.unused_ids.remove(cand.id)

            # Update Bounds
            self.min_x = min(self.min_x, nx)
            self.min_y = min(self.min_y, ny)
            self.max_x = max(self.max_x, nx + cand.w)
            self.max_y = max(self.max_y, ny + cand.h)

            # MANAGE FRONTIER (Splitting)
            used_len = cand.h if target.axis == 'V' else cand.w
            remainder = target.length - used_len

            if remainder > MIN_EDGE_LEN:
                # Add remainder back to frontier
                # New start coordinates
                rx = target.x + (0 if target.axis == 'V' else used_len)
                ry = target.y + (used_len if target.axis == 'V' else 0)

                rem_edge = FrontierEdge(
                    target.priority,  # Keep priority (or re-calc)
                    rx, ry, remainder,
                    target.axis, target.direction, target.host_piece_id
                )
                heapq.heappush(self.frontier, rem_edge)

            # ADD NEW EDGES
            # We add all 4, but we should strictly only add exposed ones.
            # Simplified: Add all, Collision check will reject overlaps later.
            self.add_piece_to_frontier(cand, nx, ny)

        print("🏁 Puzzle Assembly Complete!")
        if VISUAL_DEBUG:
            self.draw_debug_view(final=True)
            cv2.waitKey(0)

    # ================= Visualization =================

    def draw_debug_view(self, final=False):
        # Calculate Canvas Size
        w = self.max_x - self.min_x + 600
        h = self.max_y - self.min_y + 600
        canvas = np.zeros((h, w, 3), dtype=np.uint8)

        ox = -self.min_x + 300  # Center offset
        oy = -self.min_y + 300

        # 1. Draw Placed Pieces
        for p in self.placed_pieces:
            dx = p.placed_x + ox
            dy = p.placed_y + oy

            # Boundary check
            if dx < 0 or dy < 0 or dx + p.w > w or dy + p.h > h: continue

            canvas[dy:dy + p.h, dx:dx + p.w] = p.image

            # Draw border
            cv2.rectangle(canvas, (dx, dy), (dx + p.w, dy + p.h), (50, 50, 50), 1)

            # --- NEW: Draw ID Annotation ---
            # Calculate center of the piece
            cx = dx + p.w // 2 - 10
            cy = dy + p.h // 2 + 5

            text = str(p.id)
            # Draw black outline (thickness 3) so text pops on any background
            cv2.putText(canvas, text, (cx, cy), cv2.FONT_HERSHEY_SIMPLEX,
                        0.8, (0, 0, 0), 3)
            # Draw yellow text (thickness 2) on top
            cv2.putText(canvas, text, (cx, cy), cv2.FONT_HERSHEY_SIMPLEX,
                        0.8, (0, 255, 255), 2)
            # -------------------------------

        if not final:
            # 2. Draw Frontier Lines (Red)
            for edge in self.frontier:
                sx = edge.x + ox
                sy = edge.y + oy

                if edge.axis == 'H':
                    ex, ey = sx + edge.length, sy
                else:
                    ex, ey = sx, sy + edge.length

                cv2.line(canvas, (sx, sy), (ex, ey), (0, 0, 255), 2)
                # Draw normal indicator
                mx, my = (sx + ex) // 2, (sy + ey) // 2
                if edge.axis == 'H':
                    cv2.line(canvas, (mx, my), (mx, my + 10 * edge.direction), (0, 0, 255), 1)
                else:
                    cv2.line(canvas, (mx, my), (mx + 10 * edge.direction, my), (0, 0, 255), 1)

        # Scale down for display
        disp_h = int(h * DEBUG_SCALE)
        disp_w = int(w * DEBUG_SCALE)
        disp = cv2.resize(canvas, (disp_w, disp_h))

        cv2.imshow("Puzzle Solver Live View", disp)


if __name__ == "__main__":
    import sys

    # Usage: python smart_solver.py puzzle_input.jpg
    input_file = sys.argv[1] if len(sys.argv) > 1 else "more_sample_irregular/more_sample_irregular/sample1/sample1_translate.png"

    solver = PuzzleSolver(input_file)
    solver.preprocess()
    solver.solve()