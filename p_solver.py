#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
拼图求解器 v4.0（translate 版）

流程：
1. 自动检测拼图块（轮廓）
2. 用透视变换把每块图块“回正”成矩形小图
3. 自动推断 rows / cols（接近正方形的布局，不需要用户输入）
4. 用 Beam Search 只做“平移匹配”（不再尝试旋转），求出每块应该放在哪个格子
5. 输出：
   - translate_solved.png：规则网格拼好的结果
   - translate_animation.gif：空白背景上，从原始散乱位置 → 规则网格 的动画
"""

import argparse
import math
from dataclasses import dataclass
from typing import List, Tuple
from collections import Counter

import cv2
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter


# ===================== 工具函数：排序四个角点 =====================

def order_points(pts: np.ndarray) -> np.ndarray:
    """
    将四个点按：左上、右上、右下、左下 的顺序排序
    供透视变换使用
    """
    rect = np.zeros((4, 2), dtype="float32")

    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]  # top-left
    rect[2] = pts[np.argmax(s)]  # bottom-right

    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]  # top-right
    rect[3] = pts[np.argmax(diff)]  # bottom-left

    return rect


# ===================== 简单边缘特征与代价 =====================

@dataclass
class EdgeFeature:
    # 一条边的 RGB 线，形状 (L, 3)，float32
    line: np.ndarray


@dataclass
class Piece:
    id: int
    image: np.ndarray                 # 标准化后的砖图（已回正）
    edges: List[EdgeFeature]          # 4 条边：0 上，1 右，2 下，3 左
    init_center: Tuple[float, float]  # 在原图中的中心（用于动画）
    init_angle: float                 # 在原图中的角度（minAreaRect 角）


def extract_edge_simple(img: np.ndarray, side: int, band: int = 5) -> EdgeFeature:
    """
    从图像某一侧取一个竖/横向带状区域，求均值作为边缘特征。
    side: 0=top, 1=right, 2=bottom, 3=left
    band: 带宽（像素）
    """
    h, w, _ = img.shape
    band = max(1, min(band, min(h, w) // 4))

    if side == 0:        # top
        strip = img[0:band, :, :]
        line = strip.mean(axis=0)      # (w, 3)
    elif side == 2:      # bottom
        strip = img[h-band:h, :, :]
        line = strip.mean(axis=0)
    elif side == 3:      # left
        strip = img[:, 0:band, :]
        line = strip.mean(axis=1)      # (h, 3)
    else:                # right
        strip = img[:, w-band:w, :]
        line = strip.mean(axis=1)

    return EdgeFeature(line=line.astype(np.float32))


def edge_cost(a: EdgeFeature, b: EdgeFeature) -> float:
    """
    两条边之间的 MSE 代价，越小越匹配。
    """
    la = a.line.shape[0]
    lb = b.line.shape[0]
    L = min(la, lb)
    if L <= 0:
        return 1e9
    diff = a.line[:L] - b.line[:L]
    mse = (diff * diff).mean()
    return float(mse)


# ===================== 回正 + 提取所有拼图块 =====================

def straighten_and_extract_pieces(board_bgr: np.ndarray,
                                  min_area: int = 500) -> List[Piece]:
    """
    使用透视矫正思路：
    1. 灰度阈值 -> 找每个砖块的轮廓
    2. minAreaRect + order_points 得到四角
    3. 透视变换得到“回正”的砖块 warped
    4. 统一 resize，并为每块提取 4 条边特征

    返回：Piece 列表（不再做旋转搜索，只做平移匹配）
    """
    gray = cv2.cvtColor(board_bgr, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(gray, 10, 255, cv2.THRESH_BINARY)

    contours, _ = cv2.findContours(
        thresh, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
    )

    temp_tiles = []

    print(f"\n🔍 检测轮廓（候选砖块）: {len(contours)} 个")

    for cnt in contours:
        area = cv2.contourArea(cnt)
        if area < min_area:
            continue

        rect = cv2.minAreaRect(cnt)   # ((cx, cy), (w, h), angle)
        (cx, cy) = rect[0]
        angle = rect[2]

        box = cv2.boxPoints(rect)
        box = np.array(box, dtype="float32")

        rect_pts = order_points(box)
        (tl, tr, br, bl) = rect_pts

        # 计算透视后目标宽高
        widthA = np.linalg.norm(br - bl)
        widthB = np.linalg.norm(tr - tl)
        maxWidth = int(max(widthA, widthB))

        heightA = np.linalg.norm(tr - br)
        heightB = np.linalg.norm(tl - bl)
        maxHeight = int(max(heightA, heightB))

        if maxWidth < 10 or maxHeight < 10:
            continue

        dst_pts = np.array([
            [0, 0],
            [maxWidth - 1, 0],
            [maxWidth - 1, maxHeight - 1],
            [0, maxHeight - 1]
        ], dtype="float32")

        M = cv2.getPerspectiveTransform(rect_pts, dst_pts)
        warped = cv2.warpPerspective(board_bgr, M, (maxWidth, maxHeight))

        temp_tiles.append({
            "warped": warped,
            "center": (float(cx), float(cy)),
            "angle": float(angle)
        })

    if not temp_tiles:
        raise RuntimeError("没有检测到足够大的拼图块，请检查输入图或阈值。")

    # 统计最常见的宽高，统一尺寸
    sizes = [(t["warped"].shape[0], t["warped"].shape[1]) for t in temp_tiles]
    height_counter = Counter([s[0] for s in sizes])
    width_counter = Counter([s[1] for s in sizes])

    target_h = height_counter.most_common(1)[0][0]
    target_w = width_counter.most_common(1)[0][0]

    print(f"✅ 有效砖块 {len(temp_tiles)} 个，标准化到 {target_w}x{target_h}")

    pieces: List[Piece] = []

    for idx, t in enumerate(temp_tiles):
        warped = t["warped"]
        norm = cv2.resize(
            warped,
            (target_w, target_h),
            interpolation=cv2.INTER_LANCZOS4
        )

        # 提取 4 条边
        edges = [
            extract_edge_simple(norm, 0),  # top
            extract_edge_simple(norm, 1),  # right
            extract_edge_simple(norm, 2),  # bottom
            extract_edge_simple(norm, 3),  # left
        ]

        pieces.append(Piece(
            id=idx,
            image=norm,
            edges=edges,
            init_center=t["center"],
            init_angle=t["angle"]
        ))

    return pieces


# ===================== Beam Search（平移版，无旋转） =====================

@dataclass(order=True)
class BeamState:
    cost: float
    layout: tuple  # 长度 N 的元组，每项是 piece_id
    used_mask: int


class TranslateBeamSolver:
    """
    Beam Search，只考虑“哪块放在哪个格子”，不再搜索旋转。
    默认我们认为透视回正后的砖块方向已经统一。
    """

    def __init__(self,
                 pieces: List[Piece],
                 rows: int,
                 cols: int,
                 beam_width: int):
        self.pieces = pieces
        self.N = len(pieces)
        assert rows * cols == self.N
        self.rows, self.cols = rows, cols
        self.beam_width = beam_width
        self.id2piece = {p.id: p for p in pieces}

        # 为了初期搜索更充分一点
        self.adaptive_beam = min(beam_width * 3, 2400)

        # 缓存边缘代价，加速
        self._edge_cost_cache = {}

    def _pair_cost(self,
                   pid_up: int, side_up: int,
                   pid_down: int, side_down: int) -> float:
        """
        计算两块之间指定边的匹配代价，并做缓存。
        side: 0 上, 1 右, 2 下, 3 左
        """
        key = (pid_up, side_up, pid_down, side_down)
        if key in self._edge_cost_cache:
            return self._edge_cost_cache[key]

        p_up = self.id2piece[pid_up]
        p_down = self.id2piece[pid_down]
        cost = edge_cost(p_up.edges[side_up], p_down.edges[side_down])
        self._edge_cost_cache[key] = cost
        return cost

    def solve(self):
        total_cells = self.rows * self.cols
        init_layout = tuple([-1] * total_cells)  # -1 表示未填
        states = [BeamState(0.0, init_layout, 0)]

        for pos in range(total_cells):
            r = pos // self.cols
            c = pos % self.cols

            new_states: List[BeamState] = []
            current_beam = self.adaptive_beam if pos < 6 else self.beam_width

            for st in states:
                used = st.used_mask

                for p in self.pieces:
                    if used & (1 << p.id):
                        continue

                    add_cost = 0.0
                    cnt = 0

                    # 上方约束：上方的下边 vs 当前的上边
                    if r > 0:
                        up_pid = st.layout[(r - 1) * self.cols + c]
                        if up_pid >= 0:
                            add_cost += self._pair_cost(up_pid, 2, p.id, 0)
                            cnt += 1

                    # 左侧约束：左侧的右边 vs 当前的左边
                    if c > 0:
                        left_pid = st.layout[r * self.cols + (c - 1)]
                        if left_pid >= 0:
                            add_cost += self._pair_cost(left_pid, 1, p.id, 3)
                            cnt += 1

                    if cnt > 0:
                        add_cost /= cnt  # 平均一下

                    new_cost = st.cost + add_cost
                    new_layout = list(st.layout)
                    new_layout[pos] = p.id

                    new_states.append(BeamState(
                        new_cost,
                        tuple(new_layout),
                        used | (1 << p.id)
                    ))

            new_states.sort()
            states = new_states[:current_beam]

            if states:
                print(
                    f"pos {pos:2d}  beam={current_beam:4d}  "
                    f"states={len(states):4d}  best_cost={states[0].cost:.2f}"
                )

        best = min(states, key=lambda s: s.cost)
        layout_ids = best.layout

        # 转成 2D 布局
        layout = [[-1 for _ in range(self.cols)] for _ in range(self.rows)]
        for pos, pid in enumerate(layout_ids):
            r = pos // self.cols
            c = pos % self.cols
            layout[r][c] = pid

        return layout, best.cost


# ===================== 渲染拼接后的图像 =====================

def render_solution_image(pieces: List[Piece],
                          layout: List[List[int]],
                          out_path: str):
    rows = len(layout)
    cols = len(layout[0])

    cell_h, cell_w, _ = pieces[0].image.shape
    grid_h = rows * cell_h
    grid_w = cols * cell_w

    canvas = np.zeros((grid_h, grid_w, 3), dtype=np.uint8)

    for r in range(rows):
        for c in range(cols):
            pid = layout[r][c]
            piece = next(p for p in pieces if p.id == pid)
            img = piece.image
            y0 = r * cell_h
            x0 = c * cell_w
            canvas[y0:y0 + cell_h, x0:x0 + cell_w, :] = img

    cv2.imwrite(out_path, canvas)
    print(f"\n💾 拼接结果已保存: {out_path}")

    return canvas


# ===================== 动画：图块修复过程 =====================
def make_animation(board_bgr: np.ndarray,
                   pieces: List[Piece],
                   layout: List[List[int]],
                   out_path: str,
                   fps: int = 30,
                   move_duration: float = 2.0,
                   hold_duration: float = 1.0):
    """
    动画过程：
    - 画布为纯白背景（不再显示原始整图）
    - 前 move_duration 秒：每块砖从“原始散乱位置 + 旋转角”
      慢慢移动 + 旋转到“规则网格中心 + 角度 0”
    - 后 hold_duration 秒：完全静止不动，停在修复后的状态
    """

    rows = len(layout)
    cols = len(layout[0])

    H0, W0 = board_bgr.shape[:2]
    cell_h, cell_w, _ = pieces[0].image.shape
    grid_h = rows * cell_h
    grid_w = cols * cell_w

    canvas_h = max(H0, grid_h)
    canvas_w = max(W0, grid_w)

    # 网格放在画布中央
    offset_x = (canvas_w - grid_w) // 2
    offset_y = (canvas_h - grid_h) // 2

    # 原始图块的坐标也整体居中一下（只是用来算起点）
    by = (canvas_h - H0) // 2
    bx = (canvas_w - W0) // 2

    # 计算每块砖的目标中心位置
    targets = {}
    for r in range(rows):
        for c in range(cols):
            pid = layout[r][c]
            piece = next(p for p in pieces if p.id == pid)
            cx = offset_x + c * cell_w + cell_w / 2.0
            cy = offset_y + r * cell_h + cell_h / 2.0
            targets[piece.id] = (cx, cy)

    fig, ax = plt.subplots()
    ax.set_xlim(0, canvas_w)
    ax.set_ylim(canvas_h, 0)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.subplots_adjust(0, 0, 1, 1)

    # 纯白背景
    bg = np.ones((canvas_h, canvas_w, 3), dtype=np.uint8) * 255
    ax.imshow(cv2.cvtColor(bg, cv2.COLOR_BGR2RGB),
              extent=(0, canvas_w, canvas_h, 0))

    # 每块砖的 artist
    artists = {}
    for p in pieces:
        img_rgb = cv2.cvtColor(p.image, cv2.COLOR_BGR2RGB)
        im_artist = ax.imshow(img_rgb, animated=True)
        artists[p.id] = im_artist

    move_frames = max(1, int(fps * move_duration))
    hold_frames = max(1, int(fps * hold_duration))
    total_frames = move_frames + hold_frames

    def update(frame):
        # 前 move_frames 帧：从 0 → 1 插值移动
        # 后 hold_frames 帧：t 固定为 1.0（完全静止）
        if frame < move_frames:
            t = frame / (move_frames - 1) if move_frames > 1 else 1.0
        else:
            t = 1.0

        for p in pieces:
            im_artist = artists[p.id]

            # 起点（原图中的中心），整体平移到画布中心附近
            (cx0, cy0) = p.init_center
            sx = bx + cx0
            sy = by + cy0

            # 终点（网格中的中心）
            (xt, yt) = targets[p.id]

            cx = sx + (xt - sx) * t
            cy = sy + (yt - sy) * t

            # 起始角度（minAreaRect 给的角度），终点角度 0
            start_angle = p.init_angle
            angle = start_angle * (1.0 - t)

            img_rgb = cv2.cvtColor(p.image, cv2.COLOR_BGR2RGB)
            h, w = img_rgb.shape[:2]

            # 以自身中心旋转
            M = cv2.getRotationMatrix2D((w / 2, h / 2), angle, 1.0)
            rotated = cv2.warpAffine(
                img_rgb, M, (w, h), borderValue=(255, 255, 255)
            )

            im_artist.set_data(rotated)
            im_artist.set_extent(
                (cx - w / 2, cx + w / 2, cy + h / 2, cy - h / 2)
            )

        return list(artists.values())

    ani = animation.FuncAnimation(
        fig,
        update,
        frames=total_frames,
        blit=True
    )

    print(f"🎬 生成动画: {out_path} ...")
    writer = PillowWriter(fps=fps)
    ani.save(out_path, writer=writer)
    plt.close(fig)
    print("✅ 动画已保存")

# ===================== 自动检测 rows / cols =====================

def auto_infer_rows_cols(N: int) -> Tuple[int, int]:
    """
    根据块数 N 自动推 rows / cols：
    - 找 N 的所有因子对 (r, c)，r * c = N
    - 选一个 r 最接近 sqrt(N) 的，使得布局尽量接近正方形
    """
    if N <= 0:
        raise ValueError("N 必须为正整数")

    best_r = 1
    best_diff = float("inf")
    root = math.sqrt(N)

    for r in range(1, N + 1):
        if N % r == 0:
            diff = abs(r - root)
            if diff < best_diff:
                best_diff = diff
                best_r = r

    rows = best_r
    cols = N // best_r
    return rows, cols


# ===================== main =====================

def main():
    parser = argparse.ArgumentParser(
        description="先透视回正每块砖，再做平移版 Beam Search 拼图（自动 rows/cols，无需输入 beam）"
    )
    parser.add_argument("input", help="输入拼图图片（带旋转砖块）")
    parser.add_argument("--min-area", type=int, default=500,
                        help="过滤太小轮廓的最小面积（默认 500）")
    parser.add_argument("--out", type=str, default="translate_solved.png",
                        help="输出图像文件名")
    parser.add_argument("--no-anim", action="store_true",
                        help="不生成动画则加上该参数")

    args = parser.parse_args()

    board_bgr = cv2.imread(args.input)
    if board_bgr is None:
        raise RuntimeError(f"无法读取输入图像: {args.input}")

    pieces = straighten_and_extract_pieces(board_bgr, min_area=args.min_area)
    N = len(pieces)
    print(f"\n🧩 检测到 {N} 块砖")

    # 自动检测 rows / cols（规则矩形网格）
    rows, cols = auto_infer_rows_cols(N)
    print(f"📐 自动检测行列数: {rows} x {cols}  （总块数 = {N}）")

    if rows * cols != N:
        raise ValueError(f"rows*cols = {rows*cols} != {N}")

    # beam 宽度内部自动给定，用户不用输入
    if N <= 9:
        beam_width = 400
    elif N <= 20:
        beam_width = 800
    else:
        beam_width = 1200

    print(f"\n开始平移拼图求解 {rows}x{cols} ...  (beam_width = {beam_width})")
    solver = TranslateBeamSolver(
        pieces, rows=rows, cols=cols, beam_width=beam_width
    )
    layout, cost = solver.solve()
    print(f"\n✅ 完成! 总代价 = {cost:.2f}")

    _ = render_solution_image(pieces, layout, args.out)

    if not args.no_anim:
        anim_path = "translate_animation.gif"
        make_animation(board_bgr, pieces, layout, anim_path)


if __name__ == "__main__":
    main()