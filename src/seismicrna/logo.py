import io
import math

import numpy as np
import matplotlib
import matplotlib.patches as patches
import matplotlib.pyplot as plt

matplotlib.use('Agg')

def compute_arc_points(center, radius, theta1, theta2, n=100):
    if theta2 < theta1:
        angles = np.linspace(theta1, theta2 + 360, n) % 360
    else:
        angles = np.linspace(theta1, theta2, n)
    angles_rad = np.radians(angles)
    x = center[0] + radius * np.cos(angles_rad)
    y = center[1] + radius * np.sin(angles_rad)
    return np.column_stack((x, y))

def draw_seismic_logo(report: bool = False,
                      out_svg: bool | None = None,
                      dpi: int = 300):
    fig, ax = plt.subplots()

    R = 50
    r = R / 3
    large_linewidth = 12
    small_linewidth = 6
    arc_fudge_deg = 1

    blue = "#2186d9"
    grey = "#242424"

    box_color = grey if not report else "none"

    # Precalculate centers with translation and 45° rotation
    cos45 = sin45 = math.sqrt(2) / 2
    def transform_point(x, y):
        xt, yt = x - 25, y
        return (xt * cos45 - yt * sin45, xt * sin45 + yt * cos45)

    center_large_left  = transform_point(50, 0)
    center_large_right = transform_point(0, 0)
    center_outer_left  = transform_point(100, 0)
    center_outer_right = transform_point(-50, 0)

    def adjust_angles(theta1, theta2):
        return theta1 + 45, theta2 + 45

    arcs_info = [
        (center_large_left,  R, *adjust_angles(-arc_fudge_deg, 180 + arc_fudge_deg), large_linewidth, grey),
        (center_large_right, R, *adjust_angles(180 - arc_fudge_deg, 360 + arc_fudge_deg), large_linewidth, blue),
        (center_large_left,  r, *adjust_angles(-arc_fudge_deg, 180 + arc_fudge_deg), small_linewidth, grey),
        (center_large_right, r, *adjust_angles(180 - arc_fudge_deg, 360 + arc_fudge_deg), small_linewidth, blue),
        (center_outer_left,  r, *adjust_angles(180 - arc_fudge_deg, 360 + arc_fudge_deg), small_linewidth, blue),
        (center_outer_right, r, *adjust_angles(-arc_fudge_deg, 180 + arc_fudge_deg), small_linewidth, grey),
    ]

    points_list = []
    for center, radius, theta1, theta2, lw, color in arcs_info:
        arc_patch = patches.Arc(
            center, 2 * radius, 2 * radius,
            angle=0, theta1=theta1, theta2=theta2,
            linewidth=lw, edgecolor=color, facecolor="none", zorder=2
        )
        ax.add_patch(arc_patch)
        pts = compute_arc_points(center, radius, theta1, theta2, n=200)
        points_list.append(pts)
    all_points = np.concatenate(points_list, axis=0)

    min_x, min_y = np.min(all_points, axis=0)
    max_x, max_y = np.max(all_points, axis=0)
    box_center_x = (min_x + max_x) / 2
    box_center_y = (min_y + max_y) / 2
    box_side = max(max_x - min_x, max_y - min_y)
    margin = 5
    box_side += 2 * margin
    box_left = box_center_x - box_side / 2
    box_bottom = box_center_y - box_side / 2

    top_left = (box_left - 0.2, box_bottom + box_side + 0.2)
    top_right = (box_left + box_side, box_bottom + box_side + 0.2)
    bottom_left = (box_left - 0.2, box_bottom)
    bottom_right = (box_left + box_side + 0.2, box_bottom)

    box_patch = patches.FancyBboxPatch(
        (box_left, box_bottom), box_side, box_side,
        boxstyle=patches.BoxStyle("Round", pad=0.2, rounding_size=25),
        linewidth=0, edgecolor="none", facecolor=box_color, zorder=0
    )
    clip_path = box_patch.get_path().transformed(box_patch.get_transform())
    if report:
        bottom_right_triangle = patches.Polygon(
        [bottom_right, top_right, bottom_left],
        closed=True, facecolor=grey, edgecolor=None, zorder=1,
        antialiased=True)
        bottom_right_triangle.set_clip_path(clip_path, ax.transData)
        ax.add_patch(bottom_right_triangle)
        pad = 0
    else:
        upper_left_triangle = patches.Polygon(
        [top_left, top_right, bottom_left],
        closed=True, facecolor=blue, edgecolor=None, zorder=1,
        antialiased=True)
        upper_left_triangle.set_clip_path(clip_path, ax.transData)
        ax.add_patch(upper_left_triangle)
        ax.add_patch(box_patch)
        pad = 10

    ax.set_xlim(box_left - pad, box_left + box_side + pad)
    ax.set_ylim(box_bottom - pad, box_bottom + box_side + pad)
    ax.set_aspect('equal')
    plt.axis('off')

    svg_buffer = io.StringIO()
    if not out_svg:
        fig.savefig(svg_buffer, format='svg', transparent=True, bbox_inches='tight', dpi=dpi)
        plt.close(fig)
        svg_string = svg_buffer.getvalue()
        svg_buffer.close()

        return svg_string
    else:
        fig.savefig(out_svg, format='svg', transparent=True, bbox_inches='tight', dpi=dpi)
