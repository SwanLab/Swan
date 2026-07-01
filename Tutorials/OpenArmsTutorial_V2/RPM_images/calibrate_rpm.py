"""
calibrate.py  —  Interactive calibration tool for read_rpm.py
=============================================================

WHAT IT DOES:
    Lets you visually adjust the crop region, circle center, radius,
    and RPM reference angles on a sample image. When you press ESC,
    it prints the parameters ready to paste into read_rpm.py.

HOW TO USE:
    1. Place this file in your RPM_images folder alongside a sample image.
    2. Run:  python calibrate.py
    3. Adjust using the keys below, then press ESC to get the parameters.

KEYBOARD CONTROLS:
    ── Crop region (the rectangle cut out of the original image) ──
        Q / A   move crop LEFT / RIGHT
        W / S   move crop UP   / DOWN
        E / D   shrink / grow crop WIDTH
        R / F   shrink / grow crop HEIGHT

    ── Circle center (inside the crop) ──
        J / L   move center LEFT / RIGHT
        I / K   move center UP   / DOWN

    ── Circle radius ──
        +        increase radius
        -        decrease radius

    ── RPM reference angles ──
        1        select   0 RPM point  (green)
        2        select 200 RPM point  (orange)
        3        select 300 RPM point  (yellow)
        4        select 500 RPM point  (blue)
        Z / X   rotate selected point  -1° / +1°

    ── General ──
        ESC     quit and print parameters to copy into read_rpm.py

TIPS:
    - Start by getting the crop right (you should see only the gauge).
    - Then center the circle on the needle pivot.
    - Adjust the radius so the circle sits just inside the gauge rim.
    - For each RPM point, press 1-4 to select it, then Z/X to move the
      coloured dot until it lines up with the matching graduation mark.
"""

import cv2
import numpy as np
import sys
from pathlib import Path

# ── Image to use for calibration ─────────────────────────────────────────────
# The script will automatically use the first .jpg found in the folder
# if this file does not exist.
IMAGE_PATH = Path(__file__).parent / "captura_23-28-06.jpg"

# ── Starting parameters (edit these if you want a different starting point) ───
crop_x, crop_y, crop_w, crop_h = 80, 245, 200, 235
cx, cy  = 90, 140
radius  = 70

# Reference angles in trigonometric degrees (0°=right, counter-clockwise)
angles = {
      0: 85,
    200: 41,
    300: 19,
    500: 332,
}

COLORS = {
      0: (  0, 200,   0),   # green
    200: (255, 100,   0),   # orange
    300: (255, 165,   0),   # yellow
    500: (  0,   0, 255),   # blue
}

selected_rpm = 0
SCALE = 4


def draw(img_orig):
    crop = img_orig[crop_y:crop_y+crop_h, crop_x:crop_x+crop_w]
    big  = cv2.resize(crop, (crop_w*SCALE, crop_h*SCALE), interpolation=cv2.INTER_LINEAR)

    cx4, cy4 = cx*SCALE, cy*SCALE
    R = radius*SCALE

    # Main circle
    cv2.circle(big, (cx4, cy4), R, (0, 255, 0), 2)
    cv2.circle(big, (cx4, cy4), 6, (0, 0, 255), -1)

    # Degree markings every 15°
    for a_trigo in range(0, 360, 15):
        a_cv  = (-a_trigo) % 360
        rad   = np.radians(a_cv)
        ex    = int(cx4 + R * np.cos(rad))
        ey    = int(cy4 + R * np.sin(rad))
        ex2   = int(cx4 + (R + 18) * np.cos(rad))
        ey2   = int(cy4 + (R + 18) * np.sin(rad))
        cv2.line(big, (ex, ey), (ex2, ey2), (160, 160, 0), 1)
        cv2.putText(big, str(a_trigo), (ex2-10, ey2+5),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.38, (0, 210, 210), 1)

    # RPM reference dots + active radius line
    for rpm_val, a_trigo in angles.items():
        col   = COLORS[rpm_val]
        a_cv  = (-a_trigo) % 360
        rad   = np.radians(a_cv)
        px    = int(cx4 + (R - 18) * np.cos(rad))
        py    = int(cy4 + (R - 18) * np.sin(rad))
        thick = 2 if rpm_val == selected_rpm else -1
        cv2.circle(big, (px, py), 10, col, thick)
        cv2.putText(big, f"{rpm_val}={a_trigo}deg",
                    (px + 8, py - 6), cv2.FONT_HERSHEY_SIMPLEX, 0.44, col, 1)

    # Line toward the currently selected RPM
    a_sel = angles[selected_rpm]
    a_cv  = (-a_sel) % 360
    rad   = np.radians(a_cv)
    ex    = int(cx4 + R * np.cos(rad))
    ey    = int(cy4 + R * np.sin(rad))
    cv2.line(big, (cx4, cy4), (ex, ey), COLORS[selected_rpm], 2)

    # HUD
    hud = [
        f"CROP   x={crop_x}  y={crop_y}  w={crop_w}  h={crop_h}",
        f"CENTER cx={cx}  cy={cy}    RADIUS={radius}",
        f"Selected: {selected_rpm} RPM = {angles[selected_rpm]} deg",
        "Crop: Q/A=left/right  W/S=up/down  E/D=width  R/F=height",
        "Circle center: J/L=left/right  I/K=up/down     Radius: +/-",
        "RPM points: 1=0rpm  2=200rpm  3=300rpm  4=500rpm  angle: Z/X",
        "ESC = quit and print parameters",
    ]
    for idx, line in enumerate(hud):
        cv2.putText(big, line, (8, 18 + idx*19),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.42, (255, 255, 255), 1, cv2.LINE_AA)

    return big


def main():
    global crop_x, crop_y, crop_w, crop_h, cx, cy, radius, selected_rpm

    if IMAGE_PATH.exists():
        img_path = IMAGE_PATH
    else:
        candidates = sorted(
            list(Path(__file__).parent.glob("*.jpg")) +
            list(Path(__file__).parent.glob("*.jpeg"))
        )
        if not candidates:
            print("No .jpg image found in this folder. Place one here and retry.")
            sys.exit(1)
        img_path = candidates[0]
        print(f"Using image: {img_path.name}")

    img_orig = cv2.imread(str(img_path))
    if img_orig is None:
        print(f"Could not read: {img_path}")
        sys.exit(1)

    cv2.namedWindow("RPM Calibration", cv2.WINDOW_NORMAL)

    while True:
        cv2.imshow("RPM Calibration", draw(img_orig))
        key = cv2.waitKey(30) & 0xFF

        if key == 27:   # ESC
            break

        # ── Crop ──────────────────────────────────────────────────────────
        elif key == ord('q'): crop_x = max(0, crop_x - 2)
        elif key == ord('a'): crop_x = min(img_orig.shape[1] - crop_w, crop_x + 2)
        elif key == ord('w'): crop_y = max(0, crop_y - 2)
        elif key == ord('s'): crop_y = min(img_orig.shape[0] - crop_h, crop_y + 2)
        elif key == ord('e'): crop_w = max(50, crop_w - 2)
        elif key == ord('d'): crop_w = min(img_orig.shape[1] - crop_x, crop_w + 2)
        elif key == ord('r'): crop_h = max(50, crop_h - 2)
        elif key == ord('f'): crop_h = min(img_orig.shape[0] - crop_y, crop_h + 2)

        # ── Circle center  (J/L = left/right,  I/K = up/down) ────────────
        elif key == ord('j'): cx = max(0,      cx - 1)
        elif key == ord('l'): cx = min(crop_w, cx + 1)
        elif key == ord('i'): cy = max(0,      cy - 1)
        elif key == ord('k'): cy = min(crop_h, cy + 1)

        # ── Radius ────────────────────────────────────────────────────────
        elif key == ord('+'): radius = min(300, radius + 1)
        elif key == ord('-'): radius = max(10,  radius - 1)

        # ── Select RPM point ──────────────────────────────────────────────
        elif key == ord('1'): selected_rpm = 0
        elif key == ord('2'): selected_rpm = 200
        elif key == ord('3'): selected_rpm = 300
        elif key == ord('4'): selected_rpm = 500

        # ── Adjust angle of selected RPM point ────────────────────────────
        elif key == ord('z'): angles[selected_rpm] = (angles[selected_rpm] - 1) % 360
        elif key == ord('x'): angles[selected_rpm] = (angles[selected_rpm] + 1) % 360

    cv2.destroyAllWindows()

    arc = (angles[0] - angles[500]) % 360

    print("\n" + "=" * 55)
    print("  Copy these into read_rpm.py :")
    print("=" * 55)
    print(f"CROP        = ({crop_x}, {crop_y}, {crop_w}, {crop_h})")
    print(f"CENTER      = ({cx}, {cy})")
    print(f"RADIUS      = {radius}")
    print()
    print("CALIBRATION_POINTS = [")
    for rpm_val, a in sorted(angles.items()):
        print(f"    ({a:>3}, {rpm_val:>3}),")
    print("]")
    print()
    print(f"# Computed arc: {arc} degrees")
    print("=" * 55)


if __name__ == "__main__":
    main()