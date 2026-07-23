"""
read_rpm.py
-----------
Lit automatiquement la valeur RPM sur chaque image du dossier RPM_images
via OpenCV (détection d'aiguille par seuillage + PCA), puis écrit les
valeurs dans la colonne "rpm" de ../merged_output.csv.

La caméra ayant bougé plusieurs fois, le script utilise des paramètres de
calibration différents selon le timestamp de chaque image (voir SEGMENTS).
Les périodes où le cadran est illisible sont blacklistées.

Structure attendue :
    OpenArmsTutorial_V2/
    ├── merged_output.csv
    └── RPM_images/
        ├── read_rpm.py
        ├── captura_18-04-47.jpg
        └── ...

Prérequis :
    pip install opencv-python-headless pandas numpy

Usage :
    source venv/bin/activate
    python read_rpm.py
"""

import random
import re
import sys
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

# ── Chemins ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
CSV_PATH   = SCRIPT_DIR.parent / "merged_output.csv"

# ── Dates de la traversée ─────────────────────────────────────────────────────
DATE_EVENING  = "2025-05-15"   # heures >= MIDNIGHT_HOUR
DATE_MORNING  = "2025-05-16"   # heures <  MIDNIGHT_HOUR
MIDNIGHT_HOUR = 12

# ── Paramètres fixes ──────────────────────────────────────────────────────────
HUB_RADIUS = 15   # rayon du moyeu à ignorer (identique pour toutes les configs)

SUPPORTED_EXT = {".jpg", ".jpeg", ".png", ".bmp", ".tiff", ".webp"}

# ── Mode debug ────────────────────────────────────────────────────────────────
# Sauvegarde N images choisies aléatoirement, réparties sur tout le dataset.
# Mettre à 0 pour désactiver.
DEBUG_N_IMAGES = 10

# ── Calibrations ──────────────────────────────────────────────────────────────
# Chaque entrée : crop (x,y,w,h), center (cx,cy), radius, calibration_points
# calibration_points : [(angle_trigo, rpm), ...] en sens horaire visuel

CALIBRATIONS = {
    "A": {
        "crop":   (82, 207, 198, 235),
        "center": (90, 140),
        "radius": 62,
        "points": [(86, 0), (39, 200), (16, 300), (328, 500)],
    },
    "partial": {
        # Cadran à moitié visible (19:58:44 → 20:33:48) — 500 RPM estimé
        "crop":   (80, 245, 200, 235),
        "center": (84, 226),
        "radius": 66,
        "points": [(84, 0), (39, 200), (16, 300), (331, 500)],
    },
    "B": {
        "crop":   (80, 239, 200, 235),
        "center": (92, 141),
        "radius": 61,
        "points": [(86, 0), (40, 200), (16, 300), (327, 500)],
    },
    "C": {
        "crop":   (80, 245, 200, 235),
        "center": (88, 189),
        "radius": 63,
        "points": [(85, 0), (39, 200), (15, 300), (325, 500)],
    },
    "D": {
        "crop":   (82, 219, 200, 235),
        "center": (93, 123),
        "radius": 60,
        "points": [(86, 0), (39, 200), (16, 300), (328, 500)],
    },
}

# ── Segments temporels ────────────────────────────────────────────────────────
# (début_inclusif, fin_exclusif, clé_calibration_ou_None)
# None = période blacklistée → images marquées NaN directement sans traitement

SEGMENTS = [
    ("2025-05-15 18:04:47", "2025-05-15 19:58:44", "A"),
    ("2025-05-15 19:58:44", "2025-05-15 20:33:48", "partial"),
    ("2025-05-15 20:33:48", "2025-05-16 01:59:07", "B"),
    ("2025-05-16 01:59:07", "2025-05-16 02:06:30", "C"),
    ("2025-05-16 02:06:30", "2025-05-16 04:33:36", None),       # blacklist
    ("2025-05-16 04:33:36", "2025-05-16 10:52:59", "D"),
]

# Préconvertir les bornes en datetime pour comparer rapidement
_SEGMENTS_DT = [
    (pd.Timestamp(s), pd.Timestamp(e), k)
    for s, e, k in SEGMENTS
]


def get_calibration(img_dt: pd.Timestamp):
    """
    Retourne (cal_dict, key) pour le timestamp donné.
    cal_dict est None si la période est blacklistée.
    Retourne (None, None) si hors de tous les segments.
    """
    for start, end, key in _SEGMENTS_DT:
        if start <= img_dt < end:
            return (CALIBRATIONS[key] if key is not None else None), key
    return None, None


# ── Mapping angle → RPM ───────────────────────────────────────────────────────

def angle_to_rpm_piecewise(needle_trigo: float, points: list) -> float | None:
    """
    Convertit un angle (convention trigo) en RPM via interpolation par morceaux.
    L'aiguille tourne en sens HORAIRE visuel = décroissant en trigo.
    """
    angle_0 = points[0][0]

    def horaire_dist(a):
        return (angle_0 - a) % 360

    needle_dist = horaire_dist(needle_trigo)

    segments = []
    for i in range(len(points) - 1):
        a1, r1 = points[i]
        a2, r2 = points[i + 1]
        segments.append((horaire_dist(a1), horaire_dist(a2), r1, r2))

    for d1, d2, r1, r2 in segments:
        if d1 <= needle_dist <= d2:
            ratio = (needle_dist - d1) / (d2 - d1)
            return r1 + ratio * (r2 - r1)

    total_arc = horaire_dist(points[-1][0])
    margin = total_arc * 0.05

    if -margin <= needle_dist < 0:
        return float(points[0][1])
    if total_arc < needle_dist <= total_arc + margin:
        return float(points[-1][1])

    return None


# ── Image debug illisible ─────────────────────────────────────────────────────

def make_unreadable_debug_image(image_path: Path, cal: dict, reason: str = "UNREADABLE"):
    img = cv2.imread(str(image_path))
    x, y, w, h = cal["crop"]

    if img is not None:
        crop = img[y:y+h, x:x+w].copy()
    else:
        crop = np.full((h, w, 3), 60, dtype=np.uint8)

    dbg = cv2.resize(crop, (crop.shape[1]*2, crop.shape[0]*2))
    H, W = dbg.shape[:2]

    cv2.line(dbg, (0, 0), (W, H), (0, 0, 255), 3)
    cv2.line(dbg, (W, 0), (0, H), (0, 0, 255), 3)
    cv2.rectangle(dbg, (0, 0), (W, 26), (0, 0, 200), -1)
    cv2.putText(dbg, reason, (8, 19),
                cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 255), 2)
    return dbg


# ── Détection RPM ─────────────────────────────────────────────────────────────

def detect_rpm(image_path: Path, cal: dict, debug: bool = False):
    """
    Retourne (rpm, debug_img) si debug=True, sinon juste rpm.
    rpm est un float ou None si illisible. Ne lève jamais d'exception.
    """
    try:
        img = cv2.imread(str(image_path))
        if img is None:
            return (None, None) if debug else None

        x, y, w, h = cal["crop"]
        cx, cy     = cal["center"]
        radius     = cal["radius"]
        points     = cal["points"]

        crop = img[y:y+h, x:x+w]
        gray = cv2.cvtColor(crop, cv2.COLOR_BGR2GRAY)

        mask = np.zeros_like(gray)
        cv2.circle(mask, (cx, cy), radius, 255, -1)
        gray_m = cv2.bitwise_and(gray, gray, mask=mask)

        blur = cv2.GaussianBlur(gray_m, (5, 5), 0)
        _, thresh = cv2.threshold(blur, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
        thresh = cv2.bitwise_and(thresh, thresh, mask=mask)

        cv2.circle(thresh, (cx, cy), HUB_RADIUS, 0, -1)
        border = np.zeros_like(thresh)
        cv2.circle(border, (cx, cy), radius - 5, 255, -1)
        thresh = cv2.bitwise_and(thresh, border)

        pts_nonzero = np.column_stack(np.where(thresh > 0))
        if len(pts_nonzero) < 20:
            return (None, None) if debug else None

        data = pts_nonzero[:, ::-1].astype(np.float32)
        mean, eigvecs = cv2.PCACompute(data, mean=None)

        angle1 = np.degrees(np.arctan2(eigvecs[0, 1], eigvecs[0, 0]))
        angle2 = angle1 + 180
        mx, my = float(mean[0, 0]), float(mean[0, 1])
        dx, dy = mx - cx, my - cy
        dot1 = dx * np.cos(np.radians(angle1)) + dy * np.sin(np.radians(angle1))
        dot2 = dx * np.cos(np.radians(angle2)) + dy * np.sin(np.radians(angle2))
        needle = (angle1 if dot1 > dot2 else angle2) % 360

        rpm = angle_to_rpm_piecewise(needle, points)
        if rpm is None:
            return (None, None) if debug else None

        if debug:
            dbg   = crop.copy()
            cx4, cy4 = cx * 2, cy * 2
            dbg   = cv2.resize(dbg, (dbg.shape[1]*2, dbg.shape[0]*2))
            R2    = radius * 2

            cv2.circle(dbg, (cx4, cy4), R2, (0, 255, 0), 1)
            cv2.circle(dbg, (cx4, cy4), 5, (0, 0, 255), -1)

            needle_cv = (-needle) % 360
            rad = np.radians(needle_cv)
            ex  = int(cx4 + R2 * np.cos(rad))
            ey  = int(cy4 + R2 * np.sin(rad))
            cv2.arrowedLine(dbg, (cx4, cy4), (ex, ey), (0, 0, 255), 2, tipLength=0.2)

            cal_colors = [(0,200,0), (255,100,0), (255,165,0), (0,0,255)]
            for (a_trigo, rpm_ref), col in zip(points, cal_colors):
                a_cv = (-a_trigo) % 360
                r2   = np.radians(a_cv)
                px   = int(cx4 + (R2-10) * np.cos(r2))
                py   = int(cy4 + (R2-10) * np.sin(r2))
                cv2.circle(dbg, (px, py), 5, col, -1)
                cv2.putText(dbg, str(rpm_ref), (px+4, py-4),
                            cv2.FONT_HERSHEY_SIMPLEX, 0.4, col, 1)

            cv2.putText(dbg, f"{round(rpm)} RPM", (8, 20),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.6, (255, 255, 0), 2)

            th_rgb  = cv2.cvtColor(thresh, cv2.COLOR_GRAY2BGR)
            th_rgb2 = cv2.resize(th_rgb, (dbg.shape[1]//3, dbg.shape[0]//3))
            dbg[-th_rgb2.shape[0]:, -th_rgb2.shape[1]:] = th_rgb2

            return round(rpm), dbg

        return round(rpm)

    except Exception:
        return (None, None) if debug else None


# ── Interpolation ─────────────────────────────────────────────────────────────

def interpolate_missing(df: pd.DataFrame) -> tuple[pd.DataFrame, int, int]:
    """Remplace les NaN par la moyenne du voisin valide avant + après."""
    rpm    = df["rpm"].copy()
    missing_idx = rpm[rpm.isna()].index
    filled = 0

    for idx in missing_idx:
        pos      = df.index.get_loc(idx)
        prev_val = next((rpm.iloc[p] for p in range(pos-1, -1, -1)   if pd.notna(rpm.iloc[p])), None)
        next_val = next((rpm.iloc[p] for p in range(pos+1, len(rpm)) if pd.notna(rpm.iloc[p])), None)

        if prev_val is not None and next_val is not None:
            rpm.at[idx] = round((prev_val + next_val) / 2)
        elif prev_val is not None:
            rpm.at[idx] = prev_val
        elif next_val is not None:
            rpm.at[idx] = next_val

        if pd.notna(rpm.at[idx]):
            filled += 1

    df["rpm"] = rpm
    return df, len(missing_idx), filled


# ── Utilitaires ───────────────────────────────────────────────────────────────

def parse_image_time(filename: str) -> pd.Timestamp:
    m = re.search(r'(\d{2})-(\d{2})-(\d{2})', filename)
    if not m:
        raise ValueError(f"Heure introuvable dans : {filename}")
    hh, mm, ss = m.group(1), m.group(2), m.group(3)
    date = DATE_MORNING if int(hh) < MIDNIGHT_HOUR else DATE_EVENING
    return pd.Timestamp(f"{date} {hh}:{mm}:{ss}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if not CSV_PATH.exists():
        print(f"Erreur : CSV introuvable → {CSV_PATH}")
        sys.exit(1)

    images = sorted([f for f in SCRIPT_DIR.iterdir() if f.suffix.lower() in SUPPORTED_EXT])
    if not images:
        print(f"Aucune image trouvée dans {SCRIPT_DIR}")
        sys.exit(1)

    print(f"{len(images)} images trouvées.\n")

    debug_indices = set()
    if DEBUG_N_IMAGES > 0:
        n_pick = min(DEBUG_N_IMAGES, len(images))
        debug_indices = set(random.sample(range(1, len(images) + 1), n_pick))

    df = pd.read_csv(CSV_PATH)
    if "Time" not in df.columns:
        print(f"Erreur : colonne 'Time' absente. Colonnes : {list(df.columns)}")
        sys.exit(1)

    df["_time_dt"] = pd.to_datetime(df["Time"])
    if "rpm" not in df.columns:
        df["rpm"] = pd.NA

    unreadable   = 0
    blacklisted  = 0
    out_of_range = 0

    for i, img_path in enumerate(images, 1):
        try:
            img_dt = parse_image_time(img_path.name)
        except ValueError as e:
            print(f"[{i}/{len(images)}] {img_path.name} → IGNORÉ ({e})")
            continue

        cal, cal_key = get_calibration(img_dt)
        show_debug   = (i in debug_indices)
        debug_dir    = SCRIPT_DIR / "debug"

        label = f"[{i}/{len(images)}] {img_path.name}  ({img_dt.strftime('%H:%M:%S')} — cal:{cal_key})"

        # ── Période blacklistée ───────────────────────────────────────────
        if cal_key is None and cal is None:
            print(f"{label} → BLACKLIST (NaN)")
            # Nearest-match dans le CSV pour marquer NaN
            deltas = (df["_time_dt"] - img_dt).abs()
            if deltas.min().total_seconds() <= 8:
                df.loc[deltas.idxmin(), "rpm"] = float("nan")
            if show_debug:
                debug_dir.mkdir(exist_ok=True)
                placeholder_cal = list(CALIBRATIONS.values())[0]
                dbg = make_unreadable_debug_image(img_path, placeholder_cal, "BLACKLISTED")
                cv2.imwrite(str(debug_dir / f"debug_{i:04d}_{img_path.stem}.jpg"), dbg)
            blacklisted += 1
            continue

        # ── Image hors de tout segment ────────────────────────────────────
        if cal is None:
            print(f"{label} → hors segments (ignoré)")
            out_of_range += 1
            continue

        # ── Nearest-match CSV ─────────────────────────────────────────────
        deltas    = (df["_time_dt"] - img_dt).abs()
        min_delta = deltas.min()
        closest_idx = deltas.idxmin()

        if min_delta.total_seconds() > 8:
            print(f"{label} → hors plage CSV (écart {min_delta.total_seconds():.1f}s)")
            out_of_range += 1
            continue

        # ── Détection ─────────────────────────────────────────────────────
        if show_debug:
            rpm, dbg_img = detect_rpm(img_path, cal, debug=True)
        else:
            rpm    = detect_rpm(img_path, cal)
            dbg_img = None

        if show_debug:
            debug_dir.mkdir(exist_ok=True)
            out_path = debug_dir / f"debug_{i:04d}_{img_path.stem}.jpg"
            if dbg_img is not None:
                cv2.imwrite(str(out_path), dbg_img)
            else:
                cv2.imwrite(str(out_path), make_unreadable_debug_image(img_path, cal))

        if rpm is None:
            df.loc[closest_idx, "rpm"] = float("nan")
            print(f"{label} → illisible (NaN, sera interpolé)")
            unreadable += 1
        else:
            df.loc[closest_idx, "rpm"] = rpm
            print(f"{label} → {rpm} RPM  (écart {min_delta.total_seconds():.1f}s)")

    # ── Résumé ────────────────────────────────────────────────────────────
    print(f"\n── Résumé ───────────────────────────────────────────")
    print(f"Blacklistées : {blacklisted}")
    print(f"Hors plage   : {out_of_range}")
    print(f"Illisibles   : {unreadable}")

    print(f"\n── Interpolation ────────────────────────────────────")
    df, total_missing, filled = interpolate_missing(df)
    print(f"Valeurs interpolées : {filled}/{total_missing}")
    if total_missing - filled:
        print(f"Valeurs restées NaN : {total_missing - filled}")

    df.drop(columns=["_time_dt"], inplace=True)
    df.to_csv(CSV_PATH, index=False)
    print(f"\nTerminé ! CSV mis à jour : {CSV_PATH}")


if __name__ == "__main__":
    main()