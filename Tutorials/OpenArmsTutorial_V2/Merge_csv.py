"""
merge_csv.py — Fusion, nettoyage et visualisation de fichiers CSV

Structure attendue :
    projet/
    ├── merge_csv.py        ← ce script
    └── merge_csv/          ← dossier contenant tous les CSV
        ├── 01_18.csv
        ├── 02_1930.csv
        └── ...

Lancement :
    python merge_csv.py

Sorties créées à côté du script :
    projet/merged_output.csv
    projet/data_plots/speed_over_ground.png
    projet/data_plots/course_over_ground.png
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path

# ─── Configuration ────────────────────────────────────────────────────────────

SCRIPT_DIR    = Path(__file__).parent
SOURCE_FOLDER = SCRIPT_DIR / "merge_csv"
OUTPUT_FILE   = SCRIPT_DIR / "merged_output.csv"
PLOTS_FOLDER  = SCRIPT_DIR / "data_plots"

# Colonnes utilisées pour la détection de bugs GPS
GPS_COLS      = ["Latitude", "Longitude", "SpeedOverGround (m/s)"]

# Nombre de répétitions consécutives pour détecter un bug GPS
GPS_BUG_N     = 5

# Colonnes exclues du dédoublonnage
DEDUP_EXCLUDE = ["id", "Time"]

# ─── Troncatures manuelles ────────────────────────────────────────────────────
# Toutes les lignes APRÈS cette date seront ignorées dans le fichier final.
# Les fichiers sources ne sont PAS modifiés.
# Ajouter autant d'entrées que nécessaire si d'autres fichiers ont le même problème.
#
# Format : "nom_fichier.csv" -> "YYYY-MM-DD HH:MM:SS"
FILE_CUTOFFS = {
    "08_0130.csv": "2025-05-16 04:28:35",
}

# ─── Fonctions ────────────────────────────────────────────────────────────────

def load_csv(path: Path) -> pd.DataFrame | None:
    """Charge un fichier CSV et parse la colonne Time (2ème colonne)."""
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"  ❌ Impossible de lire {path.name} : {e}")
        return None

    if df.empty:
        print(f"  ⚠️  Fichier vide, ignoré : {path.name}")
        return None

    if df.shape[1] < 2:
        print(f"  ⚠️  Moins de 2 colonnes dans {path.name}, ignoré.")
        return None

    time_col = df.columns[1]
    df[time_col] = pd.to_datetime(df[time_col], errors="coerce")

    n_invalid = df[time_col].isna().sum()
    if n_invalid > 0:
        print(f"  ⚠️  {path.name} : {n_invalid} ligne(s) avec date invalide")

    return df


def remove_gps_bugs(df: pd.DataFrame, n: int = 5) -> pd.DataFrame:
    """
    Supprime les lignes où Latitude, Longitude et SpeedOverGround
    sont constantes sur n lignes consécutives ET vitesse > 0.
    (vitesse = 0 + position fixe = arrêt réel → conservé)
    """
    cols_present = [c for c in GPS_COLS if c in df.columns]
    if len(cols_present) < len(GPS_COLS):
        missing = set(GPS_COLS) - set(cols_present)
        print(f"  ⚠️  Colonnes GPS manquantes {missing}, détection de bugs ignorée.")
        return df

    to_delete = pd.Series(False, index=df.index)
    values    = df[GPS_COLS].reset_index(drop=True)
    nrows     = len(values)

    i = 0
    while i <= nrows - n:
        window = values.iloc[i:i + n]

        all_const = all(window[c].nunique(dropna=False) == 1 for c in GPS_COLS)
        ref_speed = window["SpeedOverGround (m/s)"].iloc[0]
        is_moving = pd.notna(ref_speed) and float(ref_speed) > 0.0

        if all_const and is_moving:
            # Étend la zone bug tant que les valeurs restent identiques
            j = i + n
            while j < nrows and all(
                values[c].iloc[j] == values[c].iloc[i] for c in GPS_COLS
            ):
                j += 1

            original_idx = df.index[i:j]
            to_delete[original_idx] = True
            i = j
        else:
            i += 1

    n_deleted = to_delete.sum()
    if n_deleted > 0:
        print(f"  🛰️  {n_deleted} ligne(s) supprimée(s) — bug GPS (valeurs fixes ≥{n} lignes, vitesse > 0)")
    else:
        print(f"  ✅ Aucun bug GPS détecté")

    return df[~to_delete].reset_index(drop=True)


def apply_cutoff(df: pd.DataFrame, filename: str) -> pd.DataFrame:
    """
    Applique une troncature temporelle si le fichier est dans FILE_CUTOFFS.
    Ne modifie pas le fichier source, uniquement les données en mémoire.
    """
    if filename not in FILE_CUTOFFS:
        return df

    cutoff   = pd.to_datetime(FILE_CUTOFFS[filename])
    time_col = df.columns[1]
    before   = len(df)
    df       = df[df[time_col] <= cutoff].reset_index(drop=True)
    after    = len(df)

    print(f"  ✂️  Tronqué à {FILE_CUTOFFS[filename]} : {before - after} ligne(s) ignorée(s) dans le fichier final")

    return df


def make_plots(df: pd.DataFrame, plots_dir: Path) -> None:
    """Génère les graphiques SpeedOverGround et CourseOverGround en PNG."""
    plots_dir.mkdir(parents=True, exist_ok=True)

    time_col = df.columns[1]
    valid    = df[time_col].notna()
    t        = df.loc[valid, time_col]

    def save_plot(col_name: str, ylabel: str, color: str, filename: str) -> None:
        if col_name not in df.columns:
            print(f"  ⚠️  Colonne '{col_name}' introuvable, graphique ignoré.")
            return

        fig, ax = plt.subplots(figsize=(14, 5))
        ax.scatter(t, df.loc[valid, col_name], color=color, s=1)
        ax.set_title(col_name, fontsize=13)
        ax.set_xlabel("Time")
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%d/%m"))
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        fig.autofmt_xdate(rotation=0, ha="center")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        out = plots_dir / filename
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"  ✅ Enregistré : {filename}")

    save_plot("SpeedOverGround (m/s)",  "Speed (m/s)",     "royalblue",  "speed_over_ground.png")
    save_plot("CourseOverGround (deg)", "Direction (deg)", "darkorange", "course_over_ground.png")


# ─── Programme principal ──────────────────────────────────────────────────────

def main():
    print(f"📁 Dossier du script  : {SCRIPT_DIR}")
    print(f"📂 Dossier CSV source : {SOURCE_FOLDER}")
    print(f"💾 Fichier de sortie  : {OUTPUT_FILE}")
    print(f"📊 Dossier graphiques : {PLOTS_FOLDER}")
    print()

    if not SOURCE_FOLDER.is_dir():
        print(f"❌ Dossier '{SOURCE_FOLDER.name}' introuvable à côté du script.")
        print("   Structure attendue :")
        print("       projet/")
        print("       ├── merge_csv.py")
        print("       └── merge_csv/   ← doit exister")
        return

    csv_files = sorted([
        Path(f) for f in glob.glob(str(SOURCE_FOLDER / "*.csv"))
        if Path(f).resolve() != OUTPUT_FILE.resolve()
    ])

    if not csv_files:
        print(f"❌ Aucun fichier CSV trouvé dans : {SOURCE_FOLDER}")
        return

    print(f"📄 {len(csv_files)} fichier(s) CSV trouvé(s) :")
    for f in csv_files:
        cutoff_info = f"  [tronqué à {FILE_CUTOFFS[f.name]}]" if f.name in FILE_CUTOFFS else ""
        print(f"   • {f.name}{cutoff_info}")
    print()

    # ── Chargement + nettoyage individuel ────────────────────────────────────
    frames   = []
    ref_cols = None

    for path in csv_files:
        print(f"📥 Traitement : {path.name}")

        df = load_csv(path)
        if df is None:
            continue

        print(f"   Chargé : {len(df)} lignes")

        # 1. Nettoyage des bugs GPS
        df = remove_gps_bugs(df, n=GPS_BUG_N)
        print(f"   Après nettoyage GPS : {len(df)} lignes")

        # 2. Troncature temporelle (uniquement en mémoire, fichier source intact)
        df = apply_cutoff(df, path.name)
        print(f"   Lignes conservées pour la fusion : {len(df)}")

        frames.append(df)

        if ref_cols is None:
            ref_cols = list(df.columns)

        print()

    if not frames:
        print("❌ Aucun fichier valide chargé.")
        return

    # ── Fusion ────────────────────────────────────────────────────────────────
    print("🔀 Fusion en cours...")
    merged      = pd.concat(frames, ignore_index=True)
    total_avant = len(merged)

    # ── Tri par date ──────────────────────────────────────────────────────────
    time_col = merged.columns[1]
    merged   = merged.sort_values(by=time_col, na_position="last").reset_index(drop=True)

    # ── Dédoublonnage (toutes colonnes sauf id et Time) ───────────────────────
    cols_dedup = [c for c in merged.columns if c not in DEDUP_EXCLUDE]
    merged     = merged.drop_duplicates(subset=cols_dedup).reset_index(drop=True)

    doublons = total_avant - len(merged)
    if doublons > 0:
        print(f"🗑️  {doublons} ligne(s) dupliquée(s) supprimée(s)")
    else:
        print("✅ Aucun doublon détecté")

    # Remet les colonnes dans l'ordre du premier fichier
    cols_ordered = [c for c in ref_cols if c in merged.columns]
    cols_extra   = [c for c in merged.columns if c not in ref_cols]
    merged       = merged[cols_ordered + cols_extra]

    # ── Graphiques ────────────────────────────────────────────────────────────
    print("\n📊 Génération des graphiques...")
    make_plots(merged, PLOTS_FOLDER)

    # ── Écriture CSV ──────────────────────────────────────────────────────────
    print(f"\n💾 Écriture du fichier de sortie...")
    merged.to_csv(OUTPUT_FILE, index=False)

    print(f"""
─────────────────────────────────────────
✅ Terminé !
   Fichier CSV    : {OUTPUT_FILE}
   Lignes totales : {len(merged)}
   Colonnes       : {len(merged.columns)}
   Graphiques     : {PLOTS_FOLDER}
─────────────────────────────────────────
""")


if __name__ == "__main__":
    main()