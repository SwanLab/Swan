"""
Tracé de trajectoire de navire avec coloration par vitesse
Utilise Folium pour la carte interactive + Matplotlib pour les graphes

Sorties dans ./data_plots/ :
  - ship_trajectory.html   : carte interactive
  - anomalies_map.html     : carte des anomalies (si detect_anomalies importé)
  - speed_over_time.png    : SpeedOverGround (m/s) en fonction du temps
  - course_over_time.png   : CourseOverGround (deg) en fonction du temps

Format CSV attendu :
  - Colonne index (ignorée)
  - time                    : "2025-05-15 18:16:16.203767" ou "15/05/2025  18:16:16"
  - latitude                : float
  - longitude               : float
  - speedoverground_(m/s)   : float
  - courseoverground_(deg)  : float  (optionnel)
"""

import pandas as pd
import folium
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path


# ──────────────────────────────────────────────
# 1. CONFIGURATION
# ──────────────────────────────────────────────

CSV_PATH   = "merged_output.csv"

# Dossier de sortie — créé automatiquement à côté du script
SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR    = SCRIPT_DIR / "data_plots"
OUT_DIR.mkdir(exist_ok=True)

OUTPUT_MAP = OUT_DIR / "ship_trajectory.html"


# ──────────────────────────────────────────────
# 2. LECTURE DU CSV
# ──────────────────────────────────────────────

def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # ── Vitesse SOG ───────────────────────────────
    sog_candidates = [c for c in df.columns if "speed" in c or "sog" in c]
    if not sog_candidates:
        raise ValueError(f"Aucune colonne vitesse trouvée. Colonnes : {list(df.columns)}")
    sog_col = sog_candidates[0]
    print(f"  → Colonne vitesse    : '{sog_col}'")

    # ── Cap COG ───────────────────────────────────
    cog_candidates = [c for c in df.columns if "course" in c or "cog" in c or "heading" in c]
    cog_col = cog_candidates[0] if cog_candidates else None
    if cog_col:
        print(f"  → Colonne cap        : '{cog_col}'")
    else:
        print("  → Colonne cap        : non trouvée (graphe COG ignoré)")

    # ── Temps ─────────────────────────────────────
    time_candidates = [c for c in df.columns if "time" in c or "date" in c]
    if time_candidates:
        time_col = time_candidates[0]
        parsed = pd.to_datetime(df[time_col].str.strip(),
                                format="%d/%m/%Y  %H:%M:%S", errors="coerce")
        mask = parsed.isna()
        if mask.any():
            parsed2 = pd.to_datetime(df.loc[mask, time_col].str.strip(), errors="coerce")
            parsed  = parsed.astype("datetime64[us]")
            parsed2 = parsed2.astype("datetime64[us]")
            parsed[mask] = parsed2
        df["time_parsed"] = parsed
    else:
        df["time_parsed"] = pd.NaT

    # ── Renommage & conversions ───────────────────
    df = df.rename(columns={sog_col: "sog"})
    df["sog"] = pd.to_numeric(df["sog"], errors="coerce").fillna(0.0)
    df["sog_ms"] = df["sog"].copy()          # garde les m/s bruts pour le graphe

    if "m/s" in sog_col.lower():
        df["sog"] = df["sog"] * 1.94384      # → nœuds pour la carte
        print("  → Vitesse convertie  : m/s → nœuds pour la carte")

    if cog_col:
        df = df.rename(columns={cog_col: "cog"})
        df["cog"] = pd.to_numeric(df["cog"], errors="coerce")

    df["latitude"]  = pd.to_numeric(df["latitude"],  errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    df = df.dropna(subset=["latitude", "longitude"]).reset_index(drop=True)

    print(f"  → {len(df)} points chargés")
    print(f"  → SOG : min={df['sog_ms'].min():.2f} m/s, "
          f"max={df['sog_ms'].max():.2f} m/s, "
          f"moy={df['sog_ms'].mean():.2f} m/s")
    return df


# ──────────────────────────────────────────────
# 3. GRAPHES MATPLOTLIB
# ──────────────────────────────────────────────

def _setup_time_axis(ax, df: pd.DataFrame):
    """Configure l'axe X commun en datetime."""
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")

    # Affiche la date dans le label si la trajectoire couvre plusieurs jours
    t_min = df["time_parsed"].min()
    t_max = df["time_parsed"].max()
    if t_max.date() != t_min.date():
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m %H:%M"))

    ax.set_xlabel("Temps (HH:MM)", fontsize=11)


def plot_speed(df: pd.DataFrame, out_dir: Path):
    """SpeedOverGround (m/s) en fonction du temps."""
    fig, ax = plt.subplots(figsize=(14, 4))

    times  = df["time_parsed"]
    speeds = df["sog_ms"]

    ax.plot(times, speeds, color="#1a78c2", linewidth=0.8, alpha=0.85, label="SOG (m/s)")
    ax.fill_between(times, speeds, alpha=0.15, color="#1a78c2")

    # Ligne de vitesse moyenne
    mean_v = speeds.mean()
    ax.axhline(mean_v, color="#d32f2f", linewidth=1.2, linestyle="--",
               label=f"Moyenne : {mean_v:.2f} m/s")

    ax.set_title("Speed Over Ground (m/s) en fonction du temps", fontsize=13, fontweight="bold")
    ax.set_ylabel("Vitesse (m/s)", fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    _setup_time_axis(ax, df)

    fig.tight_layout()
    out_path = out_dir / "speed_over_time.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → Graphe vitesse     : {out_path}")


def plot_course(df: pd.DataFrame, out_dir: Path):
    """CourseOverGround (deg) en fonction du temps."""
    if "cog" not in df.columns or df["cog"].isna().all():
        print("  → Graphe cap         : ignoré (colonne COG absente ou vide)")
        return

    fig, ax = plt.subplots(figsize=(14, 4))

    times  = df["time_parsed"]
    course = df["cog"]

    ax.scatter(times, course, color="#e67e22", s=1.5, alpha=0.6, label="COG (°)")

    ax.set_title("Course Over Ground (degrés) en fonction du temps", fontsize=13, fontweight="bold")
    ax.set_ylabel("Cap (°)", fontsize=11)
    ax.set_ylim(-5, 365)
    ax.set_yticks([0, 45, 90, 135, 180, 225, 270, 315, 360])
    ax.yaxis.set_tick_params(labelsize=9)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    _setup_time_axis(ax, df)

    fig.tight_layout()
    out_path = out_dir / "course_over_time.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  → Graphe cap         : {out_path}")


# ──────────────────────────────────────────────
# 4. COLORISATION PAR VITESSE
# ──────────────────────────────────────────────

def speed_to_color(speed: float, vmin: float, vmax: float) -> str:
    if vmax <= vmin:
        return "#1a78c2"
    t = max(0.0, min(1.0, (speed - vmin) / (vmax - vmin)))
    palette = [
        (0.00, (26,  120, 194)),
        (0.30, (0,   180, 180)),
        (0.60, (240, 200,   0)),
        (1.00, (211,  47,  47)),
    ]
    for i in range(len(palette) - 1):
        t0, c0 = palette[i]
        t1, c1 = palette[i + 1]
        if t0 <= t <= t1:
            alpha = (t - t0) / (t1 - t0)
            r = int(c0[0] + alpha * (c1[0] - c0[0]))
            g = int(c0[1] + alpha * (c1[1] - c0[1]))
            b = int(c0[2] + alpha * (c1[2] - c0[2]))
            return f"#{r:02x}{g:02x}{b:02x}"
    return "#d32f2f"


# ──────────────────────────────────────────────
# 5. CARTE FOLIUM
# ──────────────────────────────────────────────

def build_map(df: pd.DataFrame) -> folium.Map:
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()

    m = folium.Map(location=[center_lat, center_lon], zoom_start=13,
                   tiles="CartoDB positron", control_scale=True)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)
    folium.TileLayer(
        tiles="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
        attr="Esri", name="Satellite (Esri)",
    ).add_to(m)

    vmin, vmax = df["sog"].min(), df["sog"].max()
    coords = list(zip(df["latitude"], df["longitude"]))
    speeds = df["sog"].tolist()

    for i in range(len(coords) - 1):
        seg_speed = (speeds[i] + speeds[i + 1]) / 2.0
        color = speed_to_color(seg_speed, vmin, vmax)

        if "time_parsed" in df.columns and not pd.isna(df.iloc[i]["time_parsed"]):
            t_str = df.iloc[i]["time_parsed"].strftime("%H:%M:%S")
        else:
            t_str = f"Point {i}"

        folium.PolyLine(
            locations=[coords[i], coords[i + 1]],
            color=color, weight=4, opacity=0.85,
            tooltip=f"⏱ {t_str} | 🚢 {seg_speed:.1f} kn",
        ).add_to(m)

    folium.Marker(
        location=coords[0],
        popup=folium.Popup(f"<b>Départ</b><br>Lat: {coords[0][0]:.5f}<br>Lon: {coords[0][1]:.5f}", max_width=200),
        icon=folium.Icon(color="green", icon="play", prefix="fa"),
    ).add_to(m)
    folium.Marker(
        location=coords[-1],
        popup=folium.Popup(f"<b>Arrivée</b><br>Lat: {coords[-1][0]:.5f}<br>Lon: {coords[-1][1]:.5f}", max_width=200),
        icon=folium.Icon(color="red", icon="flag", prefix="fa"),
    ).add_to(m)

    legend_html = f"""
    <div style="position:fixed;bottom:40px;right:20px;z-index:9999;background:white;
                border:2px solid #ccc;border-radius:8px;padding:12px 16px;
                font-family:Arial,sans-serif;font-size:13px;
                box-shadow:2px 2px 8px rgba(0,0,0,0.3);min-width:160px;">
        <b>Vitesse (nœuds)</b><br><br>
        <div style="height:16px;background:linear-gradient(to right,#1a78c2,#00b4b4,#f0c800,#d32f2f);
                    border-radius:4px;margin-bottom:4px;"></div>
        <div style="display:flex;justify-content:space-between;">
            <span>{vmin:.1f} kn</span><span>{(vmin+vmax)/2:.1f} kn</span><span>{vmax:.1f} kn</span>
        </div>
        <div style="display:flex;justify-content:space-between;color:#666;font-size:11px;margin-top:2px;">
            <span>🔵 Lent</span><span>🔴 Rapide</span>
        </div>
        <hr style="margin:8px 0">
        <span style="color:green;">▶ Départ</span>&nbsp;&nbsp;
        <span style="color:red;">⚑ Arrivée</span>
    </div>"""
    m.get_root().html.add_child(folium.Element(legend_html))
    folium.LayerControl().add_to(m)
    return m


# ──────────────────────────────────────────────
# 6. MAIN
# ──────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n📂 Chargement de '{CSV_PATH}'...")
    df = load_data(CSV_PATH)

    print(f"\n📊 Génération des graphes → {OUT_DIR}/")
    plot_speed(df, OUT_DIR)
    plot_course(df, OUT_DIR)

    print(f"\n🗺️  Construction de la carte...")
    m = build_map(df)
    m.save(OUTPUT_MAP)
    print(f"  → Carte              : {OUTPUT_MAP}")

    print(f"\n✅ Tout sauvegardé dans : {OUT_DIR}\n")