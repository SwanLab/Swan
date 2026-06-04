"""
Détection d'anomalies GPS dans la trajectoire d'un navire
─────────────────────────────────────────────────────────
Principe :
  - Calcule la distance Haversine entre chaque paire de points consécutifs
  - En déduit la vitesse implicite (distance / intervalle de temps)
  - Signale les sauts anormaux (distance OU vitesse implicite au-delà du seuil)
  - Affiche un résumé + exporte les anomalies dans un CSV
  - Génère une carte Folium avec les segments normaux et les anomalies
"""

import pandas as pd
import numpy as np
import folium
from math import radians, sin, cos, sqrt, atan2
from pathlib import Path

# ──────────────────────────────────────────────────────────────
# CONFIGURATION
# ──────────────────────────────────────────────────────────────

CSV_PATH   = "merged_output.csv"

SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR    = SCRIPT_DIR / "data_plots"
OUT_DIR.mkdir(exist_ok=True)

OUTPUT_CSV = OUT_DIR / "anomalies.csv"
OUTPUT_MAP = OUT_DIR / "anomalies_map.html"

# Seuil 1 : saut de distance brut entre deux points consécutifs (mètres)
# Un navire ne peut pas téléporter de plus de ~500 m en 5 secondes
DIST_THRESHOLD_M = 15000

# Seuil 2 : vitesse implicite anormale (nœuds)
# Au-delà de 30 kn = quasi-impossible pour la plupart des navires commerciaux
SPEED_IMPL_MAX_KN = 3000.0

# ──────────────────────────────────────────────────────────────
# FONCTIONS UTILITAIRES
# ──────────────────────────────────────────────────────────────

def haversine_m(lat1, lon1, lat2, lon2) -> float:
    """Distance en mètres entre deux points GPS (formule de Haversine)."""
    R = 6_371_000  # rayon terrestre en mètres
    φ1, φ2 = radians(lat1), radians(lat2)
    Δφ = radians(lat2 - lat1)
    Δλ = radians(lon2 - lon1)
    a = sin(Δφ / 2)**2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)**2
    return 2 * R * atan2(sqrt(a), sqrt(1 - a))


def load_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip().lower().replace(" ", "_") for c in df.columns]

    # ── Temps ──────────────────────────────────────────────────
    time_col = next((c for c in df.columns if "time" in c or "date" in c), None)
    if time_col:
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
        print("⚠️  Aucune colonne temps trouvée — vitesse implicite non calculable.")

    # ── Vitesse ────────────────────────────────────────────────
    sog_col = next((c for c in df.columns if "speed" in c or "sog" in c), None)
    if sog_col:
        df["sog"] = pd.to_numeric(df[sog_col], errors="coerce").fillna(0.0)
        if "m/s" in sog_col.lower():
            df["sog"] *= 1.94384
    else:
        df["sog"] = np.nan

    df["latitude"]  = pd.to_numeric(df["latitude"],  errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    df = df.dropna(subset=["latitude", "longitude"]).reset_index(drop=True)
    print(f"  → {len(df)} points chargés")
    return df


# ──────────────────────────────────────────────────────────────
# DÉTECTION DES ANOMALIES
# ──────────────────────────────────────────────────────────────

def detect_anomalies(df: pd.DataFrame) -> pd.DataFrame:
    records = []

    for i in range(len(df) - 1):
        r0, r1 = df.iloc[i], df.iloc[i + 1]

        dist_m = haversine_m(r0.latitude, r0.longitude,
                              r1.latitude, r1.longitude)

        # Intervalle de temps
        dt_s = np.nan
        if pd.notna(r0.time_parsed) and pd.notna(r1.time_parsed):
            dt_s = (r1.time_parsed - r0.time_parsed).total_seconds()

        # Vitesse implicite (m/s → nœuds)
        speed_impl_kn = np.nan
        if pd.notna(dt_s) and dt_s > 0:
            speed_impl_kn = (dist_m / dt_s) * 1.94384

        # ΔLat / ΔLon bruts
        delta_lat = abs(r1.latitude  - r0.latitude)
        delta_lon = abs(r1.longitude - r0.longitude)

        is_anomaly = (
            dist_m > DIST_THRESHOLD_M or
            (pd.notna(speed_impl_kn) and speed_impl_kn > SPEED_IMPL_MAX_KN)
        )

        records.append({
            "idx_from"       : i,
            "idx_to"         : i + 1,
            "time_from"      : r0.time_parsed,
            "time_to"        : r1.time_parsed,
            "lat_from"       : r0.latitude,
            "lon_from"       : r0.longitude,
            "lat_to"         : r1.latitude,
            "lon_to"         : r1.longitude,
            "delta_lat"      : round(delta_lat, 6),
            "delta_lon"      : round(delta_lon, 6),
            "distance_m"     : round(dist_m, 1),
            "dt_seconds"     : round(dt_s, 2) if pd.notna(dt_s) else np.nan,
            "speed_impl_kn"  : round(speed_impl_kn, 1) if pd.notna(speed_impl_kn) else np.nan,
            "sog_kn"         : round(r0.sog, 2) if pd.notna(r0.sog) else np.nan,
            "is_anomaly"     : is_anomaly,
        })

    return pd.DataFrame(records)


# ──────────────────────────────────────────────────────────────
# RAPPORT CONSOLE
# ──────────────────────────────────────────────────────────────

def print_report(results: pd.DataFrame):
    anomalies = results[results["is_anomaly"]]
    print(f"\n{'═'*60}")
    print(f"  RÉSULTAT : {len(anomalies)} anomalie(s) sur {len(results)} transitions")
    print(f"  Seuil distance  : > {DIST_THRESHOLD_M} m")
    print(f"  Seuil vitesse   : > {SPEED_IMPL_MAX_KN} kn (implicite)")
    print(f"{'═'*60}\n")

    if anomalies.empty:
        print("  ✅ Aucune anomalie détectée.")
        return

    for _, row in anomalies.iterrows():
        t0 = row["time_from"].strftime("%Y-%m-%d %H:%M:%S") if pd.notna(row["time_from"]) else "?"
        t1 = row["time_to"].strftime("%Y-%m-%d %H:%M:%S")   if pd.notna(row["time_to"])   else "?"
        print(f"  ⚠️  Points #{int(row['idx_from'])} → #{int(row['idx_to'])}")
        print(f"     Heure        : {t0}  →  {t1}")
        print(f"     Position de  : ({row['lat_from']:.5f}, {row['lon_from']:.5f})")
        print(f"     Position à   : ({row['lat_to']:.5f},  {row['lon_to']:.5f})")
        print(f"     ΔLat={row['delta_lat']:.5f}°  ΔLon={row['delta_lon']:.5f}°")
        print(f"     Distance     : {row['distance_m']:.0f} m")
        if pd.notna(row["speed_impl_kn"]):
            print(f"     Vitesse impl.: {row['speed_impl_kn']:.1f} kn  "
                  f"(SOG enregistrée : {row['sog_kn']:.1f} kn)")
        print()


# ──────────────────────────────────────────────────────────────
# CARTE FOLIUM
# ──────────────────────────────────────────────────────────────

def build_anomaly_map(df: pd.DataFrame, results: pd.DataFrame) -> folium.Map:
    center_lat = df["latitude"].mean()
    center_lon = df["longitude"].mean()

    m = folium.Map(location=[center_lat, center_lon],
                   zoom_start=12, tiles="CartoDB positron", control_scale=True)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap").add_to(m)

    anomaly_idx = set(results[results["is_anomaly"]]["idx_from"].tolist())

    # Groupe : trajectoire normale
    grp_normal   = folium.FeatureGroup(name="Trajectoire normale",   show=True)
    # Groupe : segments anormaux
    grp_anomaly  = folium.FeatureGroup(name="⚠️ Anomalies GPS",       show=True)
    # Groupe : marqueurs d'anomalie
    grp_markers  = folium.FeatureGroup(name="Marqueurs anomalies",   show=True)

    coords = list(zip(df["latitude"], df["longitude"]))

    for i in range(len(coords) - 1):
        row = results.iloc[i]
        if i in anomaly_idx:
            folium.PolyLine(
                locations=[coords[i], coords[i + 1]],
                color="#FF0000", weight=5, opacity=0.9,
                tooltip=f"⚠️ Anomalie | {row['distance_m']:.0f} m | {row['speed_impl_kn']:.0f} kn impl.",
            ).add_to(grp_anomaly)
        else:
            folium.PolyLine(
                locations=[coords[i], coords[i + 1]],
                color="#2196F3", weight=3, opacity=0.6,
            ).add_to(grp_normal)

    # Marqueurs sur les points anormaux
    anomalies = results[results["is_anomaly"]]
    for _, row in anomalies.iterrows():
        t0 = row["time_from"].strftime("%H:%M:%S") if pd.notna(row["time_from"]) else "?"
        popup_html = (
            f"<b>⚠️ Anomalie GPS</b><br>"
            f"Heure : {t0}<br>"
            f"De : ({row['lat_from']:.5f}, {row['lon_from']:.5f})<br>"
            f"À  : ({row['lat_to']:.5f}, {row['lon_to']:.5f})<br>"
            f"Distance : <b>{row['distance_m']:.0f} m</b><br>"
            f"Vitesse impl. : <b>{row['speed_impl_kn']:.1f} kn</b>"
        )
        folium.CircleMarker(
            location=[row["lat_from"], row["lon_from"]],
            radius=8, color="#FF0000", fill=True, fill_color="#FF6600",
            fill_opacity=0.8,
            popup=folium.Popup(popup_html, max_width=250),
            tooltip=f"⚠️ {t0} | {row['distance_m']:.0f} m",
        ).add_to(grp_markers)

    grp_normal.add_to(m)
    grp_anomaly.add_to(m)
    grp_markers.add_to(m)
    folium.LayerControl().add_to(m)
    return m


# ──────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print(f"\n📂 Chargement de '{CSV_PATH}'...")
    df = load_data(CSV_PATH)

    print("🔍 Analyse des transitions...")
    results = detect_anomalies(df)

    print_report(results)

    # Export CSV des anomalies uniquement
    anomalies_df = results[results["is_anomaly"]].copy()
    anomalies_df.to_csv(OUTPUT_CSV, index=False)
    print(f"📄 Anomalies exportées → {OUTPUT_CSV}  ({len(anomalies_df)} ligne(s))")

    # Carte
    print("🗺️  Génération de la carte...")
    m = build_anomaly_map(df, results)
    m.save(OUTPUT_MAP)
    print(f"✅ Carte sauvegardée → {OUTPUT_MAP}\n")