"""
FEM Mesh Builder — crack detection images
==========================================
Expected folder structure on disk:

  crack_tutorial/
  └── archive/
      ├── fem_mesh_builder.py     <- THIS SCRIPT (place it here)
      ├── positive/               <- source images (227x227)
      │   ├── img001.jpg
      │   └── ...
      └── meshes/                 <- generated files (created automatically)
          ├── mesh_common.txt     <- geometry: nodes, neighbors, elements (one file for all images)
          └── intensities/
              ├── img001_intensity.txt   <- grayscale intensity per node
              ├── img002_intensity.txt
              └── ...

Output format (ASCII):
----------------------
mesh_common.txt
  Section [METADATA]   : image size, connectivity, units
  Section [NODES]      : node_id   x   y   row   col
  Section [NEIGHBORS]  : node_id   nb_1   nb_2   nb_3   nb_4   (-1 = no neighbor)
  Section [ELEMENTS]   : elem_id   n_top_left   n_top_right   n_bot_right   n_bot_left

  Node IDs ARE the mesh point IDs. The mesh starts at the center of the first pixel.
  Node 0 is at (x=0.5, y=0.5). Elements are defined by their 4 surrounding node IDs.
  For a pixel at (row r, col c):
    top-left     = node (r-1)*(W) + (c-1)    [or -1 if outside]
    top-right    = node (r-1)*(W) + c         [or -1 if outside]
    bottom-right = node  r   *(W) + c         [or -1 if outside]
    bottom-left  = node  r   *(W) + (c-1)     [or -1 if outside]

XXXXX_intensity.txt
  Section [METADATA]   : source image path
  Section [INTENSITY]  : node_id   intensity (float in [0.0, 1.0])
"""

import os
import sys
import glob
import numpy as np
from PIL import Image
from concurrent.futures import ProcessPoolExecutor, as_completed

IMAGE_EXTENSIONS = (".jpg", ".jpeg", ".png", ".bmp", ".tiff", ".tif")

# ── Settings ──────────────────────────────────────────────────────────────────
IMAGE_WIDTH  = 227      # expected image width  (pixels)
IMAGE_HEIGHT = 227      # expected image height (pixels)
CONNECTIVITY = 4        # neighbor connectivity: 4 (N/S/E/W) or 8 (+ diagonals)
NUM_WORKERS  = 4        # number of CPU cores used in parallel
LOG_EVERY    = 500      # print a progress line every N images
MAX_IMAGES   = 10     # max number of images to process (None = all available); for each categories
# ─────────────────────────────────────────────────────────────────────────────


def load_grayscale(image_path):
    """Load an image and convert it to grayscale float32 in [0, 1]."""
    img = Image.open(image_path).convert("L")
    return np.array(img, dtype=np.float32) / 255.0


def build_geometry(H, W, connectivity=4):
    """
    Build the mesh geometry for an H x W image.
    This is image-independent — called once for all images.

    Node numbering: origin at bottom-left (FEM convention), IDs start at 1.
    Node position : x = col + 0.5,  y = (H - row) - 0.5  (unit: pixel, y points up)

    The mesh uses node IDs directly as mesh point references.
    Elements are defined by the 4 surrounding nodes of each interior junction,
    ordered counter-clockwise from bottom-left (standard FEM convention).

    Parameters
    ----------
    H, W         : image height and width in pixels
    connectivity : 4 (S/N/W/E) or 8 (includes diagonals)

    Returns
    -------
    nodes     : dict {id, x, y, row, col}           each array of shape (N,)
    neighbors : (N, max_nb) int32 array, -1 = no neighbor
    elements  : (E, 4) int32 array of node IDs per quad element
                E = (H-1)*(W-1) interior junctions
                columns: bot-left, bot-right, top-right, top-left (counter-clockwise)
    """
    N = H * W

    # -- Node arrays -----------------------------------------------------------
    row_idx, col_idx = np.indices((H, W))
    row_idx  = row_idx.ravel().astype(np.int32)
    col_idx  = col_idx.ravel().astype(np.int32)
    # Node IDs start at 1, origin at bottom-left (FEM convention)
    # row=H-1 (bottom of image) → flipped_row=0 → node IDs 1..W
    # row=0   (top of image)    → flipped_row=H-1 → node IDs (H-1)*W+1..H*W
    flipped_row = (H - 1 - row_idx).astype(np.int32)
    node_ids    = (flipped_row * W + col_idx + 1).astype(np.int32)

    # Node position = center of pixel
    # x points right, y points up (mathematical convention)
    x = col_idx.astype(np.float32) + 0.5
    y = (H - row_idx).astype(np.float32) - 0.5

    nodes = {"id": node_ids, "x": x, "y": y, "row": row_idx, "col": col_idx}

    # -- Neighbors (fixed-width table: node_id | nb_1 | nb_2 | nb_3 | nb_4) --
    # Offsets in image space (row increases downward):
    # S = row+1 (y-1 in FEM), N = row-1 (y+1 in FEM), W = col-1, E = col+1
    offsets_4    = [(1, 0), (-1, 0), (0, -1), (0, 1)]    # S  N  W  E
    offsets_diag = [(1, -1), (1, 1), (-1, -1), (-1, 1)]  # SW SE NW NE
    offsets = offsets_4 if connectivity == 4 else offsets_4 + offsets_diag
    max_nb  = len(offsets)

    neighbors = np.full((N, max_nb), -1, dtype=np.int32)
    for k, (dr, dc) in enumerate(offsets):
        new_r = row_idx + dr
        new_c = col_idx + dc
        valid = (new_r >= 0) & (new_r < H) & (new_c >= 0) & (new_c < W)
        # Convert neighbor (new_r, new_c) to its flipped node ID
        flipped_new_r = (H - 1 - new_r[valid]).astype(np.int32)
        neighbors[valid, k] = (flipped_new_r * W + new_c[valid] + 1).astype(np.int32)

    # -- Quad elements ---------------------------------------------------------
    # One element per interior junction of the node grid.
    # Junction (r, c) with r in [1..H-1], c in [1..W-1]:
    jr = np.arange(H-1, 0, -1, dtype=np.int32)   # H-1 .. 1 : bottom to top in FEM
    jc = np.arange(1, W, dtype=np.int32)
    JR, JC = np.meshgrid(jr, jc, indexing="ij")
    JR = JR.ravel()
    JC = JC.ravel()

    # Convert image rows to flipped rows for FEM node IDs
    # JR is the image row of the junction's bottom pixel
    # JR-1 is the image row of the junction's top pixel
    f_bot = (H - 1 - JR    ).astype(np.int32)   # flipped row of bottom nodes
    f_top = (H - 1 - (JR-1)).astype(np.int32)   # flipped row of top nodes

    # Element nodes ordered counter-clockwise from bottom-left (FEM standard)
    bot_left  = (f_bot * W + (JC - 1) + 1).astype(np.int32)
    bot_right = (f_bot * W +  JC      + 1).astype(np.int32)
    top_right = (f_top * W +  JC      + 1).astype(np.int32)
    top_left  = (f_top * W + (JC - 1) + 1).astype(np.int32)

    elements = np.stack([bot_left, bot_right, top_right, top_left], axis=1)

    return nodes, neighbors, elements


def save_common(path, H, W, connectivity, nodes, neighbors, elements):
    """
    Write mesh_common.txt — geometry shared by all images.
    ASCII format with clearly labelled sections.
    """
    N = H * W
    E = elements.shape[0]
    max_nb = neighbors.shape[1]

    with open(path, "w") as f:

        # -- METADATA ----------------------------------------------------------
        f.write("[METADATA]\n")
        f.write(f"image_height   {H}\n")
        f.write(f"image_width    {W}\n")
        f.write(f"num_nodes      {N}\n")
        f.write(f"num_elements   {E}\n")
        f.write(f"connectivity   {connectivity}\n")
        f.write(f"unit           pixel\n")
        f.write(f"note           Node IDs start at 1, origin at bottom-left.\n")
        f.write(f"note           x=col+0.5, y=(H-row)-0.5. y points up.\n")
        f.write(f"note           Elements: counter-clockwise from bot-left.\n")
        f.write(f"note           Neighbor value -1 means no neighbor (border node).\n")
        f.write("\n")

        # -- NODES -------------------------------------------------------------
        f.write("[NODES]\n")
        f.write("# node_id   x        y        row   col\n")
        sort_idx = np.argsort(nodes['id'])   # tri par ID croissant (1, 2, 3, ...)
        for i in sort_idx:
            f.write(f"{nodes['id'][i]:<10d}"
                    f"{nodes['x'][i]:<10.1f}"
                    f"{nodes['y'][i]:<10.1f}"
                    f"{nodes['row'][i]:<6d}"
                    f"{nodes['col'][i]:<6d}\n")
        f.write("\n")

        # -- NEIGHBORS ---------------------------------------------------------
        nb_header = "   ".join([f"nb_{k+1}" for k in range(max_nb)])
        f.write("[NEIGHBORS]\n")
        f.write(f"# node_id   {nb_header}   (-1 = no neighbor)\n")
        for i in sort_idx:   # même ordre trié que les nœuds
            nb_vals = "   ".join([f"{neighbors[i, k]:<7d}" for k in range(max_nb)])
            f.write(f"{nodes['id'][i]:<10d}{nb_vals}\n")
        f.write("\n")

        # -- ELEMENTS ----------------------------------------------------------
        f.write("[ELEMENTS]\n")
        f.write("# elem_id   bot_left   bot_right   top_right   top_left\n")
        for i in range(E):
            e = elements[i]
            f.write(f"{i+1:<10d}{e[0]:<12d}{e[1]:<12d}{e[2]:<12d}{e[3]:<12d}\n")

    print(f"  -> mesh_common.txt written  ({N:,} nodes, {E:,} elements)")


def save_intensity(path, image_path, gray):
    """
    Write one XXXXX_intensity.txt file for a single image.
    Contains only node_id and grayscale intensity.
    """
    H, W = gray.shape
    N    = H * W

    with open(path, "w") as f:
        f.write("[METADATA]\n")
        f.write(f"source_image   {os.path.basename(image_path)}\n")
        f.write(f"num_nodes      {N}\n")
        f.write(f"intensity_range  [0.0, 1.0]\n")
        f.write("\n")
        f.write("[INTENSITY]\n")
        f.write("# node_id   intensity\n")
        intensity = gray.ravel()
        for i in range(N):
            f.write(f"{i+1:<10d}{intensity[i]:.6f}\n")


def process_one(args):
    """
    Process a single image (intensity only — geometry is already written).
    Runs inside a worker sub-process.
    Returns (image_path, True) on success, (image_path, error_message) on failure.
    """
    image_path, output_dir = args
    try:
        gray     = load_grayscale(image_path)
        basename = os.path.splitext(os.path.basename(image_path))[0]
        out_path = os.path.join(output_dir, f"{basename}_intensity.txt")
        save_intensity(out_path, image_path, gray)
        return (image_path, True)
    except Exception as e:
        return (image_path, str(e))


if __name__ == "__main__":

    # -- Folder paths ----------------------------------------------------------
    SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
    OUTPUT_DIR  = os.path.join(SCRIPT_DIR, "meshes")
    COMMON_PATH = os.path.join(OUTPUT_DIR, "mesh_common.txt")

    # MAX_IMAGES applies independently to each class:
    # e.g. MAX_IMAGES=10 → 10 positive + 10 negative = 20 total
    SOURCES = {
        "positive": os.path.join(SCRIPT_DIR, "Positive"),
        "negative": os.path.join(SCRIPT_DIR, "Negative"),
    }

    # -- Check source folders --------------------------------------------------
    for label, folder in SOURCES.items():
        if not os.path.isdir(folder):
            print(f"ERROR: source folder not found: {folder}")
            print(f"Make sure this script is placed inside the 'archive/' folder.")
            sys.exit(1)

    # -- Collect images per class ----------------------------------------------
    tasks_per_class = {}
    for label, folder in SOURCES.items():
        all_paths = sorted([
            p for p in glob.glob(os.path.join(folder, "*"))
            if os.path.splitext(p)[1].lower() in IMAGE_EXTENSIONS
        ])
        if not all_paths:
            print(f"ERROR: no images found in {folder}")
            sys.exit(1)
        selected = all_paths[:MAX_IMAGES] if MAX_IMAGES is not None else all_paths
        tasks_per_class[label] = (all_paths, selected)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # -- Print summary ---------------------------------------------------------
    print("=" * 55)
    for label, (all_paths, selected) in tasks_per_class.items():
        skipped = len(all_paths) - len(selected)
        print(f"  [{label:>8}]  available: {len(all_paths):>5,}  "
              f"to process: {len(selected):>5,}"
              + (f"  ({skipped:,} skipped)" if skipped else ""))
    total = sum(len(s) for _, s in tasks_per_class.values())
    print(f"  {'TOTAL':>8}   {total:,} images")
    print(f"  Connectivity : {CONNECTIVITY}-connected")
    print(f"  Workers      : {NUM_WORKERS} CPU cores")
    print(f"  Output       : {OUTPUT_DIR}")
    print("=" * 55)

    # -- Write mesh_common.txt (geometry, once for all images) -----------------
    print("\n  Building shared geometry...")
    nodes, neighbors, elements = build_geometry(IMAGE_HEIGHT, IMAGE_WIDTH, CONNECTIVITY)
    save_common(COMMON_PATH, IMAGE_HEIGHT, IMAGE_WIDTH, CONNECTIVITY,
                nodes, neighbors, elements)

    # -- Write intensity files per class (parallel) ----------------------------
    all_tasks = []
    for label, (_, selected) in tasks_per_class.items():
        intensity_dir = os.path.join(OUTPUT_DIR, "intensities", label)
        os.makedirs(intensity_dir, exist_ok=True)
        for p in selected:
            all_tasks.append((p, intensity_dir))

    print(f"\n  Writing {len(all_tasks):,} intensity files ...")
    ok, errors = 0, []
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as pool:
        futures = {pool.submit(process_one, a): a[0] for a in all_tasks}
        for future in as_completed(futures):
            image_path, result = future.result()
            if result is True:
                ok += 1
            else:
                errors.append((image_path, result))
            done = ok + len(errors)
            if done % LOG_EVERY == 0 or done == len(all_tasks):
                pct = done / len(all_tasks) * 100
                print(f"  {done:>6}/{len(all_tasks)}  ({pct:5.1f}%)  "
                      f"OK={ok}  errors={len(errors)}")

    # -- Final summary ---------------------------------------------------------
    print("=" * 55)
    print(f"  Done: {ok}/{len(all_tasks)} intensity files written")
    if errors:
        print(f"\n  Failures ({len(errors)}):")
        for p, msg in errors:
            print(f"    {os.path.basename(p)} -> {msg}")
    print("=" * 55)