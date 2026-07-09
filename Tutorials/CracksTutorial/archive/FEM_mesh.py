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
MAX_IMAGES   = 100   # max number of images to process (None = all available)
# ─────────────────────────────────────────────────────────────────────────────


def load_grayscale(image_path):
    """Load an image and convert it to grayscale float32 in [0, 1]."""
    img = Image.open(image_path).convert("L")
    return np.array(img, dtype=np.float32) / 255.0


def build_geometry(H, W, connectivity=4):
    """
    Build the mesh geometry for an H x W image.
    This is image-independent — called once for all images.

    Node numbering: node_id = row * W + col
    Node position : x = col + 0.5,  y = row + 0.5  (unit: pixel)

    The mesh uses node IDs directly as mesh point references.
    Elements are defined by the 4 surrounding nodes of each interior junction.
    A junction at grid point (r, c) — where r in [1..H-1], c in [1..W-1] —
    is surrounded by nodes:
        top-left  = (r-1)*W + (c-1)
        top-right = (r-1)*W +  c
        bot-right =  r   *W +  c
        bot-left  =  r   *W + (c-1)

    Parameters
    ----------
    H, W         : image height and width in pixels
    connectivity : 4 (N/S/E/W) or 8 (includes diagonals)

    Returns
    -------
    nodes     : dict {id, x, y, row, col}           each array of shape (N,)
    neighbors : (N, max_nb) int32 array, -1 = no neighbor
    elements  : (E, 4) int32 array of node IDs per quad element
                E = (H-1)*(W-1) interior junctions
                columns: top-left, top-right, bottom-right, bottom-left
    """
    N = H * W

    # -- Node arrays -----------------------------------------------------------
    row_idx, col_idx = np.indices((H, W))
    row_idx  = row_idx.ravel().astype(np.int32)
    col_idx  = col_idx.ravel().astype(np.int32)
    node_ids = np.arange(N, dtype=np.int32)

    # Node position = center of pixel
    x = col_idx.astype(np.float32) + 0.5
    y = row_idx.astype(np.float32) + 0.5

    nodes = {"id": node_ids, "x": x, "y": y, "row": row_idx, "col": col_idx}

    # -- Neighbors (fixed-width table: node_id | nb_1 | nb_2 | nb_3 | nb_4) --
    offsets_4    = [(-1, 0), (1, 0), (0, -1), (0, 1)]   # N  S  W  E
    offsets_diag = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    offsets = offsets_4 if connectivity == 4 else offsets_4 + offsets_diag
    max_nb  = len(offsets)

    neighbors = np.full((N, max_nb), -1, dtype=np.int32)
    for k, (dr, dc) in enumerate(offsets):
        new_r = row_idx + dr
        new_c = col_idx + dc
        valid = (new_r >= 0) & (new_r < H) & (new_c >= 0) & (new_c < W)
        neighbors[valid, k] = (new_r[valid] * W + new_c[valid]).astype(np.int32)

    # -- Quad elements ---------------------------------------------------------
    # One element per interior junction of the node grid.
    # Junction (r, c) with r in [1..H-1], c in [1..W-1]:
    jr = np.arange(1, H, dtype=np.int32)
    jc = np.arange(1, W, dtype=np.int32)
    JR, JC = np.meshgrid(jr, jc, indexing="ij")
    JR = JR.ravel()
    JC = JC.ravel()

    top_left  = ((JR - 1) * W + (JC - 1)).astype(np.int32)
    top_right = ((JR - 1) * W +  JC     ).astype(np.int32)
    bot_right = ( JR      * W +  JC     ).astype(np.int32)
    bot_left  = ( JR      * W + (JC - 1)).astype(np.int32)

    elements = np.stack([top_left, top_right, bot_right, bot_left], axis=1)

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
        f.write(f"note           Node_id=row*W+col. Position: x=col+0.5, y=row+0.5.\n")
        f.write(f"note           Elements defined by 4 surrounding node IDs.\n")
        f.write(f"note           Neighbor value -1 means no neighbor (border node).\n")
        f.write("\n")

        # -- NODES -------------------------------------------------------------
        f.write("[NODES]\n")
        f.write("# node_id   x        y        row   col\n")
        for i in range(N):
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
        for i in range(N):
            nb_vals = "   ".join([f"{neighbors[i, k]:<7d}" for k in range(max_nb)])
            f.write(f"{i:<10d}{nb_vals}\n")
        f.write("\n")

        # -- ELEMENTS ----------------------------------------------------------
        f.write("[ELEMENTS]\n")
        f.write("# elem_id   top_left   top_right   bot_right   bot_left\n")
        for i in range(E):
            e = elements[i]
            f.write(f"{i:<10d}{e[0]:<12d}{e[1]:<12d}{e[2]:<12d}{e[3]:<12d}\n")

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
            f.write(f"{i:<10d}{intensity[i]:.6f}\n")


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
    SCRIPT_DIR    = os.path.dirname(os.path.abspath(__file__))
    INPUT_DIR     = os.path.join(SCRIPT_DIR, "Positive")
    OUTPUT_DIR    = os.path.join(SCRIPT_DIR, "meshes")
    INTENSITY_DIR = os.path.join(OUTPUT_DIR, "intensities")
    COMMON_PATH   = os.path.join(OUTPUT_DIR, "mesh_common.txt")

    # -- Check source folder ---------------------------------------------------
    if not os.path.isdir(INPUT_DIR):
        print(f"ERROR: source folder not found: {INPUT_DIR}")
        print(f"Make sure this script is placed inside the 'archive/' folder.")
        sys.exit(1)

    # -- Collect images --------------------------------------------------------
    all_image_paths = sorted([
        p for p in glob.glob(os.path.join(INPUT_DIR, "*"))
        if os.path.splitext(p)[1].lower() in IMAGE_EXTENSIONS
    ])

    if not all_image_paths:
        print(f"ERROR: no images found in {INPUT_DIR}")
        sys.exit(1)

    if MAX_IMAGES is not None:
        image_paths = all_image_paths[:MAX_IMAGES]
        skipped     = len(all_image_paths) - len(image_paths)
    else:
        image_paths = all_image_paths
        skipped     = 0

    os.makedirs(OUTPUT_DIR,    exist_ok=True)
    os.makedirs(INTENSITY_DIR, exist_ok=True)

    # -- Print summary ---------------------------------------------------------
    print("=" * 55)
    print(f"  Available images : {len(all_image_paths):,}")
    print(f"  Images to process: {len(image_paths):,}"
          + (f"  ({skipped:,} skipped)" if skipped else ""))
    print(f"  Connectivity     : {CONNECTIVITY}-connected")
    print(f"  Output format    : ASCII")
    print(f"  Workers          : {NUM_WORKERS} CPU cores")
    print(f"  Output folder    : {OUTPUT_DIR}")
    print("=" * 55)

    # -- Write mesh_common.txt (geometry, once for all images) -----------------
    print("\n  Building shared geometry...")
    nodes, neighbors, elements = build_geometry(IMAGE_HEIGHT, IMAGE_WIDTH, CONNECTIVITY)
    save_common(COMMON_PATH, IMAGE_HEIGHT, IMAGE_WIDTH, CONNECTIVITY,
                nodes, neighbors, elements)

    # -- Write one intensity file per image (parallel) -------------------------
    print(f"\n  Writing intensity files to {INTENSITY_DIR} ...")
    task_args = [(p, INTENSITY_DIR) for p in image_paths]

    ok, errors = 0, []
    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as pool:
        futures = {pool.submit(process_one, a): a[0] for a in task_args}
        for future in as_completed(futures):
            image_path, result = future.result()
            if result is True:
                ok += 1
            else:
                errors.append((image_path, result))
            done = ok + len(errors)
            if done % LOG_EVERY == 0 or done == len(image_paths):
                pct = done / len(image_paths) * 100
                print(f"  {done:>6}/{len(image_paths)}  ({pct:5.1f}%)  "
                      f"OK={ok}  errors={len(errors)}")

    # -- Final summary ---------------------------------------------------------
    print("=" * 55)
    print(f"  Done: {ok}/{len(image_paths)} intensity files written")
    if errors:
        print(f"\n  Failures ({len(errors)}):")
        for p, msg in errors:
            print(f"    {os.path.basename(p)} -> {msg}")
    print("=" * 55)