import os
import glob
import re

code_dir = r"c:\Users\ssanz\OneDrive\Documentos\GitHub\Swan\OPTIMIZATION_CODES"
files = glob.glob(os.path.join(code_dir, "*.m"))

# The block to add for the T300 material
t300_block = """
            db.T300_914_C.type = 'ORTHOTROPIC';
            db.T300_914_C.E  = [138.0, 11.0, 11.0] * 1e9;
            db.T300_914_C.nu = [0.28, 0.28, 0.40];
            db.T300_914_C.G  = [5.5, 5.5, 3.928] * 1e9;
            db.T300_914_C.density = 1580;
"""

t300_struct = "            db.T300_914_C = struct('type','ORTHOTROPIC','E',[138.0,11.0,11.0]*1e9, 'nu',[0.28,0.28,0.40],'G',[5.5,5.5,3.928]*1e9,'density',1580);\n"

for fpath in files:
    with open(fpath, "r", encoding="utf-8") as f:
        content = f.read()
    
    modified = False

    # 1. Update obj.materialLayers
    if "'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'" in content:
        content = content.replace("'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT'", "'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'")
        modified = True
    
    # Update compositeLayers in Al7075
    if "compositeLayers = {'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT';" in content:
        content = content.replace("compositeLayers = {'EpT'; 'EpT'; 'EpT'; 'EpT'; 'EpT';", "compositeLayers = {'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C'; 'T300_914_C';")
        modified = True

    # 2. Inject db definition
    if "db.EpT.density = 1600;" in content and "db.T300_914_C" not in content:
        content = content.replace("db.EpT.density = 1600;", "db.EpT.density = 1600;" + t300_block)
        modified = True
    elif "db.EpT  = struct(" in content and "db.T300_914_C" not in content:
        content = re.sub(r"(db\.EpT\s*=\s*struct[^\n]+)", r"\1\n" + t300_struct, content)
        modified = True

    if modified:
        with open(fpath, "w", encoding="utf-8") as f:
            f.write(content)
        print(f"Updated {os.path.basename(fpath)}")
