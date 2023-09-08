import shutil
import subprocess
import zipfile
import binascii
import sys
import os

from pathlib import Path

IMAGEJ='/home/dleonpe/Fiji.app/ImageJ-linux64'
MACRO='/data/rajewsky/home/dleonpe/projects/image_restoration/scripts/image_stitching/auto_stitch_macro.ijm'
TMPDIR='/scratch/dleonpe/stitching'

INDIR = sys.argv[1]
OUTDIR = sys.argv[2]

def copytree2(source,dest):
    Path(dest).mkdir(parents=True, exist_ok=True)
    dest_dir = os.path.join(dest,os.path.basename(source))
    if os.path.exists(dest_dir) and os.path.getsize(dest_dir) == os.path.getsize(source):
        print("The directory {OUTFILE} was already copied. Skipping!")
    else:
        shutil.copytree(source,dest_dir,dirs_exist_ok=True)
    return dest_dir

# Getting grid from binary file
archive = zipfile.ZipFile(os.path.join(INDIR, "Image.bcf"), 'r')
imgdata = archive.read('GroupFileProperty/Marker/StackList')
GRIDX = int(binascii.hexlify(imgdata)[16:18],16)
GRIDY = int(binascii.hexlify(imgdata)[24:26],16)

# Configuring file names and directories
OUTFILE = os.path.join(OUTDIR, "Image_Stitched_Composite.tif")
Path(OUTDIR).mkdir(parents=True, exist_ok=True)

if os.path.exists(OUTFILE) and os.path.getsize(OUTFILE) > 0:
    raise FileExistsError("The file {OUTFILE} was already generated. Skipping!")

print(f"Copying files from {INDIR} to {TMPDIR}")
INDIR_TMP = copytree2(INDIR, TMPDIR)

print(f"Grid size: X = {GRIDX}, Y = {GRIDY}")
print(f"Running ImageJ macro {MACRO}")

rc = subprocess.run([f"{IMAGEJ}", "--headless", "--console", f"-macro", f"{MACRO}", f"{GRIDX};{GRIDY};{INDIR_TMP};{OUTFILE}"])
