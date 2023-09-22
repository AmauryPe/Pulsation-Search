import os
import argparse
import shutil
import fnmatch
from astropy.io import fits
import requests
import pandas as pd
import wget


"""
Move the interesting candidates from the processed data directory into the candidate directory
"""

parser = argparse.ArgumentParser(description='Move an interesting candidate into the candidate directory')
parser.add_argument("--obsID", help="Observation ID",type=str)
parser.add_argument("--src", help="Interesting sources (in hexadecimal), must separate them by a space",nargs='+', default = [])
args = parser.parse_args()
src_list=args.src
src_list=['8{:0>3}'.format(item) for item in src_list]
obsID=args.obsID

data_dir=f"/mnt/xmmcat/4XMM_data/DR13_incr_data/{obsID}"
process_dir=f"/mnt/data/Amaury/Processed_data/{obsID}"
outdir=f"/mnt/data/Amaury/candidates/{obsID}"

if not os.path.exists(outdir):
    os.mkdir(outdir)

image = [file for file in os.listdir(data_dir+'/product') if fnmatch.fnmatch(file, '*PN*IMAGE_8000.FTZ')][0]
pievli = [file for file in os.listdir(data_dir+'/product') if fnmatch.fnmatch(file, '*PN*PIEVLI0000.FTZ')][0]

if os.path.exists(data_dir+'/product/'+image):
    shutil.copy(data_dir+'/product/'+image,outdir)
if os.path.exists(data_dir+'/product/'+pievli):
    shutil.copy(data_dir+'/product/'+pievli,outdir)

files = os.listdir(process_dir)

for src in src_list:
    matching_files=[file for file in files if src in file]
    if not os.path.exists(outdir+'/'+src):
        os.mkdir(outdir+'/'+src)
    for file in matching_files :
        shutil.copy(process_dir+'/'+file,outdir+'/'+src)
    src_spec = [file for file in os.listdir(data_dir+'/product') if fnmatch.fnmatch(file, '*PN*SRSPEC*'+src[1:]+'.FTZ')][0]
    bkg_spec = [file for file in os.listdir(data_dir+'/product') if fnmatch.fnmatch(file, '*PN*BGSPEC*'+src[1:]+'.FTZ')][0]
    arf_mat = [file for file in os.listdir(data_dir+'/product') if fnmatch.fnmatch(file, '*PN*SRCARF*'+src[1:]+'.FTZ')][0]
    
    shutil.copy(f"{data_dir}/product/{src_spec}", f"{outdir}/{src}")
    shutil.copy(f"{data_dir}/product/{bkg_spec}", f"{outdir}/{src}")
    shutil.copy(f"{data_dir}/product/{arf_mat}", f"{outdir}/{src}")

    with fits.open(f"{data_dir}/product/{src_spec}") as hdul :
        hdr=hdul[0].header
        pos=hdr['XPROC0'].find('.rmf')
        rmf_file=hdr['XPROC0'][pos-16:pos]
    
    if not os.path.exists(f"{outdir}/{src}/{rmf_file}.rmf"):
        url=f"https://sasdev-xmm.esac.esa.int/pub/ccf/constituents/extras/responses/PN/{rmf_file}_v21.0.rmf"
        wget.download(url, out=f"{outdir}/{src}/{rmf_file}.rmf")

    with open(f"{outdir}/{src}/P{obsID}PNS003_Evts_{src}.reg", 'r') as file:
        lines = file.readlines()
        coordinates = lines[-1]
        expression = f"evselect table=/mnt/data/Amaury/candidates/{obsID}/P{obsID}PNS003PIEVLI0000.FTZ energycolumn=PI expression='#XMMEA_EP&&(PATTERN<=4)&& ((X,Y) IN {coordinates}) &&(PI in [200:12000])' withrateset=yes rateset='PN_source_lightcurve_raw.lc' timebinsize=0.01 maketimecolumn=yes makeratecolumn=yes"
        os.chdir(f"{outdir}/{src}")
        os.system(expression)

print('Done')