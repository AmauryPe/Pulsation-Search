import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from stingray.pulse.pulsar import fold_events
from stingray.pulse.search import plot_profile, z_n_search
from stingray import EventList
from astropy.io import fits

"""
Plot the pulse profile of the observation corresponding to the source in soft and hard X-ray bands
"""



parser = argparse.ArgumentParser(description='Plot the pulse profil of the best candidate')
parser.add_argument("--csv", help="Path to the csv file of interest",type=str)
parser.add_argument("--identifier", help="Name of the source from Simbad",type=str, default ='Unknown')
parser.add_argument("--ra", help="ra of the source",type=str, default='Unknown')
parser.add_argument("--dec", help="dec of the source",type=str, default='Unknown')
parser.add_argument("--otype", help="type of the source",type=str, default = 'Unknown')
parser.add_argument("--f", help="Frequency to fold", type = float, default = None)
args = parser.parse_args()


identifier = args.identifier
otype=args.otype
ra = args.ra
dec = args.dec
csv_file = args.csv
words = csv_file.split('/')
outdir=''
obsID=words[5]
src=str(int(words[-1][24:27],16))
for i in range (len(words)-1):
    outdir+=words[i]+'/'

def get_info(fitsfile):
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        frame_time=int(hdr['FRMTIME'])
        duration=float(hdr['DURATION'])
        obs_mode=hdr['SUBMODE']
        filter_ID=hdr['FILTER']
    return frame_time, duration, obs_mode, filter_ID


if 'bary' in csv_file :
    fits_file = outdir+words[-1][:27]+'_bary.fits'
else :
    fits_file = outdir+words[-1][:27]+'.fits'

frame_time, duration, obs_mode, filter_ID = get_info(fits_file)

data = pd.read_csv(csv_file)
data = data.sort_values(by = 'power')

f = data['frequency'].iat[-1]

events=EventList.read(fits_file,'ogip')

# df_min = 1/duration
# oversampling=15
# df = df_min / oversampling

# if 'LS' in csv_file :
#     frequencies = np.arange(f - 10 * df, f + 10 * df, 0.1*df)
# else :
#     frequencies = np.arange(f - 100 * df, f + 100 * df, df)

# freq_z, zstat = z_n_search(events.time, frequencies, nbin=16, nharm=1)   
# z2 = np.max(zstat)
# i = np.where(zstat== z2)
# if args.f is None :
#     fr=freq_z[i][0]
# else :
#     fr= args.f

fr =args.f
if 1/fr >10 :
    nbin=32
else :
    nbin=10

fig1,axis1=plt.subplots(figsize=(10,7))
fig2,axis2=plt.subplots(figsize=(10,7))
fig3,axis3=plt.subplots(figsize=(10,7))

fig, (ax1,ax2,ax3) = plt.subplots(nrows=3,ncols=1, figsize=(10,20))
plt.rcParams['font.size'] = 16

ph, profile, profile_err = fold_events(events.time, fr, nbin=nbin)
ax1 = plot_profile(ph, profile, err = profile_err, ax=ax1)
ax1.legend(['Period : '+str(round(1/fr,6))+' s'], fontsize = 20)
ax1.xaxis.label.set_size(16)
ax1.yaxis.label.set_size(16)
ax1.tick_params(axis='x', labelsize=16)
ax1.tick_params(axis='y', labelsize=16)
ax1.set(title=f'{obsID}-{src} 0.2-12 keV')
axis1 = plot_profile(ph, profile, err = profile_err, ax=axis1)
axis1.legend(['Period : '+str(round(1/fr,6))+' s'], fontsize = 30)
axis1.xaxis.label.set_size(30)
axis1.yaxis.label.set_size(30)
axis1.tick_params(axis='x', labelsize=30)
axis1.tick_params(axis='y', labelsize=30)
axis1.set_title(f'{obsID}-{src} 0.2-12 keV',fontsize=35)
fig1.tight_layout()

ph, profile, profile_err = fold_events(events.time[(events.energy>0.2) & (events.energy<2)], fr, nbin=nbin)
ax2 = plot_profile(ph, profile, err = profile_err, ax=ax2)
ax2.legend(['Period : '+str(round(1/fr,6))+' s'], fontsize=20)
ax2.xaxis.label.set_size(16)
ax2.yaxis.label.set_size(16)
ax2.tick_params(axis='x', labelsize=16)
ax2.tick_params(axis='y', labelsize=16)
ax2.set(title=f'{obsID}-{src} 0.2-2 keV')
axis2 = plot_profile(ph, profile, err = profile_err, ax=axis2)
axis2.legend(['Period : '+str(round(1/fr,6))+' s'], fontsize=30)
axis2.xaxis.label.set_size(30)
axis2.yaxis.label.set_size(30)
axis2.tick_params(axis='x', labelsize=30)
axis2.tick_params(axis='y', labelsize=30)
axis2.set_title(f'{obsID}-{src} 0.2-2 keV',fontsize=35)
fig2.tight_layout()


ph, profile, profile_err = fold_events(events.time[(events.energy>2) & (events.energy<12)], fr, nbin=nbin)
ax3 = plot_profile(ph, profile, err = profile_err, ax=ax3)
ax3.legend(['Period : '+str(round(1/fr,6))+' s'],fontsize = 20)
ax3.xaxis.label.set_size(16)
ax3.yaxis.label.set_size(16)
ax3.tick_params(axis='x', labelsize=16)
ax3.tick_params(axis='y', labelsize=16)
ax3.set(title=f'{obsID}-{src} 2-12 keV')
axis3 = plot_profile(ph, profile, err = profile_err, ax=axis3)
axis3.legend(['Period : '+str(round(1/fr,6))+' s'],fontsize = 30)
axis3.xaxis.label.set_size(30)
axis3.yaxis.label.set_size(30)
axis3.tick_params(axis='x', labelsize=30)
axis3.tick_params(axis='y', labelsize=30)
axis3.set_title(f'{obsID}-{src} 2-12 keV',fontsize=35)
fig3.tight_layout()


fig.suptitle(f'obsID : {obsID}, src_num : {src} : ' + identifier + ', type : '+otype+ '\n RA : ' + ra + ' , DEC : ' + dec)
fig.savefig(outdir+'profile_'+words[-1][23:27]+'.png')

fig1.savefig(f"{outdir}{obsID}_{words[-1][23:27]}_0.2-12.png")
fig2.savefig(f"{outdir}{obsID}_{words[-1][23:27]}_0.2-2.png")
fig3.savefig(f"{outdir}{obsID}_{words[-1][23:27]}_2-12.png")