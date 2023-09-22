import pandas as pd
import multiprocessing
import glob
import os
import numpy as np 
import hendrics
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt 
import hendrics.io
import hendrics.lcurve
import argparse
import csv as csv_p
from hendrics.efsearch import search_with_qffa
import stingray.pulse.pulsar
from stingray.pulse.pulsar import htest
from stingray.stats import z2_n_logprobability, equivalent_gaussian_Nsigma_from_logp, z2_n_detection_level,fold_detection_level
from stingray.pulse.search import search_best_peaks, epoch_folding_search, z_n_search
from stingray import EventList
from astropy.io import fits
from astropy import units as u
from astroquery.simbad import Simbad
import astropy.coordinates as coord
from functools import partial 


"""
Find all the csv files produced by the processing pipeline and apply the 'Quite fast folding' algorithm around the most powerful candidate
Compute the EF and Zn statistics for the best candidates and call for the pulse profile.
"""

parser = argparse.ArgumentParser(description='Analyze the candidates found by the processing pipeline')
parser.add_argument("--obsID", help="Path to observation id already processed",type=str)
args = parser.parse_args()

def get_info(fitsfile):
    # Gets the frame time, duration, the observation mode and the filter used for an observation
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        frame_time=int(hdr['FRMTIME'])
        duration=float(hdr['DURATION'])
        obs_mode=hdr['SUBMODE']
        filter_ID=hdr['FILTER']

    return frame_time, duration, obs_mode, filter_ID

def get_csv(parent_dir): 
    # List the csv and event files in the observation directory
    csv_files=[]
    ev_files=[]
    for file in os.listdir(parent_dir):
        if file.endswith('.csv') and not file.startswith('Summary'):
            df=pd.read_csv(parent_dir+'/'+file)
            if len(df) > 1: #the header is not counted as a line in the csv file
                csv_files.append(file)
        elif file.endswith('_ev.nc'):
            ev_files.append(file)

    return sorted(csv_files), ev_files

def get_source_number(csv_files):
    # List all the different source number for a list of csv files
    obs_num_list=[]
    for file in csv_files:
        obs_num_list.append(file[23:27])
    set_res = set(obs_num_list) 
    obs_num_list = (list(set_res))
    count_list=[]
    csv_list=[]
    for num in obs_num_list : 
        count=0
        num_csv_list=[]
        for csv in csv_files:
            if csv[23:27] == num and csv[13]=='S':
                count+=1
                num_csv_list.append(csv)
        count_list.append(count)
        csv_list.append(num_csv_list)
    
    return obs_num_list, count_list, csv_list

def get_ra_dec(fits_file):
    # return the RA and DEC in degree of a source from its fits file
    words = fits_file.split('/')
    sky_reg_file=''
    for i in range (len(words)-1):
        sky_reg_file+=words[i]+'/'
    sky_reg_file+=words[-1][:27]
    sky_reg_file+='_sky.reg'
    with open(sky_reg_file,'r') as f :
        lines = f.readlines()
        line = lines[-1].replace('circle(','').replace(')','')
        numbers = line.split(',')
        ra = numbers[0]
        dec = numbers[1]

        return ra, dec


def query_simbad(ra, dec):
    # Perform the query
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('otype')
    c1=coord.SkyCoord(ra=float(ra)*u.degree, dec=float(dec)*u.degree, unit='deg')
    ra = coord.Angle(float(ra)*u.degree)
    dec = coord.Angle(float(dec)*u.degree)

    if dec > 0:
        qry = (f"region(circle, {str(ra)} +{str(dec)}, 10s) &"
                " otypes in ('X','Pulsar','Neutron*','WhiteDwarf','**','Supernova','Galaxy','V*','InteractingG','PairG','BH') ")
    else : 
        qry = (f"region(circle, {str(ra)} {str(dec)}, 10s) &"
                " otypes in ('X','Pulsar','Neutron*','WhiteDwarf','**','Supernova','Galaxy','V*','InteractingG','PairG','BH') ")
    
    result_table = custom_simbad.query_criteria(qry)
    if result_table is not None:
        sep_l=[]
        for obj in result_table :
            ra = obj['RA']
            dec = obj['DEC']
            c2=coord.SkyCoord(ra=str(ra), dec=str(dec), unit=(u.hourangle,u.deg))
            sep = c1.separation(c2)
            sep_l.append(sep)
        df = result_table.to_pandas()
        df['sep']=sep_l
        df=df.sort_values(by='sep')
        obj=df.iloc[0]
        identifier = obj['MAIN_ID']
        RA = obj['RA']
        DEC = obj['DEC']
        otype = obj['OTYPE']
        sep=obj['sep']

        return identifier, RA, DEC, otype, sep
    else:

        return None

def compute_stats(obs_length, period,fdot_i,ls, events, nbin=16, nharm=1):

    # perform the EF and Zn search around a period for a given events serie
    # return the best candidate fr, the list of frequencies, the detection levels, and the EF/Zn statistics

    df_min = 1/obs_length
    oversampling=15
    df = df_min / oversampling
    frequencies = np.arange(max(1/period - 100 * df,4*df_min), 1/period + 100 * df, 0.1*df)
    freq_ef, efstat = epoch_folding_search(events.time, frequencies, nbin=16)
    freq_z, zstat = z_n_search(events.time, frequencies, nbin=16, nharm=1)   
    z2 = np.max(zstat)
    i = np.where(zstat== z2)
    fr=freq_z[i]
    fdotr=0

    ntrial = (frequencies[-1] - frequencies[0]) / df_min
    ef_detlev = fold_detection_level(nbin, epsilon=0.0013, ntrial=ntrial)
    z_detlev = z2_n_detection_level(n=1, epsilon=0.0013, ntrial=ntrial)

    cand_freqs_ef, cand_stat_ef = search_best_peaks(freq_ef, efstat, ef_detlev)
    cand_freqs_z, cand_stat_z = search_best_peaks(freq_z, zstat, z_detlev)

    if len(cand_freqs_z) :
        best_peak_power=max(cand_stat_z)
        best_peak_freq=cand_freqs_z[cand_stat_z==best_peak_power]
        fr=best_peak_freq

    return fr, frequencies, freq_z,freq_ef,z_detlev,ef_detlev,zstat,efstat ,cand_freqs_ef,cand_freqs_z

def find_proba(df,ev,ls,duration,oversample):
    # Perform the quite fast folding algorithm around the best candidate
    rslt_df = df[df['frequency'].between(df['frequency'].iat[-1]-0.0005,  df['frequency'].iat[-1]+0.0005)]
    mean_freq=np.mean(rslt_df.frequency)
    std_freq=np.std(rslt_df.frequency)
    if ls :
        fdot_i=0
    else :
        fdot_i=df['fdot'].iat[-1]
    print(f"Using {len(rslt_df)}/{len(df)} frequencies around {mean_freq}") 
    print(f"Searching for a pulsation in the interval [{np.max([4/duration,mean_freq-0.0005])},{mean_freq+0.0005}]")
    times=ev.time
    all_fgrid, all_fdotgrid,all_stats,step,fdotstep,length = search_with_qffa(times,np.max([4/duration,mean_freq-0.0005]),mean_freq+0.0005, oversample=16)
    z2 = np.max(all_stats[all_fgrid >0.0001])
    i,j = np.where(all_stats == z2)
    fr=all_fgrid[i,j]
    fdotr=all_fdotgrid[i,j]
    print(f" Found a pulsation at f={fr} Hz and f_dot={fdotr} \n Mean frequency was {mean_freq} Hz")

    ntrial=int(len(all_fgrid[0])/oversample)
    p=z2_n_logprobability(z2, n=1, ntrial=ntrial)
    sigma=equivalent_gaussian_Nsigma_from_logp(p)
    print(" The pulsation was detected at "+str(sigma)+ " sigma")
    return sigma, fr, fdotr, mean_freq, std_freq, all_fgrid, fdot_i


def plot_step(count,k, csv, df, sigma, fr, fdotr,fdot_i,ls, frame_time, all_fgrid,events,obs_length,mean_freq):
    # Plot a "row" depending on the algorithm used for the accelsearch
    # First column is a zoom on the periodogram of the second column, with the area around the best candidate
    # Third column is the EF/Zn statistic

    plt.subplot(count,3,k)
    plt.title(csv[28:-4])
    plt.scatter(df.frequency,df.power,color='blue')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power')
    plt.axvline(x=fr,label=f"Frequency found : {fr} \n sigma={sigma}",color='r',ls='--')
    if ls and 'level' in df.columns :
        plt.axhline(df['level'].iat[-1],color='r', linestyle='--')
    plt.xlim([all_fgrid[0,0]-0.002,all_fgrid[0,-1]+0.002])
    plt.fill_betweenx(np.arange(np.min(df.power),np.max(df.power),0.1),all_fgrid[0,0],all_fgrid[0,-1],alpha=0.5,color='green',label='Searching range')
    plt.legend()
    
    k+=1
    plt.subplot(count,3,k)
    plt.title(csv[28:-4])
    plt.scatter(df.frequency,df.power, color='blue')
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power')
    plt.axvline(x=fr,label=f"Frequency found : {fr} \n sigma={sigma}",color='r',ls='--')
    if ls:
        plt.xlim([0,1])
        if 'level' in df.columns:
            plt.axhline(df['level'].iat[-1],color='r', linestyle='--')
    else : 
        plt.xlim([-1,0.5/(frame_time/1000)])
    plt.legend()

    k+=1
    nbin=16
    nharm=1
    period=1/fr
    f, frequencies, freq_z,freq_ef,z_detlev,ef_detlev,zstat,efstat ,cand_freqs_ef,cand_freqs_z = compute_stats(obs_length, period, fdot_i,ls, events, nbin=nbin, nharm=nharm)
    
    ax=plt.subplot(count,3,k)
    plt.title(csv[28:-4])
    plt.axhline(z_detlev - nharm, label='$Z^2_1$ det. lev.')
    plt.axhline(ef_detlev - nbin + 1, label='EF det. lev.', color='gray')

    plt.plot(freq_z, (zstat - nharm), label='$Z^2_1$ statistics')
    plt.plot(freq_ef, efstat - nbin + 1, color='gray', label='EF statistics', alpha=0.5)

    for c in cand_freqs_ef:
        plt.axvline(c, ls='-.', zorder=10)
    for c in cand_freqs_z:
        plt.axvline(c, ls='--', zorder=10)

    plt.axvline(1/period, color='r', lw=3, alpha=0.5, label='Initial Frequency')
    plt.xlim([frequencies[0], frequencies[-1]])
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Statistics - d.o.f.')
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.6f'))
    plt.legend()

    return f, len(cand_freqs_ef),len(cand_freqs_z)

def par_analysis(index,path):
    #parralelized function 

    with open(path, 'a', newline='') as file:
        writer = csv_p.writer(file)
        num=obs_num_list[index]
        count=count_list[index]
        num_csv_list=csv_list[index]
        k=1
        fig=plt.figure(figsize=(30,int(5*count)))
        for csv in num_csv_list:
            print(csv)
            df = pd.read_csv(parent_dir+'/'+csv)
            if 'bary' in csv :
                fits_file=csv[:27]+'_bary.fits'
            else : 
                fits_file=csv[:27]+'.fits'
            ev = EventList.read(parent_dir+'/'+fits_file, 'ogip')
            frame_time, duration, obs_mode, filter_ID = get_info(parent_dir+'/'+fits_file)
            df=df.sort_values(by='power')
            if 'LS' in csv:
                ls=True
            else : 
                ls=False
            if 1/df['frequency'].iat[-1]<100:
                df=pd.DataFrame()
            if not df.empty:
                sigma, fr, fdotr, mean_freq, std_freq, all_fgrid, fdot_i =find_proba(df,ev,ls,duration,oversample=16)
                f, n_cand_ef,n_cand_z = plot_step(count,k, csv, df, sigma, fr, fdotr, fdot_i, ls,frame_time, all_fgrid, ev,duration,mean_freq)
                ra , dec = get_ra_dec(parent_dir+'/'+fits_file)
                if n_cand_z > 0 :
                    pulse = True
                    if query_simbad(ra, dec) is not None :
                        identifier, ra, dec, otype,sep = query_simbad(ra, dec)
                    else :
                        identifier = 'Unknown'
                        otype = 'Unknown'
                        sep = None
                        print("No known object found around the specified coordinates.")
                    os.system(f"python3 pulse_profil.py --f {f[0]} --csv {parent_dir}/{csv} --identifier " +'"' +f"{identifier}" +'"' +" --otype " +'"' +f"{otype}" +'"'  +" --ra " +'"' +f"{ra}" +'" --dec' +' "' +f"{dec}"+'"')
                else :
                    pulse = False
                    identifier = None
                    otype= None
                    sep=None
                writer.writerow([num, identifier, otype, ra, dec,sep, csv[28:-4], fr, sigma, int(sigma>5),pulse])  
                k+=1
            else :
                print('No candidates found')
                k+=1
            k+=2

        plt.tight_layout()
        plt.savefig(parent_dir+'/'+num)
        file.close()

parent_dir=args.obsID
obsID=parent_dir[-10:]
csv_files,ev_files=get_csv(parent_dir)
obs_num_list,count_list,csv_list=get_source_number(csv_files)
path=parent_dir+'/Summary_'+obsID+'.csv'
with open(path, 'w+', newline='') as file:
    writer = csv_p.writer(file)
    field = ["Source", "Identifier","Otype", "RA", "DEC","Distance", "Search","Freqency", "Sigma", ">5 Sigma","Pulse Profil"]
    writer.writerow(field)
    file.close()

index_liste = range(len(obs_num_list))
num_processes = multiprocessing.cpu_count()

with multiprocessing.Pool(num_processes) as pool :
    par_analysis_writer=partial(par_analysis,path=path)
    pool.map(par_analysis_writer,index_liste)


# src_list=''
# for index in index_liste:
#     src=obs_num_list[index][1:]
#     if src[0]=='0':
#         src=src[1:]
#         if src[0]=='0':
#             src=src[1:]
#     src_list+=src+' '

# os.system(f"python3 move_candidate.py --obsID {obsID} --src {src_list}")







