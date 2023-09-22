import os
import csv
import argparse
import pandas as pd
from astropy.io import fits

"""
Process and analyze some number of observations and keep track of it
"""
parser = argparse.ArgumentParser(description='Process all the observations in a directory')
parser.add_argument("--directory", type=str, help="Directory containing all observations to process")
parser.add_argument("--number", type=int, default=100, help="Number of observation to process in the drectory")
args = parser.parse_args()

directory=args.directory
print(str(len(os.listdir(directory))-1)+' observations found in ' +directory)
output_directory='/mnt/data/Amaury/Processed_data'

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def check_process(obsId, direc):
    for obs in os.listdir(direc) :
        if is_number(obs) and obs!='.ipynb_checkpoints' and int(obs)==int(obsId):
            return False
    return True

def check_obs_mode(obsID):
    for file in os.listdir(f"{obsID}/product"):
        if 'PNS' in file and 'SRSPEC' in file and file.endswith('.FTZ'):
            fitsfile=file
            with fits.open(f"{obsID}/product/{fitsfile}") as hdul:
                hdr = hdul[0].header
                obs_mode=hdr['SUBMODE']
                if obs_mode =="PrimeSmallWindow":
                    print(f'Observation mode : {obs_mode}')
                    return False
    return True 

count_p=0
count_a=0
tracker_file= f"{output_directory}/tracker.csv"
df = pd.read_csv(tracker_file, dtype ={'ObsID': str})
with open(tracker_file, 'a', newline='') as file:
    writer = csv.writer(file)
    for path in sorted(os.listdir(directory)):
        if count_p < args.number and path!='.ipynb_checkpoints' and check_process(path, output_directory) and check_obs_mode(directory+'/'+path) and path not in df['ObsID'].values:
            try :
                print('--------------------------')
                print('Processing observation : '+path)
                exit_status = os.system(f"python3 process_obsid_pipeline.py --directories {directory}/{path}/{product} --outdir {output_directory}")
                if exit_status != 0:
                    print('Observation '+path+ ' failed to process')
                    if path not in df['ObsID'] : 
                        writer.writerow([path, 'failed'])
                else :
                    print('Observation '+path+ ' processed succesfully') 
                    writer.writerow([path, 'succes'])
                    count_p+=1
            except :
                print('Observation '+path+ ' failed to process')
                if path not in df['ObsID'] : 
                    writer.writerow([path, 'failed'])
            try :   
                print('--------------------------')
                print('Analazing the candidates from observation '+path)
                exit_status = os.system("python3 analyse_candidate.py --obsID "+output_directory +'/'+path)
                if exit_status != 0:
                    print('Observation '+path+ ' failed to analyze')
                else :
                    print('Observation '+path+ ' analyzed succesfully')
                    count_a+=1
            except Exception:
                print('Observation '+path+ ' failed to analyze')

print(f"Summary : {count_p} observations processed, {count_a} observations analyzed")

            
            