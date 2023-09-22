import os
import sys
import glob
import argparse
import re
import subprocess as sp
from astropy.io import fits
from astropy import log
from astropy.wcs import WCS
from astropy.timeseries import LombScargle
from hendrics.read_events import treat_event_file
from hendrics.efsearch import main_accelsearch
import hendrics
import luigi
from regions import Regions
import tempfile
import shutil
import pandas as pd
import numpy as np

"""
Find the main event files and all the timeseries from the product data.

Get the source extraction regions from the time series of each source

Filter

That's it.
"""

parser = argparse.ArgumentParser(description='Process an ObsID : apply the accelsearch algortihm')
parser.add_argument("--outdir", help="Output directory",type=str)
parser.add_argument("--directories", help="PPS directory",type=str)
parser.add_argument("--src", help="Source to process in the Observation, give the three characters in hexadecimal, i.e 001 or 0A5,",type=str, default = None)
args = parser.parse_args()
wrtdir=args.outdir
directories=[args.directories]    

obsid_re = re.compile(r"^P([0-9]{10})")
# src_re = re.compile(r"^P([0-9]{10})((PN)|(M1)|(M2))UEEETTTTTTSXXX.ZZZ [0-9]{10}")
src_re = re.compile(r"^.*([0-9A-F]{3})\.FTZ")

def get_sum_flag(obsID, source_num):
    hdu = fits.open('/home/amaury/code/4XMM_DR13cat_v1.0.fits')
    data = hdu[1].data
    src = data['SRC_NUM']
    obsid = data['OBS_ID']
    hdu.close()
    for i in range(len(src)):
        if src[i]==source_num and obsid[i]==obsID:
            return data['SUM_FLAG'][i]


def get_region(file):
    with fits.open(file) as hdul:
        regdata = hdul["REGION_SRC"].data
        shape = regdata["SHAPE"].strip().lower()

        if len(shape) > 1 or shape[0] != "circle":
            raise ValueError("Invalid region")

        shape = shape[0]
        x = regdata["X"][0]
        y = regdata["Y"][0]
        r = regdata["R"][0]
    return shape, x, y, r


def create_region_file(ts_file, outfile):
    shape, x, y, r = get_region(ts_file)
    print(x, y, r)
    regfile_str = f"""# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image
{shape}({x},{y},{r})"""

    outdir = os.path.split(outfile)[0]    
    os.makedirs(outdir, exist_ok=True)
    with open(outfile, "w") as fobj:
        print(regfile_str, file=fobj)
    return outfile


def create_sky_region_file(pn_file, region, outfile):
    reg = Regions.read(region).regions[0]
    with fits.open(pn_file) as hdul:
        wcs = WCS(hdul[1].header, keysel=['pixel'], naxis=['celestial'])
    reg_sky = reg.to_sky(wcs)
    reg_sky.write(outfile)
    return outfile


def filter_events(pn_file, region, outfile):
    r = Regions.read(region).regions[0]

    shape = "circle"
    x = r.center.x
    y = r.center.y
    r = r.radius

    expr = f"{shape}({x},{y},{r},X,Y)"
    command = f'evselect table={pn_file}:EVENTS expression={expr} filteredset={outfile}'

    if os.path.exists(outfile):
       return outfile

    print(command)
    sp.check_call(command.split())

    return outfile


def bary_events(fitsfile, odf_dir, region, outfile):
    reg = Regions.read(region).regions[0]
    ra = reg.center.ra.to("deg").value
    dec = reg.center.dec.to("deg").value

    shutil.copy(fitsfile, outfile)

    # script = f"""heainit
    # sasinit
    # SAS_ODF={odf_dir.strip()} barycen table={outfile}:EVENTS srcra={ra} srcdec={dec} withsrccoordinates=yes ephemeris=DE405
    # """
    script = f"""
    SAS_ODF={odf_dir.strip()} barycen table={outfile}:EVENTS srcra={ra} srcdec={dec} withsrccoordinates=yes ephemeris=DE405
    """
    scriptfile = outfile + '.sh'
    with open(scriptfile, "w") as fobj:
        print(script, file=fobj)

    cmd = f"sh {scriptfile}"
    sp.check_call(cmd.split())

    os.unlink(scriptfile)
    
def get_info(fitsfile):
    with fits.open(fitsfile) as hdul:
        hdr = hdul[0].header
        frame_time=int(hdr['FRMTIME'])
        duration=float(hdr['DURATION'])
        obs_mode=hdr['SUBMODE']
        filter_ID=hdr['FILTER']
    return frame_time, duration, obs_mode, filter_ID

# def filter_events(pn_file, directory, ts_file, outdir, outfile_root=None):
#     shape, x, y, r = get_region(ts_file)

#     src = src_re.match(ts_file).group(1)

#     log.info(f"Extracting source events for {src}")

#     if outfile_root is None:
#         outfile_root = f"src_{src}"

#     outfile = os.path.join(outdir, outfile_root + ".fits")
#     regfile = os.path.join(outdir, outfile_root + ".reg")

#     regfile_str = f"""# Region file format: DS9 version 4.1
# global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
# physical
# {shape}({x},{y},{r})"""

#     os.makedirs(outdir, exist_ok=True)

#     with open(regfile, "w") as fobj:
#         print(regfile_str, file=fobj)

#     expr = f"{shape}({x},{y},{r},X,Y)"
#     command = f'evselect table={pn_file}:EVENTS expression={expr} filteredset={outfile}'

#     print(command)

#     if os.path.exists(outfile):
#        return outfile

#     sp.check_call(command.split())

#     return outfile


class WrapAll(luigi.WrapperTask):
    fname = luigi.Parameter()
    # config_file = luigi.Parameter(default=None)
    config_file = luigi.Parameter('/home/amaury/luigi.cfg')
    worker_timeout = luigi.IntParameter(default=None)

    def requires(self):
        yield LS_periodogram(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearch(self.fname, self.config_file, self.worker_timeout)
        # yield AccelSearchDetrend(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearchRedNoise(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearchInterbin(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearchBary(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearchInterbinBary(self.fname, self.config_file, self.worker_timeout)
        yield AccelSearchRedNoiseBary(self.fname, self.config_file, self.worker_timeout)
        


class LS_periodogram(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEventsBary(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_LS_periodogram.csv"))

    def run(self):
        outfile = self.output().path
        print(outfile)
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        ev = hendrics.io.load_events(ncfile)
        ls=LombScargle(ev.time, ev.energy)
        frequency, power = ls.autopower(minimum_frequency=fmin,maximum_frequency=1)
        probability= 0.0013
        level=[ls.false_alarm_level(probability, method='naive')]*len(frequency)
        data = {'frequency': frequency, 'power': power, 'level': level}
        df = pd.DataFrame(data=data)
        df.to_csv(outfile, index=False)
        
class AccelSearch(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_accelsearch.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmax=0.5/(frame_time/1000)
        fmin=1/(0.25*duration)
        main_accelsearch(f"{ncfile} --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10 --outfile {outfile}".split())

class AccelSearchDetrend(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_detrend_accelsearch.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmax=0.5/(frame_time/1000)       
        fmin=1/(0.25*duration)
        main_accelsearch(f"{ncfile} --detrend 100 --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10 --outfile {outfile}".split())
        
class AccelSearchRedNoise(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_rednoise_accelsearch.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        fmax=0.5/(frame_time/1000)
        main_accelsearch(f"{ncfile} --red-noise-filter --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10 --outfile {outfile}".split())
        
class AccelSearchRedNoiseBary(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEventsBary(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_rednoise_bary_accelsearch.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        fmax=0.5/(frame_time/1000)
        main_accelsearch(f"{ncfile} --red-noise-filter --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10 --outfile {outfile}".split())


class AccelSearchInterbin(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_accelsearch_interbin.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEvents(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        fmax=0.5/(frame_time/1000)
        main_accelsearch(f"{ncfile} --fmin {fmin} --fmax {fmax} --emin 0.2 --emax 12 --nproc 10 --interbin --outfile {outfile}".split())


class AccelSearchBary(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=None)

    def requires(self):
        return ReadEventsBary(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_bary_accelsearch.csv"))

    def run(self):
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        fmax=0.5/(frame_time/1000)
        main_accelsearch(f"{ncfile} --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10".split())


class AccelSearchInterbinBary(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ReadEventsBary(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_bary_accelsearch_interbin.csv"))

    def run(self):
        outfile = self.output().path
        ncfile = ReadEventsBary(self.fname, self.config_file, self.worker_timeout).output().path
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        frame_time, duration, obs_mode, filter_ID = get_info(fitsfile)
        fmin=1/(0.25*duration)
        fmax=0.5/(frame_time/1000)
        main_accelsearch(f"{ncfile} --emin 0.2 --emax 12 --fmin {fmin} --fmax {fmax} --nproc 10 --interbin --outfile {outfile}".split())


class ReadEvents(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return ExtractEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_ev.nc"))

    def run(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = treat_event_file(fitsfile)[0]
        print(ncfile, fitsfile)
        shutil.move(ncfile, self.output().path)


class ReadEventsBary(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return BaryEvents(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = BaryEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_ev.nc"))

    def run(self):
        fitsfile = BaryEvents(self.fname, self.config_file, self.worker_timeout).output().path
        ncfile = treat_event_file(fitsfile)[0]
        print(ncfile, fitsfile)
        shutil.move(ncfile, self.output().path)


class BaryEvents(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        yield ExtractEvents(self.fname, self.config_file, self.worker_timeout)
        yield GetSkyRegion(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        return luigi.LocalTarget(fitsfile.replace(".fits", "_bary.fits"))

    def run(self):
        fitsfile = ExtractEvents(self.fname, self.config_file, self.worker_timeout).output().path
        directory, basename = os.path.split(self.fname)
        odf_dir  = None
        for dirname in ["product", "pps"]:
            if dirname in directory:
                odf_dir = directory.replace(dirname, "odf")
                # odf_dir = directory.replace(dirname, "ODF")
                break
        else:
            raise ValueError(f"Cannot find valid way to obtain odf_dir from {directory}")

        region = GetSkyRegion(self.fname, self.config_file, self.worker_timeout).output().path
        outfile = self.output().path

        tempdir = tempfile.mkdtemp(suffix=None, prefix=None, dir=None)

        tmpfile = os.path.join(tempdir, os.path.split(outfile)[1])

        bary_events(fitsfile, odf_dir, region, tmpfile)

        shutil.move(tmpfile, outfile)


class ExtractEvents(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return GetRegion(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        """
        Returns the target output for this task.
        In this case, a successful execution of this task will create a file
        on the local filesystem.
        :return: the target output for this task.
        :rtype: object (:py:class:`luigi.target.Target`)
        """
        # config = read_config(self.config_file)
        ts_file = self.fname

        directory, basename = os.path.split(ts_file)
        obsid = obsid_re.match(basename).group(1)
        outdir = obsid

        fitsfile = basename.replace("SRCTSR", "_Evts_").replace(".FTZ", "") +  ".fits"

        return luigi.LocalTarget(os.path.join(outdir, fitsfile))

    def run(self):
        # config = read_config(self.config_file)
        fname = self.fname

        ts_file = fname
        directory, basename = os.path.split(ts_file)
        obsid = obsid_re.match(basename).group(1)
        outdir = obsid

        src = src_re.match(ts_file).group(1)

        pn_file = glob.glob(os.path.join(directory, f"P{obsid}*PIEVLI*.FTZ"))[0]
        region = GetRegion(self.fname, self.config_file, self.worker_timeout).output().path
        outfile = self.output().path

        fitsfile = filter_events(pn_file, region, outfile)


class GetSkyRegion(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def requires(self):
        return GetRegion(self.fname, self.config_file, self.worker_timeout)

    def output(self):
        """
        Returns the target output for this task.
        In this case, a successful execution of this task will create a file
        on the local filesystem.
        :return: the target output for this task.
        :rtype: object (:py:class:`luigi.target.Target`)
        """
        # config = read_config(self.config_file)
        ts_file = self.fname

        directory, basename = os.path.split(ts_file)
        obsid = obsid_re.match(basename).group(1)
        outdir = obsid

        regfile = basename.replace("SRCTSR", "_Evts_").replace(".FTZ", "") + "_sky.reg"

        return luigi.LocalTarget(os.path.join(outdir, regfile))

    def run(self):
        # config = read_config(self.config_file)
        fname = self.fname

        ts_file = fname
        directory, basename = os.path.split(ts_file)
        obsid = obsid_re.match(basename).group(1)
        outdir = obsid

        src = src_re.match(ts_file).group(1)

        pn_file = glob.glob(os.path.join(directory, f"P{obsid}*PIEVLI*.FTZ"))[0]

        outfile = self.output().path

        region = GetRegion(self.fname, self.config_file, self.worker_timeout).output().path
        regfile = create_sky_region_file(pn_file, region, outfile)


class GetRegion(luigi.Task):
    fname = luigi.Parameter()
    config_file = luigi.Parameter(default=None)
    worker_timeout = luigi.IntParameter(default=120)

    def output(self):
        """
        Returns the target output for this task.
        In this case, a successful execution of this task will create a file
        on the local filesystem.
        :return: the target output for this task.
        :rtype: object (:py:class:`luigi.target.Target`)
        """
        # config = read_config(self.config_file)
        ts_file = self.fname

        directory, basename = os.path.split(ts_file)
        obsid = obsid_re.match(basename).group(1)
        outdir = obsid

        regfile = basename.replace("SRCTSR", "_Evts_").replace(".FTZ", "") + ".reg"

        return luigi.LocalTarget(os.path.join(outdir, regfile))

    def run(self):
        # config = read_config(self.config_file)
        fname = self.fname

        ts_file = fname

        outfile = self.output().path

        regfile = create_region_file(ts_file, outfile)


def main():
    obsID=directories[0].split('/')[-2]
    for directory in directories:
        if ("product" not in directory) and ("pps" not in directory):
            log.info(f"product or pps not in {directory}")
            continue
        if args.src is None :
            time_series_files = glob.glob(os.path.join(directory, f"P*PNS*TSR*"))
        else :
            src = args.src
            time_series_files = glob.glob(os.path.join(directory, f"P*PNS*TSR*{src}*"))
        for ts_file in time_series_files:
            src_num=int(ts_file[-7:-4],16)
            print(f"Processing {ts_file}")
            flag = get_sum_flag(obsID,src_num)
            if flag == 0 :
                print(f"Data quality is good enough, flag is {flag}, processing this source")
                luigi.build([WrapAll(ts_file)], workers=100)
            else :
                print(f"Data quality is not good enough, flag is {flag}, skipping this source")
            
    if os.path.exists(wrtdir+'/'+obsID):
         shutil.rmtree(wrtdir+'/'+obsID)
    shutil.move(os.getcwd()+'/'+obsID, wrtdir)

if __name__ == "__main__":
    main()