#!/usr/bin/env python
'''

Latest update Apr 28th 2023 by David Cid (david.sanchez@ciemat.es)

This script will create a foreground objects map file using the Y3 Gold footprint 
Author: Nacho Sevilla (based Eli Rykoff's and Alex Drlica-Wagner's code)
'''

""" Need in test"""
#cd '~/local_python_packages/science_release/foreground_mask/y6/'
# cd '/scratch/davsan06/Doctorado/foreground_mask_Y6A2/y6'

import healpy as hp
import healsparse as hsp
import numpy as np
#from astropy.io import fits
# import fitsio
import sys
import os
from optparse import OptionParser
import copy
import re
# import pandas
import y6a2_mask_bits
from y6a2_mask_bits import *
from astropy.io import fits
from shapely import geometry
from shapely.geometry import Polygon, Point
# Timing and progress
import time
from tqdm import tqdm

""" Paths """
path_poly='/afs/ciemat.es/user/d/davsan06/local_python_packages/skymap/skymap/data'

""" Read polygon definition """
des19_poly=np.loadtxt(os.path.join(path_poly,'des-round19-poly.txt'))

""" Generate DES footprint polygon """
poly = geometry.Polygon([[p[0], p[1]] for p in des19_poly])

"""
 Function: check if the (RA,DEC) of a catalogue is inside DES footprint and
generate a mask
"""
def inPolygon(ra_in, dec_in, polygon):
    bool_in = []
    for r, d in zip(ra_in, dec_in):
        P = Point(r,d)
        if polygon.contains(P):
            bool_in.append(True)
        else:
            bool_in.append(False)
    return (np.array(bool_in))

def make_foremask_file(workdir,footprint_filename,foremask_filename,nside,bits,output_format):
    # nside is the nside for Healpix and the high resolution nside for Healsparse
    # output_format can be "healpix" or "healsparse"
    
    # nside = 16384
    # nside = 2048
    
    # Healsparse quantities
    nside_coverage = 32
    nside_sparse = np.copy(nside)
    
    print('Reading footprint in',footprint_filename)
    # Initialize empty foreground mask
    foremask = np.zeros(hp.nside2npix(nside),dtype=np.int32)    
    # Compute the pixel area for the given nside.
    pixarea = hp.nside2pixarea(nside,degrees=True)

    for b,d in tqdm(bits.items()):
#        print(b,d)
        # b = 128
        # d = {'name': 'Near the LMC', 'filename': None, 'mag': None, 'm1': None, 'm2': None, 'r1': None, 'r2': None, 'maxrad': None, 'minrad': None, 'cushion': None}
        # b = 1
        # d = {'name':'Gaia bright stars (G<7.)','filename':'/scratch/davsan06/Doctorado/foreground_mask_Y6A2/data/gaia_stars/gaia_stars_bright.fits','mag':'PHOT_G_MEAN_MAG','m1':min(gaia_bright_x),'m2':max(gaia_bright_x),'r1':0.03,'r2':0.008,'maxrad':max(gaia_bright_y),'minrad':min(gaia_bright_y),'cushion':0.5/3600.}
        
        # Extract the name of the sample
        mask_name = d['name']
        
        if mask_name == 'Total area': 
            continue
        print('Analyzing',mask_name)
        start_time = time.time()
        if re.search('lmc',mask_name, re.IGNORECASE):
            # Initialize empty pixel array
            hpix = np.arange(hp.nside2npix(nside))
            # Convert pixel to ra, dec
            ra,dec=hp.pix2ang(nside,hpix,lonlat=True,nest=True)
            # Acommodate RA
            ra = np.where(ra > 180, ra - 360, ra )
            # Apply rectangular cut on RA, DEC for the LMC
            badpix,=np.where((ra > 60) & (ra < 100) & (dec > -70) & (dec < -58))
        else:
            # Rest of the cases are Stars or bright objects for which we define an specific radii
            badpix = compute_radmask_badpix(workdir=workdir,
                                            footprint=footprint_filename,
                                            mask=d,
                                            nside=nside)
        # temp - specific mask per star subsample
        # foremask - final foreground mask with all foreground flags info
        # Initialize auxiliary map
        temp = np.ones(hp.nside2npix(nside))*hp.UNSEEN
        # Put bad pixels to 1
        temp[badpix]=1.
        # and put those pixels with the corresponding FOREGROUND_FLAG in the global foreground mask
        foremask[badpix] = foremask[badpix] | b
        print(mask_name,'masked',np.round(len(badpix)*pixarea,2),'square degrees')
        # Checking maps
        # hp.mollview(foremask,nest=True)
        # hp.mollview(hp.ma(foremask),nest=True)
        # hp.mollview(temp,nest=True)
        # hp.mollview(hp.ma(temp),nest=True)
        # Output filename
        if output_format == "healpix":
            # Output filename
            fname = mask_name.replace(" ", "").replace("<", "").replace(">", "").replace("(", "").replace(")", "")+'healpix.fits'
            hp.write_map(fname,hp.ma(temp),nest=True,coord='C',fits_IDL=False,partial=True,overwrite=True)
        elif output_format == "healsparse":
            # Convert healpix map to healsparse
            # Output filename
            fname = mask_name.replace(" ", "").replace("<", "").replace(">", "").replace("(", "").replace(")", "")+'healsparse.fits'
            hsp_map = hsp.HealSparseMap(nside_coverage=nside_coverage,
                                        nside_sparse=nside_sparse,
                                        healpix_map=hp.ma(temp))
            hsp_map.write(fname,clobber=True)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time -> {elapsed_time:.4f} secs")
    print("Writing map to",foremask_filename)
    #zero_,=np.where(foremask ==0)
    # Reduce the information, keep all the indices of the foreground map which are informed i.e.
    # whose pixel is marked with a flag different than zero
    notzero_,=np.where(foremask !=0)
    realfore = np.ones(hp.nside2npix(nside))*hp.UNSEEN
    realfore[notzero_]=foremask[notzero_]
    #foremask[zero_]=hp.UNSEEN
    if output_format == "healpix":
        hp.write_map(foremask_filename,hp.ma(realfore),nest=True,coord='C',fits_IDL=False,partial=True,overwrite=True)
    elif output_format == "healsparse":
        hsp_map = hsp.HealSparseMap(nside_coverage=nside_coverage,
                                    nside_sparse=nside_sparse,
                                    healpix_map=hp.ma(realfore))
        hsp_map.write(foremask_filename,clobber=True)

def compute_radmask_badpix(workdir,footprint,mask,nside):
    # nside must be the one from the initial definition and not read from the map
    # Extract the nside to work from the available footprint file
    # if footprint:
    #     nside=hp.npix2nside(footprint.size)
    # we now filter out the objects really outside the coarse footprint, and create radii for objects without them
    if mask['mag'] == None: 
        # Extract data
        hdu = fits.open(mask['filename'])
        cat = hdu[1].data
        # Homogenize RA, DEC column names
        cat.columns['_RAJ2000'].name = 'raj2000'
        cat.columns['_DEJ2000'].name = 'dej2000'
        # Extract radius
        maskrad=cat['radius']
    else:
        # Extract data
        hdu = fits.open(mask['filename'])
        cat = hdu[1].data
        # Homogenize the column names and extract the radius estimates for
        # each star sample
        if re.search('Gaia bright',mask['name'], re.IGNORECASE):
            cat.columns['RA'].name = 'raj2000'
            cat.columns['DEC'].name = 'dej2000'
            # maskrad = gaia_bright_interp(cat[mask['mag']])
            maskrad=y6a2_mask_bits_david.gaia_bright_interp(cat[mask['mag']])
        elif re.search('Gaia mid',mask['name'], re.IGNORECASE):
            cat.columns['RA'].name = 'raj2000'
            cat.columns['DEC'].name = 'dej2000'
            maskrad=y6a2_mask_bits_david.gaia_mid_interp(cat[mask['mag']])
        elif re.search('yale',mask['name'], re.IGNORECASE):
            cat.columns['RAJ2000'].name = 'raj2000'
            cat.columns['DEJ2000'].name = 'dej2000'
            maskrad=y6a2_mask_bits_david.yale_interp(cat[mask['mag']])
        elif re.search('2MASS bright',mask['name'], re.IGNORECASE):
            cat.columns['ra'].name = 'raj2000'
            cat.columns['dec'].name = 'dej2000'
            maskrad=y6a2_mask_bits_david.mass_bright_interp(cat[mask['mag']])
        elif re.search('2MASS faint',mask['name'], re.IGNORECASE):
            cat.columns['ra'].name = 'raj2000'
            cat.columns['dec'].name = 'dej2000'
            maskrad=y6a2_mask_bits_david.mass_faint_interp(cat[mask['mag']])
        elif re.search('LEDA',mask['name'], re.IGNORECASE):
            cat.columns['RAJ2000'].name = 'raj2000'
            cat.columns['DEJ2000'].name = 'dej2000'
            maskrad=y6a2_mask_bits_david.leda_interp(cat[mask['mag']])
        else:
            print('no good')
            sys.exit()
        # Define pixels where the estimated radius is higher than the maximum
        # radius and also the cases where the radius is lower than the minimum.
        hi,=np.where(maskrad > mask['maxrad'])
        lo,=np.where(maskrad < mask['minrad'])
        # Set those cases to the maximum and minimum values.
        maskrad[hi] = mask['maxrad']
        maskrad[lo] = mask['minrad']
    # Convert RA,DEC info to 3D position vector    
    vec=hp.ang2vec(cat['raj2000'],cat['dej2000'],lonlat=True)
    # Initialize empty map
    map_temp = np.zeros(hp.nside2npix(nside),dtype=np.int32)
    print('Defining foreground mask pixels...')
    for i in np.arange(maskrad.size): #for each object in coarse footprint
        # i = 10
        #check all pixels inside the avoidance radii + some cushion in high resolution
        pixint=hp.query_disc(nside,vec[i,:],(maskrad[i]+mask['cushion'])*np.pi/180.,inclusive=False,nest=True)
        map_temp[pixint] = -1#*(i+1)
    # Return bad pixels for the sample under consideration.
    #badpix,=np.where((map_temp < 0) & (footprint >= 0))
    badpix,=np.where((map_temp < 0)) #& (footprint > 0)) #pixels marked bad should in the radius AND in the footprint
    return badpix

def main():
    '''
    Run code with options
    '''
    usage = "%prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-d","--workdir",dest="syspath",help="Directory to read/store maps",default='./')
    parser.add_option("--table", action="store_true", dest="toggle_table", default=False, help="Printout coverage table")
    parser.add_option("-n","--nside",type="int",dest="nside",help="Healpix nside",default=4096)
    parser.add_option("--footprint_filename",dest="footprint_filename",help="Footprint filename (input, fits format)",default=False)
    parser.add_option("--foremask_filename",dest="foremask_filename",help="Foreground mask filename (output, fits format)",default='y6a2_foreground_mask_v1.0.fits.gz')
    parser.add_option("--output_format",dest="output_format",help="Format to save the output foreground mask. Options are healpix or healsparse",default='healsparse')
    # Parse command line
    (options, args) = parser.parse_args()
    workdir = options.syspath
    if options.footprint_filename:
        footprint_filename = workdir + options.footprint_filename
    else:
        footprint_filename = False
    foremask_filename = workdir + options.foremask_filename
    nside = options.nside
    output_format = options.output_format
    y6a2_mask_bits_david.init()
    bits = copy.deepcopy(y6a2_mask_bits_david.FOREGROUND_BITS)
    make_foremask_file(workdir,footprint_filename,foremask_filename,nside,bits,output_format)

if __name__ == '__main__':
    sys.exit(main())