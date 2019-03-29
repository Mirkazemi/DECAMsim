#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 14:05:40 2019
Unit test for descamsim.py
@author: mirkazemi
"""

import numpy as np
import pandas as pd
import os
from .decamsim import *

''' Testing DECamCCD class '''

decam_ccds_file = 'CCDs/wcs_62ccd.csv'

decam_ccds = pd.read_csv(os.path.join(os.path.dirname(__file__), decam_ccds_file))

''' Creating a DECamCCD object ''' 
ccd1 = DECamCCD(ID = decam_ccds['CCD'].values[0],
                naxis1 = decam_ccds['NAXIS1'].values[0], 
                naxis2 = decam_ccds['NAXIS2'].values[0], 
                ctype = [decam_ccds['CTYPE1'].values[0], 
                         decam_ccds['CTYPE2'].values[0]],
                crval = [decam_ccds['CRVAL1'].values[0], 
                         decam_ccds['CRVAL2'].values[0]], 
                crpix = [decam_ccds['CRPIX1'].values[0], 
                         decam_ccds['CRPIX2'].values[0]], 
                cd = [[decam_ccds['CD1_1'].values[0], decam_ccds['CD1_2'].values[0]],
                      [decam_ccds['CD2_1'].values[0], decam_ccds['CD2_2'].values[0]]])

''' Creating a random catalog of stars for testing a DECamCCD object '''
_pos = [[_ra,_dec]  for _ra in np.arange(14, 18, 0.01) for _dec in np.arange(-51, -48, 0.01) ]
catalog = pd.DataFrame(_pos, columns = ['RA', 'Dec'])

''' Test for DECamCCD.set_catalog method '''
def test_set_catalog():  
    ccd1.set_catalog(catalog = catalog, 
                     ra_column = 'RA', 
                     dec_column = 'Dec', 
                     margin = 0)
    np.testing.assert_equal(len(ccd1.catalog), 657)
    ccd1.set_catalog(catalog = catalog, 
                     ra_column = 'RA', 
                     dec_column = 'Dec', 
                     margin = 3)
    np.testing.assert_equal(len(ccd1.catalog), 646)
    
''' Test for DECamCCD.divide2chunks method '''    
def test_divide2chunks():
    
    ccd1.divide2chunks(dx1 = 32, dx2 = 16)
    
    _x1_min = [0.5, 32.5, 64.5, 96.5]
    _x2_min = [0.5, 16.5, 32.5, 48.5]
        
    _x1_max = [32.5, 64.5, 96.5, 128.5]
    _x2_max = [16.5, 32.5, 48.5, 64.5]
        
    np.testing.assert_equal( [ccd1.chunks['x1_min'].values[i] for i in [0, 256, 512, 768]],
                            _x1_min)    
    
    np.testing.assert_equal(ccd1.chunks['x2_min'].values[:4], _x2_min)    
    
    np.testing.assert_equal( [ccd1.chunks['x1_max'].values[i] for i in [0, 256, 512, 768]],
                            _x1_max)    
    
    np.testing.assert_equal(ccd1.chunks['x2_max'].values[:4], _x2_max)    
    
    np.testing.assert_equal(len(ccd1.chunks), 64*256)
    
    
''' Testing DECam class '''    
# setting arguments for craeting an instance of DECam class

ID = decam_ccds['CCD'].values
n_size = len(ID)

naxis1 = decam_ccds['NAXIS1'].values 
naxis2 = decam_ccds['NAXIS2'].values

ctype = [ [decam_ccds['CTYPE1'].values[i],  decam_ccds['CTYPE2'].values[i]] for i in range(n_size) ]

crval = [ [decam_ccds['CRVAL1'].values[i], decam_ccds['CRVAL2'].values[i]] for i in range(n_size) ] 
crpix = [ [decam_ccds['CRPIX1'].values[i], decam_ccds['CRPIX2'].values[i]] for i in range(n_size) ]  
cd = [ [[decam_ccds['CD1_1'].values[i], decam_ccds['CD1_2'].values[i]],
        [decam_ccds['CD2_1'].values[i], decam_ccds['CD2_2'].values[i]]]
      for i in range(n_size) ]
                      
# creating a DECam object
decam1 = DECam(ID = ID, naxis1 = naxis1, naxis2 = naxis2,  ctype = ctype,
               crval = crval, crpix = crpix, cd = cd)   

''' Test for DECam.set_catalog method '''
def test_DECam_set_catalog():  
    
       decam1.set_catalog(catalog = catalog, 
                          ra_column = 'RA', 
                          dec_column = 'Dec', 
                          margin = 0)
       
       np.testing.assert_equal(len(decam1.CCDs[0].catalog), 657)
       np.testing.assert_equal(len(decam1.CCDs[31].catalog), 678)
       
       decam1.set_catalog(catalog = catalog, 
                          ra_column = 'RA', 
                          dec_column = 'Dec', 
                          margin = 3)
       
       np.testing.assert_equal(len(decam1.CCDs[0].catalog), 646)       
       np.testing.assert_equal(len(decam1.CCDs[31].catalog), 673)
       
''' Test for DECam.set_catalog method '''
def test_DECam_divide2chunks():  
    
    decam1.divide2chunks(dx1 = 32, dx2 = 16)
    
    _x1_min = [0.5, 32.5, 64.5, 96.5]
    _x2_min = [0.5, 16.5, 32.5, 48.5]
        
    _x1_max = [32.5, 64.5, 96.5, 128.5]
    _x2_max = [16.5, 32.5, 48.5, 64.5]
        
    np.testing.assert_equal( [decam1.CCDs[10].chunks['x1_min'].values[i] for i in [0, 256, 512, 768]],
                            _x1_min)    
    
    np.testing.assert_equal(decam1.CCDs[20].chunks['x2_min'].values[:4], _x2_min)    
    
    np.testing.assert_equal( [decam1.CCDs[30].chunks['x1_max'].values[i] for i in [0, 256, 512, 768]],
                            _x1_max)    
    
    np.testing.assert_equal(decam1.CCDs[40].chunks['x2_max'].values[:4], _x2_max)    
    
    np.testing.assert_equal(len(decam1.CCDs[50].chunks), 64*256)    