#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:55:11 2018
Unit test for catalogsimulator.py
@author: mirkazemi
"""

#import pandas as pd
import numpy as np
from .catalogsimulator import *


''' Testing random_coordinates function '''    
def test_random_coordinates():
    _ra, _dec = random_coordinates(ra0 = 10, dec0 = -20, radius = 2, n = 100, random_seed = 1000)
    np.testing.assert_equal(len(_ra), 100)
    ra_ref = np.array([11.35838005, 10.65656531, 11.79031823, 11.37521127,  
                       9.70536068,  9.33433898, 10.03937567, 10.77478609, 
                       10.34039572,  9.81398961])    
    
    dec_ref = np.array([-19.01267501, -20.28294941, -19.02360496, -20.51404233,
                        -21.84773297, -20.67804522, -20.40183695, -18.97274992,
                        -20.91148709, -21.82667695])
    
    np.testing.assert_almost_equal(_ra[:10], ra_ref )
    np.testing.assert_almost_equal(_dec[:10], dec_ref )   
    
    
''' Testing CatalogSimulator class '''    
CatSim = CatalogSimulator()

''' Testing CatalogSimulator.set_param method ''' 
def test_set_param():
    CatSim.set_param(radius = 1, catalog_number = 10, ra0 = 10, dec0 = -20)
    np.testing.assert_equal(CatSim.radius, 1)
    np.testing.assert_equal(CatSim.catalog_number, 10)
    np.testing.assert_equal(CatSim.ra0, 10)
    np.testing.assert_equal(CatSim.dec0, -20)

''' Testing CatalogSimulator.make_stars method '''     
def test_make_stars():
    CatSim.make_stars(random_seed = 2000)
    np.testing.assert_equal(len(CatSim.catalog), 28471)
    _ra = np.array([10.77737594, 9.90689099, 10.68948922, 9.56666549,
                    10.56995034, 10.55602183, 10.27576439, 10.54612462,
                    9.48379165, 10.29458086])
    
    _dec = np.array([-19.80960171, -20.74626358, -19.73930231, -19.58741008,
                     -19.70253922, -19.49190579, -20.03113309, -20.56755261,
                     -19.9497714,  -20.2902949])
    
    _mag = np.array([14.26850413, 14.68700891, 14.91493844, 14.71178022,
                     14.44567472, 14.69650801, 14.60086598, 14.98918613,
                     14.75384273, 14.37535552])

    np.testing.assert_almost_equal(CatSim.catalog.RA.values[:10], _ra )
    np.testing.assert_almost_equal(CatSim.catalog.Dec.values[:10], _dec )
    np.testing.assert_almost_equal(CatSim.catalog.mag.values[:10], _mag )
    
''' Testing CatalogSimulator.compute_magnitudes method '''         
def test_compute_magnitudes():

    CatSim.compute_magnitudes()
    _des_g = np.array([14.33410413, 14.75260891, 14.98053844, 14.77738022,
                       14.51127472, 14.76210801, 14.66646598, 15.05478613,
                       14.81944273, 14.44095552])
    
    np.testing.assert_almost_equal(CatSim.catalog.des_g.values[:10], _des_g )
    
    

    