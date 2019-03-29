#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 14:54:26 2018

@author: mirkazemi
"""

import pandas as pd
import numpy as np
from . import besancon 
from astropy.coordinates import SkyCoord
from astropy.coordinates import SkyOffsetFrame
import copy
from scipy.optimize import least_squares
import os 

def random_coordinates(ra0 = 0, dec0 = 0, radius = 1, n = 100, random_seed = None):
    """
    This function uniformly creates random RAs and DECs within a circle around origin.
     
    Parameters:
    -----------
    r : int or float, defaul value = 1
        Radius of the circle.
        
    n : int, default value = 100
        Number of random coordinates.
        
    random_seed : int, optional
        If the random seed is set to value then the ouput will be reproducable.
        
    Returns:
    --------
    A pandas.DataFrame object including 'RA' and 'DEC' columns. 
    """
    # minimizing this funcion allow us to find required offset [delta_ra, delta_dec]
    # for changing the center of random points
    def offset(x):
        center = SkyCoord(ra = x[0], dec = x[1], unit='deg')
        points0 = SkyCoord(ra = 0, dec = 0, unit='deg')
        a = points0.transform_to(SkyOffsetFrame(origin=center)) 
        return a.lon.value - ra0, a.lat.value - dec0
        
    # creating random points [_ra, _dec] around (0, 0)
    if random_seed is not None:
        np.random.seed(random_seed)    
    r = np.random.uniform(0, radius**2, n)**0.5
    teta = np.random.uniform(0,2*np.pi,n)
    _dec = r*np.sin(teta)
    _ra = r*np.cos(teta)

    """ moving the centre to ra0, dec0 """
    # finding the offset ra and dec (res.x[0], res.x[1]) to apply the offset
    res = least_squares(offset, [-ra0, -dec0])
    coords0 = SkyCoord(ra = _ra, dec = _dec, unit='deg')
    center = SkyCoord(ra = res.x[0], dec = res.x[1], unit='deg')
    coords = coords0.transform_to(SkyOffsetFrame(origin = center)) 

    return coords.lon.value, coords.lat.value


class CatalogSimulator():
    """
    A class for creating stars catalogs using Besancon model.
    It creates a pandas.DataFrame including given number of catalogs for stars 
    within a given radius around a given center.
    
    
    """
    def __init__(self):
        self.set_param()
        
    def set_param(self, radius = 1.1, catalog_number = 1, ra0 = 15.99, dec0 = -49.33):
        """
        sets the parameter of simulation.
        
        Parameters:
        -----------
        radius : int or float
                 Radius (in degrees) of a circle inclduing the simulated stars. 
                 The Default value set to 1.
                 
        catalog_number : int
                Number of catalogs to be generated.
                
        ra0 : int or float:
              Right acsention of the center of the circle.

        dec0 : int or float:
               Declination of the center of the circle.                 
        """
        self.radius = radius
        self.catalog_number = catalog_number
        self.ra0 = ra0
        self.dec0 = dec0
        self.besoncon_obj = besancon.besancon()
        
    def make_stars(self,random_seed = None):
        """
        This method creates a given number of simulated catalogs including 'type', 
        'RA', 'Dec', 'catalob_number' and a 'mag' column (a reference magnitude 
        for used in Besancon model'. It does not add any other magnitudes to the 
        catalog. The simulated catalogs are stored in a single pandas.DataFrame
        attribute called 'catalog'.
        
        Parameters:
        -----------            
        random_seed: int, optional
            If the random seed is set to value then the ouput will be reproducable.
        """
        
        if random_seed is not None:
                self.random_seed = random_seed 
        
        # creating a list of catalogs using Besoncon model.
        _cat_list = [self.besoncon_obj.sim_stars_type(area = np.pi * self.radius**2,
                                                      random_seed = self.random_seed + i) for i in range(self.catalog_number)]
    
        # adding a column as catalog number to each catalog.
        for i, _cat in enumerate(_cat_list):
            _cat['catalog_number'] = i
         
        # concating all the catalogs to a single pandas.DataFrame.
        self.catalog = pd.concat(_cat_list)
        del _cat_list
        
        # adding random RA and Dec to the stacked catalog
        _ra, _dec = random_coordinates(ra0 = self.ra0, 
                                       dec0 = self.dec0, 
                                       radius = self.radius, 
                                       n = len(self.catalog), 
                                       random_seed = self.random_seed)
        
        self.catalog['RA'] = _ra
        self.catalog['Dec'] = _dec
        
    def compute_magnitudes(self, template_file = '',
                       template_reference_mabgnitude = 'gp',
                       magnitudes = ['des_u', 'des_g', 'des_r', 'des_i', 'des_z',
                                     'gaia_G', 'gaia_BP', 'gaia_RP', 'gaia_RVS',
                                     'up', 'rp', 'ip', 'zp']):
        """
        This method adds simulated magnitudes from a template to the
        'catalog' attribute.
        
        Parameters:
        -----------  
        template_file : string
            Address file of the templates. The file must be in CSV format with 
            a header for columns name. It must include at a column 'type' that
            labels the stellar type in the templates. Using the 'type' column 
            in the template file and 'type' column in the 'catalog' attribute 
            (built by 'create_stars' method), stars in the 'catalog' matched
            by similar type of stars in the template.
            
        template_reference_mabgnitude : string
            the column name of the magnitude in the 'template_file' file that 
            is used to compute magnitude correction. All the mgnitudes in the templates 
            are computed for unknown luminosity. So the difference between
            'mag' in the catalog attribute DataFrame and 'template_reference_mabgnitude'
            in the template computed and added to all desired magnitudes in the template.
            
        magnitudes : list of string
            List of magnitude names in the templates that are intended to be added
            to the simulated catalog DataFrame ('catalog').
        """
    
        
        if template_reference_mabgnitude not in magnitudes:
            magnitudes.append(template_reference_mabgnitude)
        columns = copy.deepcopy(magnitudes)
        columns.append('type')
        
        if template_file == '':
            self.pickles_cat = pd.read_csv(os.path.join(os.path.dirname(__file__), 
                                                        "pickles/pickles_subtype_mcolor.csv"))[columns]
            
        else:
            self.pickles_cat = pd.read_csv(template_file)[columns]

        self.catalog = self.catalog.merge(self.pickles_cat, how = 'inner', on = ['type'])
        
        # computing the magnitude offset (delta_mag) between stars and template
        self.catalog['delta_mag'] = self.catalog['mag'] - self.catalog[template_reference_mabgnitude]
        
        # computing the the magnitude by add delta_mag to the template magnitudes
        for m in magnitudes:
            self.catalog[m] = self.catalog[m] + self.catalog['delta_mag']


