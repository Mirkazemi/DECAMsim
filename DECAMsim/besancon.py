#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 23:23:32 2018

@author: mirkazemi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
class besancon():
    """
    Besoncon model of stellar density
    
    This class provides density of F, G, K, and M stars afor a given magnitude and spatial 
    area.
    Parameters:
    -----------
    model_file : string
        Address of a file in format of CSV (',' separated) including following
        columns:
        min_mag: minimum magnitude range
        max_mag: maximum magnitude range        
        type: stellar type
        number: expected number of stars with given type and with magnitude 
        between min_mag_g and max_mag_g.
        The dfault model_file is a besoncon model for F, G, K, and M stars with
        g band magnitude between 14 and 20.
    """
    def __init__(self, model_file = None):
        if model_file is not None:
            self.model_density = pd.read_csv(model_file)
            
        else:   
            model_file = 'besancon_model/density_subtype_ra180_dec-60_mags.csv'
            self.model_density = pd.read_csv(os.path.join(os.path.dirname(__file__), model_file))
            '''
            try:
                
                model_file = 'DECAMsim/besancon_model/density_subtype_ra180_dec-60_mags.csv'
                self.model_density = pd.read_csv(os.path.join(os.path.dirname(__file__), model_file))
            except:
            '''  

            
    
    def sim_number(self, area = 1, random_seed = None):
        if random_seed is not None:
            np.random.seed(random_seed)

        # simulating the number density
        self.model_density['sim_number'] = [ np.random.poisson(n*area, 1)[0] for n in self.model_density['number']]

    def sim_stars_type(self, area = 1, random_seed = None):
        """
        returns a pandas.DataFrame of stars with magnitude in g band and stellar 
        type for given solid angle.
    
        Parameters:
        -----------
        area : int or float
            solid angle of the sky
            
        random_seed: int, optional
            If the random seed is set to value then the ouput will be reproducable.
    
        Returns:
        -----------
        pandas.DataFrame
        A DataFrame of randomly generated stars including their magnitudes in
        g band and stellar type (F, G, K, and M).
        """        
        self.sim_number(area = area, random_seed = random_seed)
        types = np.repeat(self.model_density['type'].values, self.model_density['sim_number'].values) 
        _min_mags = np.repeat(self.model_density['min_mag'].values, self.model_density['sim_number'].values) 
        _max_mags = np.repeat(self.model_density['max_mag'].values, self.model_density['sim_number'].values) 
        mags = np.random.uniform( low = _min_mags, high = _max_mags)

        return pd.DataFrame({'mag' : mags, 'type' : types})
    
