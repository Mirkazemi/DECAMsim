#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 13:51:56 2018
Unit test for besancon.py
@author: mirkazemi
"""

from .besancon import *
import numpy as np

''' Testing besancon.sim_number method '''         
def test_sim_number():
    area = 2
    random_seed = 0

    besancon1 = besancon()
    besancon1.sim_number(area = area, random_seed = random_seed)
    sim_number = [23, 7, 4, 1, 13, 13, 11, 0, 11, 0, 0, 0, 0, 0, 0]
    np.testing.assert_almost_equal(besancon1.model_density.sim_number.values[:15], sim_number)

''' Testing besancon.sim_catalog method '''         
def test_sim_catalog():
    area = 2
    random_seed = 0
    besancon1 = besancon()
    df = besancon1.sim_stars_type(area = area, random_seed = random_seed)
    sim_mag = [14.65342116, 14.72634246, 14.536923, 14.11047711, 14.40503561, 14.40537358,
               14.32104299, 14.02995032, 14.73725424, 14.10978446]
    np.testing.assert_equal(len(df), 1791)
    np.testing.assert_almost_equal(df.mag[:10].values, sim_mag)

    
    
    