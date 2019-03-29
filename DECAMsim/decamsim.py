#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 14:39:55 2018

@author: mirkazemi
"""

import pandas as pd
import numpy as np
import os
import copy
from astropy import wcs

class DECamCCD():
    def __init__(self, ID, naxis1, naxis2, ctype, crval, crpix, cd):
        """
        A class for each CCD in a camera. A CCD define by its WCS infos and its
        side size in pixels * pixels.

        Parameters:
        -----------        
        ID : string or int or float
        An ID for naming the CCD
        
        naxis1 : int
            Size of the X side of CCD in pixel
        
        naxis2 : int
            size of the Y side of CCD in pixel
            
        ctype : a list of string with the size of 2.
            the coordinate types in list for example: ['RA---TAN', ,'DEC--TAN'] 
        
        crval : list of float
            [RA of reference pixel (deg), Dec of reference pixel (deg)]
            
        crpix : list of float
            [X of reference pixel, Y of reference pixel]
            
        cd : A 2D list or array of floats with the shape of (2, 2)
            [ [RA deg per column pixel, RA deg per row pixel], [Dec deg per column pixel, Dec deg per row pixel] ]
            In othaer words:
            [[CD1_1, CD1_2], [CD2_1, CD2_2]]
        """
        self.ID = ID
        self.naxis1 = naxis1
        self.naxis2 = naxis2
        self.w = wcs.WCS(naxis=2)
        self.w.wcs.crpix = crpix
        self.w.wcs.crval = crval
        self.w.wcs.ctype = ctype
        self.w.wcs.cd = cd
        self.chunks = pd.DataFrame()
        self.catalog = pd.DataFrame()
        
    def set_catalog(self, catalog, ra_column, dec_column, margin = 0):
        """
        This method get a catalog in the format of pandas.DataFrame and only stores
        the objects in the catalog that are located within the CCD. It stores them
        in 'catalog' attribute. If any chunks was defined in the CCD (after running
        divide2chunks() method) it automatically find the host chunk for each object in the catalog
        and creates a column as 'chunk_ID' in the column to show which object
        locates in which chunk.

        Parameters:
        -----------
        catalog : pandas.DataFrame
            A catalog of objects.
            
        ra_column : string
            name of right acsension column in the catalog
        
        dec_column : string
            name of declination column in the catalog
            
        margin : int or float, optional and default is 3
            It defined the edge of CCD in pixel unit. If a object locates on the
            edge it won't be included in the stored catalog (in the 'catalog'
            attribute).
        """
        _x = self.w.wcs_world2pix(catalog[ra_column].values,
                                  catalog[dec_column].values, 1) 
        catalog['x1'] = _x[0]
        catalog['x2'] = _x[1]
        self.catalog = copy.deepcopy(catalog[ (margin <= catalog['x1']) &
                                              (margin <= catalog['x2']) &
                                              (catalog['x1'] <= self.naxis1 - margin) & 
                                              (catalog['x2'] <= self.naxis2 - margin) ])
       
        # find host chunks if the chunks were aleardy defined
        if 0 < len(self.chunks):
            self.assign_chunk_to_catalog()
    
    def divide2chunks(self, dx1, dx2):
        """
        This method divides the CCD into chunks. It stored the border of chunks 
        (their minimum and maximum pixel in X and Y) in 'chunk' attribute.
        If a catalog was already given (by calling 'set_catalog' method) then it
        automatically find the host chunk for each objects (rows) in the 'catalog'.
        
        Parameters:
        -----------
        dx1 : int, float
            size of chunks in pixel unit in the X direction        
        
        dx2 : int, float
            size of chunks in pixel unit in the Y direction        
        """
        _x1_min = np.arange(0.5, self.naxis1, dx1)
        _x2_min = np.arange(0.5, self.naxis2, dx2)
        _chunks = [[_x1, _x2] for _x1 in _x1_min for _x2 in _x2_min]
        self.chunks = pd.DataFrame(_chunks, columns = ['x1_min','x2_min'])
        self.chunks['x1_max'] = self.chunks['x1_min'] + dx1
        self.chunks['x2_max'] = self.chunks['x2_min'] + dx2
        self.chunks['chunk_ID'] = self.chunks.index
        # find host chunks if a catalog were aleardy given
        if 0 < len(self.catalog):
            self.assign_chunk_to_catalog()
            
    def assign_chunk_to_catalog(self):
        """
        This method finds the host chunk for each object and adds a column 
        ('chunk_ID') into the catalog to imply the host chunks.
        """
        self.catalog['chunk_ID'] = -1
        _chunk_ID = self.catalog['chunk_ID'].values
        _x1_min = self.chunks['x1_min'].values
        _x1_max = self.chunks['x1_max'].values
        _x2_min = self.chunks['x2_min'].values
        _x2_max = self.chunks['x2_max'].values
        _n_chunks = self.chunks['chunk_ID'].values
        
        _x1 = np.array(self.catalog['x1'].values)
        _x2 = np.array(self.catalog['x2'].values)

        def numba_worker():
            for i in _n_chunks:
                _f = (_x1_min[i] <= _x1) & (_x1 < _x1_max[i]) &\
                     (_x2_min[i] <= _x2) & (_x2 < _x2_max[i])
                _chunk_ID[ _f ] = i
            return _chunk_ID
              
        self.catalog['chunk_ID'] = numba_worker()
'''        
        for index, row in self.chunks.iterrows():
            _flag = (row['x1_min'] <= self.catalog.x1) &\
                    (self.catalog.x1 <  row['x1_max']) &\
                    (row['x2_min'] <= self.catalog.x2) &\
                    (self.catalog.x2 <  row['x2_max']) 
            self.catalog.loc[_flag, 'chunk_ID'] = row['chunk_ID']
'''   

class DECam():
    def __init__(self, ID = [], naxis1 = [], naxis2 = [], ctype = [], crval = [],
                 crpix = [], cd = []):
        """
        A class for simulating DECam camera. It consists of a list of DECam CCDs
        (a list of DECamCCD objects). The methods have the same name as the 
        mathods in the DECamCCD class. When a DECam method is called, the same
        method is called for DECamCCD objects in DECam.

        Parameters:
        -----------        
        ID : a list or array of string or int or float
            IDs for naming  CCDs
        
        naxis1 : a list or array of int
            Size of the X side of CCDs in pixel
        
        naxis2 : a list or array of int
            Size of the Y side of CCDs in pixel
            
        ctype : a 2d list or array of string with the shape of (n, 2)
            The coordinate types, for example: [['RA---TAN', 'DEC--TAN'],...] 
        
        crval : a 2d list or array of floats with the shape of (n, 2)
            RA and Dec of reference pixels for each CCDs (in degrees)
            For example for reference pixels in n CCDs:
            [[crval1_1, crval2_1], ..., [crval1_n, crval2_n]]
            
        crpix : a 2d list or array of floats with the shape of (n, 2)
            X1 and X2 of reference pixels for each CCDs (in unit of pixels)
            For example for reference pixels in n CCDs:
            [[crpix1_1, crpix2_1], ..., [crpix1_n, ccrpix2_n]]
            
        cd : a 3d list or array of floats with the shape of (n, 2, 2)
            A list of roration matrixes for each CCDs .
            For example for reference pixels in n CCDs:
            [ [[CD1_1_1, CD1_2_1], [CD2_1_1, CD2_2_1]],
            ...,
            [[CD1_1_n, CD1_2_n], [CD2_1_n, CD2_2_n]] ]
        """
        assert len(ID) == len(naxis1), 'ID and naxis1 must have the same length.'
        assert len(ID) == len(naxis2), 'ID and naxis2 must have the same length.'
        assert len(ID) == len(ctype), 'ID and ctype must have the same length.'
        assert len(ID) == len(crval), 'ID and crval must have the same length.'
        assert len(ID) == len(crpix), 'ID and crpix must have the same length.'
        assert len(ID) == len(cd), 'ID and cd must have the same length.'
        
        ''' If CCDs WCS parameters are not provided, they are set to values read
        from /CCDs/wcs_62ccd.csv file. '''
        
        if len(ID) == 0:
            decam_ccds_file = 'CCDs/wcs_62ccd.csv'
            decam_ccds = pd.read_csv(os.path.join(os.path.dirname(__file__), decam_ccds_file))

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

        self.CCDs = [ DECamCCD(ID = ID[i], 
                               naxis1 = naxis1[i],
                               naxis2 = naxis2[i],
                               ctype = ctype[i],
                               crval = crval[i],
                               crpix = crpix[i],
                               cd = cd[i]) for i in range(len(ID)) ]
            
            
    def set_catalog(self, catalog, ra_column, dec_column, margin = 0):
        """
        This method get a catalog in the format of pandas.DataFrame and only stores
        the objects in the catalog that are located within the CCDs. It stores them
        in 'CCDs[i].catalog' attributes. If any chunks was defined in the CCD (after running
        divide2chunks() method) it automatically find the host chunk for each object in the catalog
        and creates a column as 'chunk_ID' in the column to show which object
        locates in which chunk.

        Parameters:
        -----------
        catalog : pandas.DataFrame
            A catalog of objects.
            
        ra_column : string
            name of right acsension column in the catalog
        
        dec_column : string
            name of declination column in the catalog
            
        margin : int or float, optional and default is 3
            It defined the edge of CCDs in pixel unit. If a object locates on the
            edge it won't be included in the stored catalog (in the 'catalog'
            attribute).
        """
        for i, CCD in enumerate(self.CCDs):
            CCD.set_catalog(catalog = catalog, ra_column = ra_column,
                            dec_column = dec_column, margin = margin)

    
    def divide2chunks(self, dx1, dx2):            
        """
        This method divides the CCDs into chunks. It stored the border of chunks 
        (their minimum and maximum pixel in X and Y) in 'CCD[i].chunk' attribute.
        If a catalog was already given (by calling 'set_catalog' method) then it
        automatically find the host chunk for each objects (rows) in the 'catalog'.
        
        Parameters:
        -----------
        dx1 : int, float
            size of chunks in pixel unit in the X direction        
        
        dx2 : int, float
            size of chunks in pixel unit in the Y direction        
        """        
        for i, CCD in enumerate(self.CCDs):
            CCD.divide2chunks(dx1 = dx1, dx2 = dx2)       


    def return_catalog(self):
        for i, ccd in enumerate(self.CCDs):
            ccd['CCD_ID'] = i
        
        return pd.concat([ccd.catalog for ccd in self.CCDs])
