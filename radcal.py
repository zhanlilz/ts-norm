#!/usr/bin/env python
#******************************************************************************
#  Name:     radcal.py
#  Purpose:  Automatic radiometric normalization
#  Usage:             
#       python radcal.py  
#
#  Copyright (c) 2013, Mort Canty
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    Modified by Nicholas Leach

import sys, os, time
import numpy as np 
from scipy import stats
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly, GDT_Float32, GDT_Int32, GDT_UInt16, GDT_Byte
import matplotlib.pyplot as plt
import json
import math
from scipy.special import betainc
 
def run_radcal(image1, image2, outfile_name, iMAD_img, full_target_scene, band_pos1=(1,2,3,4), band_pos2=(1,2,3,4),
               nochange_thresh=0.95, view_plots=True, save_invariant=True, save_residuals=True,
               datatype_out=GDT_UInt16, outdir=None):
    """

    :param image1: Image which will receive radiometric calibration (target image). Include path.
    :param image2: Image with desired radiometry (reference image). Include path.
    :param outfile_name: Name for radcal image output. DO NOT INCLUDE PATH.
    :param iMAD_img: Output image from iMAD. Include path.
    :param full_target_scene: The full-resolution, full-size image to receive radiometric corrections.
    :param band_pos1: positions of the bands to be calibrated in image 1.
    :param band_pos2: positions of the bands to be calibrated in image 2.
    :param nochange_thresh:
    :param view_plots: (bool) display the radcal fit
    :param datatype_out: GDT datatype to save radcal file
    :return:
    """

    gdal.AllRegister()
    path, img1_name = os.path.split(image1)
    if path:
        os.chdir(path)      
#  reference image    
    file1 = image2
    if file1:                  
        inDataset1 = gdal.Open(file1, GA_ReadOnly)
        cols = inDataset1.RasterXSize
        rows = inDataset1.RasterYSize    
        bands = inDataset1.RasterCount
    else:
        return
    pos1 = band_pos2
    if not pos1:
        return   
    dims = [0,0,cols,rows]
    if dims:
        x10,y10,cols1,rows1 = dims
    else:
        return 
#  target image     
    file2 = image1
    if file2:                  
        inDataset2 = gdal.Open(file2,GA_ReadOnly)     
        cols = inDataset2.RasterXSize
        rows = inDataset2.RasterYSize    
        bands = inDataset2.RasterCount
    else:
        return   
    pos2 = band_pos1
    if not pos2:
        return 
    dims=[0,0,cols,rows]
    if dims:
        x20,y20,cols2,rows2 = dims
    else:
        return  
#  match dimensions       
    bands = len(pos2)
    if (rows1 != rows2) or (cols1 != cols2) or (len(pos1) != bands):
        sys.stderr.write("Size mismatch")
        sys.exit(1)             
#  iMAD image     
    file3 = iMAD_img
    if not os.path.exists(file3):
        file3 = os.path.join(outdir, os.path.split(iMAD_img)[1])
        if not os.path.exists(file3):
            exit("Can't find iMAD image. ")
    if file3:                  
        inDataset3 = gdal.Open(file3, GA_ReadOnly)
        cols = inDataset3.RasterXSize
        rows = inDataset3.RasterYSize    
        imadbands = inDataset3.RasterCount
    else:
        return   
    dims = [0, 0, cols, rows]
    if dims:
        x30,y30,cols,rows = dims
    else:
        return     
    if (rows1 != rows) or (cols1 != cols):
        sys.stderr.write("Size mismatch")
        sys.exit(1)    
#  outfile
    if outdir:
        dir_target = outdir
    else:
        dir_target = path
    outfile, fmt = os.path.join(dir_target, outfile_name), "GTiff"
    if not outfile:
        return    
#  full scene
    fsfile = full_target_scene
#  no-change threshold    
    ncpThresh = nochange_thresh
    if ncpThresh is None:
        return                 
    chisqr = inDataset3.GetRasterBand(imadbands).ReadAsArray(x30,y30,cols,rows).ravel()  # chi2 band of iMAD image
    ncp = 1 - stats.chi2.cdf(chisqr, [imadbands-1])  # NL: no change probability. chi2 cdf w 'bands' degrees of freedom
    idx = np.where(ncp > ncpThresh)[0]  # NL: 1D array w indices of pixels above the no-change threshold
#  split train/test in ratio 2:1 
    tmp = np.asarray(range(len(idx)))
    tst = idx[np.where(np.mod(tmp,3) == 0)]
    trn = idx[np.where(np.mod(tmp,3) > 0)]
    
    print('=========================================')
    print('             RADCAL')
    print('=========================================')
    print(time.asctime())
    print('reference: '+file1)
    print('target   : '+file2)
    print('no-change probability threshold: '+str(ncpThresh))
    print('no-change pixels (train): '+str(len(trn)))
    print('no-change pixels (test): '+str(len(tst)))
    driver = gdal.GetDriverByName(fmt)    
    outDataset = driver.Create(outfile, cols, rows, bands, datatype_out)
    projection = inDataset1.GetProjection()
    geotransform = inDataset1.GetGeoTransform()
    if geotransform is not None:
        gt = list(geotransform)
        gt[0] = gt[0] + x10*gt[1]
        gt[3] = gt[3] + y10*gt[5]
        outDataset.SetGeoTransform(tuple(gt))
    if projection is not None:
        outDataset.SetProjection(projection)      
    aa = []
    bb = []  
    i = 1
    log = []
    residuals = []  # list to save residuals
    predicted = []  # list to save predicted values at the invariant pixels
    for k in pos1:
        x = inDataset1.GetRasterBand(k).ReadAsArray(x10,y10,cols,rows).astype(float).ravel()  # x=reference image
        y = inDataset2.GetRasterBand(k).ReadAsArray(x20,y20,cols,rows).astype(float).ravel()  # y=target image
        b, a, R = orthoregress(y[trn], x[trn])  # trn is the vector of training points
        mean_tgt, mean_ref, mean_nrm = np.mean(y[tst]), np.mean(x[tst]), np.mean(a+b*y[tst])
        t_test = stats.ttest_rel(x[tst], a+b*y[tst])
        var_tgt, var_ref, var_nrm = np.var(y[tst]), np.var(x[tst]), np.var(a+b*y[tst])
        F_test = fv_test(x[tst], a+b*y[tst])
        if save_residuals:
            # need to write out the predicted value and the residual at that value
            # need to do this for each band since the predicted value and residual will be different for each
            predicted_k = a+b*y[idx]
            resid_k = x[idx] - (a+b*y[idx])  # taking residuals of both training and test datasets for band k
            residuals.append(resid_k)
            predicted.append(predicted_k)
        print('--------------------')
        print('spectral band:      ', k)
        print('slope:              ', b)
        print('intercept:          ', a)
        print('correlation:        ', R)
        print('means(tgt,ref,nrm): ', mean_tgt, mean_ref, mean_nrm)
        print('t-test, p-value:    ', t_test)
        print('vars(tgt,ref,nrm)   ', var_tgt, var_ref, var_nrm)
        print('F-test, p-value:    ', F_test)
        aa.append(a)
        bb.append(b)
        fit_info = {'k':k, 'b':b, 'a':a, 'R':R, 'mean_tgt':mean_tgt, 'mean_ref':mean_ref, 'mean_nrm':mean_nrm,
                    't_test':t_test, 'var_tgt':var_tgt, 'var_ref':var_ref, 'var_nrm':var_nrm, 'F_test':F_test}
        log.append(fit_info)
        outBand = outDataset.GetRasterBand(i)
        outBand.WriteArray(np.resize(a+b*y, (rows,cols)), 0, 0)
        outBand.FlushCache()
        if i <= 10:
            if view_plots:
                plt.figure(i)
                ymax = max(y[idx])
                xmax = max(x[idx])
                plt.plot(y[idx],x[idx],'k.',[0,ymax],[a,a+b*ymax],'k-')
                plt.axis([0,ymax,0,xmax])
                plt.title('Band '+str(k))
                plt.xlabel('Target')
                plt.ylabel('Reference')
                plt.show()
        i += 1

    # NL - save an image showing the invariant pixels
    if save_invariant:
        if img1_name.endswith("downsample.tif"):
            invar_name = os.path.join(dir_target, img1_name[:-14] + "invariants.tif")
        else:
            invar_name = os.path.join(dir_target, img1_name[:-4] + "invariants.tif")
        invariant_ds = driver.Create(invar_name, cols, rows, 1, GDT_Byte)
        if geotransform is not None:
            gt = list(geotransform)
            gt[0] = gt[0] + x10*gt[1]
            gt[3] = gt[3] + y10*gt[5]
            invariant_ds.SetGeoTransform(tuple(gt))
        if projection is not None:
            invariant_ds.SetProjection(projection)
        invar_band = np.zeros(rows*cols)
        invar_band[idx] = 1
        invar_band = np.resize(invar_band, (rows, cols))
        invariant_ds.GetRasterBand(1).WriteArray(invar_band)
        invariant_ds.FlushCache()

    # write out a log with radcal fit information
    if img1_name.endswith("downsample.tif"):
        log_outpath = os.path.join(dir_target, img1_name[:-14] + "radcal_parameters.json")
        # write out an array with all the predicted invariant pixel values and their residuals
        if save_residuals:
            resid_file_contents = predicted + residuals  # first rows are the predicted pixel values
            np.savetxt(img1_name[:-14] + "residuals.csv", np.array(resid_file_contents), delimiter=',',
                       header="1st half of rows are the predicted values. 2nd half are residuals.")
    else:
        log_outpath = os.path.join(dir_target, img1_name[:-4] + 'radcal_parameters.json')
        # write out an array with all the residuals
        if save_residuals:
            resid_file_contents = predicted + residuals  # first rows are the predicted pixel values
            np.savetxt(img1_name[:-4] + "residuals.csv", np.array(resid_file_contents), delimiter=',',
                       header="1st half of rows are the predicted values. 2nd half are residuals.")
    with open(log_outpath, "w") as write_file:
        json.dump(log, write_file)

    outDataset = None
    print('result written to: '+outfile)
    if fsfile is not None:
        path = dir_target
        basename = os.path.basename(fsfile)
        root, ext = os.path.splitext(basename)
        fsoutfile = path+'/'+root+'_norm'+ext        
        print('normalizing '+fsfile+'...')
        fsDataset = gdal.Open(fsfile,GA_ReadOnly)
        cols = fsDataset.RasterXSize
        rows = fsDataset.RasterYSize    
        driver = fsDataset.GetDriver()
        outDataset = driver.Create(fsoutfile, cols, rows, bands, datatype_out)
        projection = fsDataset.GetProjection()
        geotransform = fsDataset.GetGeoTransform()
        if geotransform is not None:
            outDataset.SetGeoTransform(geotransform)
        if projection is not None:
            outDataset.SetProjection(projection) 
        j = 0
        for k in pos2:
            inBand = fsDataset.GetRasterBand(k)
            outBand = outDataset.GetRasterBand(j+1)
            for i in range(rows):
                y = inBand.ReadAsArray(0,i,cols,1)
                outBand.WriteArray(aa[j]+bb[j]*y,0,i)   # this is where the operations happen. Operates by row
            outBand.FlushCache() 
            j += 1      
        outDataset = None    
        print('result written to: '+fsoutfile)

    print('-------done-----------------------------')
    if fsfile is not None:
        return fsoutfile
    else:
        return


def orthoregress(x, y):
    Xm = np.mean(x)
    Ym = np.mean(y)
    s = np.cov(x,y)  # diagonal elements (s[0,0] and s[1,1]) are standard deviations
    R = s[0, 1]/math.sqrt(s[1, 1]*s[0, 0])  # corr coef. covar divided by sqrt of product of standard deviations
    lam, vs = np.linalg.eig(s)
    idx = np.argsort(lam)
    vs = vs[:, idx]      # increasing order, so
    b = vs[1,1]/vs[0,1] # first pc is second column
    return [b, Ym-b*Xm, R]


def fv_test(x0, x1):
    # taken from IDL library
    nx0 = len(x0)
    nx1 = len(x1)
    v0 = np.var(x0)
    v1 = np.var(x1)
    if v0 > v1:
        f = v0/v1
        df0 = nx1-1
        df1 = nx0-1
    else:
        f = v1/v0
        df0 = nx1-1
        df1 = nx0-1
    prob = 2.0*betainc(0.5*df1, 0.5*df0, df1/(df1+df0*f))
    if prob > 1:
        return f, 2.0-prob
    else:
        return f, prob
