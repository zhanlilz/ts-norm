#!/usr/bin/env python

import sys
import os
import argparse
import subprocess

import numpy as np
from osgeo import gdal
gdal.AllRegister()

from iMad import run_MAD
from radcal import run_radcal

def getCmdArgs():
    p = argparse.ArgumentParser(description='''Normalize a target image by
        matching it to a reference image.''')
    p.add_argument(dest='target_raster', metavar='TARGET', 
        help='''A target raster to be matched with a given reference
        raster.''')
    p.add_argument(dest='reference_raster', metavar='REFERENCE', 
        help='''A reference raster that needs to cover the entire target raster
        and has the same spatial reference system as the target raster.''')
    p.add_argument('--thresh_ncp', dest='thresh_ncp', metavar='THRESH_NCP',
            type=float, default=0.95, required=False, help='''A threshold of
            no-change probability. Invariant pixels to use will be those pixels
            with a no-change probability given by the Multivariate Alteration
            Detection (MAD) algorithm above this threshold. The larger
            THRESH_NCP is, the more accurate but the fewer invariant pixels
            are. Default=0.95''')
    p.add_argument('--thresh_ndvi', dest='thresh_ndvi', metavar='THRESH_NDVI',
            type=float, default=None, required=False, help='''A threshold of
            NDVI to mask out non-vegetation pixels in the normalization
            procedure.  Do not set this option if you do not like to mask out
            any non-vegetation. Default: None.''')
    p.add_argument('-R', '--registration', dest='allow_reg', action='store_true')
    p.add_argument('-O', '--out_dir', dest='out_dir', metavar='OUT_DIR',
            default=None, required=True, help='''A directory to save output
            files.''')

    cmdargs = p.parse_args()

    return cmdargs

def gdalDataTypeCode2Name(code):
    if code==gdal.GDT_Unknown:
        return 'Unknown'
    elif code==gdal.GDT_Byte:
        return 'Byte'
    elif code==gdal.GDT_UInt16:
        return 'UInt16'
    elif code==gdal.GDT_Int16:
        return 'Int16'
    elif code==gdal.GDT_UInt32:
        return 'UInt32'
    elif code==gdal.GDT_Int32:
        return 'Int32'
    elif code==gdal.GDT_Float32:
        return 'Float32'
    elif code==gdal.GDT_Float64:
        return 'Float64'
    elif code==gdal.GDT_CInt16:
        return 'CInt16'
    elif code==gdal.GDT_CInt32:
        return 'CInt32'
    elif code==gdal.GDT_CFloat32:
        return 'CFloat32'
    elif code==gdal.GDT_CFloat64:
        return 'CFloat64'
    else:
        return 'Unknown'

def main(cmdargs):
    tgt_raster = cmdargs.target_raster
    ref_raster = cmdargs.reference_raster
    thresh_ncp = cmdargs.thresh_ncp
    thresh_ndvi = cmdargs.thresh_ndvi
    allow_reg = cmdargs.allow_reg
    out_dir = cmdargs.out_dir

    if not os.path.isdir(out_dir):
        try:
            os.mkdir(out_dir)
        except:
            msg = 'Failed to create output directory {0:s}'
            msg = msg.format(out_dir)
            raise RuntimeError(msg)
   
    # Reference and target rasters have the same projections? 
    tgt_ds = gdal.Open(tgt_raster, gdal.GA_ReadOnly)
    ref_ds = gdal.Open(ref_raster, gdal.GA_ReadOnly)
    tgt_dtype_code = tgt_ds.GetRasterBand(1).DataType
    tgt_dtype_name = gdalDataTypeCode2Name(tgt_dtype_code)
    ref_dtype_name = gdalDataTypeCode2Name(ref_ds.GetRasterBand(1).DataType)

    tgt_srs = tgt_ds.GetSpatialRef()
    ref_srs = ref_ds.GetSpatialRef()
    if not tgt_srs.IsSame(ref_srs):
        msg =  'Target and Reference rasters have '
        msg += 'different spatial reference systems!'
        del tgt_ds, ref_ds
        raise RuntimeError(msg)

    # run_MAD and run_radcal use 0 as the nodata values in input rasters.
    # let's make two intermediate rasters, one as target and the other as
    # reference, of the same resolution, the same dimension and both with 0 as
    # nodata values.
    tgt_gt = tgt_ds.GetGeoTransform()
    ref_gt = ref_ds.GetGeoTransform()
    tgt_xres, tgt_yres = abs(tgt_gt[1]), abs(tgt_gt[5])
    ref_xres, ref_yres = abs(ref_gt[1]), abs(ref_gt[5])
    # common coarse resolution
    com_xres, com_yres = max(tgt_xres, ref_xres), max(tgt_yres, ref_yres)
    tgt_com_raster = os.path.join(out_dir, 
            os.path.splitext(os.path.basename(tgt_raster))[0]+'_downsample.tif')
    ref_com_raster = os.path.join(out_dir, 
            os.path.splitext(os.path.basename(ref_raster))[0]+'_downsample.tif')
    subprocess.run(['gdalwarp', '-r', 'average', 
        '-tr', str(com_xres), str(com_yres), 
        '-dstnodata', str(0), 
        tgt_raster, tgt_com_raster], check=True)
    subprocess.run(['gdalwarp', '-r', 'average', 
        '-tr', str(com_xres), str(com_yres), 
        '-dstnodata', str(0), 
        ref_raster, ref_com_raster], check=True)

    # Downsampled reference raster needs to be larger than downsampled target
    # raster.
    tgt_com_ds = gdal.Open(tgt_com_raster, gdal.GA_ReadOnly)
    ref_com_ds = gdal.Open(ref_com_raster, gdal.GA_ReadOnly)
    if (ref_com_ds.RasterXSize < tgt_com_ds.RasterXSize 
            or ref_com_ds.RasterYSize < tgt_com_ds.RasterYSize):
        msg  = 'Reference raster cannot cover the entire target raster. '
        msg += 'Either use a larger reference raster or '
        msg += 'decrease the dimension of your target raster.'
        raise RuntimeError(msg)
    # Clip downsampled reference raster to the extent of target raster and
    # align raster grids to the downsampled target raster.
    ref_clip_raster = os.path.join(out_dir, 
            os.path.splitext(os.path.basename(ref_com_raster))[0]+'_clip.tif')
    subprocess.run(['gdal_calc.py', '-A', tgt_com_raster, '--calc', 'A-A',
        '--allBands', 'A', '--NoDataValue', str(0), '--type', ref_dtype_name, 
        '--outfile', ref_clip_raster], 
        check=True)
    subprocess.run(['gdalwarp', '-r', 'average', 
        ref_com_raster, ref_clip_raster], check=True)

    # Use NDVI to filter out non-veg. pixels
    if thresh_ndvi is not None:
        ndvi_com_raster = os.path.join(out_dir, 
                (os.path.splitext(os.path.basename(tgt_com_raster))[0]
                    +'_ndvi.tif'))
        subprocess.run(['gdal_calc.py', 
            '-A', tgt_com_raster, '--A_band', str(4),
            '-B', tgt_com_raster, '--B_band', str(3), 
            '--calc', '(A-B)/(A+B).astype(float)', 
            '--NoDataValue', str(0), '--type', 'Float32', 
            '--outfile', ndvi_com_raster], check=True)

    # Now build a master mask that considers nodata values of both target and
    # reference rasters and NDVI thresholding if requested.
    mask_raster = os.path.join(out_dir, 'master_mask.tif')
    cmd2run = ['gdal_calc.py', 
        '-A', tgt_com_raster, '-B', ref_clip_raster]
    if thresh_ndvi is None:
        cmd2run += ['--calc', '(A!=0)*(B!=0)']
    else:
        cmd2run += ['-C', ndvi_com_raster, 
                '--calc', '(A!=0)*(B!=0)*(C>{0:f})'.format(thresh_ndvi)]
    cmd2run += ['--type', 'Byte', '--outfile', mask_raster]
    subprocess.run(cmd2run, check=True)

    # Now apply the master mask to both target and reference rasters.
    tgt_masked_raster = os.path.join(out_dir, 
            (os.path.splitext(os.path.basename(tgt_com_raster))[0]
                +'_masked.tif'))
    ref_masked_raster = os.path.join(out_dir, 
            (os.path.splitext(os.path.basename(ref_clip_raster))[0]
                +'_masked.tif'))
    subprocess.run(['gdal_calc.py', 
        '-A', tgt_com_raster, '-B', mask_raster,
        '--allBands', 'A', 
        '--calc', '(B==1)*A', 
        '--NoDataValue', str(0), '--type', tgt_dtype_name, 
        '--outfile', tgt_masked_raster], check=True)
    subprocess.run(['gdal_calc.py', 
        '-A', ref_clip_raster, '-B', mask_raster, 
        '--allBands', 'A', 
        '--calc', '(B==1)*A', 
        '--NoDataValue', str(0), '--type', ref_dtype_name, 
        '--outfile', ref_masked_raster], check=True)

    if allow_reg:
        tgt_reg_raster = os.path.join(out_dir, 
                (os.path.splitext(os.path.basename(tgt_masked_raster))[0] 
                    + '_coreg.tif'))
        tiepoints_shp = os.path.join(out_dir, 
                (os.path.splitext(os.path.basename(tgt_masked_raster))[0] 
                    + '_coreg_tiepoints.shp'))
        subprocess.run(['arosics_cli.py', 'local',
            '-mp', '0', 
            '-o', tgt_reg_raster, 
            '-fmt_out', 'GTiff', 
            '-nodata', '0', '0', 
            '-br', '4', 
            '-bs', '4',
            '-tiepoint_grid', tiepoints_shp, 
            ref_masked_raster, tgt_masked_raster, '200'], check=True)
        tgt_masked_raster = tgt_reg_raster

    # Now we have a target raster and a reference raster ready for run_MAD and
    # run_radcal, that is, same resolution, same dimension in aligned grids,
    # and with nodata being zero.
    out_mad_raster = (os.path.splitext(os.path.basename(tgt_masked_raster))[0]
            +'_mad.tif')
    run_MAD(tgt_masked_raster, ref_masked_raster, out_mad_raster, 
            outdir=out_dir, 
            band_pos1=(1,2,3,4), band_pos2=(1,2,3,4), penalty=0.0)
    out_rad_raster = (os.path.splitext(os.path.basename(tgt_masked_raster))[0]
            +'_rad.tif')
    tgt_norm_raster = run_radcal(tgt_masked_raster, ref_masked_raster, 
            out_rad_raster, out_mad_raster, tgt_raster, 
            band_pos1=(1,2,3,4), band_pos2=(1,2,3,4),
            nochange_thresh=thresh_ncp, view_plots=True, 
            save_invariant=True, save_residuals=True,
            datatype_out=tgt_dtype_code, 
            outdir=out_dir)

    # Reapply the nodata value to the normalized target raster.
    tgt_final_raster = os.path.join(out_dir, 
            (os.path.splitext(os.path.basename(tgt_norm_raster))[0]
                +'_final.tif'))
    subprocess.run(['gdal_translate', 
        '-a_nodata', str(0), 
        '-mo', 
        'REFERENCE_IMAGE_SOURCE={0:s}'.format(os.path.basename(ref_raster)), 
        tgt_norm_raster, tgt_final_raster], check=True)

if __name__ == '__main__':
    cmdargs = getCmdArgs()
    main(cmdargs)
