#!/usr/bin/env python

# Generate a raster template for each AAFC map that aligns with VLCE final mosaic. 

import argparse

import subprocess


def getCmdArgs():
    p = argparse.ArgumentParser(description="Use one or multiple rasters of the same or different projections to clip another single master raster of either the same or a different projection. The output subset of the master raster from the clipping will remain in the projection of the master raster, and the pixels in the output subset will remain aligned with the master raster.")

    p.add_argument("--master_raster", default=None, required=True, help="The input raster you want to clip.")
    p.add_argument("--clip_raster", nargs='+', default=None, required=True, help="The raster/s you use to clip the master raster. If more than one clip raster are provided, the output raster will be clipped from the master raster to cover all the given clip rasters. If the extent of the given clip raster/s is partially/completely outside the master raster, the outside part will be discarded and ONLY the intersect part will be output.")
    p.add_argument("--out_raster", default=None, required=True, help="The output raster from the clipping.")

    p.add_argument("--out_format", default=None, required=False, help="The output raster format. Get a list of raster format names, gdalinfo --formats. Default: GTiff")
    p.add_argument("--cache", type=int, default=None, required=False, help="The GDAL raster block cache size in MB. ")

    cmdargs=p.parse_args()

    return cmdargs


def parseGdalinfo(gdalinfo_stdout, fields=None):
    if fields is None:
        return gdalinfo_stdout

    results = []
    for fd in fields:
        idx1 = gdalinfo_stdout.find(fd)
        idx2 = gdalinfo_stdout[idx1:].find("\n")
        results.append(gdalinfo_stdout[idx1:][:idx2].lstrip(fd))

    return results


def subMaster(sub_raster_list, master_raster, out_raster, options=[]):
    # Get the four corners of the subordinate raster
    sub_corners_list = []
    for sub_raster in sub_raster_list:
        cmd_return = subprocess.run(["gdalinfo", sub_raster], stdout=subprocess.PIPE)
        tmp = parseGdalinfo(cmd_return.stdout.decode(), fields=["Upper Left", "Lower Left", "Upper Right", "Lower Right"])
        sub_corners = []
        for tstr in tmp:
            sub_corners.append( tuple([float(vstr) for vstr in tstr[tstr.find("(")+1:tstr.find(")")].split(",")]) )
        sub_corners_list.append(sub_corners)

    # Get the boundary of the master raster in master raster's projection
    cmd_return = subprocess.run(["gdalinfo", master_raster], stdout=subprocess.PIPE)
    tmp = parseGdalinfo(cmd_return.stdout.decode(), fields=["Upper Left", "Lower Left", "Upper Right", "Lower Right"])
    master_corners = []
    for tstr in tmp:
        master_corners.append( tuple([float(vstr) for vstr in tstr[tstr.find("(")+1:tstr.find(")")].split(",")]) )
    master_x, master_y = zip(*master_corners)

    # Reproject the corners into the SRS of master raster
    # Get the SRS string of sub. and master rasters
    master_srs_proj4 = subprocess.run(["gdalsrsinfo", "-o", "proj4", master_raster], stdout=subprocess.PIPE).stdout.decode().strip().strip("'")
    sub_corners_in_master_list = []
    for sub_raster, sub_corners in zip(sub_raster_list, sub_corners_list):
        sub_srs_proj4 = subprocess.run(["gdalsrsinfo", "-o", "proj4", sub_raster], stdout=subprocess.PIPE).stdout.decode().strip().strip("'")
        sub_corners_in_master = []
        for sc in sub_corners:
            with subprocess.Popen(["gdaltransform", "-s_srs", sub_srs_proj4, "-t_srs", master_srs_proj4, "-output_xy"], 
                    stdin=subprocess.PIPE, stdout=subprocess.PIPE) as proc:
                tmp, _ = proc.communicate(input="{0[0]:f} {0[1]:f}".format(sc).encode())
                sub_corners_in_master.append( tuple(float(vstr) for vstr in tmp.decode().strip().split()) )
        sub_corners_in_master_list.append(sub_corners_in_master)

    x, y = zip(*sum(sub_corners_in_master_list, []))
    min_x = min(x)
    max_x = max(x)
    min_y = min(y)
    max_y = max(y)

    # Check with the boundry of the master raster
    min_x = min_x if min_x>=min(master_x) else min(master_x)
    max_x = max_x if max_x<=max(master_x) else max(master_x)
    min_y = min_y if min_y>=min(master_y) else min(master_y)
    max_y = max_y if max_y<=max(master_y) else max(master_y)
    
    cmd_return = subprocess.run(["gdal_translate"] +  options + ["-projwin", str(min_x), str(max_y), str(max_x), str(min_y), master_raster, out_raster], stdout=subprocess.PIPE)
    
    return 0


def main(cmdargs):
    options = []
    if cmdargs.out_format is not None:
        options += ["-of", cmdargs.out_format]
    if cmdargs.cache is not None:
        options += ["--config", "GDAL_CACHEMAX", "{0:d}".format(cmdargs.cache)]
    subMaster(cmdargs.clip_raster, cmdargs.master_raster, cmdargs.out_raster, options=options)

if __name__ == "__main__":
    cmdargs = getCmdArgs()
    main(cmdargs)

