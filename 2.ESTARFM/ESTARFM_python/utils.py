import numpy as np
import gdal
import os


def read_raster(infile):
    gdal.PushErrorHandler('CPLQuietErrorHandler')
    gdal.UseExceptions()
    fp = gdal.Open(infile)
    cols = fp.RasterXSize
    rows = fp.RasterYSize
    nb = fp.RasterCount
    if nb == 1:
        band = fp.GetRasterBand(1)
        data = band.ReadAsArray()
        band.GetScale()
        band.GetOffset()
        band.GetNoDataValue()
    else:
        data = np.zeros([nb, rows, cols])
        for i in range(0, nb):
            band = fp.GetRasterBand(i+1)
            data[i, :, :] = band.ReadAsArray()
            band.GetScale()
            band.GetOffset()
            band.GetNoDataValue()
    return rows, cols, data


def writeimage(bands, path, in_ds):
    suffix = os.path.splitext(in_ds)[-1]
    in_ds = gdal.Open(in_ds)
    if bands is None or bands.__len__() == 0:
        return
    else:
        band1 = bands[0]
        img_width = band1.shape[1]
        img_height = band1.shape[0]
        num_bands = bands.__len__()

        if 'int8' in band1.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in band1.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float32

        if suffix == '.tif':
            driver = gdal.GetDriverByName("GTiff")
        elif suffix == "":
            driver = gdal.GetDriverByName("ENVI")

        dataset = driver.Create(path, img_width, img_height, num_bands, datatype)
        if dataset is not None:
            for i in range(bands.__len__()):
                dataset.GetRasterBand(i + 1).WriteArray(bands[i])
        geoTransform = in_ds.GetGeoTransform()
        dataset.SetGeoTransform(geoTransform)
        proj = in_ds.GetProjection()
        dataset.SetProjection(proj)
