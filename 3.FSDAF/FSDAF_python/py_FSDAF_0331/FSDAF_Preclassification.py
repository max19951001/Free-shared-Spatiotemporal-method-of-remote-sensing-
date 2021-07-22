from utils import read_raster, writeimage
import math
import numpy as np
import gdal
import os
import datetime
from tkinter import filedialog
import tkinter as tk
import yaml
import idlwrap
from scipy.interpolate import Rbf
import statsmodels.api as sm


def value_locate(refx, x):
    refx = np.array(refx)
    x = np.atleast_1d(x)
    loc = np.zeros(len(x), dtype='int')
    for i in range(len(x)):
        ix = x[i]
        ind = ((refx - ix) <= 0).nonzero()[0]
        if len(ind) == 0:
            loc[i] = -1
        else:
            loc[i] = ind[-1]
    return loc


# ******************************************************************************************************
#                          A new spatiotemporal data fusion model
#                          Using one pairs of fine and coarse images
#                          The program can be used for whole TM scene
#      Note: this version requires users to input pre-classified fine image at t1.
#      This version is appropriate for areas with complex land cover types and training samples
#      so users can first classify the fine image by supervised classifers like SVM.

# ******************************************************************************************************
# *******************************Set parameters and read input data*************************************

root = tk.Tk()
root.withdraw()

# please set the following parameters
f = open(filedialog.askopenfilename(title=u"Open the parameter settings file:"))
param = yaml.safe_load(f)
w = param['w']  # set the half window size, if 25, the window size is 25*2+1=51
num_similar_pixel = param['num_similar_pixel']  # set number of similar pixels
num_pure = param['num_pure']  # number of most purest coarse pixels in each class selected fro change calculation
DN_min = param['DN_min']  # set the range of DN value of the image,If byte, 0 and 255
DN_max = param['DN_max']
scale_factor = param['scale_factor']  # set the scale factor, it is integer=coarse resolution/fine resolution, e.g., 480/30=16
block_size = param['block_size']  # set the size of block, e.g., 20 means 20*20 coarse pixels, if process whole ETM scene, set 30~50
background = param['background']  # set the value of background pixels. 0 means that pixels will be considered as background if one of its bands= 0
background_band = param['background_band']  # which band with value = background indicating background pixels. Sometimes, background pixels have different values in different bands

# set path of a folder to store temporary files
temp_file = filedialog.askdirectory(title=u"Set the temporary folder")

# open the fine image of the first pair
path1 = filedialog.askopenfilename(title=u"open the fine image of the first pair:")
suffix = os.path.splitext(path1)[-1]
nl, ns, FileName1 = read_raster(path1)
orig_ns = ns
orig_nl = nl
fp = gdal.Open(path1)
nb = fp.RasterCount

patch_long = block_size*scale_factor

# divide the whole scene into blocks
n_nl = math.ceil(orig_nl / patch_long)
n_ns = math.ceil(orig_ns / patch_long)

ind_patch1 = np.zeros((n_nl * n_ns, 4), dtype=np.int)
ind_patch = np.zeros((n_nl * n_ns, 4), dtype=np.int)
location = np.zeros((n_nl * n_ns, 4), dtype=np.int)

for i_ns in range(0, n_ns):
    for i_nl in range(0, n_nl):
        ind_patch1[n_ns * i_nl + i_ns, 0] = i_ns * patch_long
        ind_patch[n_ns * i_nl + i_ns, 0] = np.max([0, ind_patch1[n_ns * i_nl + i_ns, 0] - scale_factor])
        location[n_ns * i_nl + i_ns, 0] = ind_patch1[n_ns * i_nl + i_ns, 0] - ind_patch[n_ns * i_nl + i_ns, 0]

        ind_patch1[n_ns * i_nl + i_ns, 1] = np.min([ns - 1, (i_ns + 1) * patch_long - 1])
        ind_patch[n_ns * i_nl + i_ns, 1] = np.min([ns - 1, ind_patch1[n_ns * i_nl + i_ns, 1] + scale_factor])
        location[n_ns * i_nl + i_ns, 1] = ind_patch1[n_ns * i_nl + i_ns, 1] - ind_patch1[n_ns * i_nl + i_ns, 0] + location[n_ns * i_nl + i_ns, 0]

        ind_patch1[n_ns * i_nl + i_ns, 2] = i_nl * patch_long
        ind_patch[n_ns * i_nl + i_ns, 2] = np.max([0, ind_patch1[n_ns * i_nl + i_ns, 2] - scale_factor])
        location[n_ns * i_nl + i_ns, 2] = ind_patch1[n_ns * i_nl + i_ns, 2] - ind_patch[n_ns * i_nl + i_ns, 2]

        ind_patch1[n_ns * i_nl + i_ns, 3] = np.min([nl - 1, (i_nl + 1) * patch_long - 1])
        ind_patch[n_ns * i_nl + i_ns, 3] = np.min([nl - 1, ind_patch1[n_ns * i_nl + i_ns, 3] + scale_factor])
        location[n_ns * i_nl + i_ns, 3] = ind_patch1[n_ns * i_nl + i_ns, 3] - ind_patch1[n_ns * i_nl + i_ns, 2] + location[n_ns * i_nl + i_ns, 2]

tempoutname = temp_file + '\\temp_F1'

for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName1[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# open the coarse image of the first pair
path2 = filedialog.askopenfilename(title=u"open the coarse image of the first pair:")
_, _, FileName2 = read_raster(path2)

tempoutname = temp_file + '\\temp_C1'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName2[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# open the coarse image of the prediction time
path3 = filedialog.askopenfilename(title=u"open the coarse image of the prediction time:")
_, _, FileName3 = read_raster(path3)

tempoutname = temp_file + '\\temp_C0'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName3[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# open the class image
path4 = filedialog.askopenfilename(title=u"open the class image of fine image in the 1st pair:")
_, _, FileName4 = read_raster(path4)

tempoutname = temp_file + '\\class'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName4[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# *******************************************************
# process each clock
# *******************************************************

starttime = datetime.datetime.now()  # the initial time of program running

print('there are total', n_nl*n_ns, 'blocks')

for isub in range(0, n_nl * n_ns):

    # open each block image

    FileName = temp_file + '\\temp_F1' + str(isub + 1) + suffix
    nl, ns, fine1 = read_raster(FileName)

    FileName = temp_file + '\\temp_C1' + str(isub + 1) + suffix
    _, _, coarse1 = read_raster(FileName)

    FileName = temp_file + '\\temp_C0' + str(isub + 1) + suffix
    _, _, coarse2 = read_raster(FileName)

    FileName = temp_file + '\\class' + str(isub + 1) + suffix
    _, _, L1_class0 = read_raster(FileName)

    num_class = int(np.max(L1_class0))
    # recode the classification map if the subset does not have all classes
    i_new_c = 0
    L1_class = np.zeros((nl, ns)).astype(int)
    for iclass in range(0, num_class):

        ind_ic = np.logical_and(L1_class0[0] == iclass + 1, fine1[background_band - 1, :, :] != background)
        num_ic = np.sum(ind_ic)

        if num_ic > 0:
            L1_class[ind_ic] = i_new_c + 1
            i_new_c = i_new_c + 1

    num_class = np.max(L1_class)

    if num_class > 0:  # do not process if the whole subset is background

        # correct extreme noise in fine1 because extreme values will affect the allowed data range
        for ib in range(0, nb):
            fine1_band = fine1[ib, :, :]
            fine1_band_1 = fine1_band.flatten()
            sortIndex = np.argsort(fine1_band_1, kind='mergesort')
            sortIndices = (idlwrap.findgen(float(ns) * nl + 1)) / (float(ns) * nl)
            Percentiles = [0.0001, 0.9999]
            dataIndices = value_locate(sortIndices, Percentiles)
            data_1_4 = fine1_band_1[sortIndex[dataIndices]]
            # correct too small values
            ind_small = np.logical_or(fine1[ib, :, :] <= data_1_4[0], fine1[ib, :, :] < DN_min)
            temp = fine1[ib, :, :]
            temp[ind_small] = np.min((fine1[ib, :, :])[np.logical_and(fine1[ib, :, :] > data_1_4[0], fine1[ib, :, :] >= DN_min)])
            fine1[ib, :, :] = temp
            # correct too large values
            ind_large = np.logical_or(fine1[ib, :, :] >= data_1_4[1], fine1[ib, :, :] > DN_max)
            temp = fine1[ib, :, :]
            temp[ind_large] = np.max((fine1[ib, :, :])[np.logical_and(fine1[ib, :, :] < data_1_4[1], fine1[ib, :, :] <= DN_max)])
            fine1[ib, :, :] = temp

        # get index image between coarse and fine resolutions
        ii = 0
        ns_c = int(np.floor(ns / scale_factor))
        nl_c = int(np.floor(nl / scale_factor))
        index_f = np.zeros((nl, ns)).astype(int)
        index_c = np.zeros((nl_c, ns_c)).astype(int)
        for i in range(0, ns_c):
            for j in range(0, nl_c):
                index_f[j * scale_factor:(j + 1) * scale_factor, i * scale_factor:(i + 1) * scale_factor] = ii
                index_c[j, i] = ii
                ii = ii + 1.0

        # col and row index
        row_ind = np.zeros((nl, ns)).astype(int)
        col_ind = np.zeros((nl, ns)).astype(int)
        for i in range(0, ns):
            col_ind[:, i] = i

        for i in range(0, nl):
            row_ind[i, :] = i

        # resample coarse image to coarse resolution
        fine_c1 = np.zeros((nb, nl_c, ns_c)).astype(float)
        coarse_c1 = np.zeros((nb, nl_c, ns_c)).astype(float)
        coarse_c2 = np.zeros((nb, nl_c, ns_c)).astype(float)
        row_c = np.zeros((nl_c, ns_c)).astype(float)
        col_c = np.zeros((nl_c, ns_c)).astype(float)
        for ic in range(0, ns_c):
            for jc in range(0, nl_c):
                ind_c = np.where(index_f == index_c[jc, ic])
                row_c[jc, ic] = np.mean(row_ind[ind_c])
                col_c[jc, ic] = np.mean(col_ind[ind_c])
                for ib in range(0, nb):
                    fine_c1[ib, jc, ic] = np.mean((fine1[ib, :, :])[ind_c])
                    coarse_c1[ib, jc, ic] = np.mean((coarse1[ib, :, :])[ind_c])
                    coarse_c2[ib, jc, ic] = np.mean((coarse2[ib, :, :])[ind_c])

        # step 2: get fracture of each class within each coarse pixel at t1
        Fraction1 = np.zeros((num_class, nl_c, ns_c)).astype(float)
        for ic in range(0, ns_c):
            for jc in range(0, nl_c):
                ind_c = np.where(index_f == index_c[jc, ic])
                num_c = int(int(np.size(ind_c)) / len(ind_c))
                L1_class_c = L1_class[ind_c]
                for iclass in range(0, num_class):
                    ind_ic = np.where(L1_class_c == iclass+1)
                    num_ic = int(int(np.size(ind_ic)) / len(ind_ic))
                    Fraction1[iclass, jc, ic] = num_ic / num_c

                if np.sum(Fraction1[:, jc, ic]) <= 0.999:  # avoild pixels have background fine pixels
                    Fraction1[:, jc, ic] = 0

        # get the heterogenity of each fine pixel
        het_index = np.zeros((nl, ns)).astype(float)
        scale_d = w

        for i in range(0, ns):
            for j in range(0, nl):
                # the window location
                ai = int(np.max([0, i - scale_d]))
                bi = int(np.min([ns - 1, i + scale_d]))
                aj = int(np.max([0, j - scale_d]))
                bj = int(np.min([nl - 1, j + scale_d]))
                class_t = L1_class[j, i]
                # select same-class pixels
                ind_same_class = np.where(L1_class[aj:bj+1, ai:bi+1] == class_t)
                num_sameclass = int(int(np.size(ind_same_class)) / len(ind_same_class))
                het_index[j, i] = float(num_sameclass) / ((bi-ai+1.0) * (bj-aj+1.0))

        # step 3: METHOD2:estimate average spectral change of each class using pixels without land cover change
        c_rate = np.zeros((nb, num_class)).astype(float)

        # allowed change value for each band
        min_allow = np.zeros(nb).astype(float)
        max_allow = np.zeros(nb).astype(float)
        for ib in range(0, nb):
            min_allow[ib] = np.min(coarse_c2[ib, :, :] - coarse_c1[ib, :, :]) - np.std(coarse_c2[ib, :, :] - coarse_c1[ib, :, :])
            max_allow[ib] = np.max(coarse_c2[ib, :, :] - coarse_c1[ib, :, :]) + np.std(coarse_c2[ib, :, :] - coarse_c1[ib, :, :])

        for ib in range(0, nb):
            x_matrix = np.zeros((num_pure * num_class, num_class)).astype(float)
            y_matrix = np.zeros((num_pure * num_class, 1)).astype(float)
            ii = 0
            for ic in range(0, num_class):
                order_s = np.argsort((Fraction1[ic, :, :]).flatten(), kind='mergesort')
                order = order_s[::-1]
                ind_f = np.where(Fraction1[ic, :, :] > 0.01)      # make sure all selected modis pixel contain class i
                num_f = int(int(np.size(ind_f)) / len(ind_f))
                num_pure1 = np.min([num_f, num_pure])
                change_c = (coarse_c2[ib, :, :].flatten())[order[0:num_pure1]] - (coarse_c1[ib, :, :].flatten())[order[0:num_pure1]]

                # only use 0.1-0.9 samples to exclude the land cover change pixels
                sortIndex = np.argsort(change_c, kind='mergesort')
                sortIndices = (idlwrap.findgen(float(num_pure1+1))) / num_pure1
                Percentiles = [0.1, 0.9]
                dataIndices = value_locate(sortIndices, Percentiles)
                data_1_4 = change_c[sortIndex[dataIndices]]
                ind_nonchange = np.logical_and(change_c >= data_1_4[0], change_c <= data_1_4[1])
                num_nonc = np.sum(ind_nonchange)
                if num_nonc > 0:
                    y_matrix[ii:ii+num_nonc, 0] = change_c[ind_nonchange]
                    for icc in range(0, num_class):
                        f_c = (Fraction1[icc, :, :].flatten())[order[0:num_pure1]]
                        x_matrix[ii:ii+num_nonc, icc] = f_c[ind_nonchange]
                    ii = ii + num_nonc
            x_matrix = x_matrix[0:ii, :]
            y_matrix = y_matrix[0:ii, 0]

            model = sm.OLS(y_matrix, x_matrix).fit()
            opt = model.params
            c_rate[ib, :] = opt

        # step4: predict L2 assuming no land cover change
        L2_1 = fine1.copy()
        for ic in range(1, num_class+1):
            ind_L1_class = np.where(L1_class == ic)
            for ib in range(0, nb):
                temp = L2_1[ib, :, :]
                temp[ind_L1_class] = (fine1[ib, :, :])[ind_L1_class] + c_rate[ib, ic-1]

        # resample L2_1 image to coarse resolution
        coarse_c2_p = np.zeros((nb, nl_c, ns_c)).astype(float)
        for ic in range(0, ns_c):
            for jc in range(0, nl_c):
                ind_c = np.where(index_f == index_c[jc, ic])
                for ib in range(0, nb):
                    coarse_c2_p[ib, jc, ic] = np.mean((L2_1[ib, :, :])[ind_c])

        # allowed minmum value for each band at t2
        min_allow = np.zeros(nb).astype(float)
        max_allow = np.zeros(nb).astype(float)
        for ib in range(0, nb):
            min_allow0 = np.min([np.min(coarse2[ib, :, :]), np.min(L2_1[ib, :, :])])
            min_allow[ib] = np.max([min_allow0, DN_min])
            max_allow0 = np.max([np.max(coarse2[ib, :, :]), np.max(L2_1[ib, :, :])])
            max_allow[ib] = np.min([max_allow0, DN_max])

        # step5: predict L2 using TPS
        L2_tps = np.zeros((nb, nl, ns)).astype(float)
        for ib in range(0, nb):
            rbf = Rbf(row_c.ravel(), col_c.ravel(), (coarse_c2[ib, :, :]).ravel(), function='multiquadric')
            tps = rbf(row_ind.ravel(), col_ind.ravel()).reshape([nl, ns])
            L2_tps[ib, :, :] = tps

        print('finish TPS prediction')

        # step 6: redistribute residual
        # change residual
        predict_change_c = coarse_c2_p - fine_c1     # predict change
        real_change_c = coarse_c2 - coarse_c1        # real change
        change_R = real_change_c - predict_change_c

        # redistribution residual
        change_21_c = np.zeros((nb, nl_c, ns_c)).astype(float)
        change_21 = np.zeros((nb, nl, ns)).astype(float)

        for ic in range(0, ns_c):
            for jc in range(0, nl_c):

                ind_c = np.where(index_f == index_c[jc, ic])
                num_ii = int(int(np.size(ind_c)) / len(ind_c))

                for ib in range(0, nb):
                    diff_change = change_R[ib, jc, ic]
                    w_change_tps = (L2_tps[ib, :, :])[ind_c] - (L2_1[ib, :, :])[ind_c]
                    if diff_change <= 0:
                        ind_noc = np.where(w_change_tps > 0)
                        num_noc = int(int(np.size(ind_noc)) / len(ind_noc))
                        if num_noc > 0:
                            w_change_tps[ind_noc] = 0
                    else:
                        ind_noc = np.where(w_change_tps < 0)
                        num_noc = int(int(np.size(ind_noc)) / len(ind_noc))
                        if num_noc > 0:
                            w_change_tps[ind_noc] = 0

                    w_change_tps = np.abs(w_change_tps)
                    w_unform = np.zeros(num_ii).astype(float)     # evenly distributing residuals to sub-pixels
                    w_unform[:] = np.abs(diff_change)

                    w_change = w_change_tps * het_index[ind_c] + w_unform*(1.0-het_index[ind_c]) + 0.000001  # combine these two weights
                    w_change = w_change / (np.mean(w_change))  # nomalize weight

                    # avoid extreme weights
                    ind_extrem = np.where(w_change > 10)
                    num_extrem = int(int(np.size(ind_extrem)) / len(ind_extrem))
                    if num_extrem > 0:
                        w_change[ind_extrem] = np.mean(w_change)
                    w_change = w_change / (np.mean(w_change))

                    # distribute residuals according to WEIGHT
                    temp = change_21[ib, :, :]
                    temp[ind_c] = w_change * diff_change
                    change_21[ib, :, :] = temp

        # second prediction: L1+change
        fine2_2 = L2_1 + change_21
        # correct abnormal detected change
        for ib in range(0, nb):
            temp = fine2_2[ib, :, :]
            ind_min = np.where(temp < min_allow[ib])
            num_min = int(int(np.size(ind_min)) / len(ind_min))
            if num_min > 0:
                temp[ind_min] = min_allow[ib]
            ind_max = np.where(temp > max_allow[ib])
            num_max = int(int(np.size(ind_max)) / len(ind_max))
            if num_max > 0:
                temp[ind_max] = max_allow[ib]
            fine2_2[ib, :, :] = temp

        change_21 = fine2_2 - fine1

    else:
        change_21 = fine1 - fine1

    change_21 = change_21[:, location[isub, 2]:location[isub, 3] + 1, location[isub, 0]:location[isub, 1] + 1]

    print('finish change prediction step ', isub+1, 'block')
    tempoutname1 = temp_file + '\\temp_change'
    Out_Name = tempoutname1 + str(isub + 1) + suffix
    fp = path1
    writeimage(change_21, Out_Name, fp)

# **************************mosaic all the change patch********************************
datalist = []
minx_list = []
maxX_list = []
minY_list = []
maxY_list = []

for isub in range(0, n_ns * n_nl):
    out_name = temp_file + '\\temp_change' + str(isub + 1) + suffix
    datalist.append(out_name)

    col1 = ind_patch1[isub, 0]
    col2 = ind_patch1[isub, 1]
    row1 = ind_patch1[isub, 2]
    row2 = ind_patch1[isub, 3]

    minx_list.append(col1)
    maxX_list.append(col2)
    minY_list.append(row1)
    maxY_list.append(row2)

minX = min(minx_list)
maxX = max(maxX_list)
minY = min(minY_list)
maxY = max(maxY_list)

xOffset_list = []
yOffset_list = []
i = 0
for data in datalist:
    xOffset = int(minx_list[i] - minX)
    yOffset = int(minY_list[i] - minY)
    xOffset_list.append(xOffset)
    yOffset_list.append(yOffset)
    i += 1

in_ds = gdal.Open(path1)
path = temp_file + "\\temp_change" + suffix
if suffix == '.tif':
    driver = gdal.GetDriverByName("GTiff")
elif suffix == "":
    driver = gdal.GetDriverByName("ENVI")
dataset = driver.Create(path, orig_ns, orig_nl, nb, gdal.GDT_Float32)

i = 0
for data in datalist:
    nl, ns, datavalue = read_raster(data)
    for j in range(0, nb):
        dataset.GetRasterBand(j + 1).WriteArray(datavalue[j], xOffset_list[i], yOffset_list[i])
    i += 1

geoTransform = in_ds.GetGeoTransform()
dataset.SetGeoTransform(geoTransform)
proj = in_ds.GetProjection()
dataset.SetProjection(proj)

del dataset


# *******************************step 5: final prediction*********************************

FileName6 = temp_file + "\\temp_change" + suffix
_, _, change = read_raster(FileName6)

tempoutname = temp_file + '\\temp_change'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = change[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

for isub in range(0, n_nl * n_ns):

    # open each block image

    FileName = temp_file + '\\temp_F1' + str(isub + 1) + suffix
    nl, ns, fine1 = read_raster(FileName)

    FileName = temp_file + '\\temp_C1' + str(isub + 1) + suffix
    _, _, coarse1 = read_raster(FileName)

    FileName = temp_file + '\\temp_C0' + str(isub + 1) + suffix
    _, _, coarse2 = read_raster(FileName)

    FileName = temp_file + '\\class' + str(isub + 1) + suffix
    _, _, L1_class = read_raster(FileName)

    FileName = temp_file + '\\temp_change' + str(isub + 1) + suffix
    _, _, change_21 = read_raster(FileName)

    # place the blended result
    fine2 = np.zeros([nb, location[isub, 3]-location[isub, 2]+1, location[isub, 1]-location[isub, 0]+1]).astype(float)

    # compute the distance of each pixel in the window with the target pixel (integrate window)
    D_temp1 = w - np.tile((idlwrap.indgen(w*2+1)), (int(w*2+1), 1))
    d1 = np.power(D_temp1, 2)
    D_temp2 = w - np.tile(idlwrap.indgen(1, w*2+1), (1, int(w*2+1)))
    d2 = np.power(D_temp2, 2)
    D_D_all = np.sqrt(d1 + d2)
    D_D_all = D_D_all.flatten()

    similar_th = np.zeros(nb).astype(float)
    for iband in range(0, nb):
        similar_th[iband] = np.std(fine1[iband, :, :]) * 2.0 / float(num_class)

    for i in range(location[isub, 0], location[isub, 1] + 1):    # retrieve each target pixel
        for j in range(location[isub, 2], location[isub, 3] + 1):
            if fine1[background_band - 1, j, i] != background:     # do not process the background

                ai = int(np.max([0, i - w]))
                bi = int(np.min([ns - 1, i + w]))
                aj = int(np.max([0, j - w]))
                bj = int(np.min([nl - 1, j + w]))

                ci = i - ai   # location of target pixel
                cj = j - aj

                col_wind = np.tile(idlwrap.indgen(bi-ai+1), (int(bj-aj+1), 1))
                row_wind = np.tile(idlwrap.indgen(1, bj-aj+1), (1, int(bi-ai+1)))

                # search similar pixels within window
                similar_cand = np.zeros((bi-ai+1)*(bj-aj+1)).astype(float)      # place the similarity measure between each pixel and the target pixel
                position_cand = np.zeros((bi-ai+1)*(bj-aj+1)).astype(int) + 1   # place the location of each similar pixel
                for ib in range(0, nb):
                    cand_band = np.zeros((bi-ai+1)*(bj-aj+1)).astype(int)
                    wind_fine = fine1[ib, aj:bj+1, ai:bi+1]
                    S_S = np.abs(wind_fine - wind_fine[cj, ci])
                    similar_cand = similar_cand + (S_S / (wind_fine[cj, ci] + 0.00000001)).flatten()
                    ind_cand = np.where(S_S.flatten() < similar_th[ib])
                    cand_band[ind_cand] = 1
                    position_cand = position_cand * cand_band

                indcand = np.where(position_cand != 0)
                number_cand0 = int(int(np.size(indcand)) / len(indcand))   # select similar pixel initially
                # place the spatial distance measure between each calidated pixesls and the target pixel
                if (bi-ai+1) * (bj-aj+1) < (w*2.0+1) * (w*2.0+1):   # not an integrate window
                    distance_cand = np.sqrt((ci-col_wind)**2 + (cj-row_wind)**2) + 0.00001
                else:
                    distance_cand = D_D_all    # integrate window

                # add a small weight from spatial distance to spectral distance to avoid all pixels in the window with same similarity. This happens in perfect simulated images
                combine_similar_cand = (similar_cand+0.00001)*(10.0+distance_cand/w).flatten()          # spatial distance has very small effect
                order_dis = np.argsort(combine_similar_cand[indcand], kind='mergesort')
                number_cand = np.min([number_cand0, num_similar_pixel])
                ind_same_class = (indcand[0])[order_dis[0:int(number_cand)]]           # select the N most similar samples

                # normalize these distances
                D_D_cand = (distance_cand.flatten())[ind_same_class]
                C_D = (1.0+D_D_cand/w) * (similar_cand[ind_same_class]+1.0)
                C_D = 1.0 / C_D
                weight = C_D / np.sum(C_D)

                for iband in range(0, nb):
                    # predict the value
                    change_21_win = change_21[iband, aj:bj+1, ai:bi+1]
                    change_cand = (change_21_win.flatten())[ind_same_class]
                    fine2[iband, j-location[isub, 2], i - location[isub, 0]] = fine1[iband, j, i] + np.sum(weight * change_cand)

                    # revise the abnormal prediction
                    if fine2[iband, j-location[isub, 2], i - location[isub, 0]] < DN_min:
                        another_predict = np.max([DN_min, fine1[iband, j, i]] + coarse2[iband, j, i] - coarse1[iband, j, i])
                        fine2[iband, j-location[isub, 2], i - location[isub, 0]] = np.min([DN_max, another_predict])

                    if fine2[iband, j-location[isub, 2], i - location[isub, 0]] > DN_max:
                        another_predict = np.min([DN_max, fine1[iband, j, i]] + coarse2[iband, j, i] - coarse1[iband, j, i])
                        fine2[iband, j - location[isub, 2], i - location[isub, 0]] = np.max([DN_min, another_predict])

    print('finish final prediction', str(isub+1), 'block')
    tempoutname1 = temp_file + '\\temp_blended'
    Out_Name = tempoutname1 + str(isub+1) + suffix
    fp = path1
    writeimage(fine2, Out_Name, fp)

endtime = datetime.datetime.now()
print('time used:', (endtime - starttime).seconds, 'seconds')

# # ***************************************************************
# # mosaic all the blended patch

datalist = []
minx_list = []
maxX_list = []
minY_list = []
maxY_list = []

for isub in range(0, n_ns * n_nl):
    out_name = temp_file + '\\temp_blended' + str(isub+1) + suffix
    datalist.append(out_name)

    col1 = ind_patch1[isub, 0]
    col2 = ind_patch1[isub, 1]
    row1 = ind_patch1[isub, 2]
    row2 = ind_patch1[isub, 3]

    minx_list.append(col1)
    maxX_list.append(col2)
    minY_list.append(row1)
    maxY_list.append(row2)

minX = min(minx_list)
maxX = max(maxX_list)
minY = min(minY_list)
maxY = max(maxY_list)

xOffset_list = []
yOffset_list = []
i = 0
for data in datalist:
    xOffset = int(minx_list[i] - minX)
    yOffset = int(minY_list[i] - minY)
    xOffset_list.append(xOffset)
    yOffset_list.append(yOffset)
    i += 1

in_ds = gdal.Open(path1)
path = os.path.splitext(path3)[0] + "_FSDAF_Preclassification" + suffix
if suffix == '.tif':
    driver = gdal.GetDriverByName("GTiff")
elif suffix == "":
    driver = gdal.GetDriverByName("ENVI")
dataset = driver.Create(path, orig_ns, orig_nl, nb, gdal.GDT_Float32)

i = 0
for data in datalist:
    nl, ns, datavalue = read_raster(data)
    for j in range(0, nb):
        dd = datavalue[j, :, :]
        dataset.GetRasterBand(j + 1).WriteArray(dd, xOffset_list[i], yOffset_list[i])
    i += 1

geoTransform = in_ds.GetGeoTransform()
dataset.SetGeoTransform(geoTransform)
proj = in_ds.GetProjection()
dataset.SetProjection(proj)
