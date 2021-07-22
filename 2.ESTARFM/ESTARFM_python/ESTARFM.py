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
import statsmodels.api as sm
from scipy.stats import f
from sklearn.linear_model import LinearRegression

# ******************************************************************************************************
#                            ESTARFM PROGRAM
#               Using two pairs of fine and coarse images
#         the program can be used for whole TM scene and VI index product

# ******************************************************************************************************
# *******************************Set parameters and read input data*************************************

root = tk.Tk()
root.withdraw()

# please set the following parameters
f = open(filedialog.askopenfilename(title=u"Open the parameter settings file:"))
param = yaml.safe_load(f)
w = param['w']  # set the half window size, if 25, the window size is 25*2+1=51 fine pixels
num_class = param['num_class']  # set the estimated number of classes, please set a larger value if blending images with very few bands
DN_min = param['DN_min']  # set the range of DN value of the image,If byte, 0 and 255
DN_max = param['DN_max']
background = param['background']  # set the value of background pixels. 0 means that pixels will be considered as background if one of its bands= 0
patch_long = param['patch_long']  # set the size of each block,if process whole ETM scene, set 500-1000

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

n_nl = math.ceil(orig_nl / patch_long)
n_ns = math.ceil(orig_ns / patch_long)

ind_patch = np.zeros((n_nl * n_ns, 4), dtype=np.int)

for i_ns in range(0, n_ns):
    for i_nl in range(0, n_nl):
        ind_patch[n_ns * i_nl + i_ns, 0] = i_ns * patch_long
        ind_patch[n_ns * i_nl + i_ns, 1] = np.min([ns - 1, (i_ns + 1) * patch_long - 1])
        ind_patch[n_ns * i_nl + i_ns, 2] = i_nl * patch_long
        ind_patch[n_ns * i_nl + i_ns, 3] = np.min([nl - 1, (i_nl + 1) * patch_long - 1])

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

# open the fine image of the second pair
path3 = filedialog.askopenfilename(title=u"open the fine image of the second pair:")
_, _, FileName3 = read_raster(path3)

tempoutname = temp_file + '\\temp_F2'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName3[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# open the coarse image of the second pair
path4 = filedialog.askopenfilename(title=u"open the coarse image of the second pair:")
_, _, FileName4 = read_raster(path4)

tempoutname = temp_file + '\\temp_C2'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName4[:, row1:row2 + 1, col1:col2 + 1]
    out_name = tempoutname + str(isub + 1) + suffix
    fp = path1
    writeimage(data, out_name, fp)

# open the coarse image of the prediction time
path5 = filedialog.askopenfilename(title=u"open the coarse image of the prediction time:")
_, _, FileName5 = read_raster(path5)

tempoutname = temp_file + '\\temp_C0'
for isub in range(0, n_nl * n_ns):
    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]
    data = FileName5[:, row1:row2 + 1, col1:col2 + 1]
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

    FileName = temp_file + '\\temp_F2' + str(isub + 1) + suffix
    _, _, fine2 = read_raster(FileName)

    FileName = temp_file + '\\temp_C2' + str(isub + 1) + suffix
    _, _, coarse2 = read_raster(FileName)

    FileName = temp_file + '\\temp_C0' + str(isub + 1) + suffix
    _, _, coarse0 = read_raster(FileName)

    fine0 = np.zeros((nb, nl, ns)).astype(float)    # place the blended result

    # row index of images
    row_index = np.zeros((nl, ns)).astype(int)
    for i in range(0, nl):
        row_index[i, :] = i

    # column index of images
    col_index = np.zeros((nl, ns)).astype(int)
    for i in range(0, ns):
        col_index[:, i] = i

    # compute the uncertainty,0.2% of each band is uncertain
    uncertain = (DN_max*0.002) * np.sqrt(2)

    # compute the threshold of similar pixel seeking
    similar_th = np.zeros((2, nb)).astype(float)
    for iband in range(0, nb):
        similar_th[0, iband] = np.std(fine1[iband, :, :] * 2.0 / num_class)
        similar_th[1, iband] = np.std(fine2[iband, :, :] * 2.0 / num_class)

    # compute the distance of each pixel in the window with the target pixel (integrate window)
    D_temp1 = w - np.tile((idlwrap.indgen(w*2+1)), (int(w*2+1), 1))
    d1 = np.power(D_temp1, 2)
    D_temp2 = w - np.tile(idlwrap.indgen(1, w*2+1), (1, int(w*2+1)))
    d2 = np.power(D_temp2, 2)
    D_D_all = 1.0 + np.sqrt(d1 + d2) / float(w)
    D_D_all = D_D_all.flatten()

    # find interaction of valid pixels of all input images: exclude missing pixels and background
    valid_index = np.zeros((nl, ns)).astype(int)
    ind_valid = np.where((fine1[0, :, :] != background) & (fine2[0, :, :] != background) & (coarse1[0, :, :] != background) \
        & (coarse2[0, :, :] != background) & (coarse0[0, :, :] != background))
    num_valid = int(int(np.size(ind_valid)) / len(ind_valid))
    if num_valid > 0:
        valid_index[ind_valid] = 1  # mark good pixels in all images

    for j in range(0, nl):    # retrieve each target pixel
        for i in range(0, ns):

            if valid_index[j, i] == 1:     # do not process the background

                ai = int(np.max([0, i - w]))
                bi = int(np.min([ns - 1, i + w]))
                aj = int(np.max([0, j - w]))
                bj = int(np.min([nl - 1, j + w]))

                ind_wind_valid = np.where((valid_index[aj:bj+1, ai:bi+1]).ravel() == 1)
                position_cand = idlwrap.intarr((bi-ai+1)*(bj-aj+1)) + 1    # place the location of each similar pixel
                row_wind = row_index[aj:bj+1, ai:bi+1]
                col_wind = col_index[aj:bj + 1, ai:bi + 1]

                # searching for similar pixels
                for ipair in [0, 1]:
                    for iband in range(0, nb):
                        cand_band = idlwrap.intarr((bi-ai+1)*(bj-aj+1))
                        if ipair == 0:
                            S_S = np.abs(fine1[iband, aj:bj+1, ai:bi+1] - fine1[iband, j, i])
                        elif ipair == 1:
                            S_S = np.abs(fine2[iband, aj:bj + 1, ai:bi + 1] - fine2[iband, j, i])
                        ind_cand = np.where(S_S.ravel() < similar_th[ipair, iband])
                        cand_band[ind_cand] = 1
                        position_cand = position_cand * cand_band

                cand_band = 0
                indcand = np.where((position_cand != 0) & ((valid_index[aj:bj+1, ai:bi+1]).ravel() == 1))
                number_cand = int(int(np.size(indcand)) / len(indcand))

                if number_cand > 5:    # compute the correlation
                    S_D_cand = np.zeros(number_cand).astype(float)
                    x_cand = (col_wind.ravel())[indcand]
                    y_cand = (row_wind.ravel())[indcand]
                    finecand = np.zeros((nb*2, number_cand)).astype(float)
                    coarsecand = np.zeros((nb*2, number_cand)).astype(float)

                    for ib in range(0, nb):
                        finecand[ib, :] = (fine1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        finecand[ib+nb, :] = (fine2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        coarsecand[ib, :] = (coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        coarsecand[ib+nb, :] = (coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]

                    if nb == 1:   # for images with one band, like NDVI
                        S_D_cand = 1.0 - 0.5*(nb.abs((finecand[0, :]-coarsecand[0, :]) / (finecand[0, :]+coarsecand[0, :])) +
                                              np.abs((finecand[1, :]-coarsecand[1, :]) / (finecand[1, :]+coarsecand[1, :])))
                    else:
                        # for images with multiple bands
                        sdx = np.std(finecand, axis=0, ddof=1)
                        sdy = np.std(coarsecand, axis=0, ddof=1)
                        meanx = np.mean(finecand, axis=0)
                        meany = np.mean(coarsecand, axis=0)

                        x_meanx = np.zeros((nb*2, number_cand)).astype(float)
                        y_meany = np.zeros((nb*2, number_cand)).astype(float)
                        for ib in range(0, nb*2):
                            x_meanx[ib, :] = finecand[ib, :] - meanx
                            y_meany[ib, :] = coarsecand[ib, :] - meany

                        S_D_cand = nb*2.0*np.mean(x_meanx*y_meany, axis=0) / (sdx*sdy) / (nb*2.0-1)

                    ind_nan = np.where(S_D_cand != S_D_cand)
                    num_nan = int(int(np.size(ind_nan)) / len(ind_nan))
                    if num_nan > 0:
                        S_D_cand[ind_nan] = 0.5    # correct the NaN value of correlation

                    D_D_cand = np.zeros(number_cand).astype(float)   # spatial distance
                    if (bi-ai+1)*(bj-aj+1) < (w*2.0+1)*(w*2.0+1):   # not an integrate window
                        D_D_cand = 1.0 + np.sqrt((i-x_cand)**2+(j-y_cand)**2) / w
                    else:
                        D_D_cand[0:number_cand] = D_D_all[indcand]      # integrate window

                    C_D = (1.0-S_D_cand) * D_D_cand + 0.0000001           # combined distance
                    weight = (1.0/C_D)/np.sum(1.0/C_D)

                    for ib in range(0, nb):   # compute V
                        fine_cand = np.hstack(((fine1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand], (fine2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]))
                        coarse_cand = np.hstack(((coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand], (coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]))
                        coarse_change = np.abs(np.mean((coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]) - np.mean((coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]))
                        if coarse_change >= DN_max*0.02:  # to ensure changes in coarse image large enough to obtain the conversion coefficient

                            X = coarse_cand.reshape(-1, 1)
                            Y = fine_cand.reshape(-1, 1)
                            XX = sm.add_constant(X)
                            model = sm.OLS(Y, XX).fit()
                            regress_result = model.params
                            sig = model.f_pvalue

                            # correct the result with no significancy or inconsistent change or too large value
                            if sig <= 0.05 and 0 < regress_result[1] <= 5:
                                V_cand = regress_result[1]
                            else:
                                V_cand = 1.0

                        else:
                            V_cand = 1.0

                        # compute the temporal weight
                        difc_pair1 = np.abs(np.mean((coarse0[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])-np.mean((coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid]))+0.01**5
                        difc_pair2 = np.abs(np.mean((coarse0[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])-np.mean((coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid]))+0.01**5
                        T_weight1 = (1.0/difc_pair1) / (1.0/difc_pair1+1.0/difc_pair2)
                        T_weight2 = (1.0/difc_pair2) / (1.0/difc_pair1+1.0/difc_pair2)

                        # predict from pair1
                        coase0_cand = (coarse0[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        coase1_cand = (coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        fine01 = fine1[ib, j, i] + np.sum(weight * V_cand * (coase0_cand-coase1_cand))
                        # predict from pair2
                        coase2_cand = (coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand]
                        fine02 = fine2[ib, j, i] + np.sum(weight * V_cand * (coase0_cand-coase2_cand))
                        # the final prediction
                        fine0[ib, j, i] = T_weight1 * fine01 + T_weight2 * fine02
                        # revise the abnormal prediction
                        if fine0[ib, j, i] <= DN_min or fine0[ib, j, i] >= DN_max:
                            fine01 = np.sum(weight*(fine1[ib, aj:bj+1, ai:bi+1]).ravel()[indcand])
                            fine02 = np.sum(weight*(fine2[ib, aj:bj+1, ai:bi+1]).ravel()[indcand])
                            fine0[ib, j, i] = T_weight1 * fine01 + T_weight2 * fine02

                else:  # for the case of no enough similar pixel selected

                    for ib in range(0, nb):
                        # compute the temporal weight
                        difc_pair1 = np.mean((coarse0[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])-np.mean((coarse1[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])+0.01**5
                        difc_pair1_a = np.abs(difc_pair1)
                        difc_pair2 = np.mean((coarse0[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])-np.mean((coarse2[ib, aj:bj+1, ai:bi+1]).ravel()[ind_wind_valid])+0.01**5
                        difc_pair2_a = np.abs(difc_pair2)
                        T_weight1 = (1.0/difc_pair1_a) / (1.0/difc_pair1_a+1.0/difc_pair2_a)
                        T_weight2 = (1.0/difc_pair2_a) / (1.0/difc_pair1_a+1.0/difc_pair2_a)
                        fine0[ib, j, i] = T_weight1 * (fine1[ib, j, i] + difc_pair1) + T_weight2 * (fine2[ib, j, i] + difc_pair2)

    print('finish ', str(isub + 1), 'block')
    tempoutname1 = temp_file + '\\temp_blended'
    Out_Name = tempoutname1 + str(isub + 1) + suffix
    fp = path1
    writeimage(fine0, Out_Name, fp)

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

    col1 = ind_patch[isub, 0]
    col2 = ind_patch[isub, 1]
    row1 = ind_patch[isub, 2]
    row2 = ind_patch[isub, 3]

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
path = os.path.splitext(path5)[0] + "_ESTARFM" + suffix
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
