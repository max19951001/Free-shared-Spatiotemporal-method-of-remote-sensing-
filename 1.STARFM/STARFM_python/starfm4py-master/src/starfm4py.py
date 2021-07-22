# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 15:33:49 2018

@author: Nikolina Mileva
"""

import zarr
import numpy as np
import dask.array as da
from dask.diagnostics import ProgressBar
from parameters import (windowSize, logWeight, temp, mid_idx, numberClass, spatImp, 
                        specUncertainty, tempUncertainty, path)

 


# Flatten blocks inside a dask array            
def block2row(array, row, folder, block_id=None):
    if array.shape[0] == windowSize:
        # Parameters	
        name_string = str(block_id[0] + 1)
        m,n = array.shape
        u = m + 1 - windowSize
        v = n + 1 - windowSize

    	# Get Starting block indices
        start_idx = np.arange(u)[:,None]*n + np.arange(v)

    	# Get offsetted indices across the height and width of input array
        offset_idx = np.arange(windowSize)[:,None]*n + np.arange(windowSize)

    	# Get all actual indices & index into input array for final output
        flat_array = np.take(array,start_idx.ravel()[:,None] + offset_idx.ravel())

        # Save to (dask) array in .zarr format
        file_name = path + folder + name_string + 'r' + row + '.zarr'
        zarr.save(file_name, flat_array)
    
    return array


# Divide an image in overlapping blocks   
def partition(image, folder):
    image_da = da.from_array(image, chunks = (windowSize,image.shape[1]))
    image_pad = da.pad(image_da, windowSize//2, mode='constant')
    
    for i in range(0,windowSize):
        row = str(i)
        block_i = image_pad[i:,:]
        block_i_da = da.rechunk(block_i, chunks=(windowSize,image_pad.shape[1]))
        block_i_da.map_blocks(block2row, dtype=int, row=row, folder=folder).compute()


# Create a list of all files in the folder and stack them into one dask array
def da_stack(folder, shape):
    da_list = [] 
    full_path = path + folder
    max_blocks = shape[0]//windowSize + 1 
    
    for block in range(1,max_blocks + 1):
        for row in range(0,windowSize):
            name = str(block) + 'r' + str(row)
            full_name = full_path + name + '.zarr'
            try:
                da_array = da.from_zarr(full_name)
                da_list.append(da_array) 
            except Exception:
                continue
      
    return da.rechunk(da.concatenate(da_list, axis=0), chunks = (shape[1],windowSize**2))


# Calculate the spectral distance
def spectral_distance(fine_image_t0, coarse_image_t0):
    spec_diff = fine_image_t0 - coarse_image_t0
    spec_dist = 1/(abs(spec_diff) + 1.0)
    print ("Done spectral distance!", spec_dist)
    
    return spec_diff, spec_dist


# Calculate the temporal distance    
def temporal_distance(coarse_image_t0, coarse_image_t1):
    temp_diff = coarse_image_t1 - coarse_image_t0
    temp_dist = 1/(abs(temp_diff) + 1.0)
    print ("Done temporal distance!", temp_dist)
    
    return temp_diff, temp_dist
   

# Calculate the spatial distance    
def spatial_distance(array):
    coord = np.sqrt((np.mgrid[0:windowSize,0:windowSize]-windowSize//2)**2)
    spat_dist = np.sqrt(((0-coord[0])**2+(0-coord[1])**2))
    rel_spat_dist = spat_dist/spatImp + 1.0 # relative spatial distance
    rev_spat_dist = 1/rel_spat_dist # relative spatial distance reversed
    flat_spat_dist = np.ravel(rev_spat_dist)
    spat_dist_da = da.from_array(flat_spat_dist, chunks=flat_spat_dist.shape)
    print ("Done spatial distance!", spat_dist_da)
    
    return spat_dist_da


# Define the threshold used in the dynamic classification process
def similarity_threshold(fine_image_t0):#, st_dev):
    fine_image_t0 = da.where(fine_image_t0==0, np.nan, fine_image_t0)
    st_dev = da.nanstd(fine_image_t0, axis=1)# new
    sim_threshold = st_dev*2/numberClass 
    print ("Done similarity threshold!", sim_threshold)

    return sim_threshold


# Define the spectrally similar pixels within a moving window    
def similarity_pixels(fine_image_t0):
    sim_threshold = similarity_threshold(fine_image_t0)
    # possible to implement as sparse matrix
    similar_pixels = da.where(abs(fine_image_t0 - 
                                  fine_image_t0[:,mid_idx][:,None])
        <= sim_threshold[:,None], 1, 0) #sim_threshold[:,mid_idx][:,None], 1, 0) # new
    print ("Done similarity pixels!", similar_pixels)
   
    return similar_pixels
        

# Apply filtering on similar pixels 
def filtering(fine_image_t0, spec_dist, temp_dist, spec_diff, temp_diff):
    similar_pixels = similarity_pixels(fine_image_t0) 
    max_spec_dist = abs(spec_diff)[:,mid_idx][:,None] + specUncertainty + 1
    max_temp_dist = abs(temp_diff)[:,mid_idx][:,None] + tempUncertainty + 1  
    spec_filter = da.where(spec_dist>1.0/max_spec_dist, 1, 0)
    st_filter = spec_filter
    
    if temp == True:
        temp_filter = da.where(temp_dist>1.0/max_temp_dist, 1, 0)
        st_filter = spec_filter*temp_filter  
        
    similar_pixels_filtered = similar_pixels*st_filter
    print ("Done filtering!", similar_pixels_filtered)

    return similar_pixels_filtered # sim_pixels_sparse
    

# Calculate the combined distance
def comb_distance(spec_dist, temp_dist, spat_dist):
    if logWeight == True:
        spec_dist = da.log(spec_dist + 1)
        temp_dist = da.log(temp_dist + 1)
    
    comb_dist = da.rechunk(spec_dist*temp_dist*spat_dist, 
                           chunks=spec_dist.chunksize)
    print ("Done comb distance!", comb_dist)
    
    return comb_dist
    
        
# Calculate weights
def weighting(spec_dist, temp_dist, comb_dist, similar_pixels_filtered):
    # Assign max weight (1) when the temporal or spectral distance is zero
    zero_spec_dist = da.where(spec_dist[:,mid_idx][:,None] == 1, 1, 0)
    zero_temp_dist = da.where(temp_dist[:,mid_idx][:,None] == 1, 1, 0)
    zero_dist_mid = da.where((zero_spec_dist == 1), 
                             zero_spec_dist, zero_temp_dist)
    shape = da.subtract(spec_dist.shape,(0,1))
    zero_dist = da.zeros(shape, chunks=(spec_dist.shape[0],shape[1]))
    zero_dist = da.insert(zero_dist, [mid_idx], zero_dist_mid, axis=1)
    weights = da.where((da.sum(zero_dist,1)[:,None] == 1), zero_dist, comb_dist)
    
    # Calculate weights only for the filtered spectrally similar pixels
    weights_filt = weights*similar_pixels_filtered
    
    # Normalize weights
    norm_weights = da.rechunk(weights_filt/(da.sum(weights_filt,1)[:,None]), 
                              chunks = spec_dist.chunksize)
    
    print ("Done weighting!", norm_weights)
    
    return norm_weights


# Derive fine resolution reflectance for the day of prediction 
def predict(fine_image_t0, coarse_image_t0, coarse_image_t1, shape):
    spec = spectral_distance(fine_image_t0, coarse_image_t0)
    spec_diff = spec[0]
    spec_dist = spec[1]
    temp = temporal_distance(coarse_image_t0, coarse_image_t1)
    temp_diff = temp[0] 
    temp_dist = temp[1]
    spat_dist = spatial_distance(fine_image_t0)
    comb_dist = comb_distance(spec_dist, temp_dist, spat_dist)
    similar_pixels = filtering(fine_image_t0, spec_dist, temp_dist, spec_diff, 
                               temp_diff)
    weights = weighting(spec_dist, temp_dist, comb_dist, similar_pixels)    
    pred_refl = fine_image_t0 + temp_diff
    weighted_pred_refl = da.sum(pred_refl*weights, axis=1)   
    prediction = da.reshape(weighted_pred_refl, shape)
    print ("Done prediction!")
    
    return prediction
    
 
# Compute the results (converts the dask array to a numpy array)   
def starfm(fine_image_t0, coarse_image_t0, coarse_image_t1, profile, shape):
    print ('Processing...')
    prediction_da = predict(fine_image_t0, coarse_image_t0, coarse_image_t1, shape)
    with ProgressBar():
         prediction = prediction_da.compute()
    
    return prediction



























        
        
        
