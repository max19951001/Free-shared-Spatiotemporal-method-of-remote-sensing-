import numpy as np



# Set the size of the moving window in which the search for similar pixels 
# is performed
windowSize = 31

# Set the path where the results should be stored
path = 'STARFM_demo/'

# Set to True if you want to decrease the sensitivity to the spectral distance
logWeight = False

# If more than one training pairs are used, set to True
temp = False

# The spatial impact factor is a constant defining the relative importance of 
# spatial distance (in meters)
# Take a smaller value of the spatial impact factor for heterogeneous regions 
# (e.g. A = 150 m)
spatImp = 150 

# increasing the number of classes limits the number of similar pixels
numberClass = 4 

# Set the uncertainty value for the fine resolution sensor
# https://earth.esa.int/web/sentinel/technical-guides/sentinel-2-msi/performance 
uncertaintyFineRes = 0.03

# Set the uncertainty value for the coarse resolution sensor
# https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-3-olci/validation
uncertaintyCoarseRes = 0.03

# Other global variables
mid_idx = (windowSize**2)//2
specUncertainty = np.sqrt(uncertaintyFineRes**2 + uncertaintyCoarseRes**2)
tempUncertainty = np.sqrt(2*uncertaintyCoarseRes**2)

# Set the size of the slices in which to divide the image
# This number should be multiple of the image height and not bigger than it
# Use bigger size for small images
sizeSlices = 150





