## Code to calculate normalized acf of unitvector ##
## python3 code.py ##

import numpy as np

# Read data from 'unitvector.txt'
data = np.loadtxt('data.txt')

# Calculate the dot product between each pair of unit vectors
dot_products = np.dot(data, data.T)

# Calculate the normalized autocorrelation function
num_frames = len(data)
autocorrelation = []
for lag in range(num_frames):
    acf_lag = np.mean(dot_products.diagonal(lag)) / np.mean(dot_products.diagonal(0))
    autocorrelation.append(acf_lag)

# Write autocorrelation to a new text file 'autocorrelation.txt'
np.savetxt('autocorrelation.txt', autocorrelation, fmt='%.6f')
