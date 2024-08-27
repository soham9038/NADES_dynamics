## Code to calcuate the normalized autocorrelation function ##
## python3 code.py ##

import numpy as np

# Load the data from the input text file
data = np.loadtxt('data.txt') # 'data.txt' contains a single column which contains the physical property of interest

# Subtract the mean from the data
data -= np.mean(data)

# Calculate the autocorrelation function
acf = np.correlate(data, data, mode='full')

# Calculate the variance of the data
variance = np.var(data)

# Normalize the autocorrelation function
normalized_acf = acf / (len(data) * variance)

# Select the last half of the normalized autocorrelation function
last_half_normalized_acf = normalized_acf[len(normalized_acf) // 2 :]

# Save the last half of the normalized autocorrelation function to the output text file
np.savetxt('result.txt', last_half_normalized_acf)
