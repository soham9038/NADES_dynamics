## Code to calculate fit relaxation curve with stretched exponential decay form ##
## python3 code.py ##

import numpy as np
from scipy.optimize import curve_fit

# Define the function to be fitted
def func(x, tau, beta):
    return np.exp(-(x/tau)**beta)

# Read in the data from a text file
data = np.loadtxt('data.txt') # 'data.txt' contains two columns, first column- time, second column- normalized acf
x_data = data[:, 0]
y_data = data[:, 1]

# Use curve_fit to find the optimal values of tau, and beta
initial_guess = [0.05, 0.25] # initial guess for tau, and beta
popt, pcov = curve_fit(func, x_data, y_data, p0=initial_guess)

# Write the results to a separate text file
with open('results.txt', 'w') as f:
    f.write(f"tau = {popt[0]}\n")
    f.write(f"beta = {popt[1]}\n")
