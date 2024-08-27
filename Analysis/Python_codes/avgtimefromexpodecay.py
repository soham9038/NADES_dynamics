## Code to calculate average relaxation time from fitted relaxation parameters ##
## python3 code.py ##

from scipy import integrate
import numpy as np

# Define the integrand function
def integrand(t):  # Comment out/in the below form according to the requirement
#    return np.exp(-((t/0.16783273649087313)**0.7330608741172036)) # Stretched exponential form
    return np.exp(-t/0.1860507349621272) # Single exponential form
#    return 0.2736757540615317*np.exp(-t/0.024041062359653722) + (1-0.2736757540615317)*np.exp(-t/0.26052372403268526) # Biexponential form
#    return a*np.exp(-x/tau1) + b*np.exp(-x/tau2) + (1-a-b)*np.exp(-x/tau3) # Tri exponential form

# Perform the integration
result, error = integrate.quad(integrand, 0, np.inf)

print("Result of the integration:", result)
print("Estimated error:", error)

