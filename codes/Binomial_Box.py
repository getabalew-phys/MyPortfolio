# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 14:30:47 2025

@author: getsh
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import factorial

# Define the total number of particles (trials) and the deviation range to inspect
num = 100
devt = 10

# Create an array of values from 0 to 99 for the x-axis
x = np.arange(0, num)

# Calculate the factorial for each x-value (Note: 'y' is calculated but unused later)
y = factorial(x)

def binomial(n):
    """Calculates the binomial distribution for n trials with a success probability of 0.9."""
    # Initialize an array of zeros to hold the probability distribution
    dist = np.zeros(num)
    
    # Loop through each possible number of successes (k)
    for k in range(n):
        # Apply the binomial formula: p^k * (1-p)^(n-k) * n! / ((n-k)! * k!)
        dist[k] = ((0.9)**k) * (0.1**(n-k)) * factorial(n) / (factorial(n-k) * factorial(k))
        
    return dist

# Generate the binomial distribution for 100 particles
z = binomial(num)

# Calculate a lower bound index (Note: 'b' is calculated but unused later)
b = int(num/2 - devt)

# Initialize a variable to store a partial sum
part_sum = 0

# Loop through a window of values starting from the center (index 50)
for i in range(2*devt):
    # Print the current partial sum and the corresponding probability value
    # Note: part_sum is currently never updated in this loop
    print(part_sum, z[int(num/2) + i])
    
# Calculate the ratio of the partial sum to the total sum of the distribution
p55 = part_sum / z.sum()

# Print the probability exactly at the center (index 50) and the final partial sum
print(z[int(num/2)], part_sum)

# Plot the normalized distribution (probabilities scaled so they sum to 1)
plt.plot(x, z/z.sum(), label="0.9", marker="o")

# Add descriptive labels and a title to the plot
plt.title("Binomial distribution")
plt.xlabel("Number of particle in the left box")
plt.ylabel("Probability density")

# Display the legend to show the label
plt.legend()