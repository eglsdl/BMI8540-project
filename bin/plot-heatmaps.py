#!/usr/bin/python3

# This program reads a signal intensity matrix file (generated
# with deeptools) and plots a centered heatmap (also known as a
# tornado plot). Arguments are input file (has to be an unzipped
# matrix file in tab format) and output file name (png).

# Importing relevant packages
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Check if correct amount of arguments is provided and throw an
# error if not. Assign input to variables.
if (len(sys.argv) < 3):
    print("Error: Not enough options given. The program takes " +
          "a signal intensity matrix file name and output " +
          "filename as input.")
    sys.exit()

if (len(sys.argv) > 3):
    print("Error: Too many options given. The program takes " +
          "a signal intensity matrix file name and output " +
          "filename as input.")
    sys.exit()

matrix_fn = sys.argv[1]
output_fn = sys.argv[2]

# Read input matrix file into a pandas data frame. Option
# skiprows is used to avoid rows with metadata and might need
# to be adjusted for a specific sample.
matrix = np.loadtxt(matrix_fn, delimiter='\t', skiprows=3)
#df = pd.read_csv(matrix_fn, sep='\t', comment='#', header=None, skiprows=3)
#matrix = df.values

# Through trial and error I have noticed that the data has some
# outliers that affect the scale of the heatmap. In order to
# avoid this issue, I chose to remove 2% of lowest and highest
# values from my plot by setting vmin and vmax values.
low = np.percentile(matrix, 2)
high = np.percentile(matrix, 98)

# In order to observe the nice tornado shape in the heatmap, the
# regions of interest need to be sorted based on average signal
# strength. To do that, I determine mean signal value in each
# row of my matrix, sort them in a decreasing manner, obtain
# indexes of this sorting and apply it to the whole matrix.
center_means = matrix.mean(axis=1)
sorted_idx = np.argsort(center_means)[::-1]
sorted_matrix = matrix[sorted_idx, :]

# Plotting the figue
plt.figure(figsize=(8, 10))
plt.imshow(sorted_matrix, aspect='auto', cmap='RdBu_r',
           interpolation='none', vmin=low, vmax=high)
plt.colorbar(label='Signal')
plt.xticks([0, 100, 200], ['-1000', '0', '1000'])
plt.xlabel('Distance from RNA Pol II peak center, kb')
plt.ylabel('Regions')
plt.title('Heatmap')
plt.tight_layout()
plt.savefig(output_fn, dpi=300)
plt.show()
