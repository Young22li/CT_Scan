import numpy as np
from scipy.io import savemat

# 1. First, read the text file containing your 400×400 array
array_data = np.loadtxt('HumanBrain.txt')

# 2. Reshape if necessary (if the data is flattened)
# If your data is already in 400×400 format in the file, you can skip this step
# If it's flattened, you'll need:
# array_data = array_data.reshape(400, 400)

# 3. Optional: Normalize the data if needed
# If you want to normalize from 0-255 to match a specific range like 0-1:
# array_data = array_data / 255.0

# 4. Save as .mat file with a variable name
savemat('Brain.mat', {'brain': array_data})