import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from skimage.transform import radon, iradon
from scipy.fft import fft, ifft, fftfreq

# Load the .mat file
mat_data = loadmat('Brain.mat')


image = mat_data['brain']

# Ensure the image is 2D
if image.ndim > 2:
    image = image[:, :, 0]

# Set projection angles (from 0 to 180 degrees)
theta = np.linspace(0., 180., max(image.shape), endpoint=False)

# Radon transform to generate sinogram
sinogram = radon(image, theta=theta, circle=True)


# Manually implement a ramp filter to obtain the filtered sinogram
def apply_ramp_filter(sinogram):
    filtered_sinogram = np.zeros_like(sinogram)

    # Apply ramp filter to each column (projection at each angle) of the sinogram
    nrows, ncols = sinogram.shape
    for i in range(ncols):
        # Get projection at current angle
        projection = sinogram[:, i]

        # Perform FFT on the projection
        projection_fft = fft(projection)

        # Create ramp filter
        n = len(projection)
        f = fftfreq(n)  # Frequency range

        # Construct ramp filter (|f|)
        ramp_filter = np.abs(f)

        # Apply filter in frequency domain
        filtered_projection_fft = projection_fft * ramp_filter

        # Inverse FFT to get filtered projection
        filtered_projection = np.real(ifft(filtered_projection_fft))

        # Save to result matrix
        filtered_sinogram[:, i] = filtered_projection

    return filtered_sinogram


# Get sinogram
filtered_sinogram = apply_ramp_filter(sinogram)

# Use iradon function to reconstruct image with ramp filter
reconstruction = iradon(sinogram, theta=theta, filter_name='ramp', circle=True)

# Create canvas and subplots
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(20, 5))

# original image
ax1.set_title('Original Image')
ax1.imshow(image, cmap=plt.cm.Greys_r)
ax1.set_axis_off()

# original sinogram
ax2.set_title('Original Sinogram')
ax2.set_xlabel('Angle (degrees)')
ax2.set_ylabel('Projection position')
ax2.imshow(sinogram, cmap=plt.cm.Greys_r,
           extent=(0, 180, 0, sinogram.shape[0]), aspect='auto')

#  filtered sinogram
ax3.set_title('Filtered Sinogram\n(Ramp Filter)')
ax3.set_xlabel('Angle (degrees)')
ax3.set_ylabel('Projection position')
ax3.imshow(filtered_sinogram, cmap=plt.cm.Greys_r,
           extent=(0, 180, 0, filtered_sinogram.shape[0]), aspect='auto')

#  reconstructed image
ax4.set_title('Reconstructed Image\n(Ramp Filter)')
ax4.imshow(reconstruction, cmap=plt.cm.Greys_r)
ax4.set_axis_off()



plt.tight_layout()
plt.show()
