## REFERENCE
# https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

## ART Equation
# x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon
from scipy.io import loadmat
from scipy.stats import poisson
from skimage.metrics import mean_squared_error as mse
from skimage.metrics import peak_signal_noise_ratio as psnr
from skimage.metrics import structural_similarity as ssim
# from skimage.metrics import mean_squared_error as mse
# from skimage.metrics import peak_signal_noise_ratio as psnr
# from skimage.metrics import structural_similarity as ssim
from ART import ART

## SYSTEM SETTING
N = 400
ANG = 180
VIEW = 360
THETA = np.linspace(0, ANG, VIEW + 1)
THETA = THETA[:-1]

X = loadmat('Brain.mat')
x = X['brain']

if x.ndim > 2:
    x = x[:, :, 0]  # 如果是RGB图像，只取第一个通道

# 设置投影角度（从0度到180度，间隔为1度）
#THETA = np.linspace(0., 180., max(x.shape), endpoint=False)

A = lambda x: radon(x, THETA, circle=False).astype(np.float32)
AT = lambda y: iradon(y, THETA, circle=False, filter_name=None, output_size=N).astype(np.float32)/(np.pi/(2*len(THETA)))
AINV = lambda y: iradon(y, THETA, circle=False, output_size=N).astype(np.float32)


## DATA GENERATION


p = A(x)
x_full = AINV(p)

## LOW-DOSE SINOGRAM GENERATION
# i0 = 5e4
# pn = np.exp(-p)
# pn = i0*pn
# pn = poisson.rvs(pn)
# pn[pn < 1] = 1
# pn = -np.log(pn/i0)
# pn[pn < 0] = 0

y = p

## Algebraic Reconstruction Technique (ART) INITIALIZATION
x_low = AINV(y)
x0 = np.zeros_like(x)
mu = 1e0
niter = 200
bpos = True

x_art = ART(A, AT, y, x0, mu, niter, bpos)

## CALCULATE QUANTIFICATION FACTOR
x_low[x_low < 0] = 0
x_art[x_art < 0] = 0
nor = np.amax(x)

mse_x_low = mse(x/nor, x_low/nor)
mse_x_art = mse(x/nor, x_art/nor)

psnr_x_low = psnr(x/nor, x_low/nor, data_range=1.0)
psnr_x_art = psnr(x/nor, x_art/nor, data_range=1.0)

ssim_x_low = ssim(x_low/nor, x/nor, data_range=1.0)
ssim_x_art = ssim(x_art/nor, x/nor, data_range=1.0)

## DISPLAY
# wndImg = [0, 0.03]
# wndPrj = [0, 6]
#
# plt.subplot(241)
# plt.imshow(x, cmap=plt.cm.Greys_r)
# plt.axis('off')
# plt.axis('image')
# plt.title('Ground truth')
#
# plt.subplot(242)
# plt.imshow(x_full, cmap=plt.cm.Greys_r)
# plt.axis('off')
# plt.axis('image')
# plt.title('full-dose')
#
# plt.subplot(243)
# plt.imshow(x_low, cmap=plt.cm.Greys_r)
# plt.axis('off')
# plt.axis('image')
# plt.title('low-dose\nMSE: %.4f\nPSNR: %.4f\nSSIM: %.4f' % (mse_x_low, psnr_x_low, ssim_x_low))
#
# plt.subplot(244)
# plt.imshow(x_art, cmap=plt.cm.Greys_r)
# plt.axis('off')
# plt.axis('image')
# plt.title('ART\nMSE: %.4f\nPSNR: %.4f\nSSIM: %.4f' % (mse_x_art, psnr_x_art, ssim_x_art))
#
# plt.subplot(246)
# plt.imshow(p, cmap=plt.cm.Greys_r)
# plt.title('full-dose\n(VIEW: %d)' % VIEW)
# plt.xlabel('View')
# plt.ylabel('Detector')
#
# plt.subplot(247)
# plt.imshow(y, cmap=plt.cm.Greys_r)
# plt.title('low-dose\n(VIEW: %d)' % VIEW)
# plt.xlabel('View')
# plt.ylabel('Detector')
#
# plt.subplot(248)
# plt.imshow(y - p, cmap=plt.cm.Greys_r)
# plt.title('full-dose - low-dose\n(Poisson noise)')
# plt.xlabel('View')
# plt.ylabel('Detector')
#
# plt.show()

# -- Plot ---------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(12, 4))

# Original
axes[0].imshow(x, cmap='gray')
axes[0].axis('off')
axes[0].set_title('Original')

# Sinogram
axes[1].imshow(p, cmap='gray', aspect='auto')
axes[1].set_xlabel('View')
axes[1].set_ylabel('Detector')
axes[1].set_title('Sinogram (VIEW=%d)' % VIEW)

# ART Reconstruction
axes[2].imshow(x_art, cmap='gray')
axes[2].axis('off')
axes[2].set_title('ART\nMSE: %.4f\nPSNR: %.4f\nSSIM: %.4f' % (mse_x_art, psnr_x_art, ssim_x_art))

plt.tight_layout()
plt.show()
