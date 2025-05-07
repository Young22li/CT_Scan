## REFERENCE
# https://en.wikipedia.org/wiki/Algebraic_reconstruction_technique

## ART Equation
# x^(k+1) = x^k + lambda * AT(b - A(x))/ATA

##
import numpy as np
import matplotlib.pyplot as plt


def ART(A, AT, b, x, mu=1e0, niter=1e2, bpos=True):
    ATA = AT(A(np.ones_like(x)))

    for i in range(int(niter)):

        x = x + np.divide(mu * AT(b - A(x)), ATA)

        if bpos:
            x[x < 0] = 0

        # plt.imshow(x, cmap='gray')
        # plt.title("%d / %d" % (i + 1, niter))
        # plt.pause(1)
        # plt.close()

    return x

# def ART(A, AT, b, x, mu=1e0, niter=100, bpos=True):
#     # Store original input as starting point
#     first_iteration = None
#
#     ATA = AT(A(np.ones_like(x)))
#
#     for i in range(int(niter)):
#         x = x + np.divide(mu * AT(b - A(x)), ATA)
#
#         if bpos:
#             x[x < 0] = 0
#
#         # Store the result after the first iteration
#         if i == 0:
#             first_iteration = np.copy(x)
#
#         # Optional: keep the visualization if you want to see progress
#         # plt.imshow(x, cmap='gray')
#         # plt.title("%d / %d" % (i + 1, niter))
#         # plt.pause(1)
#         # plt.close()
#
#     # Return both first iteration and final result
#     return first_iteration, x