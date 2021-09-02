"""
Test of the cancellation of terms in the calculation of
tidal energy flux.

This will follow Mofjeld's notation.

F is proportional to the energy flux of the original signal, and
FF is proportional to the sum of the energy fluxes of the incident and
reflected waves.

RESULT: The two net fluxes are only equal for zero friction.  I think this
may be because pressure work is a nonlinear term and some part of the
two waves pressure work can leak into the other.

"""

import numpy as np

A0 = 1 + 0j
U0 = 1 + 0.2j

F = A0.real*U0.real + A0.imag*U0.imag

alpha = 1 / np.sqrt(1 + 0j)
Ap = (A0 + U0/alpha)/2
Am = (A0 - U0/alpha)/2
Up = alpha * Ap
Um = alpha * Am
FF = (Ap.real*Up.real + Ap.imag*Up.imag) - (Am.real*Um.real + Am.imag*Um.imag)
print('No friction:')
print('F = %0.1f, FF = %0.1f' % (F, FF))

alpha = 1 / np.sqrt(1 + 1j)
Ap = (A0 + U0/alpha)/2
Am = (A0 - U0/alpha)/2
Up = alpha * Ap
Um = alpha * Am
FF = (Ap.real*Up.real + Ap.imag*Up.imag) - (Am.real*Um.real + Am.imag*Um.imag)
print('\nOrder-1 friction:')
print('F = %0.1f, FF = %0.1f' % (F, FF))

