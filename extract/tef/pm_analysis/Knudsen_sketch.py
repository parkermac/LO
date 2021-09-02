"""
Code to make a Knudsen diagram sketch.
"""

from lo_tools import plotting_functions as pfun
import matplotlib.pyplot as plt
import numpy as np

Sin = 30
Sout = np.linspace(0,Sin - .0001,100)

DS = Sin - Sout

Qin_Qr = Sout / DS

# PLOTTING
plt.close('all')
fs = 20
pfun.start_plot(fs=fs, figsize=(8,8))
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(Sout, DS, '-r', lw=4)
ax.plot(Sout, Qin_Qr, '-b', lw=4)
ax.axis([0, Sin, 0, Sin])
ax.set_xlabel(r'$S_{OUT}$')

ax.text(.05, .75, r'$\Delta S$', c='r', transform=ax.transAxes, size=1.5*fs)

ax.text(.9, .65, r'$Q_{IN} / Q_{R}$', c='b', transform=ax.transAxes, size=1.5*fs, ha='right')

ax.grid(True)

plt.show()

pfun.end_plot()