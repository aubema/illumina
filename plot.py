#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from pytools import load_bin

x = np.arange(511)-255
xx,yy = np.meshgrid(x,x)
r = np.sqrt(xx**2+yy**2)

#b = load_bin("Hawaii_pcl.bin")
b = load_bin("single.bin")
c = load_bin("Hawaii_pcl.bin")

plt.figure(1)
plt.plot(r.flat,(c).flat,'.')
plt.plot(r.flat,(b).flat,'.')
plt.xscale('log')
plt.yscale('log')
plt.savefig("tata.png")

plt.figure(2)
plt.imshow(np.log10(c*r**1),cmap="inferno")
plt.xlim(204.5,305.5)
plt.ylim(305.5,204.5)
plt.colorbar()
plt.savefig("toto.png")

plt.figure(3)
plt.plot(r.flat,(c/b).flat,'.')
plt.savefig("ratio.png")
plt.show()
