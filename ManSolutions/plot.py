#!/bin/python
import matplotlib.pyplot as plt
import math
import itertools

data = open("solveroutput","r")

lines = data.readlines()

d = lines.index(" \n")

x = [float(l) for l in lines[:d]]
y = [float(l) for l in lines[d+1:]]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#ax.plot(x,y, label='$\epsilon$')
ax.plot(x,[(math.exp(nx)-ny)**2 for nx,ny in itertools.izip(x,y)],label='$\epsilon$')
#ax.plot(x,[n**2 *2e-1 for n in x],label='$O(h^2)$')
#ax.plot(x,[n *1e-2 for n in x],label='$O(h)$')
ax.set_yscale('log')
#ax.set_xscale('log')
#ax.set_xlabel('$h$ [m]')
#ax.set_ylabel('Error [K]')
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$\epsilon$ [K]')
#ax.legend()
plt.savefig('fig.png')
