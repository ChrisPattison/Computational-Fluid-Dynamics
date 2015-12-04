#!/bin/python
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import math
import itertools

#Python 3.5

suffix = ""
mpl.rcParams['axes.titlesize'] = 26
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20
mpl.rcParams['legend.fontsize'] = 24

data = open("output"+suffix,"r")
spread = 1

lines = data.readlines()

x = [float(l) for l in lines[0].split()]
y = [float(l) for l in lines[1].split()]

base = 3
u = np.array([[float(f) for f in l.split()] for l in lines[base:base+int(lines[base-1])]])
base = base + int(lines[base-1]) + 1
v = np.array([[float(f) for f in l.split()] for l in lines[base:base+int(lines[base-1])]])
base = base + int(lines[base-1]) + 1
p = np.array([[float(f) for f in l.split()] for l in lines[base:base+int(lines[base-1])]])

x = np.array([x for i in range(u.shape[1])])
y = np.array([y for i in range(u.shape[0])]).transpose()

gs = gridspec.GridSpec(1,3)

fig = plt.figure(figsize=(45,10))
ax = fig.add_subplot(gs[:,0])
cm = ax.pcolormesh(x,y,np.sqrt(u**2+v**2))
sp = ax.streamplot(x[::spread, ::spread], y[::spread, ::spread], u[::spread, ::spread], v[::spread, ::spread], color='0.5')##DD0000
fig.colorbar(cm, ax=ax, cmap = plt.get_cmap('jet'))
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$y$ [m]')
ax.set_title('Velocity')
#plt.savefig('u'+suffix+'.png')

#fig = plt.figure()
ax = fig.add_subplot(gs[:,1])
cm = ax.pcolormesh(x,y,p)
#co = ax.contour(x,y,p,40)
fig.colorbar(cm, ax=ax, cmap = plt.get_cmap('jet'))
ax.set_xlabel('$x$ [m]')
ax.set_ylabel('$y$ [m]')
ax.set_title('Pressure')
#plt.savefig('p'+suffix+'.png')

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#cm = ax.pcolormesh(x,y,p)
#sp = ax.streamplot(x[::spread, ::spread], y[::spread, ::spread], u[::spread, ::spread], v[::spread, ::spread])
#ax.set_xlabel('$x$ [m]')
#ax.set_ylabel('$y$ [m]')
#plt.savefig('up.png')

data = open("residual" + suffix,"r")
lines = data.readlines()[2:-1]
k = [float(l.split()[0]) for l in lines]
uresid = [float(l.split()[1]) for l in lines]
presid = [float(l.split()[2]) for l in lines]
continuity = [float(l.split()[3]) for l in lines]
#fig = plt.figure()
ax = fig.add_subplot(gs[:,2])
ax.set_yscale('log')
ax.plot(k,uresid,label="$U$ Residual")
ax.plot(k,presid,label="$P$ Residual")
ax.plot(k,continuity,label="Continuity")
ax.set_xlabel('Iterations $k$')
ax.set_title('Residual')
ax.legend()

plt.savefig('solution'+suffix+'.png',bbox_inches='tight')
