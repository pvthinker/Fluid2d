import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.text as text
from netCDF4  import Dataset
from param import Param

plt.ion()
plt.close("all")

expname = "RB_128"

defaultparam = Param("default.xml")
fluid2d_datadir = os.path.expanduser(defaultparam.datadir)

exp_dir = f"{fluid2d_datadir}/{expname}"

diag_file = f"{expname}_diag.nc"
his_file = f"{expname}_his.nc"

with Dataset(f"{exp_dir}/{diag_file}") as nc:
    ke = nc.variables["ke"][:]
    time_diag = nc.variables["t"][:]
    pe = nc.variables["pe"][:]

with Dataset(f"{exp_dir}/{his_file}") as nc:
    attrs = nc.ncattrs()
    param ={name: nc.getncattr(name) for name in attrs}

    variables_names = [name for name in nc.variables]
    print(f"Variables in {his_file}:")
    print(f"  {variables_names}")

    x = nc.variables["x"][:]
    z = nc.variables["y"][:]
    time = nc.variables["t"][:]
    deltat = param['freq_his']
    print(f"{len(time)} snapshots from time from: 0 to {time[-1]}, every {deltat}")

    w = nc.variables["v"][:]
    b = nc.variables["buoyancy"][:]



kt = 100
fig, ax = plt.subplots(1,3,sharey=True,figsize=(12,5))
im0=ax[0].contourf(x,z,w[kt])
ax[0].set_xlabel("X")
ax[0].set_ylabel("Z")
ax[0].set_title(f"w / time={time[kt]:.1f}")
plt.colorbar(im0,ax=ax[0])

im1= ax[1].contourf(x,z,b[kt])
ax[1].set_xlabel("X")
ax[1].set_title(f"buoyancy")
plt.colorbar(im1,ax=ax[1])

ax[2].plot(np.mean(b[0],axis=-1),z,"k:")
ax[2].plot(np.mean(b[kt],axis=-1),z)
ax[2].set_xlabel("b")
ax[2].set_title(r"b / horizontal average ")


nz = len(z)
wsigma = np.zeros((nz,))
wplume = np.zeros((nz,))
for kz in range(nz):
    w_level = w[kt,kz]
    wsigma[kz] = np.std(w_level)
    threshold = 2*wsigma[kz]
    wp = w_level[np.abs(w_level)>threshold]
    wplume[kz] = np.mean(wp)
plt.figure(2)
plt.plot(wsigma,z,label="std")
plt.plot(wplume,z,label="plume")
plt.ylabel("z")
plt.xlabel("w")
plt.legend()
plt.grid()


H_plume = 0.4 # estimated from figure(2)
Lz = param["Ly"]
kz_plume = int((H_plume/Lz)*nz)

wmax = 0.3
nbins = 100
bins=np.linspace(-wmax,wmax,nbins+1)
plt.figure(3)
plt.hist(w[kt,:kz_plume].ravel(),bins,density=True)
plt.xlabel("w")
plt.ylabel("p.d.f.")
