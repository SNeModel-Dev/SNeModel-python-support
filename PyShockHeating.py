#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt

class SuperNova:
    def __init__(self, alpha=2.0, rstar=1.0e13 , mexp=2.78e34, eexp=0.5e51, filename='test', grid_sz=2000,
                 ni_mass_gm=2.0e32):
        self.filename = filename
        self.grid_sz = grid_sz
        self.ni_mass_gm = ni_mass_gm
        self.alpha = alpha
        self.rstar = rstar
        self.mexp = mexp
        self.eexp = eexp
        # These are the matrix variables.
        self.m = np.zeros(grid_sz)
        self.energy = np.zeros(grid_sz)
        self.v = np.zeros(grid_sz)
        self.rho = np.zeros(grid_sz)
        self.temp = np.zeros(grid_sz)
        self.edot = np.zeros(grid_sz)
        self.kappa = np.full(grid_sz, 0.4)
        self.vshr = np.zeros(grid_sz)
        self.r = np.zeros(grid_sz+1)
        self.rhosh0 = np.zeros(grid_sz)
        self.vsh0 = np.zeros(grid_sz)
        self.tsh0 = np.zeros(grid_sz)
        self.rhosh = np.zeros(grid_sz)
        self.vsh = np.zeros(grid_sz)
        self.tsh = np.zeros(grid_sz)
        self.dr = np.zeros(grid_sz)
        self.factor = np.full(grid_sz, 1.0)
        self.esh = np.zeros(grid_sz)
        # These are the scalar variables
        self.rhoshw = 0.0
        self.rhow = 0.0
        self.brad = 0.0
        self.vshw = 0.0
        self.tshw = 0.0
        self.accel = 0.0
        self.msh = 0.0
        self.mw = 0.0
        self.a = 7.5656e-15
        self.pi43 = 4.0/3.0 * np.pi
        self.pi4 = np.pi * 4.0
        self.sigma = 5.6704e-5 * 4.0 * np.pi
        self.difffac = self.a * 4.0 * np.pi * 3e10 / 3.0
        self.efac = 4.78e10
        self.tauni = 7.605e5
        self.rho0 = 0.0
        # core is 1/10 of the explosion mass
        self.mcore = mexp/10.0
        self.rold = 0.0
        self.v0 = 0.0
        # had previously been set to .1 sm of Ni
        self.mni = ni_mass_gm
        self.etest = 0.0
        self.gam = 4.0/3.0
        self.vwind = 4.0e8
        self.mdot = 2.26e26
        self.jedge = 0.0
        self.mtot = 0.0
        self.rw = 0.0
        self.alpha2 = 0.0
        self.dt = 0.0
        self.time = 0.0
        self.rhowi = 0.0
        self.tau = 0.0
        self.lum = 0.0
        self.diff = 0.0
        self.rphot = 0.0
        self.tphot = 0.0
        self.lumt = 0.0
        self.den = 0.0
        self.lums = 0.0
        self.p = 0.0
        self.tsh1 = 0.0
        self.rho1 = 0.0
        self.vel = 0.0
        self.trecom = 1.0 * 1.16e4
        self.precom = 0.25
        # These probably don't need to exist but were defined early because of fortran
        self.i = 0
        self.j = 0
        self.jcore = 0
        self.jphot = 0.

        self.bb_file = open(filename+'BB.dat', 'w+')
        self.prop_file = open(filename+'prop.dat', 'w+')
        self.init_file = open(filename+'initial.dat', 'w+')



        # starting the initial conditions setup
        self.rho0 = self.mexp / self.pi4 / (rstar**(3.0-self.alpha) - 1e9**(3.0 - alpha))
        print("arange")
        print(np.arange(1.0, grid_sz+1.0).shape)
        self.r[1:] = self.rstar / self.grid_sz * np.arange(1.0, grid_sz+1.0)
        self.rho[0] = self.rho0 * (0.5*(self.r[1] + 1e9))**(-self.alpha)
        print(self.r[1:grid_sz+1])
        print(self.r[0:-1])
        self.rho[1:] = self.rho0 * (0.5 * (self.r[2:grid_sz+1] + self.r[1:grid_sz]))**(-self.alpha)

        self.r[0] = 0.0
        self.m = self.pi43 * self.rho * (self.r[1:]**3 - self.r[0:grid_sz]**3)
        self.mtot = np.sum(self.m)
        self.m = self.m * self.mexp / self.mtot

        self.cumsum_m = np.cumsum(self.m)
        self.core_idx = np.where(self.cumsum_m < self.mcore) #This is probably replacing jcore
        self.kappa[self.core_idx] = 0.2
        print("Core Indices")
        print(self.core_idx)
        print("Kappa")
        print(self.kappa)
        print("m")
        print(self.m)
        print("Stuff")
        print(self.mni - np.cumsum(self.m)[3])
        cumulative_mass = np.cumsum(self.m)
        remaining_ni = np.where((self.mni - cumulative_mass) > 0.0)
        self.edot[remaining_ni] = self.m[remaining_ni] * self.efac
        print("Remaining ni last")
        print(remaining_ni)
        if(self.mni - cumulative_mass[remaining_ni[0][-1]] > 0):
            self.edot[np.argmax(self.edot == 0.0)] = (self.mni - cumulative_mass[remaining_ni[0][-1]]) * self.efac
            self.mni = 0.0
        self.v0 = self.v0 + self.r[1:]**2
        print(self.mtot)
        print(self.v0)
        print(self.m[-1]/1.9e33, (self.r[grid_sz] - self.r[grid_sz-1]), self.trecom)

        v0 = np.sqrt(2.0*self.eexp/self.v0)
        self.etest = 0.0
        self.v = self.v0*self.r[1:]
        self.rhosh0 = ((self.gam+1)/(self.gam-1)) * self.rho
        self.vsh0 = self.v * (self.rho/self.rhosh0)**(-0.19)
        self.vshr = self.v * (self.rho/self.rhosh0)
        self.tsh0 = (((3/2)*(self.gam+1)/self.a) * self.rhosh0*self.vshr**2)**0.25
        self.energy = (self.a * self.tsh0**4 * self.pi43 * (self.r[1:]**3 - self.r[0:grid_sz]**3))
        output_data = np.stack((np.arange(self.grid_sz), self.r[1:], self.rho, self.tsh0, self.vsh0, self.m), axis=-1)
        print(output_data)
        fig, axs = plt.subplots(3, 2)
        print(output_data[:,0].shape)
        print(output_data[:,0])
        axs[0,0].plot(output_data[:,0], output_data[:,1], 'o-')
        axs[0,0].set_title("Radius vs grid position")
        axs[0,1].plot(output_data[:,1], output_data[:,2], 'o-')
        axs[0,1].set_title("Density vs radius")
        axs[1,0].plot(output_data[:,1], output_data[:,3], 'o-')
        axs[1,0].set_title("Tsh0 vs radius")
        axs[1,1].plot(output_data[:,1], output_data[:,4], 'o-')
        axs[1,1].set_title("VSH0 vs radius")
        axs[2,0].plot(output_data[:,1], output_data[:,5], 'o-')
        axs[2,0].set_title("Mass vs radius")
        plt.show()
        self.time = 0.0
        self.dt = 2000.0 / grid_sz
        self.lumt = 0.0
        self.jedge = grid_sz-1 #index of last element of the grid
        self.jphot = grid_sz-1 #index of last element of the grid
        self.lums = 0.0
        #self.rhosh = 0.0 this should be a list of zeros??
        self.rhow = 0.0
        self.vshw = 0.0
        self.accel = 1.0
        self.rho1 = self.rho[-1]
        self.msh = 0.0

        # beginning the time steps

        for i in np.arange(10000000):
            #print("Timestep: " + str(i))
            self.r[1:] = self.r[1:] + self.vsh0 * self.accel * self.dt
            self.rho = self.m / self.pi43 / (self.r[1:]**3 - self.r[0:grid_sz]**3)
           # energy_idx = np.where(self.energy >= 1.0)
           # print("energy")
           # print(self.energy)
          #  print(energy_idx)
           # print(self.energy[energy_idx].shape)
           # print(self.r[energy_idx+1].shape)
           # self.temp[energy_idx] = (self.energy[energy_idx] / self.pi43 /(self.r[energy_idx+1]**3 - self.r[energy_idx]**3)/self.a)**0.25
            for j in np.arange(self.grid_sz):
                #print("First grid loop step: " + str(j))
                if(self.energy[j] < 1.0):
                    self.temp[j] = 0.0
                else:
                    self.temp[j] = (self.energy[j] / self.pi43 / (self.r[j+1]**3 - self.r[j]**3)/self.a)**0.25

                if(j > self.core_idx[0][-1]):
                    if(self.temp[j] < self.trecom):
                        self.factor[j] = (self.temp[j]/self.trecom)**self.precom
                        self.factor[j] = max(1.0e-4, self.factor[j])
                    else:
                        self.factor[j] = 1.0

                if(j<= self.jphot -1):
                    self.diff = self.r[j+1]**2 * self.difffac / 0.5 / (self.kappa[j+1] * self.factor[j+1] + self.kappa[j] * self.factor[j]) / (self.rho[j+1] + self.rho[j]) * (self.temp[j+1]**4 - self.temp[j]**4) / (self.r[j+2] - self.r[j+1])
                    self.energy[j+1]= self.energy[j+1] - self.diff * self.dt
                    self.energy[j] = self.energy[j] + self.diff*self.dt

            for j in np.arange(self.grid_sz):
                #print("Energy grid loop step: " + str(j))
                if(self.edot[j] > 0):
                    self.energy[j] = self.energy[j] + self.edot[j] * np.exp(-self.time/self.tauni)*self.dt

            self.tau = 0.0
            self.lum = 0.0

            for j in np.arange(self.grid_sz):
                #print("Shell loop step: " + str(j))
                self.rhosh[j] = ((self.gam+1)/(self.gam-1))*self.rho[j]
                self.vsh[j] = (self.rho[j] * self.vsh0[j] * self.accel) / self.rhosh[j]
                self.tsh[j] = (((3/2) * (self.gam +1)/ self.a) * self.rhosh[j] * self.vsh[j]**2)**0.25
                self.esh[j] = ((self.a * self.tsh[j]**4) * self.pi43 * (self.r[j+1]**3 - self.r[j]**3))
            for j in range(self.grid_sz, 0, -1):
                #print("Countdown loop step: " + str(j))
                self.tau = self.tau + self.rho[j-1] * self.kappa[j-1] * self.factor[j-1] * (self.r[j] - self.r[j-1])
                self.rhow = self.mdot /(self.pi4 * self.r[j]**2 * self.vwind)
                self.mw = self.pi4 * self.r[j]**2 * self.tau / self.kappa[j-1]

                if(self.tau > 0.4):
                    if((i % 10000) == 0):
                        print('Val: ' + str(self.tau) + " " + str(i) + " " + str(self.jedge) + " " + str(self.r[j]/8.0e12) + " " + str(self.temp[j-1]))
                    if(self.tau > 10.0):
                        self.lum = 0
                        self.msh = self.rho[j-1] * self.pi43 * (self.rstar)**3
                        self.accel = ((self.mw * self.kappa[j-1] * self.tau) / self.pi4 / self.rstar**2)**0.19
                        self.rold = self.r[j]
                    else:
                        if (self.jedge == self.jphot):
                            self.msh = (self.rho[j-1] * self.pi43 * (self.r[j]**3 - self.rold**3))
                            self.rhosh[j-1] = ((self.gam + 1)/(self.gam -1)) * self.rho[j-1]
                            self.accel = ((self.msh * self.kappa[j-1] * self.tau) / self.pi4 / self.r[j]**2)**0.19
                        self.lum = (self.sigma * self.r[j]**2 * self.temp[j-1]**4 * np.exp(-self.tau))

                    self.lum = min(self.lum, self.energy[j-1] * 0.5 / self.dt)
                    self.lum = max(self.lum, 0.0)
                    self.lums = self.lums + self.lum
                    self.energy[j-1] = self.energy[j-1] - self.lum * self.dt
                    self.rphot = self.r[j]
                    self.tphot = self.temp[j-1]
                    self.den = self.rho[j-1]
                    self.jphot = j-1
                    self.vel = self.vsh0[j-1] * self.accel
                    break
                else:
                    self.jedge = j-1
                    self.msh = self.rhow * self.pi43 * (self.r[self.jedge+1])**3
                    self.lum = self.sigma * self.r[j]**2 * self.temp[j-1]**4
                    self.lum = min(self.lum, self.energy[j-1] * 0.99 / self.dt)
                    self.lum = max(self.lum, 0.0)
                    self.lums = self.lums + self.lum
                    self.energy[j-1] = self.energy[j-1] - self.lum * self.dt
                    if(self.energy[j-1] < 1.0e10):
                        self.energy[j-1] = 0.0

            self.lumt = self.lumt + self.lum
            self.lums = 0.0
            self.time = self.time + self.dt
            if((i % 10000) == 0):
                print("Lum " + str(self.tau) + " " + str(self.time/3600.0/24.0) + " " + str(self.lumt/10000.0) + " " + str(self.esh[grid_sz-1]))
                self.lumt = 0.0
            if(self.time > 200.0 * 3600.0 * 24.0):
                break



def main():
    nova = SuperNova()

if __name__=='__main__':
    main()

# This program is open source under the BSD-3 License.
# Redistribution and use in source and binary forms, with or without modification, are permitted
# provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of c conditions andthe following disclaimer.
#
# 2.Redistributions in binary form must reproduce the above copyright notice, this list of !c conditions
# and the following disclaimer in the documentation and/or other materials provided with the
# distribution.
#
# 3.Neither the name of the copyright holder nor the names of its contributors may be used to !c endorse
# or promote products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# © (or copyright) 2020. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
# National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
# Department of Energy/National Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
# Security Administration. The Government is granted for itself and others acting on its behalf a
# nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
# derivative works, distribute copies to the public, perform publicly and display publicly, and to c permit
# others to do so.
