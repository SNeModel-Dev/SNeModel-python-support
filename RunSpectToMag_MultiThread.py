import concurrent.futures
import logging
import threading
import datetime
import time
import os
import subprocess

class ModelArgument:

    def __init__(self, alpha, rstar,mexp,eexp,filename,nzone):
        self.alpha = alpha
        self.rstar=rstar
        self.mexp=mexp
        self.eexp = eexp
        self.filename = filename
        self.nzone = nzone

def PingFuncSubProcWithName(a):
    logging.info("Thread %s: Started", a.filename)
    cmd = a.filename + "," + a.nzone + "\n"
    p = subprocess.Popen( ['SpectToMag.out'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True )
    p.communicate(cmd)
    logging.info("Thread %s: Completed", a.filename)


format = "%(asctime)s: %(message)s"
logging.basicConfig(format=format, level=logging.INFO)

logging.info("Staring Main Thread")

# Create an array of arguments to pass into the function
data = []
infile = open("infile.txt", "r")
values = infile.readlines()

for line in values:
    field = line.split(',')
    name = field[4]
    file=name.rstrip()
    model = open (file+"BB.dat", "r" )
    lines = model.readlines ()
    nzone=str(len(lines))
    arg = ModelArgument(field[0], field[1], field[2], field[3], field[4],nzone)
    data.append(arg)
# create an executor
executor = concurrent.futures.ThreadPoolExecutor(max_workers=10)

# this creates a thread for each of the arguments in data and passes it into the function
f = executor.map(PingFuncSubProcWithName, data)

logging.info("Waiting for threads to complete")

""" 
c This program is open source under the BSD-3 License.
c Redistribution and use in source and binary forms, with or without modification, are permitted
c provided that the following conditions are met:
c 1. Redistributions of source code must retain the above copyright notice, this list of c conditions andthe following disclaimer.
c 
c 2.Redistributions in binary form must reproduce the above copyright notice, this list of !c conditions
c and the following disclaimer in the documentation and/or other materials provided with the
c distribution.
c 
c 3.Neither the name of the copyright holder nor the names of its contributors may be used to !c endorse
c or promote products derived from this software without specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
c IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
c PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
c CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
c EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
c PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
c OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
c WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
c OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
c ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

c © (or copyright) 2020. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
c National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
c Department of Energy/National Nuclear Security Administration. All rights in the program are
c reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
c Security Administration. The Government is granted for itself and others acting on its behalf a
c nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
c derivative works, distribute copies to the public, perform publicly and display publicly, and to c permit
c others to do so.
"""