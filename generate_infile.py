#!/usr/bin/python

import numpy as np

min_alpha = 2.0
max_alpha = 2.0
step_alpha = .5

min_radius = 1e13
max_radius = 1e13
step_radius = .1e13

min_mass = 1.59e34
max_mass = 3.97e34
step_mass = 0.5e34

min_expe = 10e50
max_expe = 5e51
step_expe = .5e51
with open("infile3.txt", 'w') as outfile:
    for alpha in np.linspace(min_alpha, max_alpha, ((max_alpha-min_alpha)/step_alpha)+1):
        for radius in np.linspace(min_radius, max_radius, ((max_radius-min_radius)/step_radius)+1):
            for mass in np.linspace(min_mass, max_mass, ((max_mass-min_mass)/step_mass)+1):
                for expe in np.linspace(min_expe, max_expe, ((max_expe-min_expe)/step_expe)+1):
                    line = "{:.2e}".format(alpha) + "," + "{:.2e}".format(radius) + "," + "{:.2e}".format(mass) + "," + "{:.2e}".format(expe)
                    filename = "{:.2e}".format(alpha) + "A" + "{:.2e}".format(radius) + "R" + "{:.2e}".format(mass) + "M" + "{:.2e}".format(expe) + "E\n"
                    filename = filename.replace("e+", "d")
                    line = line.replace("e+", "d")
                    line = line + "," + filename
                    outfile.write(line)
