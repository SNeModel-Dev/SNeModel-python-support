#!/usr/bin/python
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
#files = glob.glob("2.00d00A1.00d13R*M1.00d51Espect.dat")

parser = argparse.ArgumentParser(description="Generate a plot from data directory")
parser.add_argument('--use_opacity', help="Process with numerical opacity", action="store_true")
parser.add_argument('--title', help="Graph title")
parser.add_argument('--label_file', help="Name of file containing labels")
parser.add_argument('--use_bb', help="Process BB file", action="store_true")
args = parser.parse_args()
labels_df = pd.read_csv(args.label_file, delim_whitespace=True, index_col=0)
labels_map = labels_df.to_dict();
if args.use_bb:
    files = glob.glob("*BB.dat")
else:
    files = glob.glob("*spect.dat")
    print(files)
data = []
for filename in sorted(files):
    data.append(pd.read_csv(filename, delim_whitespace=True).values)
#plt.yscale('log')
print(sorted(files))
for lightcurve, filename in zip(data, sorted(files)):
    label_text = str(labels_map["label"][filename]) + " Temp Coeff"
    print(type(lightcurve))
    if args.use_opacity:
       # plt.plot(lightcurve[:,0]/10**5, 10**lightcurve[:,-1], label=label_text)
        graph_data = [(i,j) for (i,j) in zip(lightcurve[:,1], lightcurve[:,-3]) if j > 40]
        graph_data = np.array(graph_data)
        plt.plot(graph_data[:,0]/(3600*24), graph_data[:, 1], label=label_text)
    elif args.use_bb:
        graph_data = [(i,j) for (i,j) in zip(lightcurve[:,0], lightcurve[:, 1]) if j > 42 and i < 40]
        graph_data = np.array(graph_data)
        plt.plot(graph_data[:,0], graph_data[:,1], label=label_text)
    else:
        graph_data = [(i,j) for (i,j) in zip(lightcurve[:,0], lightcurve[:,-2]) if j > 42]
        graph_data = np.array(graph_data)
        plt.plot(graph_data[:,0], graph_data[:,1], label=label_text)
    plt.title(args.title)
    plt.xlabel("Time (days)")
    plt.ylabel("Luminosity")
#plt.plot(data_np2[:,0], data_np2[:,2], "r")
#plt.plot(data_np3[:, 0], data_np3[:, 2], "g")
#plt.gca().invert_yaxis()
#plt.xlim(0.3, 10.0)
#plt.ylim(43, 47)
plt.legend()
plt.show()
#plt.savefig("nzones")
