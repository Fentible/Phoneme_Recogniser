import csv
import os
import sys

import jenkspy
import seaborn as sns
from matplotlib import pyplot as plt

rootdir = "../Jenks/"


def main(argv):
    print("Py here, processing Jenks with " + argv[0] + " breaks")
    for subdir, dirs, files in os.walk(rootdir):
        if not (subdir.endswith('/res') and subdir.endswith('/figs')):
            print("Processing :: " + subdir)
        for tmp in files:
            
            if tmp.endswith('.txt'):
                
                data = list()
                with open(subdir + "/" + tmp, 'r') as file:
                    reader = csv.reader(file)
                    for row in reader:
                        data += row

                name, ext = os.path.splitext(tmp)
                data = data[0:len(data) - 1]
                data = [float(i) for i in data]
                breaks = jenkspy.jenks_breaks(data, nb_class=int(argv[0]))
                data.sort()
                if not os.path.exists(subdir + "/res/"):
                    os.makedirs(subdir + "/res/")
                if not os.path.exists(subdir + "/figs/"):
                    os.makedirs(subdir + "/figs/")
                sns.kdeplot(data)
                with open(subdir + "/res/" + name + ".jen", 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile, delimiter='\n')
                    writer.writerow(breaks)

                for line in breaks:
                    plt.axvline(line, color='r', linestyle='--')

                plt.savefig(subdir + "/figs/" + name + ".png")
                del breaks, data
                plt.clf()


if __name__ == '__main__':
    main(sys.argv[1:])
