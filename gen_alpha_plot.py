from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import sys

usage = "python gen_alpha_plot.py <input file> <plot title> <x axis label> <y axis label>"
if len(sys.argv) != 4:
    print "Wrong number of arguments"
    print usage
    sys.exit(1)

plot_title = sys.argv[2]
x_axis_label = sys.argv[3]
y_axis_label = sys.argv[4]

genomes = []
observations = []
expected = []


#Read file
with open(sys.argv[1], 'r') as f:
    for line in f:
        a = line.split()
        genomes.append(a[0][:-1])
        observations.append(int(a[1][:-1]))
        expected.append(int(a[2]))

pdf = PdfPages('%s_plot.pdf' % sys.argv[1])

# Example data
x_pos = np.arange(len(genomes))
plt.xticks(x_pos, genomes)
plt.bar(x_pos, observations, width=0.25, color='r')
plt.bar(x_pos+0.25, expected, width=0.25, color='b')
plt.ylabel(y_axis_label)
plt.xlabel(x_axis_label)
plt.title(plot_title)
pdf.savefig()
pdf.close()
