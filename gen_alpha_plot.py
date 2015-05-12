import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import sys


def make_bar_plot(expected, computed, plot_title, filename, x_axis_label='genomes', y_axis_label='alphas'):
    """
	takes two arrays of same size
    """
    if len(expected) != len(computed):
	print "error, arrays should be same size"
	return 


    genomes = ["G"+str(i) for i in xrange(len(expected))]

    # Example data
    x_pos = np.arange(len(genomes))
    plt.xticks(x_pos, genomes)
    plt.bar(x_pos, computed, width=0.25, color='r')
    plt.bar(x_pos+0.25, expected, width=0.25, color='b')
    plt.ylabel(y_axis_label)
    plt.xlabel(x_axis_label)
    plt.title(plot_title)
    plt.savefig(filename)
    plt.clf()

def make_histogram(expected, computed, plot_title, x_axis_label, y_axis_label):
    def to_percent(y, position):
	# Ignore the passed in position. This has the effect of scaling the default
	# tick locations.
	s = str(100 * y)

	# The percent symbol needs escaping in latex
	if matplotlib.rcParams['text.usetex'] == True:
	    return s + r'$\%$'
	else:
	    return s + '%'


    if len(expected) != len(computed):
	print "expected and computed need to be the same size"
	return
    # x = | expected - computed | 
    x = [abs(expected[i]-computed[i]) for i in xrange(len(expected))]

# Make a normed histogram. It'll be multiplied by 100 later.
    plt.hist(x, bins=5, normed=True)

# Create the formatter using the function to_percent. This multiplies all the
# default labels by 100, making them all percentages
    formatter = FuncFormatter(to_percent)

# Set the formatter
    plt.gca().yaxis.set_major_formatter(formatter)
    plt.xlabel(x_axis_label)
    plt.ylabel(y_axis_label)
    plt.title(plot_title)
    pdf = PdfPages('%s_plot.pdf' % plot_title)
    pdf.savefig()
    pdf.close()

if __name__ == "__main__":
    make_histogram([np.random.randn(20)], [np.random.randn(20)], "title", "x", "y")
