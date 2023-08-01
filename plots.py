import sys
import matplotlib.pyplot as plt
import os

def parse(filename):
    T = []
    COUNTS = []
    SUM_FREQ = []
    AVG_FREQ = []
    AVG_LEN = []
    for line in open(filename):
        tokens = line.split(',')
        COUNTS.append(int(tokens[1]))
        SUM_FREQ.append(int(tokens[2]))
        AVG_FREQ.append(float(tokens[3])) # GB
        AVG_LEN.append(float(tokens[4]))
        T.append(int(tokens[5]))
    return T, COUNTS, SUM_FREQ, AVG_FREQ, AVG_LEN

def do_scatter_plot(X1, Y1, X2, Y2, X3, Y3, label1, label2, label3, xlabel, ylabel, title, filename):
    plt.figure()
    plt.scatter(X1, Y1, marker='x', label = label1, color = 'r')
    plt.scatter(X2, Y2, marker='o', facecolors='none', edgecolors='g', label = label2)
    plt.scatter(X3, Y3, marker='s', facecolors='none', edgecolors='b', label = label3)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title(title)

    print("Saving to", filename)
    plt.savefig(filename)

def plot_runs(shortest_infile, rarest_infile, dataset_tame, outfile_prefix): #superalphabet_2_infile,
    shortest_t, shortest_counts, shortest_sumfreq, shortest_avgfreq , shortest_avglen = parse(shortest_infile)
    rarest_t, rarest_counts, rarest_sumfreq, rarest_avgfreq, rarest_avglen = parse(rarest_infile)
    #sa2_t, sa2_counts, sa2_sumfreq, sa2_avgfreq = parse(superalphabet_2_infile)

    do_scatter_plot(shortest_t, shortest_counts, rarest_t, rarest_counts,
                    "Shortest", "rarest", "SBWT columns", "# finimizers",
                    "# finimizers ({})".format(dataset_tame),
                    "{}_fmin_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(shortest_t, shortest_sumfreq, rarest_t, rarest_sumfreq,
                    "Shortest", "Rarest", "SA-2", "SBWT columns", "Memory (GB)",
                    "Memory ({})".format(dataset_tame),
                    "{}_sumfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(shortest_t,shortest_avgfreq, rarest_counts, rarest_avgfreq,
                    "Shortest", "Rarest", "k", "Time (s)",
                    "Time ({})".format(dataset_tame),
                    "{}_avgfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(shortest_t,shortest_avglen, rarest_counts, rarest_avglen,
                    "Shortest", "Rarest", "k", "Time (s)",
                    "Time ({})".format(dataset_tame),
                    "{}_avglen_by_t.pdf".format(outfile_prefix))


def read_avgfreqmer_curve(filename):
    K, SIZE = [], []
    for line in open(filename):
        tokens = line.split(",")
        k, size = int(tokens[-2]), int(tokens[-1])
        K.append(k)
        SIZE.append(size)
    return K, SIZE

#def make_avgfreqmer_plot(human_file, coli_file, metagenome_file, outfile):
def make_avgfreqmer_plot( coli_file, outfile):
    #human_K, human_SIZE = read_avgfreqmer_curve(human_file)
    coli_K, coli_SIZE = read_avgfreqmer_curve(coli_file)
    #metagenome_K, metagenome_SIZE = read_avgfreqmer_curve(metagenome_file)

    plt.figure()

    #plt.plot(human_K, human_SIZE, label="Human genome", marker="x")
    plt.plot(coli_K, coli_SIZE, label="E. coli genomes", marker="x")
    #plt.plot(metagenome_K, metagenome_SIZE, label="Metagenome reads", marker="x")

    plt.xlabel("k")
    plt.ylabel("n")

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title("Number of sets in SBWT")

    print("Saving to", outfile)
    plt.savefig(outfile)


if not os.path.exists("plots"):
    os.makedirs("plots")
#make_avgfreqmer_plot("data_for_plots/shortest_coli.csv", "plots/kmers.pdf")

#make_avgfreqmer_plot("data_for_plots/shortest_human.csv", "data_for_plots/shortest_coli.csv", "data_for_plots/shortest_metagenome.csv", "plots/kmers.pdf")
#plot_runs("data_for_plots/shortest_human.csv", "data_for_plots/rarest_human.csv", "data_for_plots/superalphabet-2_human.csv", "Human genome", "plots/human")
plot_runs("data_for_plots/shortest_coli.csv", "data_for_plots/rarest_coli.csv", "E. coli genomes", "plots/coli")
#plot_runs("data_for_plots/shortest_metagenome.csv", "data_for_plots/rarest_metagenome.csv", "data_for_plots/superalphabet-2_metagenome.csv", "Metagenome reads", "plots/metagenome")