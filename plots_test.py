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
        T.append(int(tokens[0]))
    return T, COUNTS, SUM_FREQ, AVG_FREQ, AVG_LEN

def do_scatter_plot(X1, Y1, X2, Y2, label1, label2, xlabel, ylabel, title, filename):
    plt.figure()
    plt.scatter(X1, Y1, marker='x', label = label1, color = 'r')
    plt.scatter(X2, Y2, marker='o', facecolors='none', edgecolors='g', label = label2)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title(title)

    print("Saving to", filename)
    plt.savefig(filename)
    
def do_scatter_plot2(X1, Y1, X2, Y2, X3, Y3, X4, Y4, label1, label2,label3, label4, xlabel, ylabel, title, filename):
    plt.figure()
    plt.scatter(X1, Y1, marker='x', label=label1, color='r')
    plt.scatter(X2, Y2, marker='o', facecolors='none', edgecolors='g', label=label2)
    plt.scatter(X3, Y3, marker='s', label=label3, color='m')  
    plt.scatter(X4, Y4, marker='v', facecolors='none', edgecolors='b', label=label2)
    # Lines connecting the points
    plt.plot(X1, Y1, linestyle='-', color='r', linewidth=0.5)
    plt.plot(X2, Y2, linestyle='-', color='g', linewidth=0.5)
    plt.plot(X3, Y3, linestyle='-', color='m', linewidth=0.5)
    plt.plot(X4, Y4, linestyle='-', color='b', linewidth=0.5)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from the origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title(title)

    print("Saving to", filename)
    plt.savefig(filename)
    

def plot_runs(k3_infile, k2_infile, dataset_name, outfile_prefix):
    k3_t, k3_counts, k3_sumfreq, k3_avgfreq , k3_avglen = parse(k3_infile)
    k2_t, k2_counts, k2_sumfreq, k2_avgfreq, k2_avglen = parse(k2_infile)

    do_scatter_plot(k3_t, k3_counts, k2_t, k2_counts,
                    "Shortest", "rarest", "SBWT columns", "# finimizers",
                    "# finimizers ({})".format(dataset_name),
                    "{}_fmin_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(k3_t, k3_sumfreq, k2_t, k2_sumfreq,
                    "Shortest", "Rarest", "SBWT columns", "Memory (GB)",
                    "Memory ({})".format(dataset_name),
                    "{}_sumfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(k3_t,k3_avgfreq, k2_counts, k2_avgfreq,
                    "Shortest", "Rarest", "k", "Time (s)",
                    "Time ({})".format(dataset_name),
                    "{}_avgfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(k3_t,k3_avglen, k2_counts, k2_avglen,
                    "Shortest", "Rarest", "k", "Time (s)",
                    "Time ({})".format(dataset_name),
                    "{}_avglen_by_t.pdf".format(outfile_prefix))

def plot_runs_allk(k1_infile, k2_infile, k3_infile, k4_infile, dataset_name, outfile_prefix): 
    k1_t, k1_counts, k1_sumfreq, k1_avgfreq , k1_avglen = parse(k1_infile)
    k2_t, k2_counts, k2_sumfreq, k2_avgfreq, k2_avglen = parse(k2_infile)
    k3_t, k3_counts, k3_sumfreq, k3_avgfreq , k3_avglen = parse(k3_infile)
    k4_t, k4_counts, k4_sumfreq, k4_avgfreq, k4_avglen = parse(k4_infile)

    do_scatter_plot2(k1_t, k1_counts, k2_t, k2_counts,k3_t, k3_counts, k4_t, k4_counts,
                    "k=21", "k=31", "k=63", "k=127", "t", "# finimizers",
                    "Flipped Unitigs {}".format(dataset_name),
                    "{}_fmin_by_t.pdf".format(outfile_prefix))

    do_scatter_plot2(k1_t, k1_sumfreq, k2_t, k2_sumfreq,k3_t, k3_sumfreq, k4_t, k4_sumfreq,
                    "k=21", "k=31", "k=63", "k=127", "t", "Sum of finimizers frequencies",
                    "Flipped Unitigs {}".format(dataset_name),
                    "{}_sumfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot2(k1_t,k1_avgfreq, k2_t, k2_avgfreq, k3_t,k3_avgfreq, k4_t, k4_avgfreq,
                   "k=21", "k=31", "k=63", "k=127", "t", "Average finimmizer frequency",
                    "Flipped Unitigs {}".format(dataset_name),
                    "{}_avgfreq_by_t.pdf".format(outfile_prefix))

    do_scatter_plot2(k1_t,k1_avglen, k2_t, k2_avglen, k3_t,k3_avglen, k4_t, k4_avglen,
                    "k=21", "k=31", "k=63", "k=127", "t", "Average finimmizer length",
                    "Flipped Unitigs {}".format(dataset_name),
                    "{}_avglen_by_t.pdf".format(outfile_prefix))


if not os.path.exists("test_plots"):
    os.makedirs("test_plots")

#plot_runs("data_for_plots/k3_human.csv", "data_for_plots/k2_human.csv", "data_for_plots/superalphabet-2_human.csv", "Human genome", "plots/human")
#plot_runs("data_for_plots/k3_coli.csv", "data_for_plots/k2_coli.csv", "E. coli genomes", "plots/coli")
#plot_runs("data_for_plots/k3_metagenome.csv", "data_for_plots/k2_metagenome.csv", "data_for_plots/superalphabet-2_metagenome.csv", "Metagenome reads", "plots/metagenome")

plot_runs_allk("data_for_plots/f_coli3682-unitigs-k21.sbwttest.txt", "data_for_plots/f_coli3682-unitigs-k31.sbwttest.txt","data_for_plots/f_coli3682-unitigs-k63.sbwttest.txt", "data_for_plots/f_coli3682-unitigs-k127.sbwttest.txt", "E. coli", "test_plots/coli_funitigs")
plot_runs_allk("data_for_plots/f_ERR5035349-unitigs-k21.sbwttest.txt", "data_for_plots/f_ERR5035349-unitigs-k31.sbwttest.txt","data_for_plots/f_ERR5035349-unitigs-k63.sbwttest.txt", "data_for_plots/f_ERR5035349-unitigs-k127.sbwttest.txt", "Metagenome", "test_plots/metagenome_funitigs")
plot_runs_allk("data_for_plots/f_SRR25689478-unitigs-k21.sbwttest.txt", "data_for_plots/f_SRR25689478-unitigs-k31.sbwttest.txt","data_for_plots/f_SRR25689478-unitigs-k63.sbwttest.txt", "data_for_plots/f_SRR25689478-unitigs-k127.sbwttest.txt", "E. coli", "test_plots/nanopore_funitigs")
