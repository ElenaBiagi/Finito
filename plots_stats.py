import sys
import matplotlib.pyplot as plt
import os

def parse(filename):
    T = []
    F_QUERIES = []
    T_QUERIES = []
    TIME = []
    SPACE = []
    SPACE_K = []
    for line in open(filename):
        tokens = line.split(',')
        T.append(int(tokens[0]))
        F_QUERIES.append(int(tokens[1]))
        T_QUERIES.append(int(tokens[2]))
        TIME.append(float(tokens[3]))
        SPACE.append(int(tokens[4]))
        SPACE_K.append(float(tokens[5]))
        
    return T, F_QUERIES, T_QUERIES, TIME, SPACE, SPACE_K
    
def do_scatter_plot(X1, Y1, X2, Y2, X3, Y3, X4, Y4, label1, label2,label3, label4, xlabel, ylabel, title, filename):
    plt.figure()
    plt.scatter(X1, Y1, marker='x', label=label1, color='r')
    plt.scatter(X2, Y2, marker='o', facecolors='none', edgecolors='g', label=label2)
    plt.scatter(X3, Y3, marker='s', label=label3, color='m')  
    plt.scatter(X4, Y4, marker='v', facecolors='none', edgecolors='b', label=label2)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # Start axes from the origin
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)

    plt.legend()
    plt.title(title)

    print("Saving to", filename)
    plt.savefig(filename)
    plt.close("all")


def plot_runs_allk(k1_infile, k2_infile, k3_infile, k4_infile, unitigs, dataset_name, outfile_prefix): 
    k1_t, k1_fqueries, k1_tqueries, k1_time , k1_space, k1_spacek = parse(k1_infile)
    k2_t, k2_fqueries, k2_tqueries, k2_time, k2_space, k2_spacek = parse(k2_infile)
    k3_t, k3_fqueries, k3_tqueries, k3_time , k3_space, k3_spacek = parse(k3_infile)
    k4_t, k4_fqueries, k4_tqueries, k4_time, k4_space, k4_spacek = parse(k4_infile)

    """
    do_scatter_plot(k1_t, k1_fqueries, k2_t, k2_fqueries,k3_t, k3_fqueries, k4_t, k4_fqueries,
                    "k=21", "k=31", "k=63", "k=127", "t", "# finimizers",
                    unitigs+"{}".format(dataset_name),
                    "{}_fmin_by_t.pdf".format(outfile_prefix))

    do_scatter_plot(k1_t, k1_tqueries, k2_t, k2_tqueries,k3_t, k3_tqueries, k4_t, k4_tqueries,
                    "k=21", "k=31", "k=63", "k=127", "t", "Sum of finimizers frequencies",
                    unitigs+"{}".format(dataset_name),
                    "{}_tqueries_by_t.pdf".format(outfile_prefix))
    """
    do_scatter_plot(k1_time, k1_space, k2_time, k2_space, k3_time,k3_space, k4_time, k4_space,
                   "k=21", "k=31", "k=63", "k=127",   "Time(microsec/query)","Space in bytes",
                    unitigs+"{}".format(dataset_name),
                    "{}_time_by_space.pdf".format(outfile_prefix))
    
    do_scatter_plot(k1_time, k1_spacek,  k2_time, k2_spacek, k3_time, k3_spacek, k4_time, k4_spacek,
                   "k=21", "k=31", "k=63", "k=127",  "Time(microsec/query)", "Space(bits/kmer)", 
                    unitigs+"{}".format(dataset_name),
                    "{}_time_by_spacek.pdf".format(outfile_prefix))


if not os.path.exists("test_plots"):
    os.makedirs("test_plots")

# Flipped Unitigs
plot_runs_allk("data_for_plots/search/f_coli3682-unitigs-k21.sbwtstats.txt", "data_for_plots/search/f_coli3682-unitigs-k31.sbwtstats.txt",  "data_for_plots/search/f_coli3682-unitigs-k63.sbwtstats.txt", "data_for_plots/search/f_coli3682-unitigs-k127.sbwtstats.txt", "Flipped Unitigs", "E. coli", "test_plots/search/coli_funitigs")
plot_runs_allk("data_for_plots/search/f_ERR5035349-unitigs-k21.sbwtstats.txt", "data_for_plots/search/f_ERR5035349-unitigs-k31.sbwtstats.txt","data_for_plots/search/f_ERR5035349-unitigs-k63.sbwtstats.txt", "data_for_plots/search/f_ERR5035349-unitigs-k127.sbwtstats.txt", "Flipped Unitigs", "Metagenome", "test_plots/search/metagenome_funitigs")
plot_runs_allk("data_for_plots/search/f_SRR25689478-unitigs-k21.sbwtstats.txt", "data_for_plots/search/f_SRR25689478-unitigs-k31.sbwtstats.txt","data_for_plots/search/f_SRR25689478-unitigs-k63.sbwtstats.txt", "data_for_plots/search/f_SRR25689478-unitigs-k127.sbwtstats.txt", "Flipped Unitigs", "E. coli", "test_plots/search/nanopore_funitigs")

# Flipped Eulertigs
plot_runs_allk("data_for_plots/search/f_coli3682-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/f_coli3682-eulertigs-k31.sbwtstats.txt",  "data_for_plots/search/f_coli3682-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/f_coli3682-eulertigs-k127.sbwtstats.txt", "Flipped Eulertigs", "E. coli", "test_plots/search/coli_feulertigs")
plot_runs_allk("data_for_plots/search/f_ERR5035349-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/f_ERR5035349-eulertigs-k31.sbwtstats.txt","data_for_plots/search/f_ERR5035349-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/f_ERR5035349-eulertigs-k127.sbwtstats.txt", "Flipped Eulertigs", "Metagenome", "test_plots/search/metagenome_feulertigs")
# -all (127) plot_runs_allk("data_for_plots/search/f_SRR25689478-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/f_SRR25689478-eulertigs-k31.sbwtstats.txt","data_for_plots/search/f_SRR25689478-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/f_SRR25689478-eulertigs-k127.sbwtstats.txt", "Flipped Eulertigs", "E. coli", "test_plots/search/nanopore_feulertigs")

# Unitigs
plot_runs_allk("data_for_plots/search/coli3682-unitigs-k21.sbwtstats.txt", "data_for_plots/search/coli3682-unitigs-k31.sbwtstats.txt",  "data_for_plots/search/coli3682-unitigs-k63.sbwtstats.txt", "data_for_plots/search/coli3682-unitigs-k127.sbwtstats.txt", "Unitigs", "E. coli", "test_plots/search/coli_unitigs")
plot_runs_allk("data_for_plots/search/ERR5035349-unitigs-k21.sbwtstats.txt", "data_for_plots/search/ERR5035349-unitigs-k31.sbwtstats.txt","data_for_plots/search/ERR5035349-unitigs-k63.sbwtstats.txt", "data_for_plots/search/ERR5035349-unitigs-k127.sbwtstats.txt", "Unitigs", "Metagenome", "test_plots/search/metagenome_unitigs")
# - all (127)q plot_runs_allk("data_for_plots/search/SRR25689478-unitigs-k21.sbwtstats.txt", "data_for_plots/search/SRR25689478-unitigs-k31.sbwtstats.txt","data_for_plots/search/SRR25689478-unitigs-k63.sbwtstats.txt", "data_for_plots/search/SRR25689478-unitigs-k127.sbwtstats.txt", "Unitigs", "E. coli", "test_plots/search/nanopore_unitigs")

# Eulertigs
# -21 31 63 127 plot_runs_allk("data_for_plots/search/coli3682-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/coli3682-eulertigs-k31.sbwtstats.txt",  "data_for_plots/search/coli3682-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/coli3682-eulertigs-k127.sbwtstats.txt", "Eulertigs", "E. coli", "test_plots/search/coli_eulertigs")
# -21 31 63 127 plot_runs_allk("data_for_plots/search/ERR5035349-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/ERR5035349-eulertigs-k31.sbwtstats.txt","data_for_plots/search/ERR5035349-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/ERR5035349-eulertigs-k127.sbwtstats.txt", "Eulertigs", "Metagenome", "test_plots/search/metagenome_eulertigs")
# -all (127) plot_runs_allk("data_for_plots/search/SRR25689478-eulertigs-k21.sbwtstats.txt", "data_for_plots/search/SRR25689478-eulertigs-k31.sbwtstats.txt","data_for_plots/search/SRR25689478-eulertigs-k63.sbwtstats.txt", "data_for_plots/search/SRR25689478-eulertigs-k127.sbwtstats.txt", "Eulertigs", "E. coli", "test_plots/search/nanopore_eulertigs")

