"""
Module for chromosome comparison. 
This module is used for chromosome to chromosome comparison based on k-mer
counting.  The unique kmers are used to align chromosomes.  The kmer position
table and kmer alignment map are generated for further data analysis. 
"""

import os
import numpy as np 
import re 
from fasta import Fasta
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import gzip
import thread
import time


DEFAULT_KMER_SIZE = 20

def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print "%s: %s" %(time_str, mess)

def _test_dir(dirn=""):
    """Test wether a directory exists.
    Args:
        dirn: directory to be tested. 
    Returns:
        True if existes, and vice versa
    Raises:
        Error message in the shell.
    """
    if os.path.isdir(dirn):
        return True

    _message("%s doesn't exists" %dirn)
    return False

def update_kmer_position(kmer_table, seq="", len_kmer=DEFAULT_KMER_SIZE, size=1, index=0):
    """Obtain all the kmers from a sequence, and update the kmer table 
    Every kmer is added to a kmer_table {"kmer": list of [[list of positions]]}
    Args:
        seq: sequence to get all kmers
        len_kmer: length of the kmer to be checked
        size: total number of entries
        index: index of current sequence entry
    Returns:
        None, the dictionary of the kmer_table is updated
    Raises:
        None. 
    """
    seq = seq.upper()
    _message("kmer counting for the %d-th seq" %index)
    for i in xrange(len(seq) - (len_kmer - 1)):
        kmer = seq[i:i + len_kmer]
        kmer_table.setdefault(kmer, [[] for j in range(size)])
        kmer_table[kmer][index].append(i+1) # Position of kmer in the entry seq

def generate_kmer_table(seq_list, len_kmer=DEFAULT_KMER_SIZE, op_filename=""):
    """Generate the kmer table from the sequence list."""
    
    if not seq_list:
        return None

    seq_size = len(seq_list)
    kmer_table = {}
    for seq_id in range(seq_size):
        update_kmer_position(
            kmer_table, seq_list[seq_id], len_kmer,
            seq_size, seq_id)
        _message("%d kmers are counted in total" % len(kmer_table))

    if op_filename:
        with gzip.open(op_filename, "wb") as op_table:
            _message("Writing kmer-table file: op_filename")
            for kmer, value in kmer_table.iteritems():
                op_table.write("%s; " %kmer)
                for pos_list in value:
                    if not pos_list:
                        op_table.write("; ")
                    else:
                        for pos in pos_list:
                            op_table.write("%d, " % pos)
                        op_table.write("; ") # Always has an empty position entry
                op_table.write("\n")
            _message("Writing finished.")
    return kmer_table

def compare_kmer_table(table_r, table_q, op_filename=""):
    """Compare two kmer table. 

    Every kmer that appears in both table and put in a new table. 
    kmer: [group, [list of position in ref], [list of position is query]]
    Group:
        0: Unique in both groups
        1: Unique in reference, multi in query
        2: Unique in query, multi in reference
        3: Multi in both
    """

    comp_table = {}
    unalign_r = {}
    unalign_q = {}
    _message("Start comparing kmer tables")
    for kmer, poslist_r in table_r.iteritems():
        if kmer in table_q:
            poslist_q = table_q[kmer]
            is_unique_r = (len(poslist_r) == poslist_r.count([]) + 1)
            is_unique_q = (len(poslist_q) == poslist_q.count([]) + 1)
            group = -1
            if is_unique_r and is_unique_q:
                group = 0
            elif is_unique_r:
                group = 1
            elif is_unique_q:
                group = 2
            else:
                group = 3
            comp_table[kmer] = [group, poslist_r, poslist_q]
            del table_q[kmer]
        else:
            unalign_r[kmer] = poslist_r

    for kmer, poslist_q in table_q.iteritems():
        unalign_q[kmer] = poslist_q
    _message("kmer table comparison finished!")
    _message("Total number of shared kmers %d" %(len(comp_table)))
    _message("Reference only kmers: %d" %len(unalign_r))
    _message("Query only kmers: %d" %len(unalign_q))
    if op_filename:
        with gzip.open(op_filename, "wb") as out_p:
            _message("Start writing comparison table at %s" %op_filename)
            for kmer, value in comp_table.iteritems():
                out_p.write("%s; %d; " % (kmer, value[0]))
                for pos_list in value[1]:
                    if not pos_list:
                        out_p.write("; ")
                    else:
                        for pos in pos_list:
                            out_p.write("%d, " % pos)
                        out_p.write("; ") # Always has an empty position entry
                for pos_list in value[2]:
                    if not pos_list:
                        out_p.write("; ")
                    else:
                        for pos in pos_list:
                            out_p.write("%d, " % pos)
                        out_p.write("; ") # Always has an empty position entry
                out_p.write("\n")
            _message("Writing finished.")    
    return comp_table, unalign_r, unalign_q

def transform_kmer_comp_table(comp_table, unalign_r, unalign_q, op_filename=""):
    """Trasform a kmer comparison table to position tuple list.
        Position table:
        (ref_id, query_id, group, pos_ref, pos_query)
        Group Code:
        0: Unique in both groups
        1: Unique in reference, multi in query
        2: Unique in query, multi in reference
        3: Multi in both
        4: kmer only in one strain
    """
    position_list = []
    _message("Generating Position list...")
    for kmer, value in comp_table.iteritems():
        group = value[0]
        pos_ref = value[1]
        pos_query = value[2]
        for ref_id, pos_ls_ref in enumerate(pos_ref):
            if pos_ls_ref:
                for query_id, pos_ls_query in enumerate(pos_query):
                    if pos_ls_query:
                        for p_ref in pos_ls_ref:
                            for p_query in pos_ls_query:
                                position_list.append((
                                    ref_id, query_id, group,
                                    p_ref, p_query))
    for kmer, value in unalign_r.iteritems():
        group = 4
        for ref_id, pos_ls_ref in enumerate(value):
            if pos_ls_ref:
                for p_ref in pos_ls_ref:
                    position_list.append((ref_id, -1, group, p_ref, p_ref))

    for kmer, value in unalign_q.iteritems():
        group = 4
        for query_id, pos_ls_query in enumerate(value):
            if pos_ls_query:
                for p_query in pos_ls_query:
                    position_list.append((
                        -1, query_id, group, p_query, p_query))
    position_list.sort()
    _message("Position list generated")

    if op_filename:
        with gzip.open(op_filename, "wb") as out_p:
            _message("Start writing the position list in %s" %op_filename)
            for position in position_list:
                out_p.write("%d; %d; %d; %d; %d\n" %(
                    position[0], position[1], position[2],
                    position[3], position[4]))
            _message("Finished writing the position list.")
    return position_list

def mask_position_list(pos_list, window=100, op_filename=""):
    """Mask the postition table by a window for drawing

    The postion are rounded by the size of window. 
    Temporary dictionary used to accelerate the program. 
    """
    temp_dict = {}
    mask_poslist = []
    _message("Start masking data at window size = %d" %window)
    for position in pos_list:
        ref_id = position[0]
        query_id = position[1]
        group = position[2]
        p_ref = position[3] / window * window
        p_query = position[4] / window * window 
        mask_pos = (ref_id, query_id, group, p_ref, p_query)
        if mask_pos not in temp_dict:
            temp_dict[mask_pos] = ""
    for pos in temp_dict:
        mask_poslist.append(pos)
    mask_poslist.sort()
    _message("Finished Masking")

    if op_filename:
        with gzip.open(op_filename, "wb") as out_p:
            _message("Writing masked position data in %s" %op_filename)
            for position in mask_poslist:
                out_p.write("%d; %d; %d; %d; %d\n" %(
                    position[0], position[1], position[2],
                    position[3], position[4]))
            _message("Finished writing masked position data")
    return mask_poslist

def draw_chrome_comp_graph(
        ref_name_ls, query_name_ls, ref_len_ls, query_len_ls,
        position_list, dirn="", prefix="comp", ref_genome_name="",
        query_genome_name=""):
    """Draw chromosome to chromosome comparison."""

    if not os.path.isdir(dirn):
        _message("%s is not available." %dirn)
        return
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    if not os.path.isdir(dirn+prefix+"ref/"):
        os.mkdir(dirn+prefix+"ref/")
    if not os.path.isdir(dirn+prefix+"query/"):
        os.mkdir(dirn+prefix+"query/")

    group_color = {
        0: "r",
        1: "b",
        2: "g",
        3: "lightblue",
        4: "lightblue"
    }

    group_label = {
        0: "Unique in both strains",
        1: "Unique in reference strain, multi in query",
        2: "Unique in query strain, multi in reference",
        3: "Multi in both strains",
        4: "Only in one strain"
    }

    group_size = 5
    y_scale = 10
    y_drift = 40
    max_ref_len = max(ref_len_ls)
    max_query_len = max(query_len_ls)
    tick_distance = 1000


    # Draw reference to query graphs
    # Each sequence in query chromosome is aligned to all reference chromosomes

    for query_id in range(len(query_name_ls)):
        q_name = query_name_ls[query_id]
        q_len = query_len_ls[query_id]

        _message("Start drawing the picture for %s" %q_name)

        fig = plt.figure()
        fig.set_size_inches(20,20)
        axes = fig.add_subplot(111, projection='3d')

        for group_inv in range(group_size):
            group = group_size - group_inv -1
            qid_pos = [x for x in position_list if ((x[1] == query_id) and 
                                                    (x[2] == group))]
            if qid_pos:
                rid, qid, group_, p_ref, p_query = zip(*qid_pos)
                p_chrom = [rid_i * group_size*y_scale + group for rid_i in rid]
                axes.scatter(p_query,  p_chrom, p_ref, color=group_color[group],
                             marker=",", s=.3, depthshade=False)
                # each dot, plot a projection line
            _message("Group %d finished" % group)

        y_ticks = np.arange(- group_size * y_scale+y_drift, 
                            len(ref_len_ls) * group_size * y_scale+y_drift,
                            group_size * y_scale)
        axes.set_yticks(y_ticks)

        axes.set_xlabel("Position on %s" %q_name, linespacing=4.)
        axes.set_zlabel("Position on chromosome in %s" % ref_genome_name, 
                        linespacing=4.)
        axes.set_xlim([0, query_len_ls[query_id]])
        axes.set_ylim([-group_size*y_scale+y_drift, 
                       len(ref_len_ls) * group_size*y_scale+y_drift])
        axes.set_zlim([0, max_ref_len])
        axes.grid(True)
        y_tick_label = axes.get_yticks().tolist()
        y_tick_label[0] = "Only in %s" % query_genome_name
        for i in range(len(ref_name_ls)):
            y_tick_label[i+1] = ref_name_ls[i]
        axes.set_yticklabels(y_tick_label, horizontalalignment="left",
                             verticalalignment="baseline")
        plt.title("%s to %s" %(q_name, ref_genome_name))

        fig_name = dirn+prefix+"query/"+q_name+".png"
        plt.savefig(fig_name, dpi=300, format='png')
        plt.close()


        fig = plt.figure()
        fig.set_size_inches(20,20)
        axes = fig.add_subplot(111, projection='3d')

        for group in range(group_size):
            qid_pos = [x for x in position_list if ((x[1] == query_id) and 
                                                    (x[2] == group))]
            if qid_pos:
                rid, qid, group_, p_ref, p_query = zip(*qid_pos)
                p_chrom = [rid_i * group_size *y_scale + group for rid_i in rid]
                for i in range(len(rid)):
                    x = (p_query[i], p_query[i])
                    z = (p_ref[i], 0)
                    y = (p_chrom[i] , p_chrom[i] )
                    if p_ref[i] < 0:
                        print p_ref[i], p_query[i]
                    axes.plot(x, y, z, color=group_color[group], linewidth=.5,
                              linestyle="--")
            _message("Group %d finished" % group)

        y_ticks = np.arange(- group_size * y_scale +y_drift , 
                            len(ref_len_ls) * group_size * y_scale +y_drift,
                            group_size * y_scale)
        axes.set_yticks(y_ticks)
        axes.set_xlabel("Position on %s" %q_name, linespacing=4.)
        axes.set_zlabel("Position on chromosome in %s" % ref_genome_name, 
                        linespacing=4.)
        axes.set_xlim(0, query_len_ls[query_id])
        axes.set_ylim(-group_size*y_scale + y_drift, 
                      len(ref_len_ls) * group_size*y_scale +y_drift)
        axes.set_zlim(0, max_ref_len)
        axes.grid(True)
        y_tick_label = axes.get_yticks().tolist()
        y_tick_label[0] = "Only in %s" % query_genome_name
        for i in range(len(ref_name_ls)):
            y_tick_label[i+1] = ref_name_ls[i]
        axes.set_yticklabels(y_tick_label, horizontalalignment="left",
                             verticalalignment="baseline")
        plt.title("%s to %s" %(q_name, ref_genome_name))

        fig_name = dirn+prefix+"query/"+q_name+"-project.png"
        plt.savefig(fig_name, dpi=300, format='png')
        plt.close()

        _message("%s finished" %q_name)



    for ref_id in range(len(ref_name_ls)):
        r_name = ref_name_ls[ref_id]
        r_len = ref_len_ls[ref_id]

        _message("Start drawing the picture for %s" %r_name)

        fig = plt.figure()
        fig.set_size_inches(20,20)
        axes = fig.add_subplot(111, projection='3d')

        for group_inv in range(group_size):
            group = group_size - group_inv -1
            rid_pos = [x for x in position_list if ((x[0] == ref_id) and 
                                                    (x[2] == group))]
            if rid_pos:
                rid, qid, group_, p_ref, p_query = zip(*rid_pos)
                p_chrom = [qid_i * group_size*y_scale + group for qid_i in qid]
                axes.scatter(p_ref,  p_chrom, p_query, color=group_color[group],
                             marker=",", s=.3, depthshade=False)
                # each dot, plot a projection line
            _message("Group %d finished" % group)

        y_ticks = np.arange(- group_size * y_scale+y_drift, 
                            len(query_len_ls) * group_size * y_scale+y_drift,
                            group_size * y_scale)
        axes.set_yticks(y_ticks)

        axes.set_xlabel("Position on %s" %r_name, linespacing=4.)
        axes.set_zlabel("Position on chromosome in %s" % query_genome_name, 
                        linespacing=4.)
        axes.set_xlim([0, ref_len_ls[ref_id]])
        axes.set_ylim([-group_size*y_scale+y_drift, 
                       len(query_len_ls) * group_size*y_scale+y_drift])
        axes.set_zlim([0, max_query_len])
        axes.grid(True)
        y_tick_label = axes.get_yticks().tolist()
        y_tick_label[0] = "Only in %s" % ref_genome_name
        for i in range(len(query_name_ls)):
            y_tick_label[i+1] = query_name_ls[i]
        axes.set_yticklabels(y_tick_label, horizontalalignment="left",
                             verticalalignment="baseline")
        plt.title("%s to %s" %(r_name, query_genome_name))

        fig_name = dirn+prefix+"ref/"+r_name+".png"
        plt.savefig(fig_name, dpi=300, format='png')
        plt.close()

    for ref_id in range(len(ref_name_ls)):
        r_name = ref_name_ls[ref_id]
        r_len = ref_len_ls[ref_id]

        _message("Start drawing the picture for %s" %r_name)

        fig = plt.figure()
        fig.set_size_inches(20,20)
        axes = fig.add_subplot(111, projection='3d')

        for group_inv in range(group_size):
            group = group_size - group_inv -1
            rid_pos = [x for x in position_list if ((x[0] == ref_id) and 
                                                    (x[2] == group))]
            if rid_pos:
                rid, qid, group_, p_ref, p_query = zip(*rid_pos)
                p_chrom = [qid_i * group_size*y_scale + group for qid_i in qid]
                # axes.scatter(p_ref,  p_chrom, p_query, color=group_color[group],
                #              marker=",", s=.3, depthshade=False)
                # each dot, plot a projection line
                for i in range(len(qid)):
                    x = (p_ref[i], p_ref[i])
                    z = (p_query[i], 0)
                    y = (p_chrom[i] , p_chrom[i] )
                    axes.plot(x, y, z, color=group_color[group], linewidth=.5,
                              linestyle="--")                
            _message("Group %d finished" % group)

        y_ticks = np.arange(- group_size * y_scale+y_drift, 
                            len(query_len_ls) * group_size * y_scale+y_drift,
                            group_size * y_scale)
        axes.set_yticks(y_ticks)

        axes.set_xlabel("Position on %s" %r_name, linespacing=4.)
        axes.set_zlabel("Position on chromosome in %s" % query_genome_name, 
                        linespacing=4.)
        axes.set_xlim([0, ref_len_ls[ref_id]])
        axes.set_ylim([-group_size*y_scale+y_drift, 
                       len(query_len_ls) * group_size*y_scale+y_drift])
        axes.set_zlim([0, max_query_len])
        axes.grid(True)
        y_tick_label = axes.get_yticks().tolist()
        y_tick_label[0] = "Only in %s" % ref_genome_name
        for i in range(len(query_name_ls)):
            y_tick_label[i+1] = query_name_ls[i]
        axes.set_yticklabels(y_tick_label, horizontalalignment="left",
                             verticalalignment="baseline")
        plt.title("%s to %s" %(r_name, query_genome_name))

        fig_name = dirn+prefix+"ref/"+r_name+"-project.png"
        plt.savefig(fig_name, dpi=300, format='png')
        plt.close()
 
        _message("%s finished" %r_name)


    

def load_fasta(fname="", pattern=""):
    """Load Fasta file and return the sequences and names of the fasta file. 
    Args:
    fname:
        file name of the fasta file with the address
    pattern:
        reg pattern to sort file names
    Returns:
        Name_List, Sequence_list
    """

    fasta_ = Fasta(fname)

    if pattern:
        def _get_seq_tag(name):
            _match = re.match(pattern, name)
            return 
        _name_list, _seq_list = 
    return fasta_.Names, fasta_.Seqs


def rc_seq(seq=""):
    """Returns the reverse compliment sequence."""
    rc_nt_ls = []
    rc_dict = {
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C"
    }
    rc_nt_ls = [rc_dict[seq[i]] for i in range(len(seq)-1, -1, -1)]
    rc_seq = "".join(rc_nt_ls)
    return rc_seq



def adjust_direction(fn="", tag="--r", prefix="adj"):
    """Generate Reverse Compliment seq for sequence with the tag."""
    fasta_ = Fasta(fn)
    _extention = os.path.splitext(fn)[1][1:].strip()
    _path = os.path.splitext(fn)[0]
    op_name = "%s-%s.%s" %(_path, prefix, _extention)
    with open(op_name, "w") as out_p: 
        for seq_id, seq_name in enumerate(fasta_.Names):
            if tag in seq_name:
                ext_p = seq_name.find(tag)
                _name = seq_name
                out_p.write(">%s\n")
                _rc_seq = rc_seq(fasta_.Seqs[seq_id])
                out_p.write(_rc_seq)
                out_p.write("\n")
            else:
                out_p.write(">%s\n" %fasta_.Names[seq_id])
                out_p.write(fasta_.Seqs[seq_id])
                out_p.write("\n")


def test():
    """Temporary fucntion to test code."""
    fname_b = "/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa"
    fname_c = "/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa"
    pattern = r"_(\w)_(\w)"
    names_b, seqs_b = load_fasta(fname_b, pattern=pattern)
    names_c, seqs_c = load_fasta(fname_c, pattern=pattern)
    tab_b = generate_kmer_table(seq_list=seqs_b)
    tab_c = generate_kmer_table(seq_list=seqs_c)
    tab_comp, unalign_r, unalign_q = compare_kmer_table(tab_b, tab_c)
    pos_ls = transform_kmer_comp_table(tab_comp, unalign_r, unalign_q)
    mask_pos_ls = mask_position_list(pos_ls, window=200)
    len_b = []
    len_c = []
    for seq in seqs_b:
        len_b.append(len(seq))
    for seq in seqs_c:
        len_c.append(len(seq))

    draw_chrome_comp_graph(
        ref_name_ls=names_b, query_name_ls=names_c,
        ref_len_ls=len_b, query_len_ls=len_c,
        position_list=mask_pos_ls, dirn="/home/zhuwei/cglabarata/comp/kmer/", 
        prefix="bg2-cbs/",
        ref_genome_name="BG2",
        query_genome_name="CBS")

if __name__ == "__main__":
    test()