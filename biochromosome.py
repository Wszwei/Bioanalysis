"""
Module for chromosome comparison.
This module is used for chromosome to chromosome comparison based on k-mer
counting.  The unique kmers are used to align chromosomes.  The kmer position
table and kmer alignment map are generated for further data analysis.
"""

import os
import gzip
import thread
import time
import itertools
import re
import numpy as np
from fasta import Fasta
from operator import itemgetter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D


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

class KmerSTreeMRBD(object):
    """Suffix Tree of non-repeat bi-directional kmers.

    The suffix tree for genomic sequence specialized for kmer searching of DNA
    sequence. Each kmer is
    The kmer is removed from the suffix tree if it appears multiple times in
    one sequence. At the end of the sequence, the id of the sequence as well as
    the position of the kmer is stored.

    The purpose of this tree is to compare the chromosome in two strains. After
    transverse of the two trees, a kmer-chromosome existance matrix is
    established.  Each row represents the existance of one kmer along the
    sequence, and each column represents the kmer distribution on each sequence.
    Thus we can compare the genetic difference in two chromosomes.

    Attributes:
        num_sequence: Number of sequence loaded to the suffix tree
        len_kmer: length of kmer to build the tree
        stree: dictionary of suffix tree head node
        num_kmer: Number of kmers stored in the stree
    """

    def __init__(self, seq_list, len_kmer):
        """Inits KmerSTreeMR class.

        Args:
        seq_list:
            List of sequence to build the suffix tree.
        """
        self.num_sequence = len(seq_list)
        self.seq_list = seq_list
        self.seq_num_list = self.seq2num()
        self.len_kmer = len_kmer
        self.stree = self.build_stree()

    def build_tree(self):
        """Build the suffix tree"""
        stree_ = {}
        for seq_id, num_seq in enumerate(self.seq_num_list):
            _cur_kmer = num_seq[:self.len_kmer]
            self.add_kmer(_cur_kmer)
            for _id in xrange(len(num_seq) - self.len_kmer):
                _new_nt = num_seq[_id + self.len_kmer]
                _cur_kmer.pop(0)
                _cur_kmer.append(num_seq[_id])
                self.add_kmer(_cur_kmer, seq_id)

    def add_kmer(self, kmer, seq_id):
        """Add a specific kmer."""
        if len(kmer) < self.len_kmer:
            return
        # Get the smaller one
        _kmer = self.small_numseq(kmer)

        # if the first nt in the stree
        if kmer[0] in self.stree:
            pass

    def remove_kmer(self, kmer):
        """Remove kmer from the stree."""
        pass

    def find_kmer(self, kmer):
        """Find the seq_id of the kmer on the stree."""
        pass

    def trim_stree(self):
        """Trim the stree."""
        pass

    def small_numseq(self, child_list):
        """Return the smaller numbered seq in seq and its rcseq. """
        if not child_list:
            return []
        rc_child_list = [5-child_list[i] for i in range(len_kmer-1, -1, -1)]
        if child_list < rc_child_list:
            return child_list
        else:
            return rc_child_list

    def seq2num(self):
        """Turn seq string to list of numbers.

        Dictionary:
            A:1, C:2, G:3, T:4
        """
        nt2num = {
            "a": 1,
            "c": 2,
            "g": 3,
            "t": 4
        }
        _seq_num_list = []
        for seq in self.seq_list:
            _cur_num_list = []
            seq = seq.lower()
            _cur_num_list = [nt2num[_nt] for _nt in seq]
            _seq_num_list.append(_cur_num_list)

        return _seq_num_list

class KmerSTreeNode(object):
    """Class of the Node of a kmer STree.

    Attributes:
        node:
            A dict points to child node \d: {node}
            If the node is a leaf, then it only has one child: -1: [seq id]
    """

    def __init__(self, node):
        """Inits a STree Node.

        Args:
        node:
            a reference to the dict of current node.
        """
        self.node = {}

    def add_child(child_list, seq_id):
        """Add a child to the node.

        Args:
        child_list:
            list of integers for the suffix tree
        seq_id:
            The id of the sequence of the kmer
        """
        pass

    def remove_child(child_list):
        """Remove a child from the node."""
        pass



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
        # Mask rcSeq
        kmer_rc = rc_seq(kmer)
        if kmer > kmer_rc:
            kmer = kmer_rc
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
    kmer: [group, [list of position in ref], [list of position in query]]
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

def generate_kmer_table_v2(seq_list, len_kmer=DEFAULT_KMER_SIZE, op_filename="", mask_rep=True, mask_rc=True):
    """Generate the kmer table from the sequence list.
    The repetitive kmers in one sequence can be removed by mask_rep flag
    """

    if not seq_list:
        return None
    seq_size = len(seq_list)
    kmer_table = {}
    rc_dict = {
        'a': 't',
        'c': 'g',
        't': 'a',
        'g': 'c'
    }

    def add_kmer(kmer, rc_kmer, seq_id, pos, op_filename=""):
        """Add kmer to the kmer table. """
        if mask_rc:
            if kmer < rc_kmer:
                small_kmer = kmer
            else:
                small_kmer = rc_kmer
            if mask_rep:
                if small_kmer in kmer_table:
                    # Mask repeat
                    if seq_id in kmer_table[small_kmer]:
                        del kmer_table[small_kmer]
                    else:
                        kmer_table[small_kmer].append((seq_id, pos))
                else:
                    kmer_table[small_kmer] = [(seq_id, pos)]
            else:
                kmer_table.setdefault(small_kmer, [])
                kmer_table[small_kmer].append((seq_id, pos))
        else:
            if mask_rep:
                if kmer in kmer_table:
                    if seq_id in kmer_table[kmer]:
                        del kmer_table[kmer]
                    else:
                        kmer_table[kmer].append((seq_id, pos))
                else:
                    kmer_table[kmer] = [(seq_id, pos)]
            else:
                kmer_table.setdefault(kmer, [])
                kmer_table[kmer].append((seq_id, pos))

    for seq_id, seq_ in enumerate(seq_list):
        _message("Start kmer counting in seq_id: %d" %seq_id)
        _cur_kmer = seq_[:len_kmer]
        if len(_cur_kmer) < len_kmer:
            _message("seq_id: %d :Too short!" %seq_id)
            return None
        _cur_rc_kmer = rc_seq(_cur_kmer)
        add_kmer(_cur_kmer, _cur_rc_kmer, seq_id, 1)
        for _id in xrange(len(seq_) - len_kmer):
            _cur_nt = seq_[_id + len_kmer]
            _cur_kmer = "".join([_cur_kmer[1:], _cur_nt])
            _cur_rc_kmer = "".join([rc_dict[_cur_nt], _cur_rc_kmer[:-1]])
            add_kmer(_cur_kmer, _cur_rc_kmer, seq_id, _id + 2)

        _message("%d kmers are counted in total" % len(kmer_table))

    if op_filename:
        with gzip.open(op_filename, "wb") as op_table:
            _message("Writing kmer-table file: %s" % op_filename)
            for kmer, value in kmer_table.iteritems():
                op_table.write("%s; " %kmer)
                value.sort()
                for seq_id in value:
                    op_table.write("%d; " % seq_id)
                op_table.write("\n")
            _message("Writing finished.")

    return kmer_table

def generate_kmer_array(ref_seq_ls, query_seq_ls, ref_len_ls, query_len_ls, len_kmer=DEFAULT_KMER_SIZE, op_filename="", mask_rep=True, mask_rc=True):
    """Compare two kmer table, and get the kmer array list.

    Returns
    array_ref/query: [distribution of kmer on sequences]
    array_corr: array_ref * (array_query)t = ref X query correlation map
    ref_first_list: [location of the first kmer in the ref seq]
    query_first_list: [lcoation of the first kmer in the query seq]

    """
    from math import sqrt
    ref_size = len(ref_len_ls)
    query_size = len(query_len_ls)
    _ref_total_len = sum(ref_len_ls)
    _query_total_len = sum(query_len_ls)
    _message("Inits kmer arrays")
    # Caluculate the start postion of ref and query
    ref_first_list = [0]
    _pos = 0
    for ref_len in ref_len_ls:
        _pos += ref_len
        ref_first_list.append(_pos)
    query_first_list = [0]
    _pos = 0
    for query_len in query_len_ls:
        _pos += query_len
        query_first_list.append(_pos)
    # Inits the reference numpy array, default = 1j/sqrt(ref_size)
    _default_ref = 1j/sqrt(ref_size)
    _default_query = 1j/sqrt(query_size)
    array_ref = np.empty([_ref_total_len, ref_size])
    array_ref.fill(_default_ref)
    array_query = np.empty([_query_total_len, query_size])
    array_query.fill(_default_query)

    _message("Start making reference kmer table")
    ref_table = generate_kmer_table_v2(ref_seq_ls, len_kmer, mask_rep, mask_rc)
    _message("Finish making reference kmer table")
    _message("Start making query kmer table")
    query_table = generate_kmer_table_v2(query_seq_ls, len_kmer, mask_rep, mask_rc)
    _message("Finish making reference kmer table")
    _message("Start making the array list")

    for kmer, ref_seq_pos in  ref_table.iteritems():
        # Shared kmer
        # ONLY works for repeat-masked kmers
        # TODO: update to universal version
        if kmer in query_table:
            curr_dist_ref = [0.] * ref_size
            curr_dist_query = [0.] * query_size
            query_seq_pos = query_table[kmer]
            _sum_ref = len(ref_seq_pos) * 1.
            _sum_query = len(query_seq_pos) * 1.
            for seq_pos in ref_seq_pos:
                _id = seq_pos[0]
                curr_dist_ref[_id] = 1./sqrt(_sum_ref)
            for seq_pos in query_seq_pos:
                _id = seq_pos[0]
                curr_dist_query[_id] = 1./sqrt(_sum_query)

            # Update reference and query array
            for seq_pos in ref_seq_pos:
                _id = seq_pos[0]
                _pos = seq_pos[1]
                _abs_pos = ref_first_list[_id] + _pos -1
                array_ref[_abs_pos] = np.array(curr_dist_ref)

            for seq_pos in query_seq_pos:
                _id = seq_pos[0]
                _pos = seq_pos[1]
                _abs_pos = query_first_list[_id] + _pos -1
                array_query[_abs_pos] = np.array(curr_dist_query)
    ref_table.clear()
    query_table.clear()
    _message("kmer array is generated")

    # array_corr = np.dot(array_ref, array_query.transpose())

    # return array_ref, array_query, array_corr, ref_first_list, query_first_list
    return array_ref, array_query, ref_first_list, query_first_list

def compare_seq(ref_seq_ls, query_seq_ls, len_kmer=DEFAULT_KMER_SIZE, prefix="", dirn="", mask_rep=True, mask_rc=True):
    """Compare sequence in reference and query by kmer. """
    if not os.path.isdir(dirn):
        _message("%s is not available!" %dirn)
    fn_ref_array = dirn + prefix + "_ref_array.txt.gz"
    fn_query_array = dirn + prefix + "_query_array.txt.gz"
    fn_corr_array = dirn + prefix + "_corr_array.txt.gz"

    ref_len_ls = [len(seq) for seq in ref_seq_ls]
    query_len_ls = [len(seq) for seq in query_seq_ls]
    # Compare the two sequence by kmer, and generate the distribution
    # of the non-repetitive shared kmers
    array_ref, array_query, ref_first_list, query_first_list = generate_kmer_array(
        ref_seq_ls=ref_seq_ls,
        query_seq_ls=query_seq_ls,
        ref_len_ls=ref_len_ls,
        query_len_ls=query_len_ls)
    np.savetxt(fn_ref_array, array_ref, delimiter=',')
    np.savetxt(fn_query_array, array_query, delimiter=',')
    np.savetxt(fn_corr_array, array_query, delimiter=',')

def kmer_homology_group(ref_seq_ls, query_seq_ls, len_kmer=DEFAULT_KMER_SIZE, op_filename="", mask_rep=True, mask_rc=True):
    """Generate a homology map by kmer comparison. """

    _message("Start making reference kmer table")
    ref_table = generate_kmer_table_v2(ref_seq_ls, len_kmer, mask_rep, mask_rc)
    _message("Finish making reference kmer table")
    _message("Start making query kmer table")
    query_table = generate_kmer_table_v2(query_seq_ls, len_kmer, mask_rep, mask_rc)
    _message("Finish making reference kmer table")
    # TODO: generate statistics of homology groups

    homo_goup_ls = []
    # entry in the homo_group_ls:
    # ([seq_id in the group], [legnth of homologous region])

def compare_kmer_table_v2():
    pass

def homology_group(ref_seq_ls, query_seq_ls):
    pass

def draw_chrom_homo_matrix(
        ref_name_ls, query_name_ls, ref_len_ls, query_len_ls,
        position_list, dirn="", prefix="comp", ref_genome_name="",
        query_genome_name="", window=100):
    """Draw the homology, comparison map."""
    if not os.path.isdir(dirn):
        _message("%s is not available." %dirn)
        return
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    if not os.path.isdir(dirn+prefix+"ref/"):
        os.mkdir(dirn+prefix+"ref/")
    if not os.path.isdir(dirn+prefix+"query/"):
        os.mkdir(dirn+prefix+"query/")

    group_color = ['r', 'b', 'g', 'lightblue', 'grey']


    group_label = [
        "Unique in both strains",
        "Unique in %s, multi in %s" %(ref_genome_name, query_genome_name),
        "Unique in %s, multi in %s" %(query_genome_name, ref_genome_name),
        "Multi in both strains",
        "Only in one strain"]
    group_size = 5
    dpi = 300


    ref_start_pos = [0]
    query_start_pos = [0]
    # reference (x) axis tick position
    ref_label_pos = []
    # query (y) axis tick position
    query_label_pos = []
    # List of the start positions of reference and query
    for _len in ref_len_ls:
        ref_start_pos.append(_len + ref_start_pos[-1])
        ref_label_pos.append((ref_start_pos[-1] + ref_start_pos[-2])/2)
    for _len in query_len_ls:
        query_start_pos.append(_len + query_start_pos[-1])
        query_label_pos.append((query_start_pos[-1] + query_start_pos[-2])/2)
    total_ref = ref_start_pos[-1]
    total_query = query_start_pos[-1]

    # # Separate figure for each group
    # _message("Start drawing kmer alignment figure")
    # height = query_label_pos[-1] *(group_size-1)/window/dpi + 20
    # widith = ref_label_pos[-1]/window/dpi + 10
    # fig, axes = plt.subplots(group_size-1, 1,
    #                          figsize=(widith, height), dpi=dpi)
    # for group in range(group_size-1):
    #     pos_list = [x for x in position_list if x[2] == group]
    #     corr_list = [(ref_start_pos[x[0]] + x[3],
    #                   query_start_pos[x[1]] + x[4])
    #                  for x in pos_list]
    #     p_ref, p_query = zip(*corr_list)
    #     axes[group].scatter(p_ref, p_query, color=group_color[group],
    #                         marker=",", s=.2)
    # # Draw lines to separate different seqs

    # for _id in range(group_size-1):
    #     for _pos in ref_start_pos:
    #         axes[_id].plot((_pos, _pos), (0, total_query), color='black', linewidth=.1)
    #     for _pos in query_start_pos:
    #         axes[_id].plot((0, total_ref), (_pos, _pos), color='black', linewidth=.1)

    #     axes[_id].set_title(group_label[_id])
    #     axes[_id].set_xticks(ref_label_pos)
    #     axes[_id].set_xticklabels(ref_name_ls, rotation='vertical', ha='right')
    #     axes[_id].set_yticks(query_label_pos)
    #     axes[_id].set_yticklabels(query_name_ls, rotation='horizontal', va='top')

    # fig_name = "%s%s%s-%s.png" %(dirn, prefix, ref_genome_name, query_genome_name)
    # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    # plt.savefig(fig_name, bbox_inches='tight', format='png')
    # plt.close()
    # _message("Separare figure for kmer alignment is generated!")

    # Combined figure for all the groups
    # Origincal Group Color
    # group_color = {
    #     0: "r",
    #     1: "b",
    #     2: "g",
    #     3: "lightblue",
    #     4: "grey" -> alighment only in one strain, ommit
    # }
    # Combined pattern:
    # 0: 'r' == has group 1, without group 1 or 2: Paired Unique region
    # 1: 'b' == has group 1, without group 2: LOH in reference
    # 2: 'g' == has group 2, without group 1: LOH in query
    # 3: 'cyan' == has group 1 and 2: Gene Convertion
    # 4: 'lightblue' == only has 4: homologous region

    com_group_color = ['r', 'b', 'g', 'cyan', 'lightblue']

    com_group_label = [
        "Paired Unique region",
        "Loss of Homology in %s" % ref_genome_name,
        "Loss of Homology in %s" % query_genome_name,
        "LOH in both strains(Gene Conversion)",
        "Homologous region"]
    com_group_size = len(com_group_color)
    # sort by ref_id > query_id > ref_pos > query_pos
    _message("Generating Annotated bin list")
    position_list.sort(key=lambda x: (x[0], x[1], x[3], x[4]))
    com_pos_list = []
    def add_bin2_pos_list(pos_list_, bin_):
        """Add the bin to the combine position list."""
        # All the groups in the current bin
        _groups = bin_[1]
        # No LOH in both strains
        if (1 not in _groups) and (2 not in _groups):
            # Group 0: Unique paired
            if 0 in _groups:
                pos_list_.append(bin_[0] + (0,))
            # Group 4: Homologous region
            else:
                pos_list_.append(bin_[0] + (4,))
        # LOH in reference
        elif (1 in _groups) and (2 not in _groups):
            pos_list_.append(bin_[0] + (1,))
        # LOH in query
        elif (2 in _groups) and (1 not in _groups):
            pos_list_.append(bin_[0] + (2,))
        # LOH in both strains, gene convertion
        else:
            pos_list_.append(bin_[0] + (3,))

    _curr_bin = ()
    _last_pos = position_list[-1]
    for _pos in position_list:
        _curr_pos = (_pos[0], _pos[1], _pos[3], _pos[4])
        _curr_group = _pos[2]
        # Mask kmer only in one strain
        if _curr_group == 4:
            continue
        if not _curr_bin:
            _curr_bin = (_curr_pos, [_curr_group])
        else:
            if _curr_bin[0] == _curr_pos:
                _curr_bin[1].append(_curr_group)
                if  _pos == _last_pos:
                    add_bin2_pos_list(com_pos_list, _curr_bin)
            else:
                add_bin2_pos_list(com_pos_list, _curr_bin)
                _curr_bin = (_curr_pos, [_curr_group])

    _message("Finished generation of bin list")

    _message("Start drawing figure for the annotated bins")

    height = query_label_pos[-1]/window/dpi + 20
    widith = ref_label_pos[-1]/window/dpi + 10
    marker_size = (1./dpi)**2
    fig, axes = plt.subplots(figsize=(widith, height), dpi=dpi)

    # Draw the main figure
    for group in range(com_group_size):
        pos_list = [x for x in com_pos_list if x[4] == group]
        if not pos_list:
            _message("No bin is %s" %com_group_label[group])
            continue
        corr_list = [(ref_start_pos[x[0]] + x[2],
                      query_start_pos[x[1]] + x[3])
                     for x in pos_list]
        p_ref, p_query = zip(*corr_list)
        axes.plot(p_ref, p_query, color=com_group_color[group],
                  marker="s", fillstyle='full',
                  markersize=marker_size, lw=0, linestyle="")
    # Draw the lines for the end of sequences
    for _pos in ref_start_pos:
        axes.plot((_pos, _pos), (0, total_query), color='black', linewidth=.1)
    for _pos in query_start_pos:
        axes.plot((0, total_ref), (_pos, _pos), color='black', linewidth=.1)

    # Setup x, y labels
    axes.set_xticks(ref_label_pos)
    axes.set_yticks(query_label_pos)
    axes.set_xticklabels(ref_name_ls, rotation='vertical', va='top')
    axes.set_yticklabels(query_name_ls, rotation='horizontal', ha='right')
    # Setup legend
    legend_patch = [mpatches.Patch(
        color=com_group_color[group],
        label=com_group_label[group])
                    for group in range(com_group_size)]
    axes.legend(handles=legend_patch, loc='upper center',
                bbox_to_anchor=(0, 1.02), fancybox=True, ncol=1)

    axes.set_title("%d bp-bin comparison between %s and %s"
                   %(window, ref_genome_name, query_genome_name))

    fig_name = "%s%s%s-%s-bin.png" %(dirn, prefix, ref_genome_name, query_genome_name)
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(fig_name, bbox_inches='tight', format='png')
    plt.close()




def grouped_bin_list(position_list, reference_name, query_name,
                     dirn, prefix, gap=200):
    """Grouped Bin list."""
    _message("Generating Annotated bin list")
    position_list.sort(key=itemgetter(0,1,3,4))
    com_pos_list = []
    def add_bin2_pos_list(pos_list_, bin_):
        """Add the bin to the combine position list."""
        # All the groups in the current bin
        _groups = bin_[1]
        # No LOH in both strains
        if (1 not in _groups) and (2 not in _groups):
            # Group 0: Unique paired
            if 0 in _groups:
                pos_list_.append(bin_[0] + (0,))
            # Group 4: Homologous region
            else:
                pos_list_.append(bin_[0] + (4,))
        # LOH in reference
        elif (1 in _groups) and (2 not in _groups):
            pos_list_.append(bin_[0] + (1,))
        # LOH in query
        elif (2 in _groups) and (1 not in _groups):
            pos_list_.append(bin_[0] + (2,))
        # LOH in both strains, gene convertion
        else:
            pos_list_.append(bin_[0] + (3,))

    _curr_bin = ()
    _last_pos = position_list[-1]
    for _pos in position_list:
        _curr_pos = (_pos[0], _pos[1], _pos[3], _pos[4])
        _curr_group = _pos[2]
        # Mask kmer only in one strain
        if _curr_group == 4:
            continue
        if not _curr_bin:
            _curr_bin = (_curr_pos, [_curr_group])
        else:
            if _curr_bin[0] == _curr_pos:
                _curr_bin[1].append(_curr_group)
                if  _pos == _last_pos:
                    add_bin2_pos_list(com_pos_list, _curr_bin)
            else:
                add_bin2_pos_list(com_pos_list, _curr_bin)
                _curr_bin = (_curr_pos, [_curr_group])
    _message("Finished generation of bin list")

    def unqiue_filter(pos_list):
        """Filter Multi-Align Unique Bins."""
        unique_paired = [pos for pos in pos_list if pos[4] == 0]
        _message("Start Filtering the unique bins")
        single_paired = []
        mix_pairted = []
        # single dictionary (ref_id, ref_pos) -> (query_id, query_pos)
        single_ref_dict = {}
        # single dictionary (query_id, query_pos) -> (ref_id, ref_pos)
        single_query_dict = {}
        # mix ref dictionary (ref_id, ref_pos) -> None
        mix_ref_dict = {}
        mix_query_dict = {}

        for pos in unique_paired:
            ref_pos = (pos[0], pos[2])
            query_pos = (pos[1], pos[3])

            if ref_pos in mix_ref_dict or query_pos in mix_query_dict:
                mix_pairted.append(pos)
                continue

            if ref_pos in single_ref_dict:
                query_pos_p = single_ref_dict[ref_pos]
                del single_ref_dict[ref_pos]
                del single_query_dict[query_pos_p]



    homologous = [pos for pos in com_pos_list if pos[4] == 4]
    loh_r = [pos for pos in com_pos_list if pos[4] == 1]
    loh_q = [pos for pos in com_pos_list if pos[4] == 2]
    gene_con = [pos for pos in com_pos_list if pos[4] == 3]

    # Sort by ref_id -> ref_pos -> query_id -> query_pos
    unique_paired.sort(key=itemgetter(0, 2, 1, 3))
    homologous.sort(key=itemgetter(0, 2, 1, 3))
    loh_r.sort(key=itemgetter(0, 2, 1, 3))
    loh_q.sort(key=itemgetter(0, 2, 1, 3))
    gene_con.sort(key=itemgetter(0, 2, 1, 3))
    _message("%d Unique paired region" % len(unique_paired))
    _message("%d homologous region" % len(homologous))
    _message("%d Loss of homology in %s" %(len(loh_r), reference_name))
    _message("%d Loss of homology in %s" %(len(loh_q), query_name))
    _message("%d Gene Conversion" %len(gene_con))
    # _message("Writing Position file")
    # with gzip.open(dirn+prefix+"-up.txt.gz", "w") as fop:
    #     _lines_list = []
    #     for pos in unique_paired:
    #         _line = "%d\t%d\t%d\t%d\t%d\n" % pos
    #         _lines_list.append(_line)
    #     fop.write("".join(_lines_list))

    # with gzip.open(dirn+prefix+"-gc.txt.gz", "w") as fop:
    #     _lines_list = []
    #     for pos in gene_con:
    #         _line = "%d\t%d\t%d\t%d\t%d\n" % pos
    #         _lines_list.append(_line)
    #     fop.write("".join(_lines_list))

    # with gzip.open(dirn+prefix+"-hg.txt.gz", "w") as fop:
    #     _lines_list = []
    #     for pos in homologous:
    #         _line = "%d\t%d\t%d\t%d\t%d\n" % pos
    #         _lines_list.append(_line)
    #     fop.write("".join(_lines_list))
    # with gzip.open(dirn+prefix+"-lr.txt.gz", "w") as fop:
    #     _lines_list = []
    #     for pos in loh_r:
    #         _line = "%d\t%d\t%d\t%d\t%d\n" % pos
    #         _lines_list.append(_line)
    #     fop.write("".join(_lines_list))

    # with gzip.open(dirn+prefix+"-lq.txt.gz", "w") as fop:
    #     _lines_list = []
    #     for pos in loh_q:
    #         _line = "%d\t%d\t%d\t%d\t%d\n" % pos
    #         _lines_list.append(_line)
    #     fop.write("".join(_lines_list))


    def number_cluster(pos_list, start_count):
        """Numbering clusters
        cluster_dict r: {(ref_id, ref_pos) -> # cluster}
                  q: {(query_id, query_pos) -> # cluster}

        cluster list [[No. Cluster, [(ref_id, ref_pos)], [(query_id, query_pos)]]]
        """
        cluster_r = {}
        cluster_q = {}
        count = start_count
        for pos in pos_list:
            # Update dict
            ref_id, query_id, ref_pos, query_pos, group = pos
            # New cluster
            if (ref_id, ref_pos) not in cluster_r and (query_id, query_pos) not in cluster_q:
                count += 1
                cluster_r[(ref_id, ref_pos)] = count
                cluster_q[(query_id, query_pos)] = count
            # New member in the cluster
            elif (ref_id, ref_pos) in cluster_r and (query_id, query_pos) not in cluster_q:
                _cluster = cluster_r[(ref_id, ref_pos)]
                cluster_q[(query_id, query_pos)] = _cluster
            else:
                _cluster = cluster_q[(query_id, query_pos)]
                cluster_r[(ref_id, ref_pos)] = _cluster
        # Build cluster list

        cluster_size = count - start_count
        cluster_list = [[_id + start_count +1, [], []] for _id in xrange(cluster_size)]
        for pos, cluster in cluster_r.iteritems():
            _id = cluster - start_count -1
            cluster_list[_id][1].append(pos)

        for pos, cluster in cluster_q.iteritems():
            _id = cluster - start_count -1
            cluster_list[_id][2].append(pos)

        return cluster_list, count

    _message("Numbering the clusters")
    cluster_start = 0
    # Number homology group
    homo_cluster, homo_count = number_cluster(homologous, cluster_start)
    _message("%d Homologous binned clusters" %(homo_count - cluster_start))
    # Number Gene Conversion
    gc_cluster, gc_count = number_cluster(gene_con, homo_count)
    _message("%d Gene Conversion binned clusters" %(gc_count - homo_count))
    # Number LOH in ref
    lohr_cluster, lohr_count = number_cluster(loh_r, gc_count)
    _message("%d LOH ref binned clusters" %(lohr_count - gc_count))
    # Number LOH in query
    lohq_cluster, lohq_count = number_cluster(loh_q, lohr_count)
    _message("%d LOH query binned clusters" %(lohq_count - lohr_count))
    _message("Generating clusted bin list")
    # Generate position list
    # [(ref_id, ref_pos, query_id, query_pos, cluster_id)]
    # Cluster_start = [list of start cluster number]
    cluster_pos_list = []
    cluster_start_list = [0, homo_count, gc_count, lohr_count, lohq_count]


    for pos in unique_paired:
        cluster_pos_list.append(pos.itemgetter(0, 2, 1, 3, 4))
    for cluster in homo_cluster:
        _cluster = cluster[0]
        for ref in cluster[1]:
            for query in cluster[2]:
                cluster_pos_list.append(ref + query + (_cluster,))
    for cluster in gc_cluster:
        _cluster = cluster[0]
        for ref in cluster[1]:
            for query in cluster[2]:
                cluster_pos_list.append(ref + query + (_cluster,))
    for cluster in lohr_cluster:
        _cluster = cluster[0]
        for ref in cluster[1]:
            for query in cluster[2]:
                cluster_pos_list.append(ref + query + (_cluster,))
    for cluster in lohq_cluster:
        _cluster = cluster[0]
        for ref in cluster[1]:
            for query in cluster[2]:
                cluster_pos_list.append(ref + query + (_cluster,))
    _message("Writing bin list")
    with gzip.open(dirn+prefix+"-cluster.txt.gz", 'w') as fop:
        _line_list = []
        for pos in cluster_pos_list:
            # ref_id, ref_pos, query_id, query_pos
            if len(pos) != 5:
                print pos
            _line = "%d\t%d\t%d\t%d\t%d\n" % pos
            _line_list.append(_line)
        fop.write("".join(_line_list))



    return cluster_pos_list, cluster_start_list



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
            if group == 4 or group == 3:
                continue
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
            if group == 4 or group == 3:
                continue
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
    _name_list = fasta_.Names
    _seq_list = fasta_.Seqs
    # Transform to lowercase
    pattern = r"_(\w)_(\w)"
    for _id in range(len(_seq_list)):
        _seq_list[_id] = _seq_list[_id].lower()
    if pattern:
        def _get_seq_tag(name):
            """Return the tags from reg pattern"""
            _match = re.search(pattern, name)
            _tag_ls = []
            if not _match:
                return []
            if not _match.groups():
                return [_match.group()]
            else:
                for _id in range(len(_match.groups())):
                    _tag_ls.append(_match.groups()[_id])
            return _tag_ls
        # Filter the sequences with no match
        _filtered_name_list = [_get_seq_tag(_name) for _name in _name_list]
        if _filtered_name_list:
            _filtered_tuple = [_tuple for _tuple in
                               zip(_filtered_name_list, _name_list, _seq_list)
                               if _tuple[0]]
            _filtered_tuple.sort()
            _tmp, _name_list, _seq_list = zip(*_filtered_tuple)
            return _name_list, _seq_list
        return None, None



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
        "G": "C",
        "N": "N"
    }
    rc_nt_ls = [rc_dict[seq[i]] for i in range(len(seq)-1, -1, -1)]
    rc_seq_ = "".join(rc_nt_ls)
    return rc_seq_



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

def paired_aligntwo2one(pos_list, dirn, prefix="under2"):
    """Filter and graph the 2-1 and 1-1 bins.
    The masked bin from mask_position_list is filtered for:

    Args:
    pos_list: masked postion list
    Returns:
    None, generates graphs in dirn/prefix
    """

    color_group = [
        'red', # one to one alignment, unpaired
        'green', # two ref to one query, one paired
        'blue', # one ref to two query, one paired
        'grey', # one to one paired
    ]


    for pos in pos_list:
        ref_id = pos[0]
        query_id = pos[1]





def test():
    """Temporary fucntion to test code."""
    fname_b = "/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa"
    fname_c = "/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa"
    # fname_b = "/home/zhuwei/cglabarata/comp/kmer/BG2_EPA67.fa"
    # fname_c = "/home/zhuwei/cglabarata/comp/kmer/CBS_EPA67.fa"
    pattern = r"_(\w)_(\w)"
    # pattern = ""
    names_b, seqs_b = load_fasta(fname_b, pattern=pattern)
    names_c, seqs_c = load_fasta(fname_c, pattern=pattern)
    tab_b = generate_kmer_table(seq_list=seqs_b)
    tab_c = generate_kmer_table(seq_list=seqs_c)
    tab_comp, unalign_r, unalign_q = compare_kmer_table(tab_b, tab_c)
    pos_ls = transform_kmer_comp_table(tab_comp, unalign_r, unalign_q)
    mask_pos_ls = mask_position_list(pos_ls, window=50)

    len_b = []
    len_c = []
    for seq in seqs_b:
        len_b.append(len(seq))
    for seq in seqs_c:
        len_c.append(len(seq))
    # draw_chrom_homo_matrix(
    #     ref_name_ls=names_b, query_name_ls=names_c,
    #     ref_len_ls=len_b, query_len_ls=len_c,
    #     position_list=mask_pos_ls, dirn="/home/zhuwei/cglabarata/comp/kmer/",
    #     prefix="EPA67/",
    #     ref_genome_name="BG2",
    #     query_genome_name="CBS",
    #     window=200)

    grouped_bin_list(
        position_list=mask_pos_ls,
        reference_name="BG2",
        query_name="CBS",
        dirn="/home/zhuwei/cglabarata/comp/kmer/EPA67/",
        prefix="bg2-cbs")
    # draw_chrome_comp_graph(
    #     ref_name_ls=names_b, query_name_ls=names_c,
    #     ref_len_ls=len_b, query_len_ls=len_c,
    #     position_list=mask_pos_ls, dirn="/home/zhuwei/cglabarata/comp/kmer/",
    #     prefix="EPA67/",
    #     ref_genome_name="BG2",
    #     query_genome_name="CBS")
    # compare_seq(
    #     ref_seq_ls=seqs_b,
    #     query_seq_ls=seqs_c,
    #     prefix="EPA67", dirn="/home/zhuwei/cglabarata/comp/kmer/trans-table-1/")

if __name__ == "__main__":
    test()