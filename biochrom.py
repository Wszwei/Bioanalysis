import os
from fasta import Fasta
from sequence_analysis import Seq_Analyzer
from subprocess import call
import re
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import gzip
import time
import threading
from mpl_toolkits.mplot3d import Axes3D



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

def get_all_kmer_position(kmer_table, seq="", len_kmer=20, size=1, index=0):
    """Obtain all the kmers from a sequence. 
    Every kmer is added to a kmer_table {"kmer": list of [[list of positions]]}
    Args:
        kmer_table: kmer dictionary
        seq: sequence to get all kmers
        len_kmer: length of the kmer to be checked
        size: total number of entries
        index: index of current sequence entry
    Returns:
        None, kmer_table is updated within the function
    Raises:
        None. 
    """
    for i in xrange(len(seq) - (len_kmer - 1)):
        kmer = seq[i:i + len_kmer].upper()
        kmer_table.setdefault(kmer, [[] for i in range(size)])
        kmer_table[kmer][index].append(i+1) # Position of kmer in the entry seq



def kmer_align(ref_dirn="", query_dirn="", ext="fa",
	           len_kmer=20, prefix="aligned"):
    """Align the query fasta file to the referecne fasta file.

    All the fasta file in the query with the same extention is aligned to that 
    of the reference directory.

    Args:
        ref_dirn: directory of the reference fasta files
        query_dirn: directory of query fasta files
        ext: the file extention of fasta files (default = fa)
        len_kmer: length of the kmer to align the sequnces
        prefix: output of the query files: {prefix}-{query_fasta_name}

    Returns:
        1 for all right
        make a new directory under query directory
        {query_dirn}/{prefix}
        Write all the files in this directory
        Write fasta files in the query directory with names:
        {prefix}-{query} in P={prefix}.{ext} file
        Write the kmer cooridinate data in
        {prefix}-kmer-coordinates.txt.gz 
        Write the query to ref cooridinate data in
        {prefix}-align-coorindates.txt.gz
        Write the kmer only find in one dataset in
        {prefix}-ref-unaligned.txt.gz 
        {prefix}-query-unaligned.txt.gz
        Write the summary file in
        {prefix}-summary.txt

    Raises:
        -1 for wrong directory
    """
    
    # List of fasta file information
    ref_seq_ls = []
    ref_name_ls = []
    query_seq_ls = []
    query_name_ls = []

    # List of kmer table
    ref_kmer_table = {}
    query_kmer_table = {}

    # List of aligned kmer postions
    # each alignment is recorded as [(ref, query, ref_id, query_id, group)]
    # ref and query: 
    aligned = []

    # output file name
    op_fasta = "%s.%s" %(prefix, ext)
    op_coord = "%s-coordinates.txt.gz" %prefix

    # Load fasta files 
    if not (_test_dir(ref_dirn) and _test_dir(query_dirn)):
        return -1

    __len_ext = len(ext)

    for fname in os.listdir(ref_dirn):
        if fname[-__len_ext:].lower() == ext.lower():
            cur_fa = Fasta(ref_dirn + fname)
            ref_seq_ls.extend(cur_fa.Seqs)
            ref_name_ls.extend(cur_fa.Names)

    ref_size = len(ref_name_ls)
    _message("%d Reference Entries loaded" %ref_size)

    for fname in os.listdir(query_dirn):
        if fname[-__len_ext:].lower() == ext.lower():
            cur_fa = Fasta(query_dirn + fname)
            query_seq_ls.extend(cur_fa.Seqs)
            query_name_ls.extend(cur_fa.Names)

    query_size = len(query_name_ls)
    _message("%d Reference Entries loaded" %query_size)

    # Generate kmer coordinate table
    # Dictionary
    # kmer: 1D {size} array of list with all the postions
    # [] list for no alignment 
    for i in range(ref_size):
        _seq = ref_seq_ls[i]
        get_all_kmer_position(kmer_table=ref_kmer_table,
                              len_kmer=len_kmer,
                              size=ref_size,
                              index=i)

    for i in range(query_size):
        _seq = query_seq_ls[i]
        get_all_kmer_position(kmer_table=query_kmer_table,
                              len_kmer=len_kmer,
                              size=query_size,
                              index=i)


    # Generate coordinates for kmer counting
    # Write the kmer information
    # Generate ref to query position list
