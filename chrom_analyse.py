import os
import numpy
from fasta import Fasta
from sequence_analysis import Seq_Analyzer
from subprocess import call
import re
from Bio import pairwise2
# import matplotlib.pyplot as plt
# import pylab as plt
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines
import gzip
import time
import threading
from mpl_toolkits.mplot3d import Axes3D
# from fpdf import FPDF



class coord_record:
    def __init__(self, line=""):
        self.ref_s = -1
        self.ref_e = -1
        self.query_s = -1
        self.query_e = -1
        self.isReverse = False
        self.len_align_r = -1
        self.len_align_q = -1
        self.idy = -1
        self.ref_len = -1
        self.cov_r = -1
        self.cov_q = -1
        self.ref_name = ""
        self.query_name = ""
        if line:
            self.read_record(line)

    def read_record(self, line):
        pattern = r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+([0-9a-zA-Z]+)\s+([0-9a-zA-Z]+)"
        align = re.match(pattern, line)
        if align:
            self.ref_s = int(align.group(1))
            self.ref_e = int(align.group(2))
            self.query_s = int(align.group(3))
            self.query_e = int(align.group(4))
            if (self.ref_e - self.ref_s) * (self.query_e - self.query_s) < 0:
                self.isReverse = True
            self.len_align_r = int(align.group(5))
            self.len_align_q = int(align.group(6))
            self.idy = float(align.group(7))
            self.ref_len = int(align.group(8))
            self.query_len = int(align.group(9))
            self.cov_r = float(align.group(10))
            self.cov_q = float(align.group(11))
            self.ref_name = align.group(12)
            self.query_name = align.group(13)
        else:
            print("Error in loading recoords from .coords file\n")
            print("Line:\n")
            print(line)


class coords:
    def __init__(self, fn=""):
        self.fn = ""
        self.record_ls = []
        if fn:
            self.read_coords(fn)

    def read_coords(self, fn):
        pattern = r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+([0-9a-zA-Z]+)\s+([0-9a-zA-Z]+)"
        with open(fn, "r") as fp:
            for line in fp:
                _match = re.match(pattern, line)
                if _match:
                    curr_record = coord_record(line=line)
                    self.record_ls.append(curr_record)
            print ("%d records loads from %s\n" % (len(self.record_ls), fn))


def read_coord_file(fn):
    if os.path.isfile(fn):
        my_coord = coords(fn=fn)
        return my_coord
    else:
        print "Cannot read %s\n" % fn
        return None

def verify_r_q_d(ref="", query="", dirn=""):    
    if not os.path.isdir(dirn):
        print "Working directory is not available\n"
        return False
    if not os.path.isfile(dirn+ref):
        print "Reference file is not available\n"
        return False
    if not os.path.isfile(dirn+query):
        print "Query file is not available\n"
        return False
    return True

def generate_coords(ref="", query="", prefix="compare", dirn=""):

    if not verify_r_q_d(ref=ref, query=query, dirn=dirn):
        return None

    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    os.system("nucmer -p %s %s %s" %(prefix, dirn+ref, dirn+query))
    os.system("mv %s.* %s" %(prefix, dirn+prefix))
    os.system("cp %s %s" %(dirn+ref, dirn+prefix))
    os.system("cp %s %s" %(dirn+query, dirn+prefix))
    os.system("show-coords -crlT %s.delta > %s.coords" %(dirn+prefix+"/"+prefix, dirn+prefix+"/"+prefix))

    return dirn+prefix+"/", prefix+".coords"


def assign_query_to_ref(ref="", query="", dirn="", prefix="compare"):
    if not verify_r_q_d(ref=ref, query=query, dirn=dirn):
        return None
    comp_dir, comp_file = generate_coords(ref=ref, query=query, dirn=dirn, prefix=perfix)
    comp_fn = comp_dir + comp_file
    if not os.path.isfile(fn):
        print ()

def generate_dir_seq_from_labeled_seq(fn="", prefix="out"):
    if not os.path.isfile(fn):
        print "Labeled fasta file cannot be loaded\n"
        print fn
        return
    fa = Fasta(fn)
    op = open(os.path.dirname(fn)+"/"+prefix+"_dir.fa", "w")
    for i in range(len(fa.Names)):
        _name = fa.Names[i]
        if _name[:3].upper() != "CHR":
            continue
        if _name[-2:].upper() == "-R":
            op.write(">%s\n" % _name[:-2])
            _seq = Seq_Analyzer(fa.Seqs[i]).rcSeq()
            op.write(_seq)
            op.write("\n")
        else:
            op.write(">%s\n" %_name)
            op.write(fa.Seqs[i])
            op.write("\n")
    op.close()
    print "Directional fasta file generated"

def extract_ends(fn="", length=50000, prefix="ends"):
    if not os.path.isfile(fn):
        print "Labeled fasta file cannot be loaded\n"
        print fn
        return
    fa = Fasta(fn)    
    dirn = os.path.dirname(fn) + "/"
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    dirn = dirn+prefix+"/"
    for i in range(len(fa.Names)):
        _name = fa.Names[i]
        _seq = fa.Seqs[i]
        if len(_seq) > (length * 2):
            with open(dirn+_name+"_left.fa", "w") as op:
                op.write(">%s\n" % (_name+"_Left"))
                op.write(_seq[:length+1])
                op.write("\n")
            with open(dirn+_name+"_right.fa", "w") as op:
                op.write(">%s\n" % (_name+"_Right"))
                op.write(_seq[-length:])
                op.write("\n")
        else:
            with open(dirn+_name+".fa", "w") as op:
                op.write(">%s\n" % (_name))
                op.write(_seq)
                op.write("\n")


def pair_wise_alignment(seq1="", seq2=""):
    align_score = pairwise2.align.localxx(seq1, seq2, score_only=True)
    return align_score[0]

def extract_ends_2(fn="", length=50000, prefix="ends"):
    if not os.path.isfile(fn):
        print "Labeled fasta file cannot be loaded\n"
        print fn
        return
    fa = Fasta(fn)    
    dirn = os.path.dirname(fn) + "/"
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    dirn = dirn+prefix+"/"
    for i in range(len(fa.Names)):
        _name = fa.Names[i]
        _seq = fa.Seqs[i]
        if len(_seq) > (length * 2):
            # if "LEFT" in _name.upper():
            #     with open(dirn+_name+"_left.fa", "w") as op:
            #         op.write(">%s\n" % (_name+"_Left"))
            #         op.write(_seq[:length+1])
            #         op.write("\n")
            # else:
            with open(dirn+_name+"_right.fa", "w") as op:
                op.write(">%s\n" % (_name+"_Right"))
                op.write(_seq[-length:])
                op.write("\n")
        else:
            with open(dirn+_name+".fa", "w") as op:
                op.write(">%s\n" % (_name))
                op.write(_seq)
                op.write("\n")


def extract_single_fasta(prefix="single", fn=""):
    if not os.path.isfile(fn):
        print "Labeled fasta file cannot be loaded\n"
        print fn
        return
    fa = Fasta(fn)    
    dirn = os.path.dirname(fn) + "/"
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    dirn = dirn+prefix+"/"
    for i in range(len(fa.Names)):
        _name = fa.Names[i]
        _seq = fa.Seqs[i]
        with open (dirn+_name+".fa", "w") as op:
            op.write(">%s\n" %(_name))
            op.write(_seq)
            op.write("\n")

def kmer_count_dirn(dirn="", prefix="kmer", kmer=20, low=1):
    '''Kmer counting over all the fa file in one folder.
    Jellyfish is called to perform the counting'''
    if not os.path.isdir(dirn):
        print "Directory is not available\n"
        print dirn
        return
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)
    dirn_pre = dirn+prefix+"/"

    for fn in os.listdir(dirn):
        if fn[-2:].upper() == "FA":
            _name = fn[:-3]

            os.system("jellyfish count -m %d -s 100M -t 10 %s" %(kmer, dirn+fn))
            os.system("jellyfish histo mer_counts.jf > %s.histo" %(dirn+_name))
            os.system("jellyfish dump -L %d -ct mer_counts.jf > %s-dump.txt" %(low, dirn+_name))
            os.system("rm mer_counts.jf")
            os.system("mv %s.histo %s" %(dirn+_name, dirn_pre))
            os.system("mv %s-dump.txt %s" %(dirn+_name, dirn_pre))

def generate_kmer_table(dirn="", prefix="output", count_exist_only=False, write=True):
    if not os.path.isdir(dirn):
        print "Directory is not available\n"
        print dirn
        return None

    kmer_table = {}
    seq_list = []

    pattern = r"(\w+)\s+(\d+)"
    count = -1
    
    for fn in os.listdir(dirn):
        if not "dump" in fn:
            continue
        seq_list.append(fn[:-9])

    with open(dirn + "%s-seq_list.txt" %prefix, "w") as fp:
        for x in seq_list:
            fp.write("%s\n" %x)
    
    size = len(seq_list)
    
    count = -1
    for fn in os.listdir(dirn):
        if not "dump" in fn:
            continue
        count += 1

        with open(dirn+fn, "r") as fp:
            for line in fp:
                _match = re.match(pattern, line)
                _kmer = _match.group(1)
                _num = int(_match.group(2))
                kmer_table.setdefault(_kmer, [0] * size)
                if count_exist_only:
                    kmer_table[_kmer][count] += 1
                else:                    
                    kmer_table[_kmer][count] += _num

    # Write kmer table
    if write:
        with open(dirn+"%s-u-table.txt" %prefix, "w") as op_u:
            with open(dirn+"%s-m-table.txt" %prefix, "w") as op_m: 
                op_m.write("kmer  ")
                for seq in seq_list:
                    op_m.write("%s  " %seq[0])
                op_m.write("\n")

                op_u.write("kmer  ")
                for seq in seq_list:
                    op_u.write("%s  " %seq[0])
                op_u.write("\n")
                for kmer, num in kmer_table.items():
                    if num.count(0) == size -1:
                        op_u.write("%s  " %kmer)
                        for i in num:
                            op_u.write("%d  " %i)
                        op_u.write("\n")                
                    else:
                        op_m.write("%s  " %kmer)
                        for i in num:
                            op_m.write("%d  " %i)
                        op_m.write("\n")
    return seq_list, kmer_table


def compare_kmer_bg2_cbs_v2():
    # Compare kmer in BG2 with that in CBS
    # Chr-Specific alignment, multi alignment and self alignment are labeled
    # Alignment Coding: u = unique-alignment, s = self-alignment

    bg2_seq_list, bg2_kmer_table = generate_kmer_table(dirn="/home/zhuwei/cglabarata/comp/kmer/bg2/kmer/", prefix="bg2")
    cbs_seq_list, cbs_kmer_table = generate_kmer_table(dirn="/home/zhuwei/cglabarata/comp/kmer/cbs/kmer/", prefix="cbs")

    pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    size = len(bg2_seq_list)
    # Generate Map between bg2 and cbs

    bg2_seq_num_ls = []
    cbs_seq_num_ls = []

    translocate_event = [[0,""] for i in range(size**2)]

    # op = open("/home/zhuwei/cglabarata/comp/kmer/bg2-cbs-unique20mer-mask-translocation.txt", "w")
    op = open("/home/zhuwei/cglabarata/comp/kmer/bg2-cbs-unique20mer-translocation-v2.txt", "w")
    op_s = open("/home/zhuwei/cglabarata/comp/kmer/bg2-cbs-unique20mer-translocation-summary-v2.txt", "w")


    def chr2id(chr,end):
        _id = (ord(chr) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname

    # Sequence Name -> id Number
    for i in bg2_seq_list:
        _name = i[0]
        _match = re.match(pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        bg2_seq_num_ls.append(_id)


    for i in cbs_seq_list:
        _name = i[0]
        _match = re.match(pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        cbs_seq_num_ls.append(_id)

    op.write("kmer  BG2_Chr  CBS_Chr  BG2_count  CBS_count  Alignment_Code\n")


    for kmer, value in bg2_kmer_table.items():
        if value.count(0) == size -1: # Kmer specific to one fosmid in BG2
            if kmer in cbs_kmer_table:
                _bg2_pos = [i for i, e in enumerate(value) if e!=0]
                _cbs_value = cbs_kmer_table[kmer]
                _cbs_pos = [i for i, e in enumerate(_cbs_value) if e!=0]
                for p in _bg2_pos:
                    for q in _cbs_pos:
                        _bg2_id = bg2_seq_num_ls[p]
                        _cbs_id = cbs_seq_num_ls[q]
                        _bg2_chr = id2chr(_bg2_id)
                        _cbs_chr = id2chr(_cbs_id)
                        translocate_event[_cbs_id*size+_bg2_id][0]+= 1
                        _code = ""
                        if len(_bg2_pos) == 1 and len(_cbs_pos) == 1: # Chr-Specific alignment
                            translocate_event[_cbs_id*size+_bg2_id][1] += "u"
                            _code += 'u'
                        if _bg2_id == _cbs_id:
                            translocate_event[_cbs_id*size+_bg2_id][1] += "s"
                            _code += 's'                                                                        
                        op.write("%s  %s  %s  %d  %d %s\n" %(kmer, _bg2_chr, _cbs_chr, value[p], _cbs_value[q], _code))
    """Write summary file
       Alignemnt coding
       u = unique
       s = self- alignemnt""" 
    # Write specific 20 mer
    op_s.write("Chr-Specific 20 mer\n")
    for i in range(size):
        op_s.write("      BG2_%s" %id2chr(i))
    op_s.write("\n")

    for i in range(size):
        op_s.write("CBS_%s" %id2chr(i))
        for j in range(size):
            event = translocate_event[i*size+j]
            unique_event = event[1].count("u") 
            op_s.write("  %d" %unique_event)
        op_s.write("\n")

    op_s.write("\n\n\n")

    op_s.write("Multi-alignment 20 mer\n")
    for i in range(size):
        op_s.write("      BG2_%s" %id2chr(i))
    op_s.write("\n")

    for i in range(size):
        op_s.write("CBS_%s" %id2chr(i))
        for j in range(size):
            event = translocate_event[i*size+j]
            unique_event = event[1].count("u")
            multi_event = event[0] - unique_event
            op_s.write("  %d" %multi_event)
        op_s.write("\n")

    op.close()
    op_s.close()


def find_all_seq(seq="", ref=""):
    """ Return the position of all the seqs"""
    seq = seq.upper()
    ref = ref.upper()
    pos_ls = [m.start() for m in re.finditer('(?=%s)' %seq, ref)]
    return pos_ls


def compare_kmer(name_a="A", name_b="B", dirn_a="", dirn_b="", dirn="", prefix="kmer-comp", write_kmer_table=True):
    # Compare kmer in two Chromosomes
    # The kmer data are read from dump file by jellyfish
    # Group Label:
    # 0 = Chr-Specific in both strains and located in the same chromosome
    # 1 = Chr-Specific in both strains but located in different chromosomes
    # 2 = Chr-Specific in strain of chr-A, multiple in that of chr-B, still have at least one copy in the same chromosome
    # 3 = Chr-Specific in strain of chr-A, multiple in that of chr-B, but not located in the same chromosome
    # 4 = Chr-Specific in strain of chr-B, multiple in that of chr-A, still have at least one copy in the same chromosome
    # 5 = Chr-Specific in strain of chr-B, multiple in that of chr-A, but not located in the same chromosome
    # 6 = Multiple in both A and B, the located chromosomes are the same
    # 7 = Multiple in both A and B, extra alignment(s) in A
    # 8 = Multiple in both A and B, extra alignment(s) in B
    # 9 = Multiple in both A and B, extra alignemnt(s) in both
    # 10 = kmer only found in A
    # 11 = kmer only found in B
    # Specific Name pattern is used for self-alignment
    # pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    # {Strain-Name}_{Fosmid-Name=BXXX}_{Chr_Name}_{Left/Right}

    a_seq_list, a_kmer_table = generate_kmer_table(dirn=dirn_a, prefix=name_a, write=write_kmer_table)
    b_seq_list, b_kmer_table = generate_kmer_table(dirn=dirn_b, prefix=name_b, write=write_kmer_table)
    print "kmer table for %s and %s is generated!" % (name_a, name_b)
    if not os.path.isdir(dirn):
        print "Working directory not available:\n"
        print dirn

    pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    size_a = len(a_seq_list)
    size_b = len(b_seq_list)
    # Generate Map between strain a and b

    a_seq_num_ls = []
    b_seq_num_ls = []

    alignment_event = [[0,""] for i in range(size_a * size_b)]

    unalignment_a_event = [[0] for i in range(size_a)]
    unalignment_b_event = [[0] for i in range(size_b)]

    op = open(dirn+prefix+".txt", "w")
    op_s = open(dirn+prefix+"-summary.txt", "w")
    op_ua = open(dirn+prefix+"-unalign-%s.txt" %(name_a), "w")
    op_ub = open(dirn+prefix+"-unalign-%s.txt" %(name_b), "w")
    op_us = open(dirn+prefix+"-unalign-summary.txt", "w")
    op_stat = open(dirn+prefix+"-stat.txt", "w")


    def chr2id(chr,end):
        _id = (ord(chr) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname

    # Sequence Name -> id Number
    for i in a_seq_list:
        _name = i
        _match = re.match(pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        a_seq_num_ls.append(_id)


    for i in b_seq_list:
        _name = i
        _match = re.match(pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        b_seq_num_ls.append(_id)

    op.write("kmer  %s_Chr  %s_Chr  %s_count  %s_count  Group\n" %(name_a, name_b, name_a, name_b))
    op_ua.write("kmer %s_Chr %s_count\n" %(name_a, name_a))
    op_ub.write("kmer %s_Chr %s_count\n" %(name_b, name_b))
    op_stat.write("Statistics of Kmer alignment\n")

    group_count = [0] * 12

    for a_kmer, a_value in a_kmer_table.items():
        # Shared kmer between A and B
        if a_kmer in b_kmer_table:
            a_chr = [i for i, e in enumerate(a_value) if e!=0]
            b_value = b_kmer_table[a_kmer]
            b_chr = [i for i, e in enumerate(b_value) if e!=0]
            a_id_ls = []
            b_id_ls = []
            for p in a_chr:
                for q in b_chr:
                    a_id = a_seq_num_ls[p]
                    b_id = b_seq_num_ls[q]
                    alignment_event[b_id*size_a+a_id][0] += 1
                    a_id_ls.append(a_id)
                    b_id_ls.append(b_id)
            if len(a_chr) == 1 and len(b_chr) == 1: # Chr-Specific alignment
                if a_id == b_id : # Aligned to the same chromosome
                    group = 0
                else:              # Aligned to difference chromosome
                    group = 1
            elif len(b_chr) > 1 and len(a_chr) == 1: # Chr-Specific in A, Multi in B
                if a_id in b_id_ls : # kmer shared in the same chromosome
                    group = 2
                else:              # kmer not shared in the same chromosome
                    group = 3
            elif len(a_chr) > 1 and len(b_chr) == 1:  # Chr-Specific in B, Multi in A
                if b_id in a_id_ls: # kmer shared in the same chromosome
                    group = 4
                else:
                    group = 5
            else: # Multi alignment in A and B
                extra = 0
                # extra = 0 A and B have same chr alignment
                # extra = 1 A has extra
                # extra = 2 B has extra
                # extra = 3 A and B have extra
                for index in a_id_ls:
                    if not index in b_id_ls:
                        extra += 1
                        break
                for index in b_id_ls:
                    if not index in a_id_ls:
                        extra += 2
                        break
                if extra >= 4:
                    print "Error comparing Multialingnments\n"
                    print a_id_ls
                    print b_id_ls
                    print extra
                group = 6 + extra
            for a_id in a_id_ls:
                for b_id in b_id_ls:
                    alignment_event[b_id * size_a + a_id][1] += str(group)
                    kmer_a_chr = id2chr(a_id)
                    kmer_b_chr = id2chr(b_id)
                    op.write("%s  %s  %s  %d  %d  %d\n" %(a_kmer, kmer_a_chr, kmer_b_chr,\
                                                          a_value[a_seq_num_ls.index(a_id)],\
                                                          b_value[b_seq_num_ls.index(b_id)],\
                                                          group))
            group_count[group] += 1
            del b_kmer_table[a_kmer]
        # kmer only in A
        else:
            a_id_ls = [a_seq_num_ls[i] for i, e in enumerate(a_value) if e != 0]
            for a_id in a_id_ls:
                # unalignment_a_event[a_id] += 1
                kmer_a_chr = id2chr(a_id)
                op_ua.write("%s  %s  %d\n" %(a_kmer, kmer_a_chr, a_value[a_seq_num_ls.index(a_id)]))
            group_count[10] += 1
    for b_kmer, b_value in b_kmer_table.items():
        b_id_ls = [a_seq_num_ls[i] for i, e in enumerate(b_value) if e != 0]
        for b_id in b_id_ls:
            # unalignment_b_event[b_id] += 1
            kmer_b_chr = id2chr(b_id)
            op_ub.write("%s  %s  %d\n" %(b_kmer, kmer_b_chr, b_value[b_seq_num_ls.index(b_id)]))
            group_count[11] += 1


    """Write summary file
       Alignemnt coding
       u = unique
       s = self- alignemnt""" 
    # Write specific 20 mer
    op_s.write("Chr-Specific 20 mer for single-chromosome translocation events\n")
    for i in range(size_a):
        op_s.write("  %s_%s" %(id2chr(i), name_a))
    op_s.write("\n")

    for i in range(size_b):
        op_s.write("%s_%s" %(id2chr(i), name_b))
        for j in range(size_a):
            event = alignment_event[i*size_a + j]
            unique_event = event[1].count('1') 
            op_s.write("  %d" %unique_event)
        op_s.write("\n")

    op_s.write("\n\n\n")
    op_s.write("Chr-Specific 20 mer for multi-chromosome translocation events (Counted on the specific chromosome)\n\n")
    
    op_s.write("%s\n" %name_a)

    for i in range(size_a):
        op_s.write("%s  " %(id2chr(i)))
    op_s.write("\n")

    for i in range(size_a):
        count = 0
        for j in range(size_b):
            event = alignment_event[j*size_a +i]
            count += event[1].count('3')
        op_s.write("%d  " %count)
    op_s.write("\n")
    op_s.write("%s\n" %name_b)

    for i in range(size_b):
        op_s.write("%s  " %(id2chr(i)))
    op_s.write("\n")

    for i in range(size_b):
        count = 0
        for j in range(size_a):
            event = alignment_event[i*size_a +j]
            count += event[1].count('5')
        op_s.write("%d  " %count)

    op_s.write("\n")

    op_s.write("Chr-Specific 20 mer for multi-chromosome duplication events (Counted on the specific chromosome)\n\n")
    
    op_s.write("%s\n" %name_a)

    for i in range(size_a):
        op_s.write("%s  " %(id2chr(i)))
    op_s.write("\n")

    for i in range(size_a):
        count = 0
        for j in range(size_b):
            event = alignment_event[j*size_a +i]
            count += event[1].count('2')
        op_s.write("%d  " %count)

    op_s.write("\n%s\n" %name_b)

    for i in range(size_b):
        op_s.write("%s  " %(id2chr(i)))
    op_s.write("\n")
    for i in range(size_b):
        count = 0
        for j in range(size_a):
            event = alignment_event[i*size_a +j]
            count += event[1].count('4')
        op_s.write("%d  " %count)
    
            
    op_s.write("\n")

    group_label = { 0: "Chr-Specific in Strain X and Y and located in the same chromosome",
                    1: "Chr-Specific in Strain X and Y but located in different chromosomes",
                    2: "Chr-Specific in strain X, Multiple in Y, share at least one copy in the same chromosome",
                    3: "Chr-Specific in Strain X, Multiple in Y, but Y don't align to same chromosome in X",
                    4: "Chr-Specific in strain Y, Multiple in X, share at least one copy in the same chromosome",
                    5: "Chr-Specific in Strain Y, Multiple in X, but X don't align to same chromosome in Y",
                    6: "Multiple in X and Y, the located chromosomes are the same",
                    7: "Multiple in X and Y, extra chromosome alignment(s) in X",
                    8: "Multiple in X and Y, extra chromosome alignment(s) in Y",
                    9: "Mulitple in X and Y, extra chromosome alignemnt(s) in both",
                    10: "Only found in X",
                    11: "Only found in Y"
    }

    for group, count in enumerate(group_count):
        op_stat.write("%s; %d\n" %(group_label[group], count))

    op.close()
    op_s.close()
    op_ua.close()
    op_ub.close()
    op_stat.close()


def bg2_cbs_locate_kmer(fn_kmer="", fn_bg2="", fn_cbs="", dirn=""):

    chr_pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    kmer_pattern = r"(\w+)\s+Chr(\w)_(\w)\s+Chr(\w)_(\w)"

    bg2tocbs_kmer_ls = [[] for i in range(676)] #676 = 26^2
    bg2tocbs_pos_ls = [[] for i in range(676)]


    cbs = Fasta(fn_cbs)
    bg2 = Fasta(fn_bg2)

    bg2_seq_num_ls = []
    cbs_seq_num_ls = []

    
    bg2_len_ls = []
    cbs_len_ls = []


    def chr2id(chr,end):
        _id = (ord(chr) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname

    with open(fn_kmer, "r") as fkmer:
        for line in fkmer:
            _match = re.match(kmer_pattern, line)
            if _match:
                kmer = _match.group(1)

                bg2_chr = _match.group(2)
                bg2_end = _match.group(3)
                cbs_chr = _match.group(4)
                cbs_end = _match.group(5)

                bg2_id = chr2id(bg2_chr, bg2_end)
                cbs_id = chr2id(cbs_chr, cbs_end)

                bg2tocbs_kmer_ls[cbs_id*26+bg2_id].append(kmer)
    for _name in bg2.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        bg2_seq_num_ls.append(_id)


    for _name in cbs.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        cbs_seq_num_ls.append(_id)


    for _seq in cbs.Seqs:
        _len = len(_seq)
        cbs_len_ls.append(_len)

    for _seq in bg2.Seqs:
        _len = len(_seq)
        bg2_len_ls.append(_len) 


    bg2_chr_align_ls = [[] for i in range(26)]

    for m in range(len(bg2.Names)):
        i = bg2_seq_num_ls[m]
        chr_align = []
        for n in range(len(cbs.Names)):

            j = cbs_seq_num_ls[n]
            if i==j:
                continue
            bg2tocbs_kmer = bg2tocbs_kmer_ls[j*26+i]
            print "%d kmers found in BG2-%s and CBS-%s" %(len(bg2tocbs_kmer), id2chr(i), id2chr(j))
            for kmer in bg2tocbs_kmer:
                pos_bg2 = find_all_seq(seq=kmer, ref=bg2.Seqs[m])
                pos_cbs = find_all_seq(seq=kmer, ref=cbs.Seqs[n])
                for pb in pos_bg2:
                    for pc in pos_cbs:
                        bg2tocbs_pos_ls[j*26+i].append((pb,pc))
                        chr_align.append((pb,j))
            print "BG2-%s-CBS-%s is treated!" %(id2chr(i), id2chr(j))
            if len(bg2tocbs_pos_ls[j*26+i]) == 0:
                continue
            with open(dirn+"BG2-%s-CBS-%s.txt" %(id2chr(i), id2chr(j)), "w") as op:
                op.write("BG2\t\tCBS\n")
                for pos in bg2tocbs_pos_ls[j*26+i]:
                    x = pos[0]
                    y = pos[1]
                    op.write("%d\t\t%d\n" %(x, y))
            if len(bg2tocbs_pos_ls[j*26+i]) == 1:
                x = bg2tocbs_pos_ls[j*26+i][0][0]
                y = bg2tocbs_pos_ls[j*26+i][0][1]
            else:
                x, y = zip(*bg2tocbs_pos_ls[j*26+i])
            plt.plot(x,y, 'ro')
            plt.xlabel('BG2-%s' %id2chr(i))
            plt.ylabel('CBS-%s' %id2chr(j))
            plt.title('Translocation of specific 20mer of BG2 fosmid to CBS fosmid')
            plt.grid(True)
            plt.axis([0,bg2_len_ls[m]+1000,0,cbs_len_ls[n]+1000])
            plt.savefig(dirn+"BG2-%s-CBS-%s" %(id2chr(i), id2chr(j)), dpi=600, format='png', bbox_inches='tight')
            plt.close()
        if len(chr_align) == 0:
            continue
        with open(dirn+"BG2-%s-CBS-fosmid-align.txt" %id2chr(i), "w") as op:
            op.write("BG2-Position\t\tCBS fosmid\n")
            for pos in chr_align:
                x = pos[0]
                y = pos[1]
                y_chr = id2chr(y)
                op.write("%d\t\t%s\n" %(x, y_chr))

        if len(chr_align) == 1:
            x = chr_align[0][0]
            y = chr_align[0][1]
        else:
            x, y = zip(*chr_align)
        f, ax = plt.subplots()
        ax.plot(x,y,'ro')
        plt.xlabel('BG2-%s' %id2chr(i))
        y_ticks_labels = [id2chr(t) for t in range(26)]
        y_ticks = np.arange(0,25,1)

        ax.set_yticks(y_ticks)
        ax.set_yticklabels(y_ticks_labels, rotation='horizontal', fontsize=14)
        plt.ylabel('CBS fosmids')
        plt.grid(True)
        plt.axis([0,bg2_len_ls[m]+1000,-1,26])
        plt.savefig(dirn+"BG2-%s-CBS-fosmid-align" %id2chr(i), dpi=600, format='png', bbox_inches='tight')
        plt.close()

def draw_chr_to_chr_comp_graph(pos_ls, name_a, name_b, len_a, len_b, dirn="", prefix='chr-chr', min_continous_set=2, isSorted=False):
    # Draw Chromosome-to-Chromosome Comparison
    
    # Format of pos_ls:
    # pos_ls = [(x,y,group)]
    # x, y are coordinates of a specific kmer
    # Group Label:
    # 0 = Chr-Specific in both strains and located in the same chromosome
    # 1 = Chr-Specific in both strains but located in different chromosomes
    # 2 = Chr-Specific in strain of chr-A, multiple in that of chr-B, still have at least one copy in the same chromosome
    # 3 = Chr-Specific in strain of chr-A, multiple in that of chr-B, but not located in the same chromosome
    # 4 = Chr-Specific in strain of chr-B, multiple in that of chr-A, still have at least one copy in the same chromosome
    # 5 = Chr-Specific in strain of chr-B, multiple in that of chr-A, but not located in the same chromosome
    # 6 = Multiple in both A and B, the located chromosomes are the same
    # 7 = Multiple in both A and B, extra alignment(s) in A
    # 8 = Multiple in both A and B, extra alignment(s) in B
    # 9 = Multiple in both A and B, extra alignment(s) in both

    # Define color of each group

    group_color = { 0: 'gray',
                    1: 'r', # red
                    2: 'orange',
                    3: 'salmon',
                    4: 'hotpink',
                    5: 'brown',
                    6: 'lightsteelblue',
                    7: 'lightgreen',
                    8: 'lawngreen',
                    9: 'yellowgreen'
    }

    # Define data point style
    # Linewidth of continuous alignment
    linewidth = 2.0
    # Marker of single kmer alignment
    # '+': Share at least one copy in the same chromosome
    # 'o': Do not have any shared copy
    # 'x': Mutiple Aignment
    group_marker = { 0: '^',
                     1: 'o',
                     2: '^',
                     3: 'o',
                     4: '^',
                     5: 'o',
                     6: '.',
                     7: '.',
                     8: '.',
                     9: '.'
    }

    # Label for each group
    group_label = { 0: "Chr-Specific in Strain X and Y and located in the same chromosome",
                    1: "Chr-Specific in Strain X and Y but located in different chromosomes",
                    2: "Chr-Specific in strain X, Multiple in Y, share at least one copy in the same chromosome",
                    3: "Chr-Specific in Strain X, Multiple in Y, but Y don't align to same chromosome in X",
                    4: "Chr-Specific in strain Y, Multiple in X, share at least one copy in the same chromosome",
                    5: "Chr-Specific in Strain Y, Multiple in X, but X don't align to same chromosome in Y",
                    6: "Multiple in X and Y, the located chromosomes are the same",
                    7: "Multiple in X and Y, extra chromosome alignment(s) in X",
                    8: "Multiple in X and Y, extra chromosome alignment(s) in Y",
                    9: "Mulitple in X and Y, extra chromosome alignemnt(s) in both"
    }




    # Parameter assignment
    # Distance between ticks
    tick_dist = 5000

    def extract_xy(xy_tuple_ls):
        # Extract (x, y) for plotting
        # pos_ls = [(x,y)]
        if len(pos_ls) == 1:
            x = pos_ls[0][0]
            y = pos_ls[0][1]
        elif len(pos_ls) == 0:
            return [], []
        else:
            x, y, z = zip(*pos_ls)
        return x, y

    def separate_continuous(xy_tuple_ls):
        # Separate the x-y data points
        # The data points are pre sorted
        # (x+/-1, y+/-1) is considered continous with (x,y)
        # Minimum Size of continous set is set in min_continous_set
        # continuous xy ls = [(x1,y1,x2,y2)]
        # return continous_ls, separate_ls
        first_xy = ()
        curr_xy = ()
        continuous_xy_ls = []
        separate_xy_ls = []
        continuous_start = 0
        sep_start = 0
        p = -1
        len_coninous = -1
        for t in xy_tuple_ls:
            p+=1
            if first_xy == ():
                first_xy = t
                curr_xy = t
                len_coninous += 1
                continuous_start = p
                continue
            if t[0] ==curr_xy[0] + 1 and t[1] == curr_xy[1]+1 :
                curr_xy = t
                len_coninous +=1
                continue
            else:
                if len_coninous >= min_continous_set:
                    continuous_xy_ls.append((first_xy[0], first_xy[1], curr_xy[0], curr_xy[1]))
                    separate_xy_ls.extend(xy_tuple_ls[sep_start:continuous_start])
                    sep_start = p+1
                    curr_xy = ()
                    first_xy = ()
                    len_coninous = 0
                    continue
                else:
                    first_xy = ()
                    curr_xy = ()
        if first_xy == () or len_coninous < min_continous_set:
            separate_xy_ls.extend(xy_tuple_ls[sep_start:])
        else:
            continuous_xy_ls.append((first_xy[0], first_xy[1], curr_xy[0], curr_xy[1]))

        return continuous_xy_ls, separate_xy_ls



    def draw_separate_xy(ax, group, separate_xy_ls):
        if not separate_xy_ls:
            return
        x, y = extract_xy(separate_xy_ls)
        z = tuple([group] * len(x))
        ax.scatter(x, y,z, group_marker[group], markerfacecolor='None', markeredgecolor=group_color[group])

    def draw_continous_xy(ax, group, continuous_xy_ls):
        if not continuous_xy_ls:
            return
        for t in continuous_xy_ls:
            x1 = t[0]
            y1 = t[1]
            x2 = t[2]
            y2 = t[3]
            color = group_color[group]
            ax.plot([x1,x2], [y1,y2], [group, group], color=color, linestyle="-" , linewidth=linewidth)

    def set_plot_legend(ax, is_sep_ls, is_con_ls):
        # Set legend for the figure
        # is_sep_ls = [bool], True if group has separate data points
        # is_con_ls = [bool], True if group has continuous data points
        handle_ls = []
        for group, is_sep in enumerate(is_sep_ls):
            if is_sep:
                sep_handle = mlines.Line2D([], [], markeredgecolor=group_color[group],\
                                           marker = group_marker[group],\
                                           linestyle="",\
                                           markerfacecolor='None',\
                                           label = (group_label[group] + "-Separate Points")) 
                handle_ls.append(sep_handle)
        for group, is_con in enumerate(is_con_ls):
            if is_con:
                con_handle = mlines.Line2D([], [], color=group_color[group],\
                                           label = (group_label[group] + "-Coninuous Region"))
                handle_ls.append(con_handle)

        lgd = ax.legend(handles=handle_ls, loc=9, bbox_to_anchor=(0.5, -0.1), ncol=1)
        return lgd

    # Sort the position list
    # The list is sorted by group, x, y in accending order
    if not isSorted:
        pos_ls.sort(key=lambda x: (x[2], x[0], x[1]))

    if not os.path.isdir(dirn):
        print "Cannot Access Directory\n"
        print dirn
        return
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)

    dirn_op = dirn + prefix +"/"

    xy_tuple_ls = [[] for i in range(10)]

    # if the position list is empty
    if not pos_ls:
        return ""

    # Assign the data points into 9 categories

    for t in pos_ls:
        xy_tuple_ls[t[2]].append((t[0],t[1]))

    is_sep_ls = [None for i in range(10)]
    is_con_ls = [None for i in range(10)]
    # is_sep_ls = [bool], True if group has separate data points
    # is_con_ls = [bool], True if group has continuous data points    
    sep_ls = [[] for i in range(10)]
    con_ls = [[] for i in range(10)]
    # sep_ls, categorized list of separate data points
    # con_ls, categorized list of contigous data points

    # Generate separate data point list and continuous region list for each group
    for group, tuple_ls in enumerate(xy_tuple_ls):
        con_ls[group], sep_ls[group] = separate_continuous(tuple_ls)
        if con_ls[group]:
            is_con_ls[group] = True
        else:
            is_con_ls[group] = False
        if sep_ls[group]:
            is_sep_ls[group] = True
        else:
            is_sep_ls[group] = False


    # Start Plotting
    f = plt.figure()
    ax = plt.add_subplot(111, projection='3d')
    groups = [i for i in range(10)]
    y_ticks = np.arange(0, len_b, tick_dist)
    x_ticks = np.arange(0, len_a, tick_dist)
    ax.set_yticks(y_ticks)
    ax.set_xticks(x_ticks)
    ax.axis([0, len_a+1000, 0, len_b+1000])
    ax.set_z_label("Alignment Group")

    lgd = set_plot_legend(ax, is_sep_ls, is_con_ls)
    art = []
    art.append(lgd)
    plt.xlabel(name_a)
    plt.ylabel(name_b)
    plt.grid(True)
    plt.title("%s-%s" %(name_a, name_b))

    # Draw the figure by group

    for group, tuple_ls in enumerate(con_ls):
        draw_continous_xy(ax, group, tuple_ls)

    for group, tuple_ls in enumerate(sep_ls):
        draw_separate_xy(ax, group, tuple_ls)

    # Save the figure
    fn = dirn_op+"%s-%s.png" %(name_a, name_b)
    plt.savefig(fn, dpi=600, format='png', bbox_inches='tight', additional_artist=art)
    plt.close()
    return fn

def bg2_cbs_locate_kmer_v2(fn_kmer="", fn_bg2="", fn_cbs="", dirn=""):

    chr_pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    kmer_pattern = r"(\w+)\s+Chr(\w)_(\w)\s+Chr(\w)_(\w)\s+\d+\s+\d+\s+(\w*)\n"

    bg2tocbs_kmer_ls = [[] for i in range(676)] 


    bg2tocbs_pos_ut_ls = [[] for i in range(676)] #Specific and Translocated kmer
    bg2tocbs_pos_us_ls = [[] for i in range(26)]  #Specific and Self-alignemnt kmer
    bg2tocbs_pos_mt_ls = [[] for i in range(676)] #Multi-fosmid alignment and translocated kmer
    bg2tocbs_pos_ms_ls = [[] for i in range(26)]  #Multi-fosmid alignment and Self-alignemnt kmer



    cbs = Fasta(fn_cbs)
    bg2 = Fasta(fn_bg2)

    bg2_seq_num_ls = []
    cbs_seq_num_ls = []

    
    bg2_len_ls = []
    cbs_len_ls = []


    def chr2id(chr,end):
        _id = (ord(chr) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname

    def extract_xy(pos_ls):
        # pos_ls = [(x,y)]
        if len(pos_ls) == 1:
            x = pos_ls[0][0]
            y = pos_ls[0][1]
        elif len(pos_ls) == 0:
            return [], []
        else:
            x, y = zip(*pos_ls)
        return x, y

    with open(fn_kmer, "r") as fkmer:
        for line in fkmer:
            _match = re.match(kmer_pattern, line)
            if _match:
                kmer = _match.group(1)

                bg2_chr = _match.group(2)
                bg2_end = _match.group(3)
                cbs_chr = _match.group(4)
                cbs_end = _match.group(5)
                code = _match.group(6)

                bg2_id = chr2id(bg2_chr, bg2_end)
                cbs_id = chr2id(cbs_chr, cbs_end)

                bg2tocbs_kmer_ls[cbs_id*26+bg2_id].append([kmer,code])
    for _name in bg2.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        bg2_seq_num_ls.append(_id)


    for _name in cbs.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        cbs_seq_num_ls.append(_id)


    for _seq in cbs.Seqs:
        _len = len(_seq)
        cbs_len_ls.append(_len)

    for _seq in bg2.Seqs:
        _len = len(_seq)
        bg2_len_ls.append(_len) 


    bg2_chr_align_ls = [[] for i in range(26)]

    for m in range(len(bg2.Names)):
        i = bg2_seq_num_ls[m]
        chr_align_us = []
        chr_align_ut = []
        chr_align_ms = []
        chr_align_mt = []
        for n in range(len(cbs.Names)):
            j = cbs_seq_num_ls[n]
            bg2tocbs_kmer = bg2tocbs_kmer_ls[j*26+i]
            print "%d kmers found in BG2-%s and CBS-%s" %(len(bg2tocbs_kmer), id2chr(i), id2chr(j))
            if i==j:
                for kmer in bg2tocbs_kmer:
                    pos_bg2 = find_all_seq(seq=kmer[0], ref=bg2.Seqs[m])
                    pos_cbs = find_all_seq(seq=kmer[0], ref=cbs.Seqs[n])
                    for pb in pos_bg2:
                        for pc in pos_cbs:
                            if 'u' in kmer[1]:
                                bg2tocbs_pos_us_ls[i].append((pb,pc))
                                chr_align_us.append((pb,j))
                            else:
                                bg2tocbs_pos_ms_ls[i].append((pb,pc))
                                chr_align_ms.append((pb,j))
                print "BG2-%s-CBS-%s is treated!" %(id2chr(i), id2chr(j))
                if len(bg2tocbs_pos_us_ls[i]) + len(bg2tocbs_pos_ms_ls[i]) == 0:
                    continue
                with open(dirn+"BG2-%s-CBS-%s.txt" %(id2chr(i), id2chr(j)), "w") as op:
                    op.write("BG2\t\tCBS\t\tAlign\n")
                    for pos in bg2tocbs_pos_us_ls[i]:
                        x = pos[0]
                        y = pos[1]
                        op.write("%d\t\t%d\t\tSpecific\n" %(x, y))
                    for pos in bg2tocbs_pos_ms_ls[i]:
                        x = pos[0]
                        y = pos[1]
                        op.write("%d\t\t%d\t\tMulti\n" %(x, y))
                ux, uy = extract_xy(bg2tocbs_pos_us_ls[i])
                mx, my = extract_xy(bg2tocbs_pos_ms_ls[i])

                plt.plot(ux,uy, '^', label='Specific', markerfacecolor='None', markeredgecolor="blue")
                plt.plot(mx, my, 'p', label="Multiple", markerfacecolor='None', markeredgecolor="green")
                plt.xlabel('BG2-%s' %id2chr(i))
                plt.ylabel('CBS-%s' %id2chr(j))
                plt.title('Self-alignment of specific 20mer of BG2 fosmid to CBS fosmid')
                plt.grid(True)
                lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
                art = []
                art.append(lgd)
                plt.axis([0,bg2_len_ls[m]+1000,0,cbs_len_ls[n]+1000])
                plt.savefig(dirn+"BG2-%s-CBS-%s.png" %(id2chr(i), id2chr(j)), dpi=600, format='png', bbox_inches='tight', additional_artist=art)
                plt.close()                
                continue

def draw_strain_to_chr_comp_graph(pos_ls, name_ch, name_st, len_ch, len_max_st, 
                                  dirn="", prefix='chr-st', 
                                  min_continous_set=2, isSorted=False):
    # Draw Chromosome-to-Chromosome Comparison
    
    # Format of pos_ls:
    # pos_ls = list of [(x,y,group)]
    # x, y are coordinates of a specific kmer
    # Group Label:
    # 0 = Chr-Specific in both strains and located in the same chromosome
    # 1 = Chr-Specific in both strains but located in different chromosomes
    # 2 = Chr-Specific in strain of chr-A, multiple in that of chr-B, still have at least one copy in the same chromosome
    # 3 = Chr-Specific in strain of chr-A, multiple in that of chr-B, but not located in the same chromosome
    # 4 = Chr-Specific in strain of chr-B, multiple in that of chr-A, still have at least one copy in the same chromosome
    # 5 = Chr-Specific in strain of chr-B, multiple in that of chr-A, but not located in the same chromosome
    # 6 = Multiple in both A and B, the located chromosomes are the same
    # 7 = Multiple in both A and B, extra alignment(s) in A
    # 8 = Multiple in both A and B, extra alignment(s) in B
    # 9 = Multiple in both A and B, extra alignment(s) in both

    # Drawing grous:
    # 0 = specific in both strains
    # 1 = specific in strain of the chromosome, multiple in the compared strain
    # 2 = Multiple in the strain of the chromosome, specfific in the compaired strain
    # 3 = Multiple in both strains

    # Define color of each group

    group_color = { 0: 'red',
                    1: 'blue',
                    2: 'green',
                    3: 'cyan'

    }

    # Define data point style
    # Linewidth of continuous alignment
    linewidth = 2.0
    # Marker of single kmer alignment
    # '+': Share at least one copy in the same chromosome
    # 'o': Do not have any shared copy
    # 'x': Mutiple Aignment
    group_marker = { 0: 'o',
                     1: 's',
                     2: 's',
                     3: '+'

    }

    # Label for each group
    group_label = { 0: "Chr-Specific in Strain X and Y",
                    1: "Chr-Specific in Strain X, Multi-chromosome alignment in Y",
                    2: "Multi-chromosome alignment in Strain X, chr-specific in Y",
                    3: "Multi-chromosome alignment in X and Y",

    }


    group_size = 4


    # Parameter assignment
    # Distance between ticks
    # tick_dist = 5000

    def extract_xy(xy_tuple_ls):
        # Extract (x, y) for plotting
        # pos_ls = [(x,y)]
        if not xy_tuple_ls:
            return [], []
        else:
            x, y = zip(*xy_tuple_ls)
        return x, y

    def separate_continuous(xy_tuple_ls):
        # Separate the x-y data points
        # The data points are pre sorted
        # (x+/-1, y+/-1) is considered continous with (x,y)
        # Minimum Size of continous set is set in min_continous_set
        # continuous xy ls = [(x1,y1,x2,y2)]
        # return continous_ls, separate_ls
        first_xy = ()
        curr_xy = ()
        continuous_xy_ls = []
        separate_xy_ls = []
        continuous_start = 0
        sep_start = 0
        p = -1
        len_coninous = -1
        for t in xy_tuple_ls:
            p += 1
            if first_xy == ():
                first_xy = t
                curr_xy = t
                len_coninous += 1
                continuous_start = p
                continue
            if t[0] == curr_xy[0] + 1 and t[1] == curr_xy[1]+1:
                curr_xy = t
                len_coninous += 1
                continue
            else:
                if len_coninous >= min_continous_set:
                    continuous_xy_ls.append((first_xy[0], first_xy[1], curr_xy[0], curr_xy[1]))
                    separate_xy_ls.extend(xy_tuple_ls[sep_start:continuous_start])
                    sep_start = p+1
                    curr_xy = ()
                    first_xy = ()
                    len_coninous = 0
                    continue
                else:
                    first_xy = ()
                    curr_xy = ()
        if first_xy == () or len_coninous < min_continous_set:
            separate_xy_ls.extend(xy_tuple_ls[sep_start:])
        else:
            continuous_xy_ls.append((first_xy[0], first_xy[1], curr_xy[0], curr_xy[1]))

        return continuous_xy_ls, separate_xy_ls

    def draw_separate_xy(ax, group, separate_xy_ls, y_id):
        if not separate_xy_ls:
            return
        x, y = extract_xy(separate_xy_ls)
        z = tuple([y_id * group_size + group + 0.5 for i in range(len(x))])
        ax.plot(x, y, z, group_marker[group], markerfacecolor='None', 
        	    markeredgecolor=group_color[group], linestyle="")

    def draw_continous_xy(ax, group, continuous_xy_ls, y_id):
        if not continuous_xy_ls:
            return
        for t in continuous_xy_ls:
            x1 = t[0]
            y1 = t[1]
            x2 = t[2]
            y2 = t[3]
            color = group_color[group]
            z = y_id * group_size
            ax.plot([x1, x2], [y1, y2], [z, z], color=color, linestyle="-", linewidth=linewidth)

    # def set_plot_legend(ax, is_sep_ls, is_con_ls):
    #     # Set legend for the figure
    #     # is_sep_ls = [bool], True if group has separate data points
    #     # is_con_ls = [bool], True if group has continuous data points
    #     handle_ls = []
    #     for group, is_sep in enumerate(is_sep_ls):
    #         if is_sep:
    #             sep_handle = mlines.Line2D([], [], markeredgecolor=group_color[group],\
    #                                        marker = group_marker[group],\
    #                                        linestyle="",\
    #                                        markerfacecolor='None',\
    #                                        label = (group_label[group] + "-Separate Points")) 
    #             handle_ls.append(sep_handle)
    #     for group, is_con in enumerate(is_con_ls):
    #         if is_con:
    #             con_handle = mlines.Line2D([], [], color=group_color[group],\
    #                                        label = (group_label[group] + "-Coninuous Region"))
    #             handle_ls.append(con_handle)
    # 
    #     lgd = ax.legend(handles=handle_ls, loc=9, bbox_to_anchor=(0.5, -0.1), ncol=1)
    #     return lgd
    
    # if the position list is empty
    if not pos_ls:
        return ""
    # Sort the position list
    # The list is sorted by group, x, y in accending order
    if not isSorted:
        for b_ls in pos_ls:
            b_ls.sort(key=lambda x: (x[2], x[0], x[1]))

    if not os.path.isdir(dirn):
        print "Cannot Access Directory\n"
        print dirn
        return
    if not os.path.isdir(dirn+prefix):
        os.mkdir(dirn+prefix)

    dirn_op = dirn + prefix +"/"

    xy_tuple_ls = [[[] for i in range(group_size)] for j in range(26)]

    # Assign the data points into 4 categories

    for b_id, b_ls in enumerate(pos_ls):
        for t in b_ls:
            if t[2] < 6:
                group = t[2]/2
            else:
                group = 3
            xy_tuple_ls[b_id][group].append((t[1], t[2]))

    is_sep_ls = [[None for i in range(group_size)] for i in range(26)]
    is_con_ls = [[None for i in range(group_size)] for i in range(26)]
    # is_sep_ls = [bool], True if group has separate data points
    # is_con_ls = [bool], True if group has continuous data points    
    sep_ls = [[[] for i in range(group_size)] for i in range(26)]
    con_ls = [[[] for i in range(group_size)] for i in range(26)]
    # sep_ls, categorized list of separate data points
    # con_ls, categorized list of contigous data points

    # Generate separate data point list and continuous region list for each group
    for b_id, b_tuple_ls in enumerate(xy_tuple_ls):
        for group, tuple_ls in enumerate(b_tuple_ls): 
            con_ls[b_id][group], sep_ls[b_id][group] = separate_continuous(tuple_ls)
        is_con_ls[b_id][group] = bool(con_ls[b_id][group])
        is_sep_ls[b_id][group] = bool(sep_ls[b_id][group])
    # Save point locations
    fn_c = dirn_op+"%s-%s-continuous-sorted.txt.gz" % (name_ch, name_st)
    fn_s = dirn_op+"%s-%s-separate-sorted.txt.gz" % (name_ch, name_st)
    with gzip.open(fn_c, "wb") as fcp, gzip.open(fn_s, "wb") as fsp:
        fcp.write("Y-Chr-id\tX1\tY1\tX2\tY2\tGroup\n")
        fsp.write("Y-Chr-id\tX\tY\tGroup\n")
        for b_id, tuple_ls in enumerate(con_ls):
            for group, t_ls in enumerate(tuple_ls):
            	for t in t_ls:
                    fcp.write("%d\t%d\t%d\t%d\t%d\t%d\n" 
    				                      %(b_id, t[0], t[1], t[2], t[3], group))
        for b_id, tuple_ls in enumerate(sep_ls):
            for group, t_ls in enumerate(tuple_ls):
                for t in t_ls:
                    fsp.write("%d\t%d\t%d\t%d\n" 
                              %(b_id, t[0], t[1], group))
    print "Coordinates of %s-%s is written" %(name_ch, name_st)


    # Start Plotting
    f = plt.figure()
    ax = f.add_subplot(111, projection='3d')
    # y_ticks = np.arange(0, len_b, tick_dist)
    # x_ticks = np.arange(0, len_a, tick_dist)
    # ax.set_yticks(y_ticks)
    # ax.set_xticks(x_ticks)

    ax.axis([0, len_ch+1000, 0, len_max_st+1000])
    ax.set_zlabel("Alignment on chromosomes")

    # lgd = set_plot_legend(ax, is_sep_ls, is_con_ls)
    # art = []
    # art.append(lgd)
    plt.xlabel(name_ch)
    plt.ylabel(name_st)
    plt.grid(True)
    plt.title("%s-%s" %(name_ch, name_st))

    # Draw the figure by group

    for b_id in range(26):
        for group, tuple_ls in enumerate(con_ls[b_id]):
            draw_continous_xy(ax, group, tuple_ls, b_id)
        print "%s-%s-%d Continuous Done" % (name_ch, name_st, b_id)
        for group, tuple_ls in enumerate(sep_ls[b_id]):
            draw_separate_xy(ax, group, tuple_ls, b_id)
        print "%s-%s-%d Separate Done" % (name_ch, name_st, b_id)
    # Save the figure
    fn = dirn_op+"%s-%s.png" %(name_ch, name_st)
    # plt.savefig(fn, dpi=600, format='png', bbox_inches='tight', additional_artist=art)
    plt.savefig(fn, dpi=600, format='png', bbox_inches='tight')
    plt.close()
    # return fn


def locate_kmer(fn_kmer="", fn_a="", fn_b="", dirn=""):
    # Locate kmer alignment
    # kmer-storage file format
    # kmer Chr-A Chr-B  A-Count B-Count Group
    # Group Label:
    # 0 = Chr-Specific in both strains and located in the same chromosome
    # 1 = Chr-Specific in both strains but located in different chromosomes
    # 2 = Chr-Specific in strain of chr-A, multiple in that of chr-B, still have at least one copy in the same chromosome
    # 3 = Chr-Specific in strain of chr-A, multiple in that of chr-B, but not located in the same chromosome
    # 4 = Chr-Specific in strain of chr-B, multiple in that of chr-A, still have at least one copy in the same chromosome
    # 5 = Chr-Specific in strain of chr-B, multiple in that of chr-A, but not located in the same chromosome
    # 6 = Multiple in both A and B, the located chromosomes are the same
    # 7 = Multiple in both A and B, extra alignment(s) in A
    # 8 = Multiple in both A and B, extra alignment(s) in B
    # 9 = Multiple in both A and B, extra alignemnt(s) in both

    chr_pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"
    kmer_pattern = r"(\w+)\s+Chr(\w)_(\w)\s+Chr(\w)_(\w)\s+\d+\s+\d+\s+(\d)\n"

    atob_kmer_ls = [[] for i in range(676)] # 
    atob_tuple_ls = [[] for i in range(26)] # Tuple (x, y, group)


    b = Fasta(fn_b)
    a = Fasta(fn_a)

    a_seq_num_ls = []
    b_seq_num_ls = []

    
    a_len_ls = []
    b_len_ls = []

    f_pos = gzip.open(dirn + "bg2-cbs-postions.txt.gz", "wb")
    # list of file names of chr-chr comparison
    chr_pic_fn_ls = ["" for i in range(676)]


    class kmer_align_Thread (threading.Thread):
        def __init__(self, threadID, name, counter):
            threading.Thread.__init__(self)
            self.threadID = threadID
            self.name = name
            self.counter = counter

        def run(self):
            start_time = time.time()
            print "Starting " + self.name
            kmer_align(counter=self.counter)
            print "Exiting " + self.name
            print "Taking %s sec" % (time.time() - start_time)

    def kmer_align(name_a, name_b, atob_kmer, seq_a, seq_b, fp, counter):
        pass
 


    def chr2id(chr,end):
        _id = (ord(chr) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname

    # Read kmer file
    s_time = time.time()
    with gzip.open(fn_kmer, "rb") as fkmer:
        for line in fkmer:
            _match = re.match(kmer_pattern, line)
            if _match:
                kmer = _match.group(1)

                a_chr = _match.group(2)
                a_end = _match.group(3)
                b_chr = _match.group(4)
                b_end = _match.group(5)
                group = int(_match.group(6))

                a_id = chr2id(a_chr, a_end)
                b_id = chr2id(b_chr, b_end)

                atob_kmer_ls[b_id*26+a_id].append([kmer,group])

    print "kmer file loaded! using %s seconds" % (time.time() - s_time)

    for _name in a.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        a_seq_num_ls.append(_id)

    for _name in b.Names:
        _match = re.match(chr_pattern, _name)
        _chr = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(chr=_chr,end=_end)
        b_seq_num_ls.append(_id)


    for _seq in b.Seqs:
        _len = len(_seq)
        b_len_ls.append(_len)

    for _seq in a.Seqs:
        _len = len(_seq)
        a_len_ls.append(_len) 


    # Find postion of the kmers
    # Header of position file
    f_pos.write("BG2-Chr  CBS-Chr  BG2-pos  CBS-pos  group\n")

    for m in range(len(a.Names)):
        i = a_seq_num_ls[m]
        for n in range(len(b.Names)):
            s_time = time.time()
            j = b_seq_num_ls[n]
            atob_kmer = atob_kmer_ls[j*26+i]
            print "%d kmers found in BG2-%s and CBS-%s" %(len(atob_kmer), id2chr(i), id2chr(j))
            for kmer in atob_kmer:
                pos_a = find_all_seq(seq=kmer[0], ref=a.Seqs[m])
                pos_b = find_all_seq(seq=kmer[0], ref=b.Seqs[n])
                for pb in pos_a:
                    for pc in pos_b:
                        atob_tuple_ls[j].append((pb, pc, kmer[1]))
                        f_pos.write("%s  %s  %d  %d  %d\n" %(id2chr(i), id2chr(j), pb, pc, kmer[1]))
            print "BG2-%s-CBS-%s is treated! Taking %s sec" %(id2chr(i), id2chr(j), (time.time()-s_time))

            # Draw chr-to-chr comparison
            # fn = draw_chr_to_chr_comp_graph(atob_tuple_ls[j],\
            #                                   name_a="BG2-%s" %id2chr(i),\
            #                                    name_b="CBS-%s" %id2chr(j),\
            #                                    len_a=a_len_ls[i],\
            #                                   len_b=b_len_ls[j], \
            #                                    dirn=dirn,\
            #                                    prefix='chr',\
            #                                    min_continous_set=2,\
            #                                    isSorted=False)
            # chr_pic_fn_ls[j*26+i] = fn

    f_pos.close()
    # List of pdf with page number
    # page_ls = []
    # pdf = FPDF()
    # pdf.set_font('Arial', 'B', 15)
    # pdf.set_title("20mer alignment of BG2 to CBS")
    # for k, fn in enumerate(chr_pic_fn_ls):
    #     if not fn:
    #         continue
    #     j = k % 26
    #     i = (k-i) / 26
    #     a_chr = id2chr(i)
    #     b_chr = id2chr(j)
    #     pdf.addpage()
    #     pdf.cell(20,10, "BG@-%s-CBS-%s" % (a_chr, b_chr), 1, 1, 'C')
    #     pdf.image(fn, type='PNG')


    # pdf.output(dirn+"bg2-cbs.pdf", 'F')

def load_and_draw(kmer_fn="", fn_a="", fn_b="", dirn="", prefix="output"):

    pattern = r"Chr(\w)_(\w)\s+Chr(\w)_(\w)\s+(\d+)\s+(\d+)\s+(\d+)"
    def chr2id(ch,end):
        _id = (ord(ch) - 65) * 2
        if end == "R":
            _id += 1
        return _id

    def id2chr(num):
        _end = num % 2
        _chr = num/2
        _chrname = "Chr"+chr(_chr+65)+"_"
        if _end == 0:
            _chrname += "L"
        else:
            _chrname += "R"
        return _chrname
    fa_a = Fasta(fn_a)
    fa_b = Fasta(fn_b)
    
    a_len_ls = [0 for i in range(len(fa_a.Names))]
    b_len_ls = [0 for i in range(len(fa_b.Names))]
    
    chr_pattern = r"[A-Z0-9]+_B[A-Z0-9]+_([A-Z])_([A-Z])"

    for i, _seq in enumerate(fa_a.Seqs):
        _len = len(_seq)
        _name = fa_a.Names[i]
        _match = re.match(chr_pattern, _name)
        _ch = _match.group(1)
        _end = _match.group(2)
        a_id = chr2id(ch=_ch, end=_end)        
        a_len_ls[a_id] = _len
    a_max_len = max(a_len_ls)

    for i, _seq in enumerate(fa_b.Seqs):
        _len = len(_seq)
        _name = fa_b.Names[i]
        _match = re.match(chr_pattern, _name)
        _ch = _match.group(1)
        _end = _match.group(2)
        _id = chr2id(ch=_ch, end=_end)    
        b_len_ls[_id] = _len
    b_max_len = max(b_len_ls)

    a_tuple_ls = [[[] for i in range(26)] for j in range(26)]
    # b_tuple_ls = [[] for i in range(26)]
    with gzip.open(kmer_fn, "rb") as fp:
        for line in fp:
            _match = re.match(pattern, line)
            if not _match:
            	print line
                continue
            a_ch = _match.group(1)
            a_end = _match.group(2)
            b_ch = _match.group(3)
            b_end = _match.group(4)
            a_id = chr2id(a_ch, a_end)
            b_id = chr2id(b_ch, b_end)
            a_tuple_ls[a_id][b_id].append((int(_match.group(5)),
                                           int(_match.group(6)),
                                           int(_match.group(7))))
            # b_tuple_ls[b_id][a_id].append((int(_match.group(6)),
            #                                    _match.group(5)),
            #                                    _match.group(7))
    print "kmer postions loaded!"
    for a_id, tuple_ls in enumerate(a_tuple_ls):
        draw_strain_to_chr_comp_graph(pos_ls=tuple_ls,
                                      name_ch="BG2-" + id2chr(a_id),
                                      name_st="CBS",
                                      len_ch=a_len_ls[a_id],
                                      len_max_st=b_max_len, 
                                      dirn=dirn,
                                      prefix=prefix)
        print "Figure of BG2-%s is generated" % id2chr(a_id)





def test():
    # generate_coords(dirn="/home/zhuwei/cglabarata/comp/", ref="CBS_ref.fasta", query="CBS_canu.fasta", prefix="cbs")
    # generate_dir_seq_from_labeled_seq(fn="/home/zhuwei/cglabarata/comp/cbs/cbs_dir.fa", prefix="cbs_dir2")
    # extract_ends_2(fn = "/home/zhuwei/cglabarata/comp/cbs_ends/cbs-fosmid.fasta", length=4000, prefix = "fos-end" )
    
    # fn = "/home/zhuwei/cglabarata/comp/cbs_ends/cbsfosmids/fos-end/cbs-fos-right-4kb.fa"
    # fa = Fasta(fn)
    # n = len(fa.Names)
    # score_ls = []
    # name_ls = []
    # score = pair_wise_alignment(fa.Seqs[0], fa.Seqs[1])
    # print score

    # extract_single_fasta(prefix="cbs", fn="/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa")
    # kmer_count_dirn(dirn="/home/zhuwei/cglabarata/comp/kmer/cbs/")
    # generate_kmer_table(dirn="/home/zhuwei/cglabarata/comp/kmer/bg2/kmer/", prefix="bg2")
    # compare_kmer_bg2_cbs()
    # bg2_cbs_locate_kmer(fn_kmer="/home/zhuwei/cglabarata/comp/kmer/bg2-cbs-unique20mer-translocation.txt",\
    #                             fn_bg2="/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa",\
    #                             fn_cbs="/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa",\
    #                             dirn="/home/zhuwei/cglabarata/comp/kmer/trans-map-2/")
    # compare_kmer_bg2_cbs_v2()
    # bg2_cbs_locate_kmer_v2(fn_kmer="/home/zhuwei/cglabarata/comp/kmer/bg2-cbs-unique20mer-translocation-v2.txt",\
    #                                fn_bg2="/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa",\
    #                                fn_cbs="/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa",\
    #                                dirn="/home/zhuwei/cglabarata/comp/kmer/trans-map-3/")
    # compare_kmer(name_a="BG2", name_b="CBS",\
    #              dirn_a="/home/zhuwei/cglabarata/comp/kmer/bg2/kmer/",\
    #              dirn_b="/home/zhuwei/cglabarata/comp/kmer/cbs/kmer/",\
    #              dirn="/home/zhuwei/cglabarata/comp/kmer/",\
    #              prefix="bg2-cbs",\
    #              write_kmer_table=False)
    # locate_kmer(fn_kmer="/home/zhuwei/cglabarata/comp/kmer/bg2-cbs.txt.gz",\
    #             fn_a="/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa",\
    #             fn_b="/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa",\
    #             dirn="/home/zhuwei/cglabarata/comp/kmer/trans-map-5/")
    load_and_draw(kmer_fn="/home/zhuwei/cglabarata/comp/kmer/trans-map-5/bg2-cbs-positions.txt.gz",
                  fn_a="/home/zhuwei/cglabarata/comp/kmer/bg2-fosmid.fa",
                  fn_b="/home/zhuwei/cglabarata/comp/kmer/cbs-fosmid.fa",
                  dirn="/home/zhuwei/cglabarata/comp/kmer/trans-map-5/",
                  prefix="bg2-cbs")

if __name__ == "__main__":
    test()

