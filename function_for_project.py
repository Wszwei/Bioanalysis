import re
import time
import os
import gzip
from operator import itemgetter
import matplotlib.pyplot as plt
import numpy as np
from collections import namedtuple

def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print "%s: %s" %(time_str, mess)

def simple_fastq_load(fname):
    """Load Fastq file with Phred Score in 32-ASCII code."""
    names = []
    seqs = []
    scores = []
    if not os.path.isfile(fname):
        _message("%s is not available!" % fname)
        return [], [], []
    with open(fname) as filep:
        _cur_name = ""
        _cur_seq = []
        _cur_score = []
        for count, line in enumerate(filep):
            line = line.strip()
            if not line:
                continue
            # New entry
            if count %4 == 0:
                if line[0] != '@':
                    _message("Fastq file reading error in reading %s" %line[:20])
                    return [], [], []
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append(_cur_seq)
                    scores.append(_cur_score)
                _cur_name = line[1:]
                _cur_seq = []
                _cur_score = []
            elif count %4 == 1:
                _cur_seq = line
            # load Phred Score
            elif count %4 == 3:
                _cur_score = [ord(phred) - 33 for phred in line]
                if len(_cur_score) != len(_cur_seq):
                    _message("Length of sequence and score doesnt match for %s" %_cur_name)
                    return [], [], []
        names.append(_cur_name)
        seqs.append(_cur_seq)
        scores.append(_cur_score)
    return names, seqs, scores

def simple_fastq_write(fname, names, seqs, scores):
    """Write the simple fastq file."""
    with open(fname, 'w') as filep:
        lines = []
        for id_, name in enumerate(names):
            name = '@' + name
            seq = seqs[id_]
            score = [chr(sco_ + 33)  for sco_ in  scores[id_]]
            score = "".join(score)
            lines.extend([name, seq, '+', score])
        filep.write("\n".join(lines))
        filep.write("\n")


def simple_fasta_load(fname):
    """Load fasta file.
    Args:
    fname: name of the fasta file
    Returns:
    names, seqs: name list and seqlist
    """

    names = []
    seqs = []
    if not os.path.isfile(fname):
        _message("%s is not available!" %fname)
        return [], []

    with open(fname, 'r') as fasta_p:
        _cur_name = ""
        _cur_seq = []
        for line in fasta_p:
            line_ = line.strip()
            # Skip empty lines
            if not line_:
                continue
            # New entry
            if line_[0] == ">":
                # Add previous entry
                if _cur_name:
                    names.append(_cur_name)
                    seqs.append("".join(_cur_seq))
                _cur_name = line_[1:] # Omit >
                _cur_seq = []
            else:
                # Update sequence of the entry
                if not _cur_name:
                    _message("One seq without entry")
                    _message(line_)
                    return [], []
                _cur_seq.append(line_)

        # Update the final entry
        if _cur_name:
            names.append(_cur_name)
            seqs.append("".join(_cur_seq))

    return names, seqs


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

    _name_list, _seq_list = simple_fasta_load(fname)
    # Transform to lowercase
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


def simple_fasta_write(fname, names, seqs, linewidth=80):
    """Write the fasta file."""
    with open(fname, "w") as fasta_p:
        total = len(names)
        for _id in xrange(total):
            fasta_p.write(">%s\n" %names[_id])
            _seq = seqs[_id]
            # seq with linebreak
            _seq_wb = "\n".join([_seq[i:linewidth+i]
                                 for i in xrange(0, len(_seq), linewidth)])
            fasta_p.write(_seq_wb)
            fasta_p.write("\n")

def simple_gff_load(fname):
    """Load Gff file without sequence region.
    Returns:
    comment, entry: Comment region and list of entries
    """
    if not os.path.isfile(fname):
        _message("%s is not available!" %fname)
        return "", []
    with open(fname, 'r') as gff_p:
        comment = []
        entry = []
        for line in gff_p:
            line_ = line.strip()
            if not line_:
                continue
            if line_[0] == "#":
                comment.append(line_)
            else:
                entry.append(line_.split())

    comment = "\n".join(comment)
    comment = "".join([comment, "\n"])

    return comment, entry

def simple_gff3_load(fname, return_fasta=False):
    """Load gff3 files."""
    entries = []
    with open(fname, "r") as gff3:
        for line in gff3:
            line = line.strip()
            if line[0] == "#":
                if line == "##FASTA":
                    break
                continue
            entry = line.split("\t")
            if len(entry) < 9:
                _message("Error loading %s: Less than 9 items in the entry\n%s" %(fname, line))
                continue
            notes = entry.pop().split(";")
            entry.extend(notes)
            entry[3] = int(entry[3])
            entry[4] = int(entry[4])
            entries.append(entry)

    if not return_fasta:
        return entries

    else:
        names = []
        seqs = []
        with open(fname, "r") as gff3:
            line = gff3.readline()
            line = line.strip()
            while not line == "##FASTA":
                line = gff3.readline()
                line = line.strip()
            _cur_name = ""
            _cur_seq = []
            line = gff3.readline()
            while line:
                line_ = line.strip()
                line = gff3.readline()
                # Skip empty lines
                if not line_:
                    continue
                # New entry
                if line_[0] == ">":
                    # Add previous entry
                    if _cur_name:
                        names.append(_cur_name)
                        seqs.append("".join(_cur_seq))
                    _cur_name = line_[1:] # Omit >
                    _cur_seq = []
                else:
                # Update sequence of the entry
                    if not _cur_name:
                        _message("One seq without entry")
                        _message(line_)
                        return entries, [], []
                    _cur_seq.append(line_)
            # Update the final entry
            if _cur_name:
                names.append(_cur_name)
                seqs.append("".join(_cur_seq))
        return entries, names, seqs


def rc_seq(seq=""):
    """Returns the reverse compliment sequence."""
    rc_nt_ls = []
    rc_dict = {
        "a": "t",
        "c": "g",
        "t": "a",
        "g": "c",
        "n": "n",
        "A": "T",
        "C": "G",
        "T": "A",
        "G": "C",
        "N": "N"
    }
    rc_nt_ls = [rc_dict[seq[i]] for i in range(len(seq)-1, -1, -1)]
    rc_seq_ = "".join(rc_nt_ls)
    return rc_seq_

def extract_part(seq="", is_rc=False, pos=None):
    """Extract part of the sequence.
    Args:
    seq: sequence
    isrc: whether to return reverse compliment seq
    pos: [start_pos, end_pos] postion to extract
    Returns:
    extract_seq: DNA sequence
    """

    extract_seq = seq[pos[0] - 1: pos[1]] # Include nt and end_pos

    if is_rc:
        extract_seq = rc_seq(extract_seq)

    return extract_seq

def translate_exon(seq):
    """Translate exon by normal codon."""

    codon = {
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "CTT": "L",
        "CTC": "L",
        "CTG": "L",
        "CTA": "L",
        "TTA": "L",
        "TTG": "L",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "TTT": "F",
        "TTC": "F",
        "ATG": "M",
        "TGT": "C",
        "TGC": "C",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACT": "T",
        "ACC": "T",
        "ACG": "T",
        "ACA": "T",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "AGT": "S",
        "AGC": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TGG": "W",
        "CAA": "Q",
        "CAG": "Q",
        "AAT": "N",
        "AAC": "N",
        "CAT": "H",
        "CAC": "H",
        "GAA": "E",
        "GAG": "E",
        "GAT": "D",
        "GAC": "D",
        "AAA": "K",
        "AAG": "K",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGA": "R",
        "AGG": "R",
        "TAA": "*",
        "TAG": "*",
        "TGA": "*"
    }

    seq = seq.upper()
    size = len(seq)
    if size % 3 != 0:
        _message("Length is %d, invalid !" % size)
    protein = [codon[seq[id_: id_+3]] for id_ in range(0, size, 3)]
    protein = "".join(protein)
    return protein

# Aim: Filer out the features in cbs ref gff file in the fosmid region

def id2chr(num):
    """Transform 0-25 to ChrA-Z."""
    return 'Chr' + chr(num+65)

def chr2id(chrname):
    """Transform ChrA-Z to 0-25."""
    pattern = re.compile('Chr([A-Z])')
    match_ = pattern.match(chrname)
    if match_:
        return ord(match_.group(1)) - 65

    _message("%s is not in the Chr[X] format" %chrname)
    return -1


def tmp_filter_name_by_pos_cbs(names):
    """Filter names in the fasta file without fosmid region."""

    # Name pattern for cbs and mrna file
    pattern = re.compile("COORDS:Chr(\w)_C_glabrata_CBS138:(\d+)-(\d+)")
    # Region of non-fosmid cbs region
    cbs_region = [
        [20487, 459937],
        [13803, 483636],
        [30470, 540559],
        [24693, 631859],
        [35344, 650244],
        [34830, 901377],
        [15591, 971009],
        [20281, 1031248],
        [27900, 1081474],
        [22048, 1154829],
        [21762, 1280027],
        [43659, 1409236],
        [27473, 1385304]]

    filter_list = [True for i in xrange(len(names))]

    for _id, name in enumerate(names):
        _match = pattern.search(name)
        if not _match:
            filter_list[_id] = False
            continue
        _chr = _match.group(1)
        _chr_id = ord(_chr) - 65 # A-Z => 0-25
        _chr_s = cbs_region[_chr_id][0]
        _chr_e = cbs_region[_chr_id][1]

        _entry_region = [int(_match.group(2)), int(_match.group(3))]
        _entry_region.sort()
        # Entry in / overlaps the fosmid region
        if _entry_region[0] < _chr_s or _entry_region[1] > _chr_e:
            filter_list[_id] = False

    return filter_list

def filter_cbs_fos():
    """Filter out the features of the cbs cDNA and protein fasta file in the fosmid region."""

    dirn = "/home/zhuwei/cglabarata/anno/filter_fos/"
    cds_fn = dirn + "cbs-cds.fasta"
    mrna_fn = dirn + "cbs-mRNA.fasta"

    cds_nf = dirn + "CBS-cds-nofos.fa"
    cds_filtered = dirn + "CBS-cds-filterout.fa"

    mrna_nf = dirn + "CBS-mrna-nofos.fa"
    mrna_filtered = dirn + "CBS-mrna-filterout.fa"


    # Filter the cds file
    cds_name, cds_seq = simple_fasta_load(cds_fn)
    cds_filter_list = tmp_filter_name_by_pos_cbs(cds_name)

    cds_name_nf = [cds_name[_id] for _id in xrange(len(cds_name))
                   if cds_filter_list[_id]]

    cds_seq_nf = [cds_seq[_id] for _id in xrange(len(cds_seq))
                  if cds_filter_list[_id]]

    cds_name_filtered = [cds_name[_id] for _id in xrange(len(cds_name))
                         if not cds_filter_list[_id]]
    cds_seq_filtered = [cds_seq[_id] for _id in xrange(len(cds_seq))
                        if not cds_filter_list[_id]]
    simple_fasta_write(cds_nf, cds_name_nf, cds_seq_nf)
    simple_fasta_write(cds_filtered, cds_name_filtered, cds_seq_filtered)

    mrna_name, mrna_seq = simple_fasta_load(mrna_fn)
    mrna_filter_list = tmp_filter_name_by_pos_cbs(mrna_name)

    mrna_name_nf = [mrna_name[_id] for _id in xrange(len(mrna_name))
                    if mrna_filter_list[_id]]

    mrna_seq_nf = [mrna_seq[_id] for _id in xrange(len(mrna_seq))
                   if mrna_filter_list[_id]]
    mrna_name_filtered = [mrna_name[_id] for _id in xrange(len(mrna_name))
                          if not mrna_filter_list[_id]]
    mrna_seq_filtered = [mrna_seq[_id] for _id in xrange(len(mrna_seq))
                         if not mrna_filter_list[_id]]

    simple_fasta_write(mrna_nf, mrna_name_nf, mrna_seq_nf)
    simple_fasta_write(mrna_filtered, mrna_name_filtered, mrna_seq_filtered)

# Aim: Adjust direction

def adjust_direction(fname="", tag="--r", prefix="adj"):
    """Generate Reverse Compliment seq for sequence with the tag."""
    _extention = os.path.splitext(fname)[1][1:].strip()
    _path = os.path.splitext(fname)[0]
    op_name = "%s-%s.%s" %(_path, prefix, _extention)
    with open(op_name, "w") as out_p:
        names, seqs = simple_fasta_load(fname)
        lines = []
        for seq_id, seq_name in enumerate(names):
            if tag in seq_name:
                ext_p = seq_name.find(tag)
                _name = seq_name[:ext_p]
                lines.append(">%s\n" %_name)
                _rc_seq = rc_seq(seqs[seq_id])
                lines.append(_rc_seq)
                lines.append("\n")
            else:
                lines.append(">%s\n" %names[seq_id])
                lines.append(seqs[seq_id])
                lines.append("\n")
        out_p.write("".join(lines))

def adjust_direction_fastq(fname="", tag="--r", prefix="adj"):
    """Generate Reverse Compliment seq for sequence with the tag for fastq files."""
    _extention = os.path.splitext(fname)[1][1:].strip()
    _path = os.path.splitext(fname)[0]
    op_name = "%s-%s.%s" %(_path, prefix, _extention)
    with open(fname) as orif, open(op_name, "w") as adjf:
        lines = []
        count = 0
        rc_flag = False
        for line in orif:
            line = line.strip()
            if not line:
                continue
            if count %4 == 0:
                if tag in line:
                    rc_flag = True
                    tag_pos = line.find(tag)
                    line = line[:tag_pos]
                    lines.append(line)
                else:
                    lines.append(line)
                    rc_flag = False
            elif count % 4 == 1:
                if rc_flag:
                    seq = rc_seq(line)
                    lines.append(seq)
            elif count % 4 == 2:
                lines.append(line)
            elif count % 4 == 3:
                if rc_flag:
                    score = line[::-1]
                else:
                    score = line
                lines.append(score)
            count += 1
        adjf.write("\n".join(lines))
    _message("%s: Direction adjusted!" %fname)


def extract_fos():
    """Extract cbs fosmid from the canu assembly by the cbs-fos sequence."""

    extract_dict = {
    'ChrA': [1, 35134, 474584, 527797],
    'ChrB': [1, 16812, 486647, 512654],
    'ChrC-1': [1, 38194, -1, -1],
    'ChrC-2': [-1, -1, 472588, 492862],
    'ChrD': [1, 34422, 641586, 667896],
    'ChrE': [1, 50253, 667551, 707419],
    'ChrF': [1, 57905, 924449, 957511],
    'ChrG': [1, 29183, 984590, 1010279],
    'ChrH': [1, 23241, 1035470, 1057384],
    'ChrI': [1, 29585, 1115489, 1142513],
    'ChrJ': [1, 38521, 1199488, 1247560],
    'ChrK': [1, 23488, 1281747, 1307083],
    'ChrL-1': [1, 51341, 1422011, 1486314],
    'ChrM': [1, 37270, 1395090, 1420936]
    }

    fname = "/home/zhuwei/cglabarata/cbs/fos/cbs-2.fasta"
    dirn = "/home/zhuwei/cglabarata/cbs/fos/"

    fos_name = dirn + "cbs-canu-fos.fa"
    fos_names = []
    fos_seqs = []

    names, seqs = simple_fasta_load(fname)
    for id_, name in enumerate(names):
        if name in extract_dict:
            _message("Extracting fosmid in %s" %name)
            _message("Length: %d" %len(seqs[id_]))

            pos = extract_dict[name]
            _message(pos)
            chr_name = name[:4]
            left = pos[:2]
            right = pos[2:]
            if left[0] > 0:
                _seq = extract_part(seq=seqs[id_], is_rc=True, pos=left)
                fos_names.append(chr_name+"_Left")
                fos_seqs.append(_seq)
            if right[0] > 0:
                _seq = extract_part(seq=seqs[id_], is_rc=False, pos=right)
                fos_names.append(chr_name+"_Right")
                fos_seqs.append(_seq)

    simple_fasta_write(fos_name, fos_names, fos_seqs)


# Aim: Kmer linkage map

class LinkNode(object):
    """Node for linkage map."""
    def __init__(self, link_map, kmer, pos, kmer_p, kmer_n, pos_p, pos_n):
        self.link_map = link_map
        self.kmer = kmer
        self.kmer_p = kmer_p
        self.kmer_n = kmer_n
        self.pos = pos
        self.pos_p = pos_p
        self.pos_n = pos_n

    def pre_node(self):
        node_ = self.link_map[self.kmer_p]
        return node_

    def next_node(self):

        node_ = self.link_map[self.kmer_n]
        return node_

    def del_node(self):
        node_p = self.pre_node()
        node_n = self.next_node()

        node_p.kmer_n = node_n.kmer
        node_p.pos_n = node_n.pos
        node_n.kmer_p = node_p.kmer
        node_n.pos_p = node_p.pos
        del self.link_map[self.kmer]

class ChrLinkMap(object):
    """Kmer Linkage map."""
    def __init__(self, kmer=20, seq=""):
        self.link_map = {}
        self.kmer_size = kmer
        if seq:
            self.seq = seq.upper()
            self.generate_map()

    def generate_map(self):

        # Generate Head Node

        self.link_map["$"] = LinkNode(self.link_map, "$", -1, "$", "$", -1, -1)
        kmer_pre = ""
        pos_p = -1
        kmer = self.get_kmer(0)
        pos = 0
        kmer_next = self.get_kmer(1)
        pos_n = 1
        self.link_map[kmer] = LinkNode(self.link_map, kmer, pos, kmer_pre, kmer_next, pos_p, pos_n)
        self.link_map["$"] = LinkNode(self.link_map, "$", -1, "$", "$", -1, -1)
        self.link_map[""] = LinkNode(self.link_map, "", -1, "", kmer, -2, 0)
        rep_kmer = []
        cur_node = self.link_map[kmer]

        for pos_ in xrange(1, len(self.seq) - self.kmer_size + 1):
            # mask repeat
            if kmer_next in self.link_map:
                node_ = self.link_map[kmer_next]
                if node_ == cur_node:
                    cur_node = cur_node.pre_node()
                node_.del_node()
                rep_kmer.append(kmer_next)
                kmer_next = self.get_kmer(pos_ + 1)
                continue
            elif kmer_next in rep_kmer:
                kmer_next = self.get_kmer(pos_ + 1)
                continue
            cur_node.kmer_n = kmer_next
            cur_node.pos_n = pos_
            kmer_pre = cur_node.kmer
            kmer = kmer_next
            kmer_next = self.get_kmer(pos_ + 1)
            self.link_map[kmer] = LinkNode(self.link_map, kmer, pos_,
                                           cur_node.kmer, "$", cur_node.pos,
                                           len(self.seq) - self.kmer_size + 1)
            cur_node = self.link_map[kmer]
        # Update Tail Node
        self.link_map["$"] = LinkNode(self.link_map, "$", len(self.seq) - self.kmer_size + 1,
                                      kmer, "$", pos_, len(self.seq) - self.kmer_size + 1)
        rep_kmer = []

    def get_kmer(self, position):

        if position == len(self.seq) - self.kmer_size + 1:
            return "$"
        kmer_ = self.seq[position: position + self.kmer_size]
        kmer_rc = rc_seq(kmer_)
        if kmer_rc < kmer_:
            kmer_ = kmer_rc
        return kmer_

    def del_kmer(self, kmer):
        if kmer == "":
            print "Cannot Delete HEAD"
            return
        if kmer == "$":
            print "Cannot Delete Tail"
            return
        node_ = self.link_map[kmer]
        node_.del_node()

def generate_shared_kmer_map(fname_a, fname_b, pattern, kmer=20):
    # Generate kmer linkage map for two strains
    name_a, seq_a = load_fasta(fname_a, pattern)
    name_b, seq_b = load_fasta(fname_b, pattern)
    size_a = len(name_a)
    size_b = len(name_b)
    _message("Building reference kmer linkage map")
    map_a = [ChrLinkMap(seq=seq_a[id_], kmer=kmer) for id_ in range(size_a)]
    _message("Building query kmer linkage map")
    map_b = [ChrLinkMap(seq=seq_b[id_], kmer=kmer) for id_ in range(size_b)]

    kmer_a = {}
    kmer_share = {}
    _message("Generating shared kmer list")
    for map_ in map_a:
        for kmer_ in map_.link_map:
            if kmer_ not in kmer_a:
                kmer_a[kmer_] = None
    for map_ in map_b:
        for kmer_ in map_.link_map:
            if kmer_ in kmer_a:
                kmer_share[kmer_] = None
    kmer_a.clear()
    _message("Filtering shared kmer in reference")
    if "$" not in kmer_share:
        print "No tail"
        kmer_share["$"] = None
    if "" not in kmer_share:
        print "No head"
        kmer_share[""]=None
    for map_ in map_a:
        for kmer_ in map_.link_map.keys():
            if kmer_ not in kmer_share:
                map_.del_kmer(kmer_)
    _message("Filtering shared kmer in query")
    for map_ in map_b:
        for kmer_ in map_.link_map.keys():
            if kmer_ not in kmer_share:
                map_.del_kmer(kmer_)

    return map_a, map_b

class LinkageNode(object):
    """Node for the strain tp strain linkage map."""
    def __init__(self, id_, id_align, pos, pos_n, align_count):
        self.id_ = id_
        self.id_align = id_align
        self.align_count = align_count
        self.pos = pos
        self.pos_next = pos_n

def build_linkage_map(map_list_r, map_list_q):
    """ Build the kmer linkage map"""

    ref_linkage = []
    query_linkage = []

    _message("Start making the linkage map")
    # Linkage mapping
    _message("Mapping the query strain to the reference strain")
    for id_, map_q in enumerate(map_list_q):
        kmer_map_q = map_q.link_map
        node_ = kmer_map_q[""]
        linkage_map = []

        while node_.kmer != "$":
            node_ = node_.next_node()
            kmer = node_.kmer
            id_align, kmer_next, align_count = align_kmer(kmer, map_q, map_list_r)
            linkage_map.append(LinkageNode(id_, id_align,
                                           node_.pos, kmer_map_q[kmer_next].pos, align_count))
            node_ = kmer_map_q[kmer_next].next_node()
        query_linkage.append(linkage_map)
        _message("%d Finished" %id_)

    _message("Mapping the reference strain to the query strain")
    for id_, map_r in enumerate(map_list_r):
        kmer_map_r = map_r.link_map
        node_ = kmer_map_r[""]
        linkage_map = []

        while node_.kmer != "$":
            node_ = node_.next_node()
            kmer = node_.kmer
            id_align, kmer_next, align_count = align_kmer(kmer, map_r, map_list_q)
            linkage_map.append(LinkageNode(id_, id_align,
                                           node_.pos, kmer_map_r[kmer_next].pos, align_count))
            node_ = kmer_map_r[kmer_next].next_node()
        ref_linkage.append(linkage_map)
    _message("Linkage Map Finished")
    return ref_linkage, query_linkage

def align_kmer(kmer, map_kmer, map_list):
    """Align one kmer to the other strain."""
    align_count = 0
    kmer_next = ""
    pos_next = -1
    id_align = -1
    for id_, map_ in enumerate(map_list):
        link_map = map_.link_map
        if kmer in link_map:
            align_count += 1
            kmer_next_, pos_r_ = search_branch(map_kmer, map_, kmer)
            if pos_r_ > pos_next:
                pos_next = pos_r_
                kmer_next = kmer_next_
                id_align = id_
    return id_align, kmer_next, align_count

def search_branch(map_r, map_q, kmer):
    # Return the shared branch between map_a and map_b from shared kmer
    link_map_r = map_r.link_map
    link_map_q = map_q.link_map
    node_a = link_map_r[kmer]
    node_b = link_map_q[kmer]
    kmer_ = kmer
    pos_r = node_a.pos
    node_n_a = node_a.next_node()
    node_n_b = node_b.next_node()
    while node_n_a.kmer == node_n_b.kmer:
        kmer_ = node_n_a.kmer
        pos_r = node_n_a.pos
        node_n_a = node_n_a.next_node()
        node_n_b = node_n_b.next_node()
        if node_n_a.kmer == "$":
            return kmer_, pos_r
    return kmer_, pos_r


def id2fosmid(num):
    """Converts number to fosmid.
    Odd Number -> Left, Even numbers -> Right
    Chr = Number/2 -> A - Z
    """
    name = "Chr" + chr(65 + num/2)
    if num%2 == 0:
        name += "_Left"
    else:
        name += "_Right"
    return name


def draw_linkage_alignment(ref_linkage, query_linkage, ref_name, query_name,
                           ref_len_list, query_len_list, namefunction, dirn,
                           prefix, dpi=300):
    """Draw the kmer linkage alignment."""
    # Group Color:
    # 0 = Paired Single Alignment, grey
    # 1 = Unpaired Single Alignment, red
    # 2 = Double Alignment, blue
    # 3 = >3 Alignment, cyan

    group_color = ["grey", "red", "blue", "cyan"]
    ref_ticks = range(0, len(ref_len_list))
    ref_labels = ["%s_%s" %(ref_name, namefunction(id_)) for id_ in ref_ticks]
    query_ticks = range(0, len(query_len_list))
    query_labels = ["%s_%s" %(query_name, namefunction(id_)) for id_ in query_ticks]

    if not os.path.isdir(dirn + "/" + prefix):
        os.mkdir(dirn + "/" + prefix)

    _message("Start drawing the Reference linkage map")

    # Reference Linkage map

    for chr_link in ref_linkage:
        width = ref_len_list[id_]/dpi/5
        height = width * 0.618
        fig, axes = plt.subplots(figsize=(width, height), dpi=dpi)
        for linkage in chr_link:
            id_ = linkage.id_
            id_align = linkage.id_align
            pos = linkage.pos
            pos_next = linkage.pos_next
            align_count = linkage.align_count
            if align_count == 1 and id_ == id_align:
                group = 0
            elif align_count == 1 and id_ != id_align:
               group = 1
            elif align_count == 2:
                group = 2
            else:
                group = 3
            # Tmp
            if pos_next - pos < 500:
                axes.plot((pos), (id_align), color=group_color[group], marker="o")
            else:
                axes.plot((pos, pos_next), (id_align, id_align), color=group_color[group], linestyle="-", linewidth=2)
        # End of the sequence
        axes.plot((ref_len_list[id_], (ref_len_list[id_])), (0, len(ref_len_list)), color="grey", linewidth=0.5)
        # Set up y ticks
        axes.set_yticks(query_ticks)
        axes.set_yticklabels(query_labels, rotation='horizontal', ha='right')
        # Set title

        axes.set_title("Linkage Map of %s to %s" %(ref_labels[id_], query_name))
        fig_name = "%s%s/%s-%s-linkage.png" %(dirn, prefix, ref_labels[id_], query_name)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(fig_name, bbox_inches="tight", format="png")
        plt.close()


    _message("Start drawing the Query linkage map")
    for chr_link in query_linkage:
        width = query_len_list[id_]/dpi/5
        height = width * 0.618
        fig, axes = plt.subplots(figsize=(width, height), dpi=dpi)
        for linkage in chr_link:
            id_ = linkage.id_
            id_align = linkage.id_align
            pos = linkage.pos
            pos_next = linkage.pos_next
            align_count = linkage.align_count
            if align_count == 1 and id_ == id_align:
                group = 0
            elif align_count == 1 and id_ != id_align:
               group = 1
            elif align_count == 2:
                group = 2
            else:
                group = 3
            # Tmp
            if pos_next - pos < 500:
                axes.plot((pos), (id_align), color=group_color[group], marker="o")
            else:
                axes.plot((pos, pos_next), (id_align, id_align), color=group_color[group], linestyle="-", linewidth=2)
        # End of the sequence
        axes.plot((query_len_list[id_], (query_len_list[id_])), (0, len(query_len_list)), color="grey", linewidth=0.5)
        # Set up y ticks
        axes.set_yticks(ref_ticks)
        axes.set_yticklabels(ref_labels, rotation='horizontal', ha='right')
        # Set title
        axes.set_title("Linkage Map of %s to %s" %(query_labels[id_], ref_name))
        fig_name = "%s%s/%s-%s-linkage.png" %(dirn, prefix, query_labels[id_], ref_name)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(fig_name, bbox_inches="tight", format="png")
        plt.close()




def compare_strain(ref_name, query_name, dirn, ref_file, query_file, pattern, prefix="out"):
    """Compare two strains by kmer linkage map."""
    map_ref, map_query = generate_shared_kmer_map(dirn + ref_file, dirn + query_file, pattern)
    ref_linkage, query_linkage = build_linkage_map(map_ref, map_query)

    _message("Writing Linkage File")

    with open(dirn + prefix + "-" + ref_name + ".txt", "w") as ref_p:
        lines = []
        for chr_link in ref_linkage:
            for link in chr_link:
                line = "%d, %d, %d, %d, %d\n" %(link.id_, link.id_align,
                                                link.pos, link.pos_next, link.align_count)
                lines.append(line)
        ref_p.write("".join(lines))

    with open(dirn + prefix + "-" + query_name + ".txt", "w") as query_p:
        lines = []
        for chr_link in query_linkage:
            for link in chr_link:
                line = "%d, %d, %d, %d, %d\n" %(link.id_, link.id_align,
                                                link.pos, link.pos_next, link.align_count)
                lines.append(line)
        query_p.write("".join(lines))

    names, seqs = load_fasta(dirn + ref_file, pattern)
    ref_len_list = [len(seq) for seq in seqs]
    names, seqs = load_fasta(dirn + query_file, pattern)
    query_len_list = [len(seq) for seq in seqs]

    draw_linkage_alignment(ref_linkage, query_linkage, ref_name, query_name,
                           ref_len_list, query_len_list, namefunction=id2fosmid, dirn=dirn,
                           prefix=prefix)
    _message("Comparison Finished")


def one2one_nucmer():
    """Perform one 2 one nucmer comparison."""
    dirn = "/home/zhuwei/cglabarata/cbs/fos/single/"
    if not os.path.isdir(dirn):
        os.mkdir(dirn)
    os.chdir(dirn)
    names, seqs = load_fasta(dirn + "../cbs-fosmid.fa", pattern=r"_(\w)_(\w)")
    na_c, seq_c = load_fasta(dirn + "../cbs-canu-fos.fa", pattern=r"(\w)_(\w)")
    # Generate single fastas
    for id_, name_ in enumerate(names):
        fn_cbs_ = dirn + "cbs_" + id2fosmid(id_) + ".fa"
        simple_fasta_write(fn_cbs_, [name_], [seqs[id_]])
        fn_canu_ = dirn + "canu_" + id2fosmid(id_) + ".fa"
        simple_fasta_write(fn_canu_, [na_c[id_]], [seq_c[id_]])
        os.system("mummer -maxmatch -b %s %s > %s.mums" %(fn_cbs_, fn_canu_, id2fosmid(id_), ))
        os.system('''mummerplot --color -p %s --xrange "[0,%d]" --yrange "[0,%d]"  -png --large %s.mums'''
                  %(id2fosmid(id_), len(seqs[id_]), len(seq_c[id_]), id2fosmid(id_)))
        # os.system("show-coords -clrT %s.delta > %s.coords" %(id2fosmid(id_), id2fosmid(id_)))
        # os.system("show-snps -CrT %s.delta > %s.snps" %(id2fosmid(id_), id2fosmid(id_)))

def pairwise_nucmer():
    """Pairwise nucmer comparison."""
    dirn = "/home/zhuwei/cglabarata/comp/58x/pairwise/"
    dirn_fa = "/home/zhuwei/cglabarata/comp/58x/"
    os.chdir(dirn)
    con_5x = ["53-con-2.fa"]
    con_58x = ["580-con.fa", "581-con.fa", "582-con.fa", "583-con.fa", "585-con.fa"]
    for fn_5x in con_5x:
        for fn_58x in con_58x:
            name_ = fn_5x[:2] + "." + fn_58x[:3]
            os.system("nucmer -l 5000 -mum -p %s %s %s" %(name_, dirn_fa + fn_5x, dirn_fa + fn_58x))
            os.system("show-coords -L 5000 -clrT %s.delta > %s.coords" %(name_, name_))
            os.system("show-snps -CrT %s.delta > %s.snps" %(name_, name_))
            os.system("mummerplot -p %s -png --medium --filter %s.delta" %(name_, name_))

def maker_gff_extract(fname, dirn, prefix="out"):
    """Extact the information from maker gff3 file.
    The augustus_masked match extracted for augustus call genes
    The maker gene extracted for augusuted called genes with blastx overlap
    the protein2genome protein_match is extracted for blastx called protein."""

    # 0 if stop codon exists in protein lib, 3 if not
    STOP_SKIP = 3
    os.chdir(dirn)
    gff_entry = simple_gff3_load(fname)
    augustus_masked = []
    maker_gene = {}
    maker_pos2name = {}
    protein_blastx = []

    augustus_maksed_name = []
    augustus_abinitio = []
    maker_protein_region = []

    maker_blastx = {}

    for entry in gff_entry:
        if entry[1] == "maker" and entry[2] == "gene":
            maker_gene[entry[8][3:]] = entry
            augustus_maksed_name.append(entry[8][3:])
            pos_ = (entry[3], entry[4])
            maker_protein_region.append(pos_)
            maker_pos2name[pos_] = entry[8][3:]
            continue
        elif entry[1] == "augustus_masked" and entry[2] == "match":
            augustus_masked.append(entry)
            continue
        elif entry[1] == "protein2genome" and entry[2] == "protein_match":
            protein_blastx.append(entry)
    # Ab initio genes without any
    for entry in augustus_masked:
        name_ = entry[-1][5:]
        name_ = name_.replace("abinit", "processed")
        name_ = name_.replace("-mRNA-1", "")
        print name_
        if name_ not in augustus_maksed_name:
            augustus_abinitio.append(entry)

    # Check overlap of augustus call genes with blastx
    for entry in protein_blastx:
        if STOP_SKIP == 0:
            pos_ = (entry[3], entry[4])
        else:
            if entry[6] == "+":
                pos_ = (entry[3], entry[4] + STOP_SKIP)
            else:
                pos_ = (entry[3] - STOP_SKIP, entry[4])
        if pos_ in maker_protein_region:
            maker_blastx.setdefault(pos_, [])
            maker_blastx[pos_].append(entry[-1][5:])
    # Write exact overlapping genes
    with open(dirn + prefix + "-maker-exact.txt", "w") as overp:
        lines = []
        for pos_, name_ in maker_blastx.iteritems():
            m_name = maker_pos2name[pos_]
            m_entry = maker_gene[m_name]
            del maker_gene[m_name]
            if m_entry[6] == "+":
                sense_ = "sense"
            else:
                sense_ = "antisense"
            fosmid = m_entry[0]
            line = "%s\t%d\t%d\t%s\t%s\n" %(name_, pos_[0], pos_[1], sense_, fosmid)
            lines.append(line)
        overp.write("".join(lines))
    # Write Maker genes with out perfect overlapping
    with open(dirn + prefix + "-maker-non-extract.txt", "w") as overp:
        lines = []
        for name_, entry_ in maker_gene.iteritems():
            if entry_[6] == "+":
                sense_ = "sense"
            else:
                sense_ = "antisense"
            fosmid = entry_[0]
            line = "%s\t%d\t%d\t%s\t%s\n" %(name_, entry_[3], entry_[4], sense_, fosmid)
            lines.append(line)
        overp.write("".join(lines))
    # Write Ab initio genes
    with open(dirn + prefix + "-augustus-abinitio.txt", "w") as overp:
        lines = []
        for entry_ in augustus_abinitio:
            fosmid = entry_[0]
            if entry_[6] == "+":
                sense_ = "sense"
            else:
                sense_ = "antisense"
            name_ = entry_[-1][5:]
            line = "%s\t%d\t%d\t%s\t%s\n" %(name_, entry_[3], entry_[4], sense_, fosmid)
            lines.append(line)
        overp.write("".join(lines))

def compare2fos():
    """Compare the canu assemblies to cbs_I_Left fosmid."""
    dirn = "/home/zhuwei/cglabarata/comp/58x/MR/"
    if not os.path.isdir(dirn):
        os.mkdir(dirn)
    os.chdir(dirn)
    for fname in os.listdir(dirn):
        if not fname.endswith("fa"):
            continue
        if fname == "cbs_MR.fa":
            continue
        os.system("nucmer -mum -p %s cbs_MR.fa %s" %(fname[:3], fname))
        os.system("show-coords -clrT -L 5000 %s.delta > %s.coords" %(fname[:3], fname[:3]))

def pairwise_dnadiff():
    """Pairwise Comparison by nucmer."""
    dirn = "/home/zhuwei/cglabarata/comp/chrcomp/170913/58x_self/"
    os.chdir(dirn)
    fn_58 = ["580-con.fa", "582-con.fa", "583-con.fa", "585-con.fa"]
    # fn_5 = ["51-con.fa", "52-con.fa", "53-con.fa", "54-con.fa"]
    for x in range(3):
        fn_r = fn_58[x]
        for y in range(x + 1, 4):
            fn_q = fn_58[y]
            prefix = fn_r[:3] + "." + fn_q[:3]
            os.system("dnadiff -p %s %s %s" %(prefix, fn_r, fn_q))

def extract_gpi():
    """Extract gpi proteins."""
    dirn = "/home/zhuwei/cglabarata/anno/canu-fosmid/gpi/"
    os.chdir(dirn)
    gpi_nf = "gpi-nf.txt"
    gpi_fos = "gpi-fos.txt"
    pr_nf = "cbs-cds-nofos.fa"
    pr_fos = "cbs.fos.protein.fasta"
    out = "gpi.fa"
    gpi_names = []
    gpi_seqs = []
    names_nf, seqs_nf = simple_fasta_load(pr_nf)
    names_fos, seqs_fos = simple_fasta_load(pr_fos)
    names_nf = [item.split()[0] for item in names_nf]

    with open(gpi_nf, "r") as nfp:
        for line in nfp:
            name = line.strip()
            if name in names_nf:
                id_ = names_nf.index(name)
                print "None-fos: %s" %name
                gpi_names.append(name)
                gpi_seqs.append(seqs_nf[id_])
    with open(gpi_fos, "r") as fosp:
        for line in fosp:
            name = line.strip()
            id_ = names_fos.index(name)
            gpi_names.append(name)
            gpi_seqs.append(seqs_fos[id_])

    simple_fasta_write(out, gpi_names, gpi_seqs)

def snp_cluster(dirn, fname, distance=10000, prefix="out"):
    """cluster the snps for mummer output"""
    os.chdir(dirn)
    cluster = {}
    snps = []
    with open(fname, "r") as snpp:
        for line in snpp:
            line = line.lstrip()
            if not line:
                continue
            entry = line.split()
            snps.append((entry[-2], entry[-1], int(entry[0]), int(entry[3])))
    snps.sort()
    cur_pair = None
    cur_left_r = -1
    cur_right_r = -1
    cur_left_q = -1
    cur_right_q = -1
    # load snps
    for snp in snps:
        pair = (snp[0], snp[1])
        if cur_pair is None:
            cur_pair = pair
            cur_left_r = snp[2]
            cur_right_r = snp[2]
            cur_left_q = snp[3]
            cur_right_q = snp[3]
            cluster[pair] = []
        if  pair == cur_pair:
            if abs(snp[2] - cur_right_r) < distance:
                cur_right_r = snp[2]
                cur_right_q = snp[3]
            else:
                cluster[cur_pair].append((cur_left_r, cur_right_r, cur_left_q, cur_right_q))
                cur_left_r = snp[2]
                cur_right_r = snp[2]
                cur_left_q = snp[3]
                cur_right_q = snp[3]
        else:
            cluster[cur_pair].append((cur_left_r, cur_right_r, cur_left_q, cur_right_q))
            cur_pair = pair
            cur_left_r = snp[2]
            cur_right_r = snp[2]
            cur_left_q = snp[3]
            cur_right_q = snp[3]
            cluster[pair] = []
    cluster[cur_pair].append((cur_left_r, cur_right_r, cur_left_q, cur_right_q))
    # Write output file
    with open(prefix + ".snp.cluster.txt", "w") as outp:
        lines = ["Ref_tig\tQuery_tig\tS1\tE1\tS2\tE2\n"]
        for pair, cluster_list in cluster.iteritems():
            for clu in cluster_list:
                line = "%s\t%s\t%d\t%d\t%d\t%d\n" % (pair + clu)
                lines.append(line)
        outp.write("".join(lines))
    _message("%d clusters for %s" %(len(lines)-1, prefix))

def simple_protein_align(seq1, seq2, matrix="BLOSUM62"):
    """Return the protein alingment score by dynamic programming."""
    # Load Score Matrix
    score_matrix = load_protein_score_matrix(matrix=matrix)

    # Start Dynamic programming
    matrix = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype=int)
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    for i in range(len(seq1)):
        matrix[i + 1, 0] = matrix[i, 0] + score_matrix[(seq1[i], "X")]
    for i in range(len(seq2)):
        matrix[0, i + 1] = matrix[i, 0] + score_matrix[(seq2[i], "X")]
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            match = matrix[i, j] + score_matrix[(seq1[i], seq2[j])]
            delelte = matrix[i, j+1] + score_matrix[(seq1[i], "X")]
            insert = matrix[i+1, j] + score_matrix[(seq2[j], "X")]
            matrix[i+1, j+1] = max(match, delelte, insert)
    return matrix[len(seq1), len(seq2)]

def load_protein_score_matrix(matrix="BLOSUM62"):
    """Load protein score matrix"""
    dirn = "/home/zhuwei/seq_tools/my_python_scripts/analysis_lib/"
    score_matrix = {}
    amino_list = []
    with open(dirn + matrix + ".dat", "r") as matp:
        for line in matp:
            line = line.strip()
            if line[0] == "#" or not line:
                continue
            if not amino_list:
                amino_list = line.split()
                continue
            entry = line.split()
            amino1 = entry[0]
            for id_, score in enumerate(entry[1:]):
                score = int(score)
                amino2 = amino_list[id_]
                score_matrix[(amino1, amino2)] = score
    return score_matrix

def protein_match_score(seq1, seq2, matrix="BLOSUM62"):
    """Return the exact match score from the shared part from N-term"""
    len_ = min([len(seq1), len(seq2)])
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    score_matrix = load_protein_score_matrix(matrix=matrix)
    score = 0
    seq_pair = zip(list(seq1)[:len_], list(seq2)[:len_])
    for pair in seq_pair:
        score += score_matrix[pair]
    return score

def protein_self_align(dirn, fname, window=10, threshold=.62, maxhit=1, matrix="BLOSUM62", prefix="out"):
    """Protein self alignment to discover repeat region.
    Only the non-repetitive N terminus is pooled out.

    """
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fname)
    n_unique_end = []
    for id_, seq in enumerate(seqs):
        find_multi = False
        p = 0
        while p < len(seq) - window + 1:
            seq1 = seq[p: p+window]
            count = 0
            # score_threshold = int(protein_match_score(seq1, seq1, matrix=matrix) * threshold)
            score_threshold = 23
            for q in xrange(p+1, len(seq) - window + 1):
                seq2 = seq[q: q+window]
                score = protein_match_score(seq1, seq2, matrix=matrix)
                if score > score_threshold:
                    count += 1
                    if count > maxhit:
                        find_multi = True
                        break
            if find_multi:
                break
            p += 1
        n_unique_end.append(p)
        _message("%s finished! N-unqiue until %d" %(names[id_], p))
    len_list = [len(seq) for seq in seqs]
    unique_seq_list = [seqs[id_][: n_unique_end[id_]] for id_ in xrange(len(names))]

    simple_fasta_write(prefix + ".unique.fasta", names, unique_seq_list)
    with open(prefix + ".unique.region", "w") as fp:
        lines = []
        for id_, pos in enumerate(n_unique_end):
            len_ = len_list[id_] - window + 1
            name = names[id_]
            coverage = 1. * (pos + 1) / (len_ + 1)
            line = "%s, %d, %d, %2.2f\n" %(name, pos + 1, len_, coverage)
            lines.append(line)
        fp.write("".join(lines))

def protein_self_dotmatcher(dirn, fname, prefix="out"):
    """Use dotmatcher from emboss for protein self dotmatrix"""
    os.chdir(dirn)
    if not os.path.isdir(dirn + prefix):
        os.mkdir(dirn + prefix)
    names, seqs = simple_fasta_load(fname)
    os.system("cp %s ./%s/" %(fname, prefix))
    os.chdir(dirn + prefix)
    for name in names:
        os.system("dotmatcher %s:%s %s:%s -graph png -goutfile %s" %(fname, name, fname, name, name))


def compare_snp_cluster_5x_58x():
    """Compare snp cluster between 5x and 58x"""
    dirn = "/home/zhuwei/cglabarata/comp/chrcomp/170912/snp/"
    name_5 = [51, 52, 53, 54]
    name_58 = [580, 582, 583, 585]
    for n_5 in name_5:
        for n_58 in name_58:
            fname = "%d.%d" %(n_5, n_58)
            snp_cluster(dirn=dirn, fname=fname + ".snps", prefix=fname)

def find_share_snp_5x_58x():
    """Find the conserved snp."""
    dirn = "/home/zhuwei/cglabarata/comp/chrcomp/170913/pairwise/"
    name_5 = [51, 52, 53, 54]
    name_58 = [580, 582, 583, 585]
    os.chdir(dirn)

    shared_snps = {}

    # Make pairwise Comparison
    snp_list = {}
    for i in range(len(name_5)):
        for j in range(len(name_58)):
            fn_5 = "%d-con.fa" % name_5[i]
            fn_58 = "%d-con.fa" % name_58[j]
            if not os.path.isfile("%s.%s.snps" %(name_5[i], name_58[j])):
                _message("dnadiff -p %s.%s %s %s" %(name_5[i], name_58[j], fn_5, fn_58))
                os.system("dnadiff -p %s.%s %s %s" %(name_5[i], name_58[j], fn_5, fn_58))
            snps = load_snps_file("%s.%s.snps" %(name_5[i], name_58[j]))
            snp_list[(i, None, j, None)] = snps

    for i in range(len(name_5)):
        for j in range(i):
            fn_1 = "%d-con.fa" % name_5[i]
            fn_2 = "%d-con.fa" % name_5[j]
            if not os.path.isfile("%s.%s.snps" %(name_5[i], name_5[j])):
                _message("dnadiff -p %s.%s %s %s" %(name_5[i], name_5[j], fn_1, fn_2))
                os.system("dnadiff -p %s.%s %s %s" %(name_5[i], name_5[j], fn_1, fn_2))
            snps = load_snps_file("%s.%s.snps" %(name_5[i], name_5[j]))
            snp_list[(i, j, None, None)] = snps
            snps_r = {}
            for pos_r, snp in snps.iteritems():
                pos_q = snp[1]
                type_ = snp[0][::-1]
                snps_r[pos_q] = (type_, pos_r)
            snp_list[(j, i, None, None)] = snps_r

    for i in range(len(name_58)):
        for j in range(i):
            fn_1 = "%d-con.fa" % name_58[i]
            fn_2 = "%d-con.fa" % name_58[j]
            if not os.path.isfile("%s.%s.snps" %(name_58[i], name_58[j])):
                _message("dnadiff -p %s.%s %s %s" %(name_58[i], name_58[j], fn_1, fn_2))
                os.system("dnadiff -p %s.%s %s %s" %(name_58[i], name_58[j], fn_1, fn_2))
            snps = load_snps_file("%s.%s.snps" %(name_58[i], name_58[j]))
            snp_list[(None, None, i, j)] = snps
            snps_r = {}
            for pos_r, snp in snps.iteritems():
                pos_q = snp[1]
                type_ = snp[0][::-1]
                snps_r[pos_q] = (type_, pos_r)
            snp_list[(None, None, j, i)] = snps_r

    comp_id_list = []
    for i in range(len(name_5)):
        for j in range(i):
            for m in range(len(name_58)):
                for n in range(m):
                    comp_id_list.append((i, j, m, n))

    # # SWAP 5x and 58x
    # for i in range(4):
    #     for j in range(4):
    #         if i == j:
    #             continue
    #         tmp = snp_list[(i, j, None, None)]
    #         snp_list[(i, j, None, None)] = snp_list[(None, None, i, j)]
    #         snp_list[(None, None, i, j)] = tmp
    # for i in range(4):
    #     for j in range(4):
    #         snp_list_ = snp_list[(i, None, j, None)]
    #         tmp = {}
    #         for pos_r, snp_ in snp_list_.items():
    #             pos_q = snp_[1]
    #             type_ = snp_[0]
    #             tmp[pos_q] = (type_[::-1], pos_r)
    #         snp_list[(i, None, j, None)] = tmp
    #         snp_list_.clear()

    # # i, j, m, n
    # # shared kmer between i, j, k, n
    # for id_ in comp_id_list:
    #     i = id_[0]
    #     j = id_[1]
    #     m = id_[2]
    #     n = id_[3]
    #     dif_5 = [t for t in range(len(name_5)) if t not in [i, j]]
    #     dif_58 = [t for t in range(len(name_58)) if t not in [m, n]]
    #     # snp i -> dif_5 -> j
    #     snps_5_1 = snp_list[(i, dif_5[0], None, None)]
    #     snps_5_2 = snp_list[(i, dif_5[1], None, None)]
    #     snps_5_1_3 = snp_list[(dif_5[0], j, None, None)]
    #     snps_5_2_3 = snp_list[(dif_5[1], j, None, None)]
    #     # snp i -> dif_58
    #     snps_5_58_1 = snp_list[(i, None, dif_58[0], None)]
    #     snps_5_58_2 = snp_list[(i, None, dif_58[1], None)]

    #     # snp dif_58 -> [m, n]
    #     snps_58_1_0 = snp_list[(None, None, dif_58[0], m)]
    #     snps_58_1_3 = snp_list[(None, None, dif_58[0], n)]
    #     snps_58_2_0 = snp_list[(None, None, dif_58[1], m)]
    #     snps_58_2_3 = snp_list[(None, None, dif_58[1], n)]
    #     # Shared SNPS
    #     # i -> dif_5[0] = dif5[1] -> j (reverse)
    #     # i -> dif_58[0] = dif_58[1] -> m = n (reverse)
    #     for pos_5_0, snp_501 in snps_5_1.iteritems():
    #         if ((pos_5_0 not in snps_5_2)
    #             or (pos_5_0 not in snps_5_58_1)
    #             or (pos_5_0 not in snps_5_58_2)):
    #             continue
    #         type_501 = snp_501[0]
    #         pos_5_1 = snp_501[1]
    #         snp_502 = snps_5_2[pos_5_0]
    #         type_502 = snp_502[0]
    #         pos_5_2 = snp_502[1]
    #         # if type_501 != type_502:
    #         #     continue
    #         if (pos_5_1 not in snps_5_1_3) or (pos_5_2 not in snps_5_2_3):
    #             continue
    #         if snps_5_1_3[pos_5_1] != snps_5_2_3[pos_5_2]:
    #             continue
    #         type_513 = snps_5_1_3[pos_5_1][0]
    #         # if type_513[::-1] != type_501:
    #         #     continue
    #         names = (name_5[id_[0]], name_5[id_[1]])
    #         pos_5_3 = snps_5_1_3[pos_5_1][1]
    #         line = "%s\t%s\t%s\t%d\t%s\t%d\t%s" %(names + pos_5_1 + pos_5_3 + (type_501, ))
    #         print line
    #         # The SNP is shared between #i and #j in 5x
    #         # Check whether it is shared with #m and #n in 58x
    #         snp_55801 = snps_5_58_1[pos_5_0]
    #         snp_55802 = snps_5_58_2[pos_5_0]
    #         # if snp_55801[0] != snp_55802[0] or snp_55801[0] != type_501:
    #         #     continue

    #         pos_58_1 = snp_55801[1]
    #         pos_58_2 = snp_55802[1]
    #         if ((pos_58_1 not in snps_58_1_0)
    #             or (pos_58_1 not in snps_58_1_3)
    #             or (pos_58_2 not in snps_58_2_0)
    #             or (pos_58_2 not in snps_58_2_3)):
    #             continue

    #         if ((snps_58_1_0[pos_58_1] != snps_58_2_0[pos_58_2])
    #             or (snps_58_1_3[pos_58_1] != snps_58_2_3[pos_58_2])):
    #             continue
    #         type_5810 = snps_58_1_0[pos_58_1][0]
    #         type_5823 = snps_58_2_3[pos_58_2][0]
    #         # if (type_5810 != type_5823) or (type_5810 != type_513):
    #         #     continue
    #         # Shared SNP obtained
    #         shared_snps.setdefault((i, j, m, n), [])
    #         pos_5_3 = snps_5_1_3[pos_5_1][1]
    #         pos_58_0 = snps_58_1_0[pos_58_1][1]
    #         pos_58_3 = snps_58_2_3[pos_58_2][1]
    #         shared_snps[(i, j, m, n)].append(pos_5_0 + pos_5_3 + pos_58_0 + pos_58_3 + (type_501, ))
    # # Save shared snps
    # with open("shared.snps", "w") as filep:
    #     lines = []
    #     for id_, snp_list in shared_snps.iteritems():
    #         for snp in snp_list:
    #             names = (name_5[id_[0]], name_5[id_[1]], name_58[id_[2]], name_58[id_[3]])
    #             line = "%d\t%d\t%d\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%s\n" %(names + snp)
    #             lines.append(line)
    #     filep.write("".join(lines))
    # # Save summary file
    # with open("shared.snps.summary", "w") as filep:
    #     lines = []
    #     for id_, snp_list in shared_snps.iteritems():
    #         names = (name_5[id_[0]], name_5[id_[1]], name_58[id_[2]], name_58[id_[3]])
    #         line = "%d\t%d\t%d\t%d\t%d\n" %(names + (len(snp_list), ))
    #         lines.append(line)
    #     filep.write("".join(lines))

    # Find unqiue snp
    shared_single_snps = [[] for i in range(16)]
    for i_5x in range(4):
        for i_58x in range(4):
            o_5x = [id_ for id_ in range(4) if id_ != i_5x]
            o_58x = [id_ for id_ in range(4) if id_ != i_58x]
            for pos_501 in snp_list[(i_5x, o_5x[0], None, None)].keys():
                # snp with other 5x strains
                check = False
                for j in o_5x[1:]:
                    if pos_501 not in snp_list[(i_5x, j, None, None)]:
                        check = True
                        break
                if check:
                    continue


                for j in o_5x:
                    snp = snp_list[(i_5x, j, None, None)][pos_501]
                    pos_5o = snp[1]
                    type_ = snp[0]
                    if pos_5o not in snp_list[(j, None, i_58x, None)]:
                        check = True
                        break
                    snp_2 = snp_list[(j, None, i_58x, None)][pos_5o]
                    # if snp_2[0][::-1] != type_:
                    #     check = True
                    #     break
                if check:
                    continue

                print name_5[i_5x], name_58[i_58x], pos_501

                for j in o_58x:
                    if pos_501 not in snp_list[(i_5x, None, j, None)]:
                        check = True
                        # print "No Hit"
                        break
                    snp = snp_list[(i_5x, None, j, None)][pos_501]
                    type_ = snp[0]
                    pos_58o = snp[1]
                    if pos_58o not in snp_list[(None, None, j, i_58x)]:
                        check = True
                        # print "No SecHit"

                        break
                    snp_2 = snp_list[(None, None, j, i_58x)][pos_58o]
                    # if snp_2[0][::-1] != type_:
                    #     check = True

                    #     break
                if check:
                    continue
                pos_58 = snp_2[1]
                # print name_5[i_5x], name_58[i_58x], name_58[j], pos_501, snp_2
                shared_single_snps[i_5x * 4 + i_58x].append(pos_501 + pos_58)
    # for i in range(4):
    #     for j in range(4):
    #         _message("%d shared snps between %s and %s" %(len(shared_single_snps[i*4+j]), name_5[i], name_58[j]))
    #         for pos in shared_single_snps[i*4 + j]:
    #             print  pos




def load_snps_file(snps_file):
    """Load the snps file by Mummer."""
    snps = {}
    with open(snps_file, "r") as filep:
        for line in filep:
            entry = line.split()
            pos_r = (entry[-2], int(entry[0]))
            type_ = entry[1] + entry[2]
            pos_q = (entry[-1], int(entry[3]))
            snps[pos_r] = (type_, pos_q)
    return snps

def extrac_share_snp_5x_58x():
    """Extract shared snp pair."""
    dirn = "/home/zhuwei/cglabarata/comp/chrcomp/170912/snp/"
    name_5 = [51, 52, 53, 54]
    name_58 = [580, 582, 583, 585]
    os.chdir(dirn)
    share_snp_58x = {}
    for n_58 in name_58:
        snps = {}
        shared_snps = {}
        for n_5 in name_5:
            fname = "%d.%d" %(n_5, n_58)
            with open(fname + ".snps", "r") as filep:
                for line in filep:
                    line.strip()
                    entry = line.split()

def split_fasta(maxcount=500, fname="", dirn="", prefix="out"):
    """Split fasta file by 500."""
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fname)
    names = [name.split()[0] for name in names]
    names_split = [names[id_ :id_+maxcount] for id_ in range(0, len(names), maxcount)]
    seqs_split =  [seqs[id_ :id_+maxcount] for id_ in range(0, len(names), maxcount)]
    for id_, names_s in enumerate(names_split):
        seqs_s = seqs_split[id_]
        simple_fasta_write(prefix + "-%02d.fa" % id_, names_s, seqs_s)

def _extract_blast_tab(fname, dirn, cutoff=0.1, min_pospercent=0.8):
    """Extract tabular blast result for annotation."""
    os.chdir(dirn)
    evidence = []
    alignment = {}
    with open(fname, "r") as refp:
        cur_query = None
        cur_ref = None
        cur_percent = 0.
        for line in refp:
            line.strip()
            if line.startswith("#"):
                continue
            if not line:
                continue
            entry = line.split()
            # New query
            if entry[0] != cur_query:
                evidence.append(line)
                cur_query = entry[0]
                cur_ref = entry[1]
                cur_percent = float(entry[-1])/int(entry[7])
                if cur_percent > min_pospercent:
                    alignment[cur_query] = [cur_ref]
                else:
                    alignment[cur_query] = []
            # Add alignment
            else:
                # Same ref, filter out
                if entry[1] == cur_ref:
                    continue
                # Different ref
                else:
                    new_percent = float(entry[-1])/int(entry[7])
                    if cur_percent - new_percent < cutoff and new_percent > min_pospercent:
                        evidence.append(line)
                        cur_ref = entry[1]
                        cur_percent = new_percent
                        alignment[cur_query].append(cur_ref)
    return alignment, evidence

def reciprocal_assign(ref_align, query_align):
    """Reciprocal Assignment of blastp result."""
    assigned_pair = {}
    unassigned_ref = []
    unassigned_query = {}
    unpaired_query = []
    unpaired_ref = []
    paired_ref = []
    assigned_ref = []
    # Assign by best alignment
    for q_prot, q_prot_list in query_align.iteritems():
        # Unpaired Query
        if not q_prot_list:
            unpaired_query.append(q_prot)
            continue
        r_prot = q_prot_list[0]
        # Reciprocal best
        # Partial Align
        if not ref_align[r_prot]:
            unpaired_query.append(q_prot)
            continue
        if q_prot == ref_align[r_prot][0]:
            assigned_pair[q_prot] = r_prot
            assigned_ref.append(r_prot)
            paired_ref.append(r_prot)
        else:
            unassigned_query[q_prot] = q_prot_list
            for r_prot in q_prot_list:
                if r_prot not in paired_ref:
                    paired_ref.append(r_prot)

    unpaired_ref = [key for key in ref_align.keys() if key not in paired_ref]
    unassigned_ref = [key for key in ref_align.keys() if key not in assigned_ref
                      and key not in unpaired_ref]
    _message("Reciprocal Assignment Finished")
    _message("%d Pairs Assigned" % len(assigned_pair))
    _message("%d Reference and %d Query are unassigned"
             %(len(unassigned_ref), len(unassigned_query)))
    _message("%d Reference and %d Query are unpaired"
             %(len(unpaired_ref), len(unpaired_query)))
    return assigned_pair, unassigned_ref, unassigned_query, unpaired_ref, unpaired_query

def linkage_map_by_name(namelist, pattern):
    """Build linkage map by nomenclature.
    Pattern Requires 2 items, 1st for chr assigment, 2nd for sequence info."""
    name_pattern = re.compile(pattern)
    linkage_map = {}
    chr_list = []
    name_extract = []
    for name in namelist:
        match_ = name_pattern.match(name)
        chr_ = match_.group(1)
        id_ = int(match_.group(2))
        if chr_ not in chr_list:
            chr_list.append(chr_)
        name_extract.append((chr_, id_))
    tmp = zip(name_extract, namelist)
    tmp.sort()
    name_extract, namelist = zip(*tmp)
    pre = None
    for no_, extract in enumerate(name_extract):
        chr_ = extract[0]
        if chr_ not in linkage_map:
            linkage_map[chr] = namelist[no_]
        else:
            linkage_map[pre] = namelist[no_]
        pre = namelist[no_]
    return chr_list, linkage_map

def linkage_map_gff(gff_entry):
    """Build Linkage Map by gff file."""
    chr_list = []
    link_map = {}
    id_pattern = re.compile(r"ID=([\w\d_.-])*")
    for entry in gff_entry:
        type_ = entry[2]
        if type_ != "gene":
            continue
        chr_ = entry[0]
        start = entry[3]
        note = entry[8]
        match_ = id_pattern.match(note)
        id_ = match_.group(1)


def fix_maker_blastx(maker_gff, maker_protein, dirn, prefix):
    """Fix signle exon maker protein late start codon problem by Blastx result"""
    os.chdir(dirn)
    maker_entry, maker_names, maker_seqs = simple_gff3_load(maker_gff, return_fasta=True)
    maker_seqs = [seq.upper() for seq in maker_seqs]

    stop_skip = 3
    # Extract Maker Gene
    augustus_masked = []
    maker_gene = {}
    maker_pos2name = {}
    protein_blastx = []
    other_entry = []

    augustus_masked_name = []
    augustus_abinitio = []
    maker_protein_region = []

    maker_blastx = {}
    fix_names = []
    fix_seqs = []
    fix_maker_entry = []

    for entry in maker_entry:
        if entry[1] == "maker" and entry[2] == "gene":
            gene_name = entry[8][3:]
            maker_gene[gene_name] = entry
            augustus_masked_name.append(gene_name)
            pos_ = (entry[0], entry[3], entry[4])
            maker_protein_region.append(pos_)
            maker_pos2name[pos_] = gene_name
            continue
        elif entry[1] == "augustus_masked" and entry[2] == "match":
            augustus_masked.append(entry)
            continue
        elif entry[1] == "protein2genome" and entry[2] == "protein_match":
            protein_blastx.append(entry)
        else:
            other_entry.append(entry)
    # Ab initio genes without any
    for entry in augustus_masked:
        name_ = entry[-1][5:]
        name_ = name_.replace("abinit", "processed")
        name_ = name_.replace("-mRNA-1", "")
        if name_ not in augustus_masked_name:
            augustus_abinitio.append(entry)
    # Check overlap of augustus call genes with blastx and fix late start codon
    fix_count = 0
    for entry in protein_blastx:
        direction_ = entry[6]
        if stop_skip == 0:
            pos_ = (entry[0], entry[3], entry[4])
        else:
            if entry[6] == "+":
                pos_ = (entry[0], entry[3], entry[4] + stop_skip)
            else:
                pos_ = (entry[0], entry[3] - stop_skip, entry[4])
        # Exact Overlap
        if pos_ in maker_protein_region:
            maker_blastx.setdefault(pos_, [])
            maker_blastx[pos_].append(entry[-1][5:])
        # Find End Alingment
        else:
            if (pos_[2] - pos_[1]) % 3 != 2:
                continue
            for pos_m in maker_protein_region:
                if (direction_ == "+" and pos_m[2] == pos_[2] and pos_m[0] == pos_[0] and pos_[1] < pos_m[1]) or (
                    direction_ == "-" and pos_m[1] == pos_[1] and pos_m[0] == pos_[0] and pos_[2] > pos_m[2]) :
                    # Test Codon
                    id_ = maker_names.index(pos_[0])
                    maker_name_ = maker_pos2name[pos_m]
                    seq_ = maker_seqs[id_][pos_[1] - 1 : pos_[2]]
                    if direction_ == "-":
                        seq_ = rc_seq(seq_)
                    if seq_.startswith("ATG"):
                        trans_ = translate_exon(seq_)
                        if "*" not in trans_[:-1]:
                            if maker_name_ in maker_gene:
                                fix_count += 1
                                fix_names.append(maker_name_ + "-fix")
                                maker_entry_ = maker_gene[maker_name_]
                                maker_entry_ = maker_entry_[:9]
                                maker_entry_[8] = "ID=" + maker_name_ + "-fix"
                                maker_entry_[3] = pos_[1]
                                maker_entry_[4] = pos_[2]
                                fix_maker_entry.append(maker_entry_)
                                fix_seqs.append(trans_)
                                del maker_gene[maker_name_]
                                _message("%s removed!" % maker_name_)
                            # Multi match
                            else:
                                id_ = fix_names.index(maker_name_ + "-fix")
                                len_old = len(fix_seqs[id_])
                                len_new = len(trans_)
                                if len_new > len_old:
                                    fix_seqs[id_] = trans_
                                fix_maker_entry[id_][3] = pos_[1]
                                fix_maker_entry[id_][4] = pos_[2]
                                _message("%s updated for longer version" % maker_name_)
                    break
    _message("%d Late stop codon fixed!" %(fix_count))
    # Write new gff file
    with open(prefix + ".fix.gff", "w") as filep:
        lines = ["##gff-version 3\n"]
        for entry in maker_gene.itervalues():
            entry = [str(ent) for ent in entry]
            line = "\t".join(entry[:8])
            line += "\t" + ";".join(entry[8:]) + "\n"
            lines.append(line)
        for entry in fix_maker_entry:
            entry = [str(ent) for ent in entry]
            line = "\t".join(entry[:8])
            line += "\t" + ";".join(entry[8:]) + "\n"
            lines.append(line)
        for entry in augustus_masked:
            entry = [str(ent) for ent in entry]
            line = "\t".join(entry[:8])
            line += "\t" + ";".join(entry[8:]) + "\n"
            lines.append(line)
        for entry in protein_blastx:
            entry = [str(ent) for ent in entry]
            line = "\t".join(entry[:8])
            line += "\t" + ";".join(entry[8:]) + "\n"
            lines.append(line)
        for entry in other_entry:
            entry = [str(ent) for ent in entry]
            line = "\t".join(entry[:8])
            line += "\t" + ";".join(entry[8:]) + "\n"
            lines.append(line)
        lines.append("##FASTA\n")
        for id_ in range(len(maker_names)):
            lines.append(">" + maker_names[id_] + "\n")
            seq_ = maker_seqs[id_]
            lines.extend(["%s\n" %seq_[pos: pos + 80] for pos in range(0, len(seq_), 80)])
        filep.write("".join(lines))
    # Load protein file
    pr_names, pr_seqs = simple_fasta_load(maker_protein)
    pr_names = [name.split()[0] for name in pr_names]
    for name_ in fix_names:
        name_ = name_[:-4]
        id_ = pr_names.index(name_ + "-mRNA-1")
        del pr_names[id_]
        del pr_seqs[id_]
    pr_names.extend(fix_names)
    pr_seqs.extend(fix_seqs)
    simple_fasta_write(prefix + ".fix.protein.fa", pr_names, pr_seqs)


def synteny_pair(paired_genes, assigned_ref, assigned_query, ref_linkage, query_linkage):
    """Pair the gene by synteny evidence."""
    pass

def auto_maker_annotation(query_prot="",
                          ref_prot="",
                          query_pattern="",
                          ref_pattern=r"CAGL0(\w)(\d*)",
                          query_blast_file="",
                          ref_blast_file="",
                          query_gff_file="",
                          ref_gff_file="",
                          dirn="", prefix="out",
                          cutoff=0.1):
    """Annotate maker generated maker protein file with referenece protein database."""
    os.chdir(dirn)
    # Generate blast file, the query/ref_blast_file will be overwritten

    if not (query_blast_file and ref_blast_file):
        _message("Generating Blastp database")
        os.system("makeblastdb -in %s -parse_seqids -dbtype prot" % query_prot)
        os.system("makeblastdb -in %s -parse_seqids -dbtype prot" % ref_prot)
        # Alignment file for top5 alignments
        os.system("blastp -query %s -db %s -num_threads 12 -num_alignments 5 -out %s.r2q.align"
                  %(ref_prot, query_prot, prefix))
        os.system("blastp -query %s -db %s -num_threads 12 -num_alignments 5 -out %s.q2r.align"
                  %(query_prot, ref_prot, prefix))
        os.system('''blastp -query %s -db %s \
-num_alignments 5 -num_threads 12 -out %s.r2q.tab \
-outfmt "7 qseqid sseqid qstart qend sstart send length  qlen slen evalue bitscore nident pident positive"'''
                  % (ref_prot, query_prot, prefix))
        os.system('''blastp -query %s -db %s \
-num_alignments 5 -num_threads 12 -out %s.q2r.tab \
-outfmt "7 qseqid sseqid qstart qend sstart send length  qlen slen evalue bitscore nident pident positive"'''
                  % (query_prot, ref_prot, prefix))
        ref_blast_file = "%s.r2q.tab" % prefix
        query_blast_file = "%s.q2r.tab" % prefix

    # Load tabular blast result
    # For each pair, only the first alignment is saved
    # Positive/Total is used to threshold the results
    # Alignment with Difference in Pos/Total > cutoff is filtered out
    _message("Loading Blast result")

    query2ref_blast, query2ref_evidence = _extract_blast_tab(query_blast_file, dirn, cutoff)
    ref2query_blast, ref2query_evidence = _extract_blast_tab(ref_blast_file, dirn, cutoff)

    _message("Writing Alignemnt Proof for annotation")
    with open(prefix + ".r2q.align.proof", "w") as filep:
        filep.write("".join(ref2query_evidence))
    with open(prefix + ".q2r.align.proof", "w") as filep:
        filep.write("".join(query2ref_evidence))
    # Reciprocal check for blast result
    _message("Reciprocal Assignment Initiated")
    assigned_pair, ref_unassign, query_unassign, ref_unpaired, query_unpaired = reciprocal_assign(
        ref2query_blast, query2ref_blast)

    # TODO: Test Reciprocal Assignment
    _message("Writing Reciprocal Result")
    with open(prefix + ".reci.pair.txt" ,"w") as filep:
        lines = []
        for query_p, ref_p in assigned_pair.iteritems():
            line = "%s\t%s\n" %(query_p, ref_p)
            lines.append(line)
        filep.write("".join(lines))
    with open(prefix + ".unassigned.ref.txt", "w") as filep:
        lines = []
        for ref_p in ref_unassign:
            lines.append(ref_p + "\n")
        filep.write("".join(lines))
    with open(prefix + ".unassigned.query.txt", "w") as filep:
        lines = []
        for ref_p in query_unassign:
            lines.append(ref_p + "\n")
        filep.write("".join(lines))
    with open(prefix + ".unpaired.ref.txt", "w") as filep:
        lines = []
        for ref_p in ref_unpaired:
            lines.append(ref_p + "\n")
        filep.write("".join(lines))
    with open(prefix + ".unpaired.query.txt", "w") as filep:
        lines = []
        for ref_p in query_unpaired:
            lines.append(ref_p + "\n")
        filep.write("".join(lines))


    # Load Genomic Positions
    ref_maker_genes = []
    if ref_gff_file:
        ref_gff = simple_gff3_load(ref_gff_file)
        ref_pr_names, ref_pr_seqs = simple_fasta_load(ref_prot)
        ref_pr_names = [name.split()[0] for name in ref_pr_names]
        for entry in ref_gff:
            if entry[1] == "CGD" and entry[2] == "gene":
                pr_name = entry[8][3:]
                pr_pos = (entry[0][:4], entry[3], entry[4])
                ref_maker_genes.append((pr_pos, pr_name, entry))
        ref_maker_genes.sort()
    else:
        pr_names, pr_seqs = simple_fasta_load(ref_prot)
        pr_names = [name.split()[0] for name in pr_names]
        pr_seqs = None
        pr_names.sort()
        for id_, name in enumerate(pr_names):
            ma_ = re.match(ref_pattern, name)
            chr_ = ma_.group(1)
            ref_maker_genes.append((("Chr%s" % chr_, id_*2, id_*2+1), name, "ID=%s" %name))

    _message("Reference linkage map loaded")

    query_gff = simple_gff3_load(query_gff_file)

    query_pr_names, query_pr_seqs = simple_fasta_load(query_prot)
    query_pr_names = [name.split()[0] for name in query_pr_names]
    query_maker_genes = []
    for entry in query_gff:
        if entry[1] == "maker" and entry[2] == "gene":
            pr_name = entry[8][3:]
            pr_pos = (entry[0], entry[3], entry[4])
            query_maker_genes.append((pr_pos, pr_name, entry))
    query_maker_genes.sort()
    _message("Query linkage map loaded")
    # Build Linkage Map
    cur_chr = ""
    cur_ref_id_ = ""
    query2reference_map = {}
    ref_linkage_map = {}
    query_linkage_map = {}
    ref_name2pos = {}

    for ref_pr in ref_maker_genes:
        r_pr_name = ref_pr[1]
        r_chr = ref_pr[0][0]
        if r_chr != cur_chr:
            ref_linkage_map[r_chr] = [r_pr_name]
            cur_chr = r_chr
        else:
            ref_linkage_map[r_chr].append(r_pr_name)
    for chr_, r_prot_list in ref_linkage_map.iteritems():
        for id_, r_prot in enumerate(r_prot_list):
            ref_name2pos[r_prot] = (chr_, id_)

    ref_linkage_map.clear()

    for query_pr in query_maker_genes:
        q_pr_name = query_pr[1] + "-mRNA-1"
        q_chr = query_pr[0][0]
        if q_chr != cur_chr:
            query_linkage_map[q_chr] = [q_pr_name]
            cur_chr = q_chr
        else:
            query_linkage_map[q_chr].append(q_pr_name)

    for chr_, q_prot_list in query_linkage_map.iteritems():
        query2reference_map[chr_] = []
        for id_, q_prot in enumerate(q_prot_list):
            # TODO: fix fosmid region
            if q_prot in assigned_pair:
                r_prot = assigned_pair[q_prot]
                r_prot_pos = ref_name2pos[r_prot]
                query2reference_map[chr_].append(r_prot_pos)
            elif q_prot in query_unpaired:
                query2reference_map[chr_].append(("unpaired", 0))
            else:
                query2reference_map[chr_].append(("unassigned", 0))
    # Save query2reference map
    with open(prefix + ".q2r.map", "w") as filep:
        lines = []
        for chr_, linkage in query2reference_map.iteritems():
            for id_, r_pos in enumerate(linkage):
                line = "%s\t%04d\t%s\t%04d\n" %(chr_, id_, r_pos[0], r_pos[1])
                lines.append(line)
        filep.write("".join(lines))
    _message("Linkage Map saved")
    # Draw synteny Map
    dpi = 300
    query_chr_list = []
    for chr_, linkage in query2reference_map.iteritems():
        linkage_map_y = []
        ref_chr_list = []
        ref_chr_max_min = {}
        ref_chr_top = []
        for r_pos in linkage:
            ref_chr = r_pos[0]
            ref_id = r_pos[1]
            if ref_chr not in ref_chr_list:
                ref_chr_list.append(ref_chr)
                ref_chr_max_min[ref_chr] = [ref_id, ref_id]
            else:
                if ref_chr_max_min[ref_chr][0] < ref_id:
                    ref_chr_max_min[ref_chr][0] = ref_id
                elif ref_chr_max_min[ref_chr][1] > ref_id:
                    ref_chr_max_min[ref_chr][1] = ref_id
        total = 0
        ref_chr_list.sort()
        for ref_chr in ref_chr_list:
            min_ = ref_chr_max_min[ref_chr][1]
            max_ = ref_chr_max_min[ref_chr][0]
            size_ = max_ - min_ + 20
            ref_chr_max_min[ref_chr][1] = min_ - total
            ref_chr_max_min[ref_chr][0] = total + size_
            total += size_
            ref_chr_top.append(total)
        for r_pos in linkage:
            max_min = ref_chr_max_min[r_pos[0]]
            linkage_map_y.append(r_pos[1] - max_min[1] + 10)

        fig, axes = plt.subplots(dpi=dpi)
        linkage_map_x = range(len(linkage))

        axes.scatter(tuple(linkage_map_x), tuple(linkage_map_y), marker="," , s=.2)

        # Draw Separate line
        for pos in ref_chr_top:
            axes.plot((0, len(linkage)), (pos, pos), linestyle="-", color="grey", linewidth=0.2)

        fig_name = chr_ + ".linkage.png"

        linkage_map_yticks = []
        pre = 0
        for pos in ref_chr_top:
            linkage_map_yticks.append((pos + pre) / 2)
            pre = pos
        axes.set_yticks(linkage_map_yticks)
        axes.set_yticklabels(ref_chr_list, rotation='horizontal', ha='right')

        # plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        plt.savefig(fig_name, bbox_inches="tight", format="png")
        plt.close()






def kmer_fingerprint():
    """Find the specific kmers for one sequence with a group a sequence."""
    dirn = "/home/zhuwei/cglabarata/comp/58x/MR/"
    kmer_size = 20
    os.chdir(dirn)
    fn_58 = ["580-M.fa", "582-M.fa", "583-M.fa", "585-M.fa"]
    fn_5 = ["51-M.fa", "52-M.fa", "53-M.fa", "54-M.fa"]
    specific_kmer_5x = [[] for id_ in range(4)]
    specific_kmer_58x = [[] for id_ in range(4)]
    seq_5x = []
    seq_58x = []
    name_5x = []
    name_58x = []
    match_count = [0] * 16
    for fname in fn_5:
        _name, _seq = simple_fasta_load(dirn + fname)
        _name = fname[:2]
        seq_5x.append(_seq[0].upper())
        name_5x.append(_name)
    for fname in fn_58:
        _name, _seq = simple_fasta_load(dirn + fname)
        _name = fname[:3]
        seq_58x.append(_seq[0].upper())
        name_58x.append(_name)
    # Specific kmers in reference
    for id_, seq in enumerate(seq_5x):
        _seq_ls = [seq_5x[i] for i in range(4) if i != id_]
        _kmer_dict = {}
        for pos in range(len(seq) - kmer_size):
            kmer = seq[pos: pos + kmer_size]
            kmer_rc = rc_seq(kmer)
            if kmer_rc < kmer:
                kmer = kmer_rc
            if kmer not in _kmer_dict:
                _kmer_dict[kmer] = ""
        _message("%d kmers loaded for %s" %(len(_kmer_dict), name_5x[id_]))
        for _seq in _seq_ls:
            for pos in range(len(_seq) - kmer_size):
                kmer = _seq[pos: pos + kmer_size]
                kmer_rc = rc_seq(kmer)
                if kmer_rc < kmer:
                    kmer = kmer_rc
                if kmer in _kmer_dict:
                    del _kmer_dict[kmer]
        specific_kmer_5x[id_] = _kmer_dict
        _message("%d specific kmers identified for %s" %(len(_kmer_dict), name_5x[id_]))
    # Specific kmers in query
    for id_, seq in enumerate(seq_58x):
        _seq_ls = [seq_58x[i] for i in range(4) if i != id_]
        _kmer_dict = {}
        for pos in range(len(seq) - kmer_size):
            kmer = seq[pos: pos + kmer_size]
            kmer_rc = rc_seq(kmer)
            if kmer_rc < kmer:
                kmer = kmer_rc
            if kmer not in _kmer_dict:
                _kmer_dict[kmer] = ""
        _message("%d kmers loaded for %s" %(len(_kmer_dict), name_58x[id_]))
        for _seq in _seq_ls:
            for pos in range(len(_seq) - kmer_size):
                kmer = _seq[pos: pos + kmer_size]
                kmer_rc = rc_seq(kmer)
                if kmer_rc < kmer:
                    kmer = kmer_rc
                if kmer in _kmer_dict:
                    del _kmer_dict[kmer]
        _message("%d specific kmers identified for %s" %(len(_kmer_dict), name_58x[id_]))
        specific_kmer_58x[id_] = _kmer_dict
    # Write specific kmers
    for id_, name in enumerate(name_5x):
        with gzip.open("%s.spec.kmer.txt.gz" % name, "w") as spep:
            lines = "\n".join(specific_kmer_5x[id_])
            spep.write(lines)

    for id_, name in enumerate(name_58x):
        with gzip.open("%s.spec.kmer.txt.gz" % name, "w") as spep:
            lines = "\n".join(specific_kmer_58x[id_])
            spep.write(lines)
    # Match chromosome
    _message("Start Pairing by kmer fingerprint")
    for id_, kmer_ls in enumerate(specific_kmer_58x):
        _message("Pairing %s" %(name_58x[id_]))
        for kmer in kmer_ls:
            for _id_, _kmer_ls_ in enumerate(specific_kmer_5x):
                if kmer in _kmer_ls_:
                    match_count[_id_ * 4 + id_] += 1
                    continue
    _message("Writing the comparison matrix")
    with open("match_count.txt", "w") as match_p:
        lines = []
        for row in range(4):
            line = "%d\t%d\t%d\t%d\n" % tuple(match_count[row*4: row*4 + 4])
            lines.append(line)
        line = "%d\t%d\t%d\t%d\n" % tuple([len(_kmer) for _kmer in specific_kmer_58x])
        lines.append(line)
        match_p.write("".join(lines))

def exract_cds(linewidth=60):
    """Extact Protein sequence from pre-defined txt file."""
    dirn = "/home/zhuwei/cglabarata/anno/cbs-fosmid/cbs-fosmid-fos-augustus/sys-name"
    os.chdir(dirn)
    txt_fn = "cbs.sys.name.txt"
    out_fn = "cbs.fos.protein.fasta"
    fos_fn = "cbs-fos.fa"

    names, seqs = simple_fasta_load(fos_fn)

    fosmid = {names[id_]: seqs[id_] for id_ in range(len(names))}

    protein_pos = []
    protein_fa = []

    with open(txt_fn, "r") as txt_f:
        for line in txt_f:
            entry = line.split()
            name = entry[0]
            fos = entry[4]
            if entry[3] == "sense":
                pos = (int(entry[1]) -1, int(entry[2]) -1)
            else:
                pos = (int(entry[2]) -1, int(entry[1]) -1)
            protein_pos.append((fos, name, pos))
    for pro_ in protein_pos:
        seq = fosmid[pro_[0]]
        pos = pro_[2]
        name = pro_[1]
        if pos[0] < pos[1]:
            pro_nt = seq[pos[0]: pos[1] + 1]
        else:
            pro_nt = seq[pos[1]: pos[0] + 1]
            pro_nt = rc_seq(pro_nt)
        protein = translate_exon(pro_nt)
        protein_fa.append((name, protein))

    with open(out_fn, "w") as out_f:
        lines = []
        for protein in protein_fa:
            lines.append(">%s\n" %protein[0])
            pro_seq = [protein[1][id_: id_ + linewidth] + "\n"
                       for id_ in range(0, len(protein[1]), linewidth)]
            lines.extend(pro_seq)
        out_f.write("".join(lines))


def sort_write_fasta():
    """Sort and write fasta file."""
    p_cbs = r"\w*"
    dirn = "/home/zhuwei/cglabarata/anno/bg2-canu/bg2-comp"
    os.chdir(dirn)
    names, seqs = load_fasta("cbs-canu.fa", p_cbs)
    simple_fasta_write("cbs-canu-sort.fa", names, seqs)

def get_junction(dirn, fn, query_fn, distance=10000, flank=5000, prefix="out"):
    """Obtain junction point of two chromosomes"""

    os.chdir(dirn)
    align = {}
    disjoint = []
    # Load alignment
    with open(fn, "r") as alp:
        for line in alp:
            line.lstrip()
            if not line:
                continue
            entry = line.split()
            pos_r = (int(entry[0]), int(entry[1]))
            pos_q = (int(entry[2]), int(entry[3]))
            # Ensure the Query is in W direction
            if pos_q[0] > pos_q[1]:
                pos_q = (pos_q[1], pos_q[0])
                pos_r = (pos_r[1], pos_r[0])
            is_same_direction = (pos_r[0] < pos_r[1])

            name_r = entry[4]
            name_q = entry[5]

            align_ = (name_r, pos_q, pos_r, is_same_direction)
            if name_q not in align:
                align[name_q] = [align_]
            else:
                align[name_q].append(align_)

    for name_q, aligns in align.iteritems():
        # Sort alignment
        aligns.sort()
        pre = aligns[0]
        count = 0
        for aln in aligns[1:]:
            # Different Chromosome
            if aln[0] != pre[0]:
                count += 1
                point = (name_q, pre[1][1], pre[0], "Change_Chr", pre[1][0], pre[1][1],
                         pre[2][0], pre[2][1], "%03d-1" % count)
                disjoint.append(point)
                point = (name_q, aln[1][0], aln[0], "Change_Chr", aln[1][0], aln[1][1],
                         aln[2][0], aln[2][1], "%03d-2" % count)
                disjoint.append(point)
            # Change Direction
            elif not pre[3] is aln[3]:
                count += 1
                point = (name_q, pre[1][1], pre[0], "Change_Direct", pre[1][0], pre[1][1],
                         pre[2][0], pre[2][1], "%03d-1" % count)
                disjoint.append(point)
                point = (name_q, aln[1][0], aln[0], "Change_Direct", aln[1][0], aln[1][1],
                         aln[2][0], aln[2][1], "%03d-2" % count)
                disjoint.append(point)
            # Large Distance
            elif abs(aln[1][0] - pre[1][1]) > distance:
                count += 1
                point = (name_q, pre[1][1], pre[0], "Large_Dist", pre[1][0], pre[1][1],
                         pre[2][0], pre[2][1], "%03d-1" % count)
                disjoint.append(point)
                point = (name_q, aln[1][0], aln[0], "Large_Dist", aln[1][0], aln[1][1],
                         aln[2][0], aln[2][1], "%03d-2" % count)
                disjoint.append(point)
            pre = aln
    # Write Junction point
    with open(prefix + "-junct.txt", "w") as junp:
        lines = []
        lines.append("Name_Q, Junct Pos, Name_R, Type, Q1, Q2, R1, R2, Junction_No\n")
        for junct in disjoint:
            line = "%s, %d, %s, %s, %d, %d, %d, %d, %s\n" % junct
            lines.append(line)
        junp.write("".join(lines))
    # Extract Junction sequence
    names, seqs = simple_fasta_load(query_fn)
    print names
    if not os.path.isdir(dirn + prefix):
        os.mkdir(dirn + prefix)
    os.chdir(dirn + prefix)
    for junct in disjoint:
        id_ = names.index(junct[0])
        seq = seqs[id_]
        pos = junct[1]
        pos_b =max([0, pos - flank -1])
        pos_e = min([pos + flank, len(seq)])
        flank_seq = [seq[pos_b: pos_e]]
        fname = junct[0] + "_%s.fa" % junct[8]
        name = [fname[:-3]]
        simple_fasta_write(fname, name, flank_seq)

def extract_fasta_namelist(dirn, fasta, namelist, prefix="out"):
    """Extract the fasta file by name list.
    Only the first part of the name is used
    """
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fasta)
    names_f = [name.split()[0] for name in names]
    seqs_filter = []
    names_filter = []
    with open(namelist, "r") as filep:
        for line in filep:
            line.strip()
            line = line[:-1]
            if not line or line=="\n":
                continue
            id_ = names_f.index(line)
            seqs_filter.append(seqs[id_])
            names_filter.append(names[id_])
    simple_fasta_write(prefix + ".fa", names_filter, seqs_filter)

def bridge_map(map1, map2, dirn, prefix="out"):
    """bridge mapped names."""
    os.chdir(dirn)
    map1_list = {}
    map2_list = {}
    map_bridge = {}
    with open(map1, "r") as filep:
        for line in filep:
            pair = line.split()
            map1_list[pair[0]] = pair[1]
    with open(map2, "r") as filep:
        for line in filep:
            pair = line.split()
            map2_list[pair[0]] = pair[1]
    for ori_name, map2_name in map2_list.iteritems():
        map1_name = map1_list[ori_name]
        map_bridge[map1_name] = map1_name + "-" + map2_name
    with open(prefix + ".map", "w") as filep:
        lines = []
        for map1_name, brdige_name in map_bridge.iteritems():
            line = map1_name + "\t" + brdige_name
            lines.append(line)
        filep.write("\n".join(lines))


def load_blastp_result(fname, dirn, prefix="out"):
    """Load default blastp result."""
    os.chdir(dirn)
    query_list = {}
    query_pattern = re.compile(r"Query= ([\w\d]*)")
    entry_pattern = re.compile(r">([\w\d]*)")
    socre_pattern = re.compile(r"Score = ([\d.]*) bits \(\d*\),  Expect = ([e\.\d-]*)")
    ide_pattern = re.compile(r"Identities = ([\d]*)/[\d]* \(([\d]*)\%\), Positives = ([\d]*)/[\d]* \(([\d]*)\%\), Gaps = ([\d]*)/([\d]*) \(([\d]*)\%\)")
    # Load alignment
    with open(fname, "r") as filep:
        name = None
        align = None
        for line in filep:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Query="):
                query = query_pattern.match(line)
                name = query.group(1)
                query_list.setdefault(name, [])
                continue
            if line.startswith(">"):
                entry = entry_pattern.match(line)
                align = [entry.group(1)]
                continue
            if line.startswith("Score ="):
                score = socre_pattern.match(line)
                align.extend(list(score.groups()))
                continue
            if line.startswith("Identities ="):
                ide = ide_pattern.match(line)
                align.extend(list(ide.groups()))
                tmp = align[-1]
                align[-1] = align[-2]
                align[-2] = tmp
                query_list[name].append(align)
                align = [align[0]]
    # Write alignment
    with open(prefix + ".score.txt", "w") as filep:
        lines = []
        lines.append("Query\tSubject\tScore\tExpect\tIdentity\tIde_Percent\tPositive\tPos_Percent\tGap\tGap_Percent\tTotal_Align\n")
        for name, align_list in query_list.iteritems():
            for align in align_list:
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %((name,) + tuple(align))
                lines.append(line)
        filep.write("".join(lines))
    # Write high e-value results
    with open(prefix + ".score.lowide.txt", "w") as filep:
        lines = []
        lines.append("Query\tSubject\tScore\tExpect\tIdentity\tIde_Percent\tPositive\tPos_Percent\tGap\tGap_Percent\tTotal_Align\n")
        for name, align_list in query_list.iteritems():
            align = align_list[0]
            ide = float(align[4])
            if ide < 80.:
                line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %((name,) + tuple(align))
                lines.append(line)
        filep.write("".join(lines))

def match_strain(ref, query, dirn, pattern=r"Chr(\w)", prefix="out"):
    chr_pattern = re.compile(pattern)
    os.chdir(dirn)
    if not os.path.isdir(dirn + prefix):
        os.mkdir(dirn + prefix)
    ref_names, ref_seqs = simple_fasta_load(ref)
    query_names, query_seqs = simple_fasta_load(query)
    os.chdir(dirn + prefix)

    for id_, name in enumerate(ref_names):
        simple_fasta_write("r-" + name+".fa", [name], [ref_seqs[id_]])
    for id_, name in enumerate(query_names):
        simple_fasta_write("q-" + name + ".fa", [name], [query_seqs[id_]])
    ref_id_list = []
    query_id_list = []
    for name in ref_names:
        match = chr_pattern.match(name)
        if match:
            id_ = match.group(1)
            ref_id_list.append(id_)
        else:
            ref_id_list.append(None)
    for name in query_names:
        match = chr_pattern.match(name)
        if match:
            id_ = match.group(1)
            query_id_list.append(id_)
        else:
            query_id_list.append(None)
    # Compare Chromosomes
    for id_, query_name in enumerate(query_names):
        query_id = query_id_list[id_]
        if not query_id:
            continue
        ref_index = ref_id_list.index(query_id)
        ref_name = ref_names[ref_index]
        pre_ = prefix+"-"+query_name
        # os.system("nucmer --mum -p %s %s %s"
        #           %(pre_, "r-" + ref_name+".fa", "q-" + query_name+".fa"))
        # os.system("show-coords -clrT %s.delta > %s.coords" %(pre_, pre_))
        # os.system("show-snps -ClrT -x 10 %s.delta > %s.snps" %(pre_, pre_))
        # os.system("nucmer --maxmatch -p max-%s %s %s"
        #           %(pre_, "r-" + ref_name+".fa", "q-" + query_name+".fa"))
        os.system("dnadiff -p dna-%s %s %s"
                  %(pre_, "r-" + ref_name+".fa", "q-" + query_name+".fa"))
        # os.system('''mummerplot --color -p max-%s --xrange "[0,%d]" --yrange "[0,%d]"  -png --large max-%s.delta'''
        #           %(pre_, len(ref_seqs[ref_index]), len(query_seqs[id_]), pre_))
        # os.system('''mummerplot --color -p %s --SNP --xrange "[0,%d]" --yrange "[0,%d]"  -png --large %s.delta'''
        #           %(pre_, len(ref_seqs[ref_index]), len(query_seqs[id_]), pre_))

def extract_gff3_region(entries, column, tag):
    """Extract regions from gff3 entries.
    Args:
    entry: gff3 entry from simple_gff3_load
    column: column to search for tag
    Returns:
    region = {contig_name: [(start, end, direction (1 = +, -1 = -), cds_name]}
    """
    region = {}
    for entry in entries:
        tig = entry[0]
        id_ = entry[8][3:]
        if tag in entry[column]:
            if entry[6] == "+":
                direction = 1
            else:
                direction = -1
            pos = (int(entry[3]), int(entry[4]), direction, id_)
            region.setdefault(tig, [])
            region[tig].append(pos)
    for tig, pos_list in region.items():
        pos_list.sort()
        region[tig] = pos_list
    return region

def load_coords(fname):
    """Load the coords file.
    Args: fname: name of the coords file
    Returns:
    coords = {(tig_query, tig_ref): (query_s, query_e, ref_s, ref_e, idy, cov_q, cov_r, direction)}
    """
    with open(fname) as filep:
        coords = {}
        for num, line in enumerate(filep):
            if num < 4:
                continue
            if not line:
                continue
            entry = line.split()
            tig = (entry[-1], entry[-2])
            entry = [float(item) for item in entry[:-2]]
            if (entry[0] - entry[1]) * (entry[2] - entry[3]) < 0:
                direction = -1
            else:
                direction = 1
            coords.setdefault(tig, [])
            # Ensure
            if entry[2] < entry[3]:
                coords[tig].append((entry[2], entry[3], entry[0], entry[1], entry[6], entry[10], entry[9], direction))
            else:
                coords[tig].append((entry[3], entry[2], entry[0], entry[1], entry[6], entry[10], entry[9], direction))
        return coords

def assigned_nucmer_compare(ref, query, query_ref_map):
    """Compare the assigned query sequence to the reference sequence.
    Args:
    ref: reference sequence in fasta format
    query: query sequence in fastq format
    query_ref_map:
        map query to reference sequence. {query_name: (ref_name, typeid)}
        typeid = 0: Optimal Paired; 1: Suboptimal Assigned (Multi assignment to one ref); 2: Unassigned
    ref_cds_region:
        region to filter SNPs in protein coding region
    """
    # make single fasta files for nucmer compare
    ref_names, ref_seqs = simple_fasta_load(ref)
    query_names, query_seqs = simple_fasta_load(query)
    ref_names_fa = [name + ".fa" for name in ref_names]
    query_names_fa = [name + ".fa" for name in query_names]
    snp_files = []
    coords_files = []
    figure_files = []
    pair_list = []
    for id_, ref_fa in enumerate(ref_names_fa):
        ref_seq = [ref_seqs[id_]]
        ref_name = [ref_names[id_]]
        simple_fasta_write(ref_fa, ref_name, ref_seq)
    for id_, query_fa in enumerate(query_names_fa):
        query_seq = [query_seqs[id_]]
        query_name = [query_names[id_]]
        simple_fasta_write(query_fa, query_name, query_seq)
    # Make snps and coords
    for id_, query_name in enumerate(query_names):
        ref_name = query_ref_map[query_name][0]
        type_ = query_ref_map[query_name][1]
        query_fa = query_names_fa[id_]
        # Optimal Paired
        print ref_name
        if type_ == 0:
            ref_fa = ref_name +".fa"
        # Suboptimal Paired
        elif type == 1:
            ref_fa = ref_name +".fa"
        # Unassigned
        else:
            ref_fa = ref
            ref_name = "ref"
        prefix = query_name + "-" + ref_name
        pair_list.append((query_name, ref_name))
        os.system("nucmer -mum -p %s %s %s" %(prefix, ref_fa, query_fa))
        os.system("show-coords -clqT %s.delta > %s.coords" %(prefix, prefix))
        os.system("show-snps -ClqT -x 10 %s.delta > %s.snps" %(prefix, prefix))
        os.system("mummerplot --medium -png -p %s --color %s.delta" %(prefix, prefix))
        os.system("nucmer -maxmatch -p rep-%s %s %s" %(prefix, ref_fa, query_fa))
        os.system("mummerplot --medium -png -p rep-%s --color rep-%s.delta" %(prefix, prefix))
        os.system("show-coords -clqT %s.delta > rep-%s.coords" %(prefix, prefix))
        os.system("nucmer -p rev-%s %s %s" %(prefix, query_fa, ref))
        os.system("mummerplot -medium -png -p rev-%s --color rev-%s.delta" %(prefix, prefix))
        coords_files.append(prefix + ".coords")
        snp_files.append(prefix + ".snps")
        figure_files.append(prefix + '.png')
        figure_files.append("rep-" + prefix + ".png")
        figure_files.append("rev-" + prefix + ".png")
    return pair_list, coords_files, snp_files, figure_files

def load_snp_file(fname):
    """Load snp file by numcer.
    Param to generate the snp file: -ClqT -x 10
    The micro indels are integrated together
    Args: fname: name of the snps file
    Returns:
    List of snps [(type, query_tig, (query_pos), ref_tig, (ref_pos), query_nt, ref_nt, query_x, ref_x)]
    for microIndels, the query_x/ref_x => (query/ref_x_first,query/ref_x_last) / etc

    """
    snps = []
    with open(fname, 'r') as filep:
        is_clustered = False
        ref_pos = [-1, -1]
        query_pos = [-1, -1]
        q_nt = ""
        r_nt = ""
        r_tig = ""
        q_tig = ""
        r_x = ""
        q_x = ""
        type_ = ""
        dist = -1
        for lin, line in enumerate(filep):
            if lin < 4:
                continue
            entry = line.split()
            if len(entry) != 14:
                _message("Error loading SNP file %s" % fname)
                return []
            cur_dist = int(entry[4])
            cur_r = int(entry[0])
            cur_q = int(entry[3])
            cur_r_nt = entry[1]
            cur_q_nt = entry[2]
            cur_r_x = entry[8]
            cur_q_x = entry[9]
            cur_r_tig = entry[12]
            cur_q_tig = entry[13]
            cur_q_tig = entry[-1]
            cur_r_tig = entry[-2]
            # New Alignment
            if cur_r_tig != r_tig or cur_q_tig != q_tig:
                if q_nt:
                    snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
                r_tig = cur_r_tig
                q_tig = cur_q_tig
                q_nt = ""
                r_nt = ""
                type_ = ""
                ref_pos = [-1, -1]
                query_pos = [-1, -1]
                is_clustered = False
            # clustered region
            if cur_dist == dist and (dist == 1 or dist == 0):
                # Elongation of a cluster
                if is_clustered:
                    # Micro Deletion
                    if query_pos[1] == cur_q or query_pos[0] == cur_q:
                        if entry[2] != ".":
                            _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                            _message(line)
                            return []
                        if ref_pos[1] == cur_r - 1:
                            ref_pos[1] = cur_r
                            r_nt = r_nt + cur_r_nt
                            r_x = (r_x[0], cur_r_x)
                        else:
                            ref_pos[0] = cur_r
                            r_nt = cur_r_nt + r_nt
                            r_x = (cur_r_x, r_x[1])
                        type_ = "mDel"
                    # Micro Insertion
                    elif ref_pos[1] == cur_r or ref_pos[0] == cur_r:
                        type_ = "mIn"
                        if entry[1] != ".":
                            _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                            _message(line)
                            return []
                        if query_pos[1] == cur_q - 1:
                            query_pos[1] = cur_q
                            q_nt = cur_q_nt + q_nt
                            q_x = (q_x[0], cur_q_x)
                        else:
                            query_pos[0] = cur_q
                            q_nt = cur_q_nt + q_nt
                            q_x = (cur_q_x, q_x[1])
                    # Do Not Record MicroSub for now
                    else:
                        is_clustered = False
                        if q_nt:
                            snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
                        q_nt = cur_q_nt
                        r_nt = cur_r_nt
                        ref_pos = [cur_r, cur_r]
                        query_pos = [cur_q, cur_q]
                        r_x = (cur_r_x, cur_r_x)
                        q_x = (cur_r_x, cur_q_x)
                        if entry[1] == ".":
                            type_ = "Del"
                        elif entry[2] == ".":
                            type_ = "Ins"
                        else:
                            type_ = "Sub"
                        # TODO Special Case to be solved
                        # TTGCCACCAC A T   CCTTCAACC
                        # TTGCCATCAC C T C CCTTCAAC
                        # Type AC and .C
                        # r_direction = int(entry[10])
                        # q_direction = int(entry[11])
                        # type_ = "mSub"
                        # if r_direction == 1 and q_direction == 1:
                        #     # if not (ref_pos[1] == cur_r -1 and query_pos[1] == cur_q -1):
                        #     #     _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                        #     #     _message(line)
                        #     #     return []
                        #     ref_pos[1] = cur_r
                        #     query_pos[1] = cur_q
                        #     r_nt = r_nt + cur_r_nt
                        #     q_nt = q_nt + cur_q_nt
                        #     r_x = (r_x[0], cur_r_x)
                        #     q_x = (q_x[0], cur_q_x)
                        # elif r_direction == 1 and q_direction == -1:
                        #     # if not (ref_pos[1] == cur_r -1 and query_pos[0] == cur_q + 1):
                        #     #     _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                        #     #     _message(line)
                        #     #     return []
                        #     ref_pos[1] = cur_r
                        #     query_pos[0] = cur_q
                        #     r_nt = r_nt + cur_r_nt
                        #     q_nt = cur_q_nt + q_nt
                        #     r_x = (r_x[0], cur_r_x)
                        #     q_x = (cur_q_x, q_x[1])
                        # elif r_direction == -1 and q_direction == 1:
                        #     # if not (ref_pos[0] == cur_r + 1 and query_pos[1] == cur_q - 1):
                        #     #     _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                        #     #     _message(line)
                        #     #     return []
                        #     ref_pos[0] = cur_r
                        #     query_pos[1] = cur_q
                        #     r_nt = cur_r_nt + r_nt
                        #     q_nt = q_nt + cur_q_nt
                        #     r_x = (cur_r_x, r_x[1])
                        #     q_x = (q_x[0], cur_q_x)
                        # else:
                        #     _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                        #     _message(line)
                        #     return []
                    # else:
                    #     _message("Error Loading SNP file %s at line %d: Unknown type" %(fname, lin))
                    #     _message(line)
                    #     return []
                # Init a cluster
                else:
                    dist = cur_dist
                    is_clustered = True
                    if q_nt:
                        snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
                    q_nt = cur_q_nt
                    r_nt = cur_r_nt
                    ref_pos = [cur_r, cur_r]
                    query_pos = [cur_q, cur_q]
                    r_x = (cur_r_x, cur_r_x)
                    q_x = (cur_r_x, cur_q_x)
                    if entry[1] == ".":
                        type_ = "Del"
                    elif entry[2] == ".":
                        type_ = "Ins"
                    else:
                        type_ = "Sub"
            elif cur_dist == 0 or cur_dist == 1:
                dist = cur_dist
                is_clustered = True
                if q_nt:
                    snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
                q_nt = cur_q_nt
                r_nt = cur_r_nt
                ref_pos = [cur_r, cur_r]
                query_pos = [cur_q, cur_q]
                r_x = (cur_r_x, cur_r_x)
                q_x = (cur_r_x, cur_q_x)
                if entry[1] == ".":
                    type_ = "Del"
                elif entry[2] == ".":
                    type_ = "Ins"
                else:
                    type_ = "Sub"
            else:
                # Normal SNP
                # Add previous feature
                is_clustered = False
                if q_nt:
                    snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
                q_nt = cur_q_nt
                r_nt = cur_r_nt
                ref_pos = [cur_r, cur_r]
                query_pos = [cur_q, cur_q]
                r_x = (cur_r_x, cur_r_x)
                q_x = (cur_r_x, cur_q_x)
                # Del
                if entry[1] == ".":
                    type_ = "Del"
                elif entry[2] == ".":
                    type_ = "Ins"
                else:
                    type_ = "Sub"
        if q_nt:
            snps.append((type_, q_tig, query_pos, r_tig, ref_pos, q_nt, r_nt, q_x, r_x))
    return snps



def classify_homopolymer_snp(string, homo=6, x=10):
    """Determine homopolymer variant by context string and inserted variants."""
    str_a = string[:x][::-1].upper()
    str_b = string[-x:].upper()

    def _count_first_nt(string):
        """Count the occurance of the first nt in the string."""
        first_nt = string[0]
        count = 0
        for nt_ in string:
            if  nt_ == first_nt:
                count += 1
            else:
                break
        return (first_nt, count)

    nt_a, count_a = _count_first_nt(str_a)
    nt_b, count_b = _count_first_nt(str_b)
    nt_ = ""
    if nt_a == nt_b:
        count = count_a + count_b
        nt_ = nt_a
    else:
        if count_a > count_b:
            nt_ = nt_a
            count = count_a
        else:
            nt_ = nt_b
            count = count_b
    if count > homo:
        return True, count, nt_
    else:
        return False, -1, ""

def generate_histo(number_list, bin_size=1):
    """Generate the numbers for the histogram.
    Args:
    numberlist: list of numbers to make the histogram.
    bin: size of the bin
    Returns:
    histo = [(bin_start, bin_count)]
    """
    if not number_list:
        return []
    min_ = min(number_list)
    max_ = max(number_list)
    count = [0] * (int(max_ - min_ + 1)/bin_size)
    histo = []
    for num in number_list:
        bin_id = int(num - min_)/bin_size
        count[bin_id] += 1
    for bin_id, bin_count in enumerate(count):
        bin_start = min_ + bin_id * bin_size
        histo.append((bin_start, bin_count))
    return histo

def classify_snp(snp_list, cds_dict, score_dict, x=10, homo=6):
    """Classify the snps
    The snps are generated by MUMMER with show-snps -ClqT -x
    Classsification of the snps:
    type = (is_coding, is_homopolymer, is_micro)
    Args:
    snp_list: list of snps: (type, query_tig, (query_pos), ref_tig, (ref_pos), query_nt, ref_nt, query_x, ref_x)
    cds_dict: {ref_contig_name: [(start, end, direction (1 = +, -1 = -))]}
    score_dict: {query_tig_name: Phred Score List}
    x: the up/downstream flanking region of the snp
    Returns:
    snps =
     [(query_tig, query_pos[0], query_pos[1], query_nt, ref_nt, type_, is_code, is_homo, is_m,
       score, avg_score, ref_tig, ref_pos[0], ref_pos[1], query_x[0], query_x[1], ref_x[0], ref_x[1])]
    is_code = (True, ((start, end, direction, name))) / (False, ())
    is_homo = "True; #count#homo_nt" / "False;"
    is_m = True/False
    score = "#score; #score..."
    histo_snp = [(Phred_score, count_of_snps)]
    histo_general = [(Phred_score, count_of_bases)]
    """
    snps = []
    snp_score = []
    general_score = []
    is_genetic = []
    is_homopolymer = []
    is_micro = []
    histo_snp = []
    histo_general = []

    for snp in snp_list:
        type_ = snp[0]
        query_tig = snp[1]
        query_pos = snp[2]
        ref_tig = snp[3]
        ref_pos = snp[4]
        query_nt = snp[5]
        ref_nt = snp[6]
        query_x = snp[7]
        ref_x = snp[8]

        # Check micro:
        if type_[0] == "m":
            is_micro.append(True)
        else:
            is_micro.append(False)
        # Check homopolymer
        if type_ == "mIn":
            context = ref_x[0]
            homo, number, nt_ = classify_homopolymer_snp(context)
            if homo:
                is_homopolymer.append("True; %d%s" %(number, nt_))
            else:
                is_homopolymer.append("False")
        elif type_ == "mDel":
            context = query_x[0]
            homo, number, nt_ = classify_homopolymer_snp(context)
            if homo:
                is_homopolymer.append("True; %d%s" %(number, nt_))
            else:
                is_homopolymer.append("False")
        else:
            context = ref_x[0]
            homo, number, nt_ = classify_homopolymer_snp(context)
            if homo:
                is_homopolymer.append("True; %d%s" %(number, nt_))
            else:
                is_homopolymer.append("False")
        # Check genetic
        cds_ = cds_dict[ref_tig]
        ref_pos_a = ref_pos[0]
        ref_pos_b = ref_pos[1]
        mut_exon = []
        if ref_pos_a == ref_pos_b:
            for exon in cds_:
                start, end, direction, name = exon
                if ref_pos_a > start and ref_pos_a < end:
                    mut_exon.append(exon)
        else:
            for exon in cds_:
                start, end, direction, name = exon
                if ref_pos_a > start and ref_pos_a < end:
                    mut_exon.append(exon)
                elif ref_pos_b > start and ref_pos_b < end:
                    mut_exon.append(exon)
        if mut_exon:
            is_genetic.append((True, tuple(mut_exon)))
        else:
            is_genetic.append((False, ()))

        # Load Phred Score
        score_list = [score_dict[query_tig][pos] for pos in range(query_pos[0]-1, query_pos[1])]

        snp_score.append(score_list)
    for q_tig, score_list in score_dict.iteritems():
        general_score.extend(score_list)

    histo_general = generate_histo(general_score)
    _snp_score = []
    for score in snp_score:
        _snp_score.extend(score)
    histo_snp = generate_histo(_snp_score)
    print len(histo_snp)
    # snps =
    # [(query_tig, query_pos_s, query_pos_e, query_nt, ref_nt, mut_type,is_coding_region, str score, average_score, ref_tig, (ref_pos))]
    for snp_id, snp in enumerate(snp_list):
        type_ = snp[0]
        query_tig = snp[1]
        query_pos = snp[2]
        ref_tig = snp[3]
        ref_pos = snp[4]
        query_nt = snp[5]
        ref_nt = snp[6]
        query_x = snp[7]
        ref_x = snp[8]
        is_homo = is_homopolymer[snp_id]
        is_code = is_genetic[snp_id]
        is_m = is_micro[snp_id]
        score_ = [str(sco) for sco in snp_score[snp_id]]
        score_ = "; ".join(score_)
        avg_score = sum(snp_score[snp_id]) * 1. / len(snp_score[snp_id])
        snps.append((query_tig, query_pos[0], query_pos[1], query_nt, ref_nt, type_, is_code,
                     is_homo, is_m, score_, avg_score, ref_tig, ref_pos[0], ref_pos[1],
                     query_x[0], query_x[1], ref_x[0], ref_x[1]))
    return snps, histo_snp, histo_general


def merge_histogram(list_of_histogram):
    """Merge histogram with same bin size, and no-overlapping bins."""
    histo = []
    _histos = []
    for histo_ in list_of_histogram:
        _histos.extend(histo_)
    _histos.sort()
    _bin = None
    _count = 0

    for bin_ in _histos:
        cur_bin = bin_[0]
        cur_count = bin_[1]
        if cur_bin != _bin:
            if _bin:
                histo.append((_bin, _count))
            _bin = cur_bin
            _count = cur_count
        else:
            _count += cur_count
    histo.append((_bin, _count))
    return histo

def compare_strain_by_tig_snp(ref, query, dirn, snp_files, ref_cds_gff, homo=6, x=10, prefix="out", ref_strain=""):
    """Compare strain by the snp on paired tig.
    snp format
     [(query_tig, query_pos[0], query_pos[1], query_nt, ref_nt, type_, is_code, is_homo, is_m,
       score, avg_score, ref_tig, ref_pos[0], ref_pos[1], query_x[0], query_x[1], ref_x[0], ref_x[1])]
    is_code = (True, ((start, end, direction, name))) / (False, ())
    is_homo = "True; #count#homo_nt" / "False;"
    is_m = True/False
    score = "#score; #score..."
    histo_snp = [(Phred_score, count_of_snps)]
    histo_general = [(Phred_score, count_of_bases)]
    """
    os.chdir(dirn)
    # ref_tigs, ref_seqs = simple_fasta_load(ref)
    query_tigs, query_seqs, query_scores = simple_fastq_load(query)
    score_dict = {}
    for id_, query_tig in enumerate(query_tigs):
        score_dict[query_tig] = query_scores[id_]
    snp_master_list = []
    histo_snp_master = []
    histo_general_master = []
    ref_gff_entries = simple_gff3_load(ref_cds_gff)
    # filter cds region
    cds_region = extract_gff3_region(ref_gff_entries, 2, 'CDS')
    if ref_strain:
        for tig in cds_region.keys():
            cds_region[ref_strain + "-" + tig[:-18]] = cds_region.pop(tig)

    for snp_file in snp_files:
        snp_list = load_snp_file(snp_file)
        snps, histo_snp, histo_general = classify_snp(snp_list, cds_region, score_dict)
        snp_master_list.extend(snps)
        histo_snp_master.append(histo_snp)

    general_score_list = []
    for score in query_scores:
        general_score_list.extend(score)
    histo_general_master = generate_histo(general_score_list)
    histo_snp_master = merge_histogram(histo_snp_master)

    _message("Writing snp file")
    # snp format
    #  [(query_tig, query_pos[0], query_pos[1], query_nt, ref_nt, type_, is_code, is_homo, is_m,
    #    score, avg_score, ref_tig, ref_pos[0], ref_pos[1], query_x[0], query_x[1], ref_x[0], ref_x[1])]
    # is_code = (True, ((start, end, direction, name))) / (False, ())
    # is_homo = "True; #count#homo_nt" / "False;"
    # is_m = True/False
    # score = "#score; #score..."
    with open(prefix + ".snp.csv", "w") as filep:
        lines = ["Query_tig, Query_Pos_S, Query_Pos_E, Query_nt, Ref_nt, "
                 "Type, is_exon, is_homo, is_m, score, average_score, Ref_tig, Ref_Pos_S, Ref_Pos_E, "
                 "Query_Context_1, Query_Context_2, Ref_Context_1, Ref_Context_2\n"]
        for snp in snp_master_list:
            line = "%s, %d, %d, %s, %s, %s, " % snp[:6]
            is_exon = snp[6]
            is_homo = snp[7]
            is_m = snp[8]
            if is_exon[0]:
                tmp = ["True"]
                for exon in is_exon[1]:
                    tmp.append("%s:%d-%d(%d)" %(exon[3], exon[0], exon[1], exon[2]))
                tmp = "; ".join(tmp)
                is_exon = tmp
            else:
                is_exon = "False"
            line += "%s, %s, %s, %s, %2.2f, %s, %d, %d, %s, %s, %s, %s\n" %((is_exon, is_homo, is_m) + snp[9:])
            lines.append(line)
        filep.write("".join(lines))
    _message("Writing Histograms")
    with open(prefix + ".general.histo", "w") as filep:
        lines = ["%d\t%d\n" % bin_ for bin_ in histo_general_master]
        filep.write("".join(lines))
    with open(prefix + ".snps.histo", "w") as filep:
        lines = ["%d\t%d\n" % bin_ for bin_ in histo_snp_master]
        filep.write("".join(lines))


def align_contig_to_ref(ref, query, prefix, dirn,
                        query_strain="", ref_strain="", ref_cds_gff="", phred_threshold=40, homomer=6,
                        remove_polish_tag=""):
    """Align the query contigs to reference genome.
    The query contigs from canu assembly is aligned to reference genome by NUCMER.
    The query is taken in fastq file, which Phred score is given by ASCII-32 code (Phred 0 => "!" - 32)
    The temporary fasta file is generated for nucmer alignment. The -mum option is used to minimize the
    effect of repetitive region.  The contigs with > 70% of its region assigned to one reference chromosome is assigned
    the corresponding name, where the other contigs are assigned as Chr%03d
    The direction of the chromosome is also adjusted to be in accordance with the reference genome.
    Mummerplot is used to generate 1) nucmer with maxmatch to show change in repetitive regions. 2) nucmer with mum to
    to distribution of the snps. The snps and coords are generated by nucmer -mum.
    The categorization of the snps:
    1) Phred Score from fastq file, filtered by phred_threshold
    2) Coding Region if ref_cds_gff is provided
    3) Homopolymers
    4) single snps: Sub(stitutions), Indel, microIndel are identified"""
    # Load files
    os.chdir(dirn)
    ref_names, ref_seqs = simple_fasta_load(ref)
    query_names, query_seqs, query_scores = simple_fastq_load(query)
    if remove_polish_tag:
        query_names = [name.replace(remove_polish_tag, "") for name in query_names]
    cds_region = {}
    # Load coding region
    if ref_cds_gff:
        ref_gff_entries = simple_gff3_load(ref_cds_gff)
        # filter cds region
        cds_region = extract_gff3_region(ref_gff_entries, 2, 'CDS')
    # make directory to same results
    if not os.path.isdir(dirn + prefix):
        os.mkdir(dirn + prefix)
    if not os.path.isdir(dirn + prefix + "/snps"):
        os.mkdir(dirn + prefix + "/snps")
    if not os.path.isdir(dirn + prefix + "/png"):
        os.mkdir(dirn + prefix + "/png")
    if not os.path.isdir(dirn + prefix + "/coords"):
        os.mkdir(dirn + prefix + "/coords")
    # make tmp files in tmp directory
    if not os.path.isdir(dirn + "tmp/"):
        os.mkdir(dirn + "tmp")
    os.chdir(dirn + "tmp")
    if ref_strain:
        ref_names = [ref_strain + "-" + name for name in ref_names]
        if ref_cds_gff:
            for tig in ref_cds_gff.keys():
                ref_cds_gff[ref_strain + "-" + tig] = ref_cds_gff.pop(tig)
    if query_strain:
        query_names = [query_strain + "-" + name for name in query_names]
    query_fa = query + ".fa"
    simple_fasta_write(query_fa, query_names, query_seqs)
    simple_fasta_write(ref, ref_names, ref_seqs)
    # Make nucmer comparison for chromosome name adjustment
    os.system("nucmer -p strain %s %s" %(ref, query_fa))
    os.system("show-coords -clqT strain.delta > strain.coords")
    # Load and analyze the coords file for contig assignemnt
    contig_assign_list = load_coords("strain.coords")
    query_ref_tig = {}
    query_ref_coverage = {}
    ref_query_coverage = {}
    for tig, pos_list in contig_assign_list.iteritems():
        q_plus_cov = 0
        q_minus_cov = 0
        r_plus_cov = 0
        r_minus_cov = 0
        for pos in pos_list:
            q_coverage = pos[-3]
            r_coverage = pos[-2]
            direction = pos[-1]
            if direction == 1:
                q_plus_cov += q_coverage
                r_plus_cov += r_coverage
            else:
                q_minus_cov += q_coverage
                r_plus_cov += r_coverage
        query_ref_coverage.setdefault(tig[0], [])
        query_ref_coverage[tig[0]].append((q_plus_cov + q_minus_cov, tig[1], q_plus_cov, q_minus_cov))
        ref_query_coverage.setdefault(tig[1], [])
        ref_query_coverage[tig[1]].append((r_plus_cov + r_minus_cov, tig[0], r_plus_cov,r_minus_cov))
    for tig, pos_list in ref_query_coverage.items():
        pos_list.sort(reverse=True)
        ref_query_coverage[tig] = pos_list
    for tig, pos_list in query_ref_coverage.items():
        pos_list.sort(reverse=True)
        query_ref_coverage[tig] = pos_list

    # Write the coverage information
    with open('coverage.q.txt', 'w') as filep:
        lines = []
        lines.append("Query_tig\tRef_tig\tPos_cov\tNeg_cov\tTotal\n")
        for query_tig, cov_list in query_ref_coverage.iteritems():
            for q_cov in cov_list:
                line = "%s\t%s\t%2.2f\t%2.2f\t%2.2f\n" %(query_tig, q_cov[1], q_cov[2], q_cov[3], q_cov[0])
                lines.append(line)
        filep.write("".join(lines))
        _message("Query Coverage Information saved")
    with open('coverage.r.txt', 'w') as filep:
        lines = []
        lines.append("Ref_tig\tQuery_tig\tPos_cov\tNeg_cov\tTotal\n")
        for ref_tig, cov_list in ref_query_coverage.iteritems():
            for r_cov in cov_list:
                line = "%s\t%s\t%2.2f\t%2.2f\t%2.2f\n" %(ref_tig, r_cov[1], r_cov[2], r_cov[3], r_cov[0])
                lines.append(line)
        filep.write("".join(lines))
        _message("Ref Coverage Information saved")
    # Reciprocal Best Assignment
    un_assigned = 0
    for query_tig, query_cov_list in query_ref_coverage.iteritems():
        best_ref = query_cov_list[0][1]
        best_query = ref_query_coverage[best_ref][0][1]
        # Reciprocol Best
        if best_query == query_tig:
            query_ref_tig[query_tig] = best_ref
        else:
            best_choice = False
            # Not Reciprocol Best, the best one is overwhelming good
            if len(query_cov_list) == 1:
                best_choice = True
            else:
                best_ref_coverage = query_cov_list[0][0]
                second_ref_coverage = query_cov_list[1][0]
                best_choice = (best_ref_coverage > 70 ) and (second_ref_coverage < 30)
            if best_choice:
                best_ref_query_list = ref_query_coverage[best_ref]
                for id_, query_ in enumerate(best_ref_query_list, 1):
                    if query_[1] == query_tig:
                        break
                query_assign = "%s-%02d" %(best_ref, id_)
            else:
                un_assigned += 1
                if query_strain:
                    query_assign = "%s-Chr%03d" % (query_strain, un_assigned)
                else:
                    query_assign = "Chr%03d" % (un_assigned)
            query_ref_tig[query_tig] = query_assign

    # Dict = [ori_query_name : {assgined_name, 1=Same Direction/ -1=RC seq)}
    query_update_list = {}
    for q_tig, assigned_q_tig in query_ref_tig.iteritems():
        if assigned_q_tig.startswith(ref_strain + "-"):
            p = len(ref_strain + "-")
            if query_strain:
                assigned_q_tig = query_strain + "-" + assigned_q_tig[p:]
            else:
                assigned_q_tig = assigned_q_tig[p:]
        plus = query_ref_coverage[q_tig][0][-2]
        minus = query_ref_coverage[q_tig][0][-1]
        if plus > minus:
            direction = 1
        else:
            direction = -1
        query_update_list[q_tig] = (assigned_q_tig, direction)
    # Saved Assigned Information
    for q_name in query_names:
        if q_name not in query_ref_tig:
            un_assigned += 1
            if query_strain:
                query_assign = "%s-Chr%03d" % (query_strain, un_assigned)
            else:
                query_assign = "Chr%03d" % (un_assigned)
            query_ref_tig[q_name] = query_assign
            query_update_list[q_name] = (query_assign, 1)
    with open("query_tig_assign.txt", 'w') as filep:
        lines = ["Query_tig\tAssigned_tig\tIs_Assigned2Ref\tQuery_Cov\tRef_Coverage\n"]
        for q_tig, r_tig in query_ref_tig.iteritems():

            if re.match(r'\d\d\d', r_tig[-3:]):
                lines.append("%s\t%s\tUnassigned\n" %(q_tig, r_tig))
            elif re.match(r'-\d\d', r_tig[-3:]):
                r_tig = r_tig[:-3]
                q_coverage = query_ref_coverage[q_tig][0]
                r_coverage = ref_query_coverage[r_tig][0]
                r_tig = query_update_list[q_tig][0]
                lines.append("%s\t%s\tAssigned\t%2.2f\t%2.2f\n" %(q_tig, r_tig, q_coverage[0], r_coverage[0]))
            else:
                q_coverage = query_ref_coverage[q_tig][0]
                r_coverage = ref_query_coverage[r_tig][0]
                r_tig = query_update_list[q_tig][0]
                lines.append("%s\t%s\tAssigned\t%2.2f\t%2.2f\n" %(q_tig, r_tig, q_coverage[0], r_coverage[0]))
        filep.write("".join(lines))
        _message('Contig Assignment information stored')
    # Make Updated Query Fasta and Fastq file

    # Write the query fasta file
    query_names_adj = [query_update_list[query_name][0] for query_name in query_names]
    print query_names_adj
    for id_, seq in enumerate(query_seqs):
        direction = query_update_list[query_names[id_]][1]
        if direction == -1:
            seq = rc_seq(seq)
            query_seqs[id_] = seq
            query_scores[id_] = query_scores[id_][::-1]

    q_file_seperate = query.split('.')
    q_name_major = ".".join(q_file_seperate[:-1])
    query_fa = q_name_major + ".adj.fa"
    simple_fasta_write(query_fa, query_names_adj, query_seqs)
    query = q_name_major + ".adj.fq"
    simple_fastq_write(query, query_names_adj, query_seqs, query_scores)
    _message("Adjusted Query Fasta and Fastq file are generated")
    os.system("cp %s ../%s/" % (query, prefix))
    os.system("cp %s ../%s/" % (query_fa, prefix))

    # Make Paired Comparison by nucmer
    # Make the map of assigned query to ref
    query_ref_map = {}
    for q_tig, r_tig in query_ref_tig.iteritems():
        q_assigned = query_update_list[q_tig][0]
        if re.match(r"\d\d\d", q_assigned[-3:]):
            type_ = 2
        elif re.match(r"-\d\d", q_assigned[-3:]):
            if r_tig.startswith(query_strain + "-"):
                p = len(query_strain + "-")
                if ref_strain:
                    r_tig = ref_strain + "-" + r_tig[p:][:-3]
                else:
                    r_tig = r_tig[p:]
            type_ = 1
        else:
            if ref_strain:
                r_tig = ref_strain + "-" + r_tig[p:]
            else:
                r_tig = r_tig[p:]
            type_ = 0
        query_ref_map[q_assigned] = (r_tig, type_)
    pair_list, coords, snps, figures = assigned_nucmer_compare(ref, query_fa, query_ref_map)
    os.system("mv *.snps ../%s/snps" % prefix)
    os.system("mv *.coords ../%s/coords" % prefix)
    os.system("mv *.png ../%s/png" % prefix)
    os.system("mv *.txt ../%s/" % prefix)
    # # Analyze the snps

def align_contig_to_ref_fasta(ref, query, prefix, dirn,
                              query_strain="", ref_strain="", ref_cds_gff="", phred_threshold=40, homomer=6,
                              remove_polish_tag=""):
    """Align the query contigs to reference genome. (Query is fasta file)
    The query contigs from canu assembly is aligned to reference genome by NUCMER.
    The query is taken in fastq file, which Phred score is given by ASCII-32 code (Phred 0 => "!" - 32)
    The temporary fasta file is generated for nucmer alignment. The -mum option is used to minimize the
    effect of repetitive region.  The contigs with > 70% of its region assigned to one reference chromosome is assigned
    the corresponding name, where the other contigs are assigned as Chr%03d
    The direction of the chromosome is also adjusted to be in accordance with the reference genome.
    Mummerplot is used to generate 1) nucmer with maxmatch to show change in repetitive regions. 2) nucmer with mum to
    to distribution of the snps. The snps and coords are generated by nucmer -mum.
    The categorization of the snps:
    1) Phred Score from fastq file, filtered by phred_threshold
    2) Coding Region if ref_cds_gff is provided
    3) Homopolymers
    4) single snps: Sub(stitutions), Indel, microIndel are identified"""
    # Load files
    os.chdir(dirn)
    ref_names, ref_seqs = simple_fasta_load(ref)
    query_names, query_seqs = simple_fasta_load(query)
    if remove_polish_tag:
        query_names = [name.replace(remove_polish_tag, "") for name in query_names]
    cds_region = {}
    # Load coding region
    if ref_cds_gff:
        ref_gff_entries = simple_gff3_load(ref_cds_gff)
        # filter cds region
        cds_region = extract_gff3_region(ref_gff_entries, 2, 'CDS')
    # make directory to same results
    if not os.path.isdir(dirn + prefix):
        os.mkdir(dirn + prefix)
    if not os.path.isdir(dirn + prefix + "/snps"):
        os.mkdir(dirn + prefix + "/snps")
    if not os.path.isdir(dirn + prefix + "/png"):
        os.mkdir(dirn + prefix + "/png")
    if not os.path.isdir(dirn + prefix + "/coords"):
        os.mkdir(dirn + prefix + "/coords")
    # make tmp files in tmp directory
    if not os.path.isdir(dirn + "tmp/"):
        os.mkdir(dirn + "tmp")
    os.chdir(dirn + "tmp")
    if ref_strain:
        ref_names = [ref_strain + "-" + name for name in ref_names]
        if ref_cds_gff:
            for tig in ref_cds_gff.keys():
                ref_cds_gff[ref_strain + "-" + tig] = ref_cds_gff.pop(tig)
    if query_strain:
        query_names = [query_strain + "-" + name for name in query_names]
    query_fa = query + ".fa"
    simple_fasta_write(query_fa, query_names, query_seqs)
    simple_fasta_write(ref, ref_names, ref_seqs)
    # Make nucmer comparison for chromosome name adjustment
    os.system("nucmer -p strain %s %s" %(ref, query_fa))
    os.system("show-coords -clqT strain.delta > strain.coords")
    # Load and analyze the coords file for contig assignemnt
    contig_assign_list = load_coords("strain.coords")
    query_ref_tig = {}
    query_ref_coverage = {}
    ref_query_coverage = {}
    for tig, pos_list in contig_assign_list.iteritems():
        q_plus_cov = 0
        q_minus_cov = 0
        r_plus_cov = 0
        r_minus_cov = 0
        for pos in pos_list:
            q_coverage = pos[-3]
            r_coverage = pos[-2]
            direction = pos[-1]
            if direction == 1:
                q_plus_cov += q_coverage
                r_plus_cov += r_coverage
            else:
                q_minus_cov += q_coverage
                r_plus_cov += r_coverage
        query_ref_coverage.setdefault(tig[0], [])
        query_ref_coverage[tig[0]].append((q_plus_cov + q_minus_cov, tig[1], q_plus_cov, q_minus_cov))
        ref_query_coverage.setdefault(tig[1], [])
        ref_query_coverage[tig[1]].append((r_plus_cov + r_minus_cov, tig[0], r_plus_cov,r_minus_cov))
    for tig, pos_list in ref_query_coverage.items():
        pos_list.sort(reverse=True)
        ref_query_coverage[tig] = pos_list
    for tig, pos_list in query_ref_coverage.items():
        pos_list.sort(reverse=True)
        query_ref_coverage[tig] = pos_list

    # Write the coverage information
    with open('coverage.q.txt', 'w') as filep:
        lines = []
        lines.append("Query_tig\tRef_tig\tPos_cov\tNeg_cov\tTotal\n")
        for query_tig, cov_list in query_ref_coverage.iteritems():
            for q_cov in cov_list:
                line = "%s\t%s\t%2.2f\t%2.2f\t%2.2f\n" %(query_tig, q_cov[1], q_cov[2], q_cov[3], q_cov[0])
                lines.append(line)
        filep.write("".join(lines))
        _message("Query Coverage Information saved")
    with open('coverage.r.txt', 'w') as filep:
        lines = []
        lines.append("Ref_tig\tQuery_tig\tPos_cov\tNeg_cov\tTotal\n")
        for ref_tig, cov_list in ref_query_coverage.iteritems():
            for r_cov in cov_list:
                line = "%s\t%s\t%2.2f\t%2.2f\t%2.2f\n" %(ref_tig, r_cov[1], r_cov[2], r_cov[3], r_cov[0])
                lines.append(line)
        filep.write("".join(lines))
        _message("Ref Coverage Information saved")
    # Reciprocal Best Assignment
    un_assigned = 0
    for query_tig, query_cov_list in query_ref_coverage.iteritems():
        best_ref = query_cov_list[0][1]
        best_query = ref_query_coverage[best_ref][0][1]
        # Reciprocol Best
        if best_query == query_tig:
            query_ref_tig[query_tig] = best_ref
        else:
            best_choice = False
            # Not Reciprocol Best, the best one is overwhelming good
            if len(query_cov_list) == 1:
                best_choice = True
            else:
                best_ref_coverage = query_cov_list[0][0]
                second_ref_coverage = query_cov_list[1][0]
                best_choice = (best_ref_coverage > 70 ) and (second_ref_coverage < 30)
            if best_choice:
                best_ref_query_list = ref_query_coverage[best_ref]
                for id_, query_ in enumerate(best_ref_query_list, 1):
                    if query_[1] == query_tig:
                        break
                query_assign = "%s-%02d" %(best_ref, id_)
            else:
                un_assigned += 1
                if query_strain:
                    query_assign = "%s-Chr%03d" % (query_strain, un_assigned)
                else:
                    query_assign = "Chr%03d" % (un_assigned)
            query_ref_tig[query_tig] = query_assign

    # Dict = [ori_query_name : {assgined_name, 1=Same Direction/ -1=RC seq)}
    query_update_list = {}
    for q_tig, assigned_q_tig in query_ref_tig.iteritems():
        if assigned_q_tig.startswith(ref_strain + "-"):
            p = len(ref_strain + "-")
            if query_strain:
                assigned_q_tig = query_strain + "-" + assigned_q_tig[p:]
            else:
                assigned_q_tig = assigned_q_tig[p:]
        plus = query_ref_coverage[q_tig][0][-2]
        minus = query_ref_coverage[q_tig][0][-1]
        if plus > minus:
            direction = 1
        else:
            direction = -1
        query_update_list[q_tig] = (assigned_q_tig, direction)
    # Saved Assigned Information
    for q_name in query_names:
        if q_name not in query_ref_tig:
            un_assigned += 1
            if query_strain:
                query_assign = "%s-Chr%03d" % (query_strain, un_assigned)
            else:
                query_assign = "Chr%03d" % (un_assigned)
            query_ref_tig[q_name] = query_assign
            query_update_list[q_name] = (query_assign, 1)
    with open("query_tig_assign.txt", 'w') as filep:
        lines = ["Query_tig\tAssigned_tig\tIs_Assigned2Ref\tQuery_Cov\tRef_Coverage\n"]
        for q_tig, r_tig in query_ref_tig.iteritems():

            if re.match(r'\d\d\d', r_tig[-3:]):
                lines.append("%s\t%s\tUnassigned\n" %(q_tig, r_tig))
            elif re.match(r'-\d\d', r_tig[-3:]):
                r_tig = r_tig[:-3]
                q_coverage = query_ref_coverage[q_tig][0]
                r_coverage = ref_query_coverage[r_tig][0]
                r_tig = query_update_list[q_tig][0]
                lines.append("%s\t%s\tAssigned\t%2.2f\t%2.2f\n" %(q_tig, r_tig, q_coverage[0], r_coverage[0]))
            else:
                q_coverage = query_ref_coverage[q_tig][0]
                r_coverage = ref_query_coverage[r_tig][0]
                r_tig = query_update_list[q_tig][0]
                lines.append("%s\t%s\tAssigned\t%2.2f\t%2.2f\n" %(q_tig, r_tig, q_coverage[0], r_coverage[0]))
        filep.write("".join(lines))
        _message('Contig Assignment information stored')
    # Make Updated Query Fasta and Fastq file

    # Write the query fasta file
    query_names_adj = [query_update_list[query_name][0] for query_name in query_names]
    print query_names_adj
    for id_, seq in enumerate(query_seqs):
        direction = query_update_list[query_names[id_]][1]
        if direction == -1:
            seq = rc_seq(seq)
            query_seqs[id_] = seq

    q_file_seperate = query.split('.')
    q_name_major = ".".join(q_file_seperate[:-1])
    query_fa = q_name_major + ".adj.fa"
    simple_fasta_write(query_fa, query_names_adj, query_seqs)
    query = q_name_major + ".adj.fq"
    simple_fasta_write(query, query_names_adj, query_seqs)
    _message("Adjusted Query Fasta and Fastq file are generated")
    os.system("cp %s ../%s/" % (query, prefix))
    os.system("cp %s ../%s/" % (query_fa, prefix))

    # Make Paired Comparison by nucmer
    # Make the map of assigned query to ref
    query_ref_map = {}
    for q_tig, r_tig in query_ref_tig.iteritems():
        q_assigned = query_update_list[q_tig][0]
        if re.match(r"\d\d\d", q_assigned[-3:]):
            type_ = 2
        elif re.match(r"-\d\d", q_assigned[-3:]):
            if r_tig.startswith(query_strain + "-"):
                p = len(query_strain + "-")
                if ref_strain:
                    r_tig = ref_strain + "-" + r_tig[p:][:-3]
                else:
                    r_tig = r_tig[p:]
            type_ = 1
        else:
            if ref_strain:
                r_tig = ref_strain + "-" + r_tig[p:]
            else:
                r_tig = r_tig[p:]
            type_ = 0
        query_ref_map[q_assigned] = (r_tig, type_)
    pair_list, coords, snps, figures = assigned_nucmer_compare(ref, query_fa, query_ref_map)
    os.system("mv *.snps ../%s/snps" % prefix)
    os.system("mv *.coords ../%s/coords" % prefix)
    os.system("mv *.png ../%s/png" % prefix)
    os.system("mv *.txt ../%s/" % prefix)



def batch_assign():
    ref = "cbs-ref.fasta"
    # query_ls = ["51-con.fq", "52-con.fq", "53-con.fq", "54-con.fq", "580.con.fq", "582.con.fq", "583.con.fq", "585.con.fq"]
    # pre_ls = ["51o", "52o", "53o", "54o", "51n", "52n", "53n", "54n"]
    query_ls = ["cbs-2-con.fq"]
    pre_ls = ["cbs-canu"]
    for id_, query in enumerate(query_ls):
        pre = pre_ls[id_]
        align_contig_to_ref(ref=ref, query=query, prefix=pre, dirn="/home/zhuwei/cglabarata/anno/canu_ref_cbs/canu/",
                            query_strain=pre, ref_strain="cbs138", ref_cds_gff="", phred_threshold=40, homomer=6,
                            remove_polish_tag="")

def filter_cbs_snp_fosmid():
    dirn = "/home/zhuwei/cglabarata/anno/canu_ref_cbs/canu/cbs-canu/"
    snp = "out.snp.csv"
    os.chdir(dirn)

    cbs_region = [
        [20487, 459937],
        [13803, 483636],
        [30470, 540559],
        [24693, 631859],
        [35344, 650244],
        [34830, 901377],
        [15591, 971009],
        [20281, 1031248],
        [27900, 1081474],
        [22048, 1154829],
        [21762, 1280027],
        [43659, 1409236],
        [27473, 1385304]]


    filter_lines = []

    with open(snp, 'r') as filep:
        for lin, line in enumerate(filep):
            if lin == 0:
                filter_lines.append(line)
                continue
            if not line.strip():
                continue
            pos = line.find("cbs138-Chr")
            chr_ = line[pos + 10]
            chr_id = ord(chr_) - 65
            region = cbs_region[chr_id]
            entry = line.split(',')
            ref_pos_1 = int(entry[-6])
            ref_pos_2 = int(entry[-5])
            if ref_pos_1 <region[0] or ref_pos_2 > region[1]:
                continue
            filter_lines.append(line)

    with open("out.snp.filtered.csv", 'w') as filep:
        filep.write("".join(filter_lines))

def filter_cbs_canu_fosmid(fname, dirn):
    """Filter the fosmid region of the assembled cbs sequence by canu."""
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fname)
    fos_names = []
    fos_seqs = []
    # Fosmid region by nucmer alignment result
    fosmid_region = {
        "cbs-canu-ChrA": (35134, 474584, "A"),
        "cbs-canu-ChrB": (16812, 486647, "B"),
        "cbs-canu-ChrC-02": (38194, None, "C"),
        "cbs-canu-ChrC": (None, 472588, "C"),
        "cbs-canu-ChrD": (34422, 641586, "D"),
        "cbs-canu-ChrE": (50253, 667551, "E"),
        "cbs-canu-ChrF": (57905, 924449, "F"),
        "cbs-canu-ChrG": (29183, 984590, "G"),
        "cbs-canu-ChrH": (23241, 1035470, "H"),
        "cbs-canu-ChrI": (29585, 1115489, "I"),
        "cbs-canu-ChrJ": (38521, 1199488, "J"),
        "cbs-canu-ChrK": (23488, 1281747, "K"),
        "cbs-canu-ChrL": (51341, 1422011, "L"),
        "cbs-canu-ChrM": (37270, 1395090, "M")
    }
    for id_, name in enumerate(names):
        if name not in fosmid_region:
            continue
        pos = fosmid_region[name]
        pos_l, pos_r, fos_chr = pos
        if pos_l:
            f_name = "CBS_" + fos_chr + "_L"
            f_seq = seqs[id_][:pos_l]
            fos_names.append(f_name)
            fos_seqs.append(rc_seq(f_seq))
        if pos_r:
            f_name = "CBS_" + fos_chr + "_R"
            f_seq = seqs[id_][pos_r:]
            fos_names.append(f_name)
            fos_seqs.append(f_seq)
    simple_fasta_write("cbs.canu.fosmid.fa", fos_names, fos_seqs)


def comp_fos_in_pair():
    """Compare the canu assembled fosmid to the paired fosmid."""
    dirn = "/home/zhuwei/cglabarata/cbs/canu.v2/fos/"
    if not os.path.isdir(dirn + "fos_figure/"):
        os.mkdir(dirn + "fos_figure")
    if not os.path.isdir(dirn + "fos_coords/"):
        os.mkdir(dirn + "fos_coords")
    if not os.path.isdir(dirn + "single_fos"):
        os.mkdir(dirn + "single_fos/")
    os.chdir(dirn)
    fos_names, fos_seqs = simple_fasta_load("cbs-fosmid.fasta")
    canu_fos_names, canu_fos_seqs = simple_fasta_load("cbs.canu.fosmid.fa")

    # compair by nucmer in pairs
    os.chdir(dirn + "single_fos")
    for id_, cf_name in enumerate(canu_fos_names):
        cf_seq = canu_fos_seqs[id_]
        cf_fn = cf_name + ".canu.fa"
        fos_chr = cf_name[-3:]
        fos_name = ""
        fos_seq = ""
        f_fn = ""
        for id_f, f_name in enumerate(fos_names):
            f_chr = f_name[-3:]
            if f_chr == fos_chr:
                f_fn = f_name + ".fa"
                fos_name = f_name
                fos_seq = fos_seqs[id_f]
                break
        simple_fasta_write(cf_fn, [cf_name], [cf_seq])
        simple_fasta_write(f_fn, [fos_name], [fos_seq])
        # os.system("nucmer --noextend -b 1 -g 0 -c 20 -p %s -maxmatch %s %s" %(fos_chr, f_fn, cf_fn))
        # os.system("show-coords -clrT %s.delta > %s.coords" %(fos_chr, fos_chr))
        # os.system("mummerplot -p %s -png -medium -color %s.delta" %(fos_chr, fos_chr))
        os.system("mummer -maxmatch -b %s %s > %s.mums" %(f_fn, cf_fn, fos_chr))
        os.system('''mummerplot --color -p %s --xrange "[0,%d]" --yrange "[0,%d]"  -png --medium %s.mums'''
                  %(fos_chr, len(fos_seq), len(cf_seq), fos_chr))
    os.system("mv *.coords ../fos_coords/")
    os.system("mv *.png ../fos_figure/ ")

def score_to_value(score):
    """Return the numeric value of a quality score."""
    return ord(score) - 33

def _mpileup_match_string_to_list(match_string):
    """Translate the mpileup string to readable list."""
    special_match = None
    match_list = []
    mis_match = ["A", "T", "C", "G", "a", "t", "c", "g"]

    str_len = len(match_string)

    pos = 0

    while pos < str_len:
        ch_ = match_string[pos]
        if ch_ == "." or ch_ == ",":
            match_list.append("M")
            pos += 1
            continue
        elif ch_ == "*":
            match_list.append("Del")
            pos += 1
            continue
        elif ch_ == "^":
            pos += 2
            continue
        elif ch_ == "$":
            pos += 1
            continue
        elif ch_ == "+":
            pos += 1
            num_match = re.match(r"(\d+)", match_string[pos: pos+5])
            num = int(num_match.group(1))
            if num < 10:
                ins = 1
            elif num <100:
                ins = 2
            elif num < 1000:
                ins = 3
            else:
                ins = 4
            ch_ = match_string[pos + ins: pos + num + ins]
            match_list[-1] = match_list[-1] + ", +%d%s" % (num, ch_)
            pos += num + ins
            continue
        elif ch_ == "-":
            pos += 1
            num_match = re.match(r"(\d+)", match_string[pos: pos+5])
            num = int(num_match.group(1))
            if num < 10:
                ins = 1
            elif num <100:
                ins = 2
            elif num < 1000:
                ins = 3
            else:
                ins = 4
            ch_ = match_string[pos + ins: pos + num + ins]
            match_list[-1] = match_list[-1] + ", -%d%s" % (num, ch_)
            pos += num + ins
            continue
        elif ch_ in mis_match:
            ch_ = ch_.upper()
            match_list.append(ch_)
            pos += 1
            continue
        elif ch_ == "<":
            ch_ = match_string[pos: pos + 3]
            if ch_ == "</>":
                match_list.append("Skip")
                pos += 3
                continue
            else:
                _message("Error Loading %s at %d match" %(match_string, pos))
                return []
        else:
            _message("Error Loading %s at %d match" %(match_string, pos))
            return []
    return match_list


def _load_mpileup_entry(entry):
    """Load the entry of mpileup file for filter_mpileup_by_fos function.

    Returns:
    entry = [tig_name, pos, match_list]
    match_list = [[match_type, [bq_list], [mq_list]]]
    """

    # Change Match String to readable list
    match_string = entry[4]
    match_list = _mpileup_match_string_to_list(match_string)

    read_count = int(entry[3])
    if len(match_list) != read_count:
        _message("Error reading line: %s.  %d Matches analyzed, %d annotated"
            % ("\t".join(entry), len(match_list), read_count))
        print match_list
        return []

    tig = entry[0]
    pos = int(entry[1])
    ref_nt = entry[2]
    count = int(entry[3])
    bq_list = [score_to_value(sco) for sco in entry[5]]
    mq_list = [score_to_value(sco) for sco in entry[6]]

    # Merge same matches
    match = zip(match_list, bq_list, mq_list)
    match.sort()
    match_merge = []
    cur_match = None
    cur_bq = []
    cur_mq = []
    for mat in match:
        if cur_match:
            if mat[0] == cur_match:
                cur_bq.append(mat[1])
                cur_mq.append(mat[2])
            else:
                match_merge.append([cur_match, cur_bq, cur_mq])
                cur_match = mat[0]
                cur_bq = [mat[1]]
                cur_mq = [mat[2]]
        else:
            cur_match = mat[0]
            cur_bq = [mat[1]]
            cur_mq = [mat[2]]
    match_merge.append([cur_match, cur_bq, cur_mq])
    entry = [tig, pos, ref_nt, count]
    entry.append(match_merge)
    return entry

def analyze_mpileup(fname, dirn, prefix="out", filter_bq=0, filter_mq=0):
    """Anaylze the mpileup file
    Args:
    fname: name of the mpileup file by samtools mpileup -f $REF -s $BAMFILE > $fname
    poslist: ('Chr_name': (pos1, pos2, ...), )
    prefix: Output file %prefix.filter.txt %prefix.filter.bq.histo %prefix.filter.mq.histo
    filter_bq: quality score to filter the read
    filter_mapq: quality score to filter the read
    Returns:
    Void
    The output file is generated.
    Format of the output:
    $prefix.parsedpileup.txt
    $Chr_Name\t$pos\t$Ref_nt\t$Read_Count\t$Type1:count;$Type2:count;$Type3:count\t$Av_Type1_BQ;...\t$Av_Type1_MQ;...\t$Type1_BQ1,BG2...;...\t$Type1_MQ1,MQ2,...
    $prefix.filter.txt_f

    """
    os.chdir(dirn)
    m_entries = []
    with open(fname, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            if len(entry) != 7:
                continue
            m_entry = _load_mpileup_entry(entry)
            m_entries.append(m_entry)

    with open(prefix + ".parsedpileup.txt", "w") as filep:
        lines = []
        for entry in m_entries:
            entry_h = entry[:3]
            entry_m = entry[4]
            count = 0
            filter_match = []
            total_bq = 0
            total_mq = 0
            for match in entry_m:
                match_bq = match[1]
                match_mq = match[2]
                match_q = zip(match_bq, match_mq)
                match_q = [score for score in match_q if score[0] > filter_bq and score[1] > filter_mq]
                if not match_q:
                    continue
                match_bq, match_mq = zip(*match_q)

                count += len(match_bq)
                total_bq += sum(match_bq)
                total_mq += sum(match_mq)
                filter_match.append([match[0], len(match_bq), match_bq, match_mq])
            av_bq = float(total_bq) / count
            av_mq = float(total_mq) / count
            # Match Type format
            # $Header\tAll_Filtered_Reads\t%Read_Count\tAve_BQ\tAve_MQ\n
            # $Header\t$Type\t$Count\t$Ave_BQ\t$Ave_MQ\t$BQ1,$BQ2..\t$MQ1, $MQ2..
            line = "%s\t%d\t%s\tAll_Filtered_Reads\t%d\t%2.2f\t%2.2f\n" %(
                tuple(entry_h) + (count, av_bq, av_mq))
            lines.append(line)
            for match in filter_match:
                m_count = match[1]
                m_av_bq = float(sum(match[2])) / m_count
                m_av_mq = float(sum(match[3])) / m_count
                bq_line = ", ".join(map(str, match[2]))
                mq_line = ", ".join(map(str, match[3]))
                line = "%s\t%d\t%s\t%s\t%d\t%2.2f\t%2.2f\t%s\t%s\n" %(
                    tuple(entry_h) + (match[0], m_count, m_av_bq, m_av_mq, bq_line, mq_line))
                lines.append(line)
        filep.write("".join(lines))
    # Generate Statistics of Matched reads
    with open(prefix + ".parsedpileup.stats", "w") as filep:
        lines = []
        for entry in m_entries:
            entry_h = entry[:3]
            entry_m = entry[4]
            count = 0
            filter_match = []
            total_bq = 0
            total_mq = 0
            for match in entry_m:
                match_bq = match[1]
                match_mq = match[2]
                match_q = zip(match_bq, match_mq)
                match_q = [score for score in match_q if score[0] > filter_bq and score[1] > filter_mq]
                if not match_q:
                    continue
                match_bq, match_mq = zip(*match_q)

                count += len(match_bq)
                total_bq += sum(match_bq)
                total_mq += sum(match_mq)
                filter_match.append([match[0], len(match_bq), match_bq, match_mq])
            av_bq = float(total_bq) / count
            av_mq = float(total_mq) / count
            for match in filter_match:
                if match[0] != "M":
                    continue
                m_count = match[1]
                m_av_bq = float(sum(match[2])) / m_count
                m_av_mq = float(sum(match[3])) / m_count
                m_coverage = float(m_count)/count * 100
                break

            line = "%s\t%d\t%s\t%d\t%d\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n" %(
                tuple(entry_h) + (count, m_count, m_coverage, av_bq, av_mq, m_av_bq, m_av_mq))
            lines.append(line)
        filep.write("".join(lines))





def filter_mpileup_by_pos(fname, dirn, pos_file, prefix="out"):
    """Export the filtered mpileupfile by provided positions
    Args:
    fname: name of the mpileup file by samtools mpileup -f $REF -s $BAMFILE > $fname
    poslist: ('Chr_name': (pos1, pos2, ...), )
    Returns:
    Void
    The filtered mpileup file is written in $prefix.filter.mpileup
    """
    os.chdir(dirn)
    filter_lines = []
    pos_dict = {}
    with open(pos_file, 'r') as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            tig = entry[0]
            pos = int(entry[1])
            pos_dict.setdefault(tig, [])
            pos_dict[tig].append(pos)


    with open(fname, 'r') as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            if len(entry) < 2:
                continue
            tig = entry[0]
            pos = int(entry[1])
            if pos not in pos_dict[tig]:
                continue
            filter_lines.append(line)
    with open(prefix + ".filter.mpileup", "w") as filep:
        filep.write("\n".join(filter_lines))

def filter_mpileup_by_range(fname, dirn, prefix="out"):
    """Export the filtered mpileupfile by predefined regions
    Args:
    fname: name of the mpileup file by samtools mpileup -f $REF -s $BAMFILE > $fname
    poslist: ('Chr_name': (pos1, pos2, ...), )
    Returns:
    Void
    The filtered mpileup file is written in $prefix.filter.mpileup
    """
    os.chdir(dirn)
    filter_lines = []

    non_fosmid_region = {
        "cbs-canu-ChrA": (35134, 474584, "A"),
        "cbs-canu-ChrB": (16812, 486647, "B"),
        "cbs-canu-ChrC-02": (38194, 12813068, "C"),
        "cbs-canu-ChrC": (1, 472588, "C"),
        "cbs-canu-ChrD": (34422, 641586, "D"),
        "cbs-canu-ChrE": (50253, 667551, "E"),
        "cbs-canu-ChrF": (57905, 924449, "F"),
        "cbs-canu-ChrG": (29183, 984590, "G"),
        "cbs-canu-ChrH": (23241, 1035470, "H"),
        "cbs-canu-ChrI": (29585, 1115489, "I"),
        "cbs-canu-ChrJ": (38521, 1199488, "J"),
        "cbs-canu-ChrK": (23488, 1281747, "K"),
        "cbs-canu-ChrL": (51341, 1422011, "L"),
        "cbs-canu-ChrM": (37270, 1395090, "M"),
        "cbs-canu-ChrL-02": (58766, 58766,"L02")
    }

    with open(fname, 'r') as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            if len(entry) < 2:
                continue
            tig = entry[0]
            pos = int(entry[1])
            if tig not in non_fosmid_region:
                continue
            non_region = non_fosmid_region[tig]
            non_l = non_region[0]
            non_r = non_region[1]
            if pos < non_l or pos > non_r:
                filter_lines.append(line)
    with open(prefix + ".filter.mpileup", "w") as filep:
        filep.write("\n".join(filter_lines))


def load_snp_position(fname, dirn, ref_filter=None, query_filter=None):
    """Export the position of all the snps.
    Args:
    fname: Name of the snp file
    dirn: working directory
    ref_filter: {ref_tig: (left, right)}, only the snps within the filtered region is loaded
    query_filter: {query_tig: (left, right)}, only the snps within the filteded region is loaded
    Returns:
    snp_position = [(query_tig, postition, ref_tig, position)]
    raw_snps = [lines of the snp file]
    """

    os.chdir(dirn)
    raw_snps = []
    snp_position = []
    with open(fname, "r") as filep:
        for lin, line in enumerate(filep):
            if lin < 4:
                continue
            line = line.strip()
            entry = line.split()
            ref_pos = int(entry[0])
            query_pos = int(entry[3])
            ref_tig = entry[-2]
            query_tig = entry[-1]
            if ref_filter:
                if ref_tig not in ref_filter:
                    continue
                else:
                    if ref_pos < ref_filter[ref_tig][0] or ref_pos > ref_filter[ref_tig][1]:
                        continue
            if query_filter:
                if query_tig not in query_filter:
                    continue
                else:
                    if query_pos < query_filter[query_tig][0] or query_pos > query_filter[query_tig][1]:
                        continue
            snp_position.append((query_tig, query_pos, ref_tig, ref_pos))
            if ref_filter or query_filter:
                raw_snps.append(line)
    return snp_position, raw_snps

def filter_cbs_canu_snps():
    """Filter the snps in the fosmid region of the cbs-canu assmembly compared with published CBS138"""
    non_fosmid_region = {
        "cbs-canu-ChrA": (35134, 474584, "A"),
        "cbs-canu-ChrB": (16812, 486647, "B"),
        "cbs-canu-ChrC-02": (38194, 12813068, "C"),
        "cbs-canu-ChrC": (1, 472588, "C"),
        "cbs-canu-ChrD": (34422, 641586, "D"),
        "cbs-canu-ChrE": (50253, 667551, "E"),
        "cbs-canu-ChrF": (57905, 924449, "F"),
        "cbs-canu-ChrG": (29183, 984590, "G"),
        "cbs-canu-ChrH": (23241, 1035470, "H"),
        "cbs-canu-ChrI": (29585, 1115489, "I"),
        "cbs-canu-ChrJ": (38521, 1199488, "J"),
        "cbs-canu-ChrK": (23488, 1281747, "K"),
        "cbs-canu-ChrL": (51341, 1422011, "L"),
        "cbs-canu-ChrM": (37270, 1395090, "M")
    }
    dirn = "/home/zhuwei/cglabarata/cbs/canu.v2/canu/cbs-canu/snps"
    os.chdir(dirn)
    snp_position = []
    snp_lines = []
    for fname in os.listdir(dirn):
        if fname.endswith("snps"):
            snps, lines = load_snp_position(fname=fname, dirn=dirn, query_filter=non_fosmid_region)
            if snps:
                snp_position.extend(snps)
                snp_lines.extend(lines)
    with open("non_fosmid.snps", "w") as filep:
        filep.write("\n".join(snp_lines))
        filep.write("\n")

    with open("nf_snp_pos.txt", "w") as filep:
        lines = []
        for snp in snp_position:
            lines.append("%s\t%d\t%s\t%d\n" %(snp[0], snp[1], snp[2], snp[3]))
        filep.write("".join(lines))



def comp_cbs_canu():
    """Compare Canu assembled sequecne with cbs138 published version."""
    ref = "cbs-ref.fasta"
    query = "cbs-2-con.adj.fq"
    pre = "cbs-canu"
    ref_cds_gff = "cbs138.gff"
    dirn = "/home/zhuwei/cglabarata/anno/canu_ref_cbs/canu/cbs-canu/"
    dirn_snp = dirn + "snps/"
    snp_files = []
    for file in os.listdir(dirn_snp):
        if file.endswith("snps"):
            snp_files.append(dirn_snp + file)
    compare_strain_by_tig_snp(ref=ref, query=query, dirn=dirn, snp_files=snp_files,
                              ref_cds_gff=ref_cds_gff, homo=6, x=10, prefix="out", ref_strain="cbs138")
    # align_contig_to_ref(ref=ref, query=query, prefix=pre, dirn="/home/zhuwei/cglabarata/anno/canu_ref_cbs/genome/",
    #                     query_strain="cbs-canu", ref_strain="cbs138", ref_cds_gff="", phred_threshold=40, homomer=6,
    #                     remove_polish_tag="")

def mummer_pair(dirn, min_repeat=3, merge_diff=200, min_length=20,
    min_cluster=100, overlap=500, ref_name="cbs.fos",
    query_name="cbs.canu", ref_tag="CBS_B", query_tag="canu.fa"):
    """Run mummer in pair to identify repetitive region."""
    os.chdir(dirn)

    if not os.path.isdir(dirn + "../fos_mums/"):
        os.mkdir(dirn + "../fos_mums/")
    if not os.path.isdir(dirn + "../fos_mums_png/"):
        os.mkdir(dirn + "../fos_mums_png/")

    fosmid_names = []
    canu_names = []

    for fname in os.listdir(dirn):
        if query_tag in fname:
            canu_names.append(fname)
        if ref_tag in fname:
            fosmid_names.append(fname)
        continue

    # run mummer to pair compare the two sequence
    mums_list = []
    for canu in canu_names:
        tig = canu.replace(query_tag, "")
        fos = ""
        for fosmid in fosmid_names:
            tig_f = fosmid.replace(ref_tag, "")
            if tig_f == tig:
                fos = fosmid
                break
        if fos:
            os.system("mummer -maxmatch -b -L -c %s %s > %s.mums" %(fos, canu, tig))
            os.system("mummerplot -png -color -medium -p %s %s.mums" %(tig, tig))
            mums_list.append("%s.mums" %tig)
        else:
            _message("NO counterpart for %s" %canu)

    os.system("mv *.mums ../fos_mums/")
    os.system("mv *.png ../fos_mums_png/")
    os.chdir(dirn + "../fos_mums")

    repeat_list = {}
    for mum in mums_list:
        print mum, repeat_list
        with open(mum, "r") as filep:
            mums = []
            tig = mum[:-5]
            for line in filep:
                line = line.strip()
                if not line:
                    continue
                if line[0] == ">":
                    continue
                entry = line.split()
                ref_s = int(entry[0])
                query_s = int(entry[1])
                len_ = int(entry[2])
                ref_e = ref_s + len_
                query_e = query_s + len_
                mums.append((ref_s, ref_e, query_s, query_e))
            ref_starts, ref_ends, query_starts, query_ends = zip(*mums)
            ref_s = min(ref_starts)
            ref_e = max(ref_ends)
            query_s = min(query_starts)
            query_e = max(query_ends)
            ref_starts = None
            ref_ends = None
            query_starts = None
            query_ends = None
            ref_ends = None
            ref_coords = np.zeros(ref_e - ref_s + 1, dtype=int)
            query_coords = np.zeros(query_e - query_s + 1, dtype=int)
            # count bins
            for match in mums:
                ref_coords[match[0] - ref_s: match[1] - ref_s +1] += 1
                query_coords[match[2] - query_s: match[3] - query_s +1] += 1
            # Export Projected match counts
            ref_regions = []
            left = None
            right = None
            cur_count = 0
            for pos, count in enumerate(np.nditer(ref_coords), start=ref_s):
                if count <= min_repeat:
                    if not left:
                        continue
                    elif right - left + 1 > min_length and cur_count > min_repeat:
                        ref_regions.append((left, right, cur_count))
                    left = None
                    right = None
                    cur_count = 0
                    continue
                if not left:
                    left = pos
                    right = pos
                    cur_count = float(count)
                else:
                    if cur_count > count - merge_diff -1 and cur_count < count + merge_diff + 1:
                        cur_count = (cur_count * (right - left +1) + count)/(pos - left +1)
                        right = pos
                    else:
                        if right - left + 1 > min_length:
                            ref_regions.append((left, right, cur_count))
                        left = pos
                        right = pos
                        cur_count = float(count)
            if left:
                if right - left + 1 > min_length:
                    ref_regions.append((left, right, cur_count))

            query_regions = []
            left = None
            right = None
            cur_count = 0
            for pos, count in enumerate(np.nditer(query_coords), start=ref_s):
                if count <= min_repeat:
                    if not left:
                        continue
                    elif right - left + 1 > min_length:
                        query_regions.append((left, right, cur_count))
                    left = None
                    right = None
                    cur_count = 0.
                    continue
                if not left:
                    left = pos
                    right = pos
                    cur_count = float(count)
                else:
                    if cur_count > count - merge_diff -1 and cur_count < count + merge_diff + 1:
                        cur_count = (cur_count * (right - left +1) + count)/(pos - left +1)
                        right = pos
                    else:
                        if right - left + 1 > min_length:
                            query_regions.append((left, right, cur_count))
                        left = pos
                        right = pos
                        cur_count = float(count)
            if left:
                if right - left + 1 > min_length:
                    query_regions.append((left, right, cur_count))
            ref_coords = None
            query_coords = None
            repeat_list[tig] = (ref_regions, query_regions)
    with open(ref_name + ".project.count", "w") as ref_p,  \
         open(query_name + ".project.count", "w") as query_p:
        ref_lines = []
        query_lines = []
        for tig, region in repeat_list.iteritems():
            ref_regions = region[0]
            query_regions = region[1]
            ref_tig = ref_name + "." + tig
            query_tig = query_name + "." + tig
            ref_line = ["%s\t%d\t%d\t%2.1f\n"
                %((ref_tig,) + ref_region) for ref_region in ref_regions]
            query_line = ["%s\t%d\t%d\t%2.1f\n"
                %((query_tig,) + query_region) for query_region in query_regions]
            ref_lines.extend(ref_line)
            query_lines.extend(query_line)
        ref_p.write("".join(ref_lines))
        query_p.write("".join(query_lines))

    # Merge regions
    ref_merge_regions = {}
    query_merge_regions = {}

    for tig, regions in repeat_list.iteritems():
        ref_regions = regions[0]
        query_regions = regions[1]
        ref_merge_regions[tig] = []
        query_merge_regions[tig] = []
        left = None
        right = None
        count = 0.
        for region in ref_regions:
            if not left:
                left = region[0]
                right = region[1]
                count = region[2]
                continue
            if region[0] - right < overlap and abs(region[2] - count) < merge_diff + 1:
                count = (count * (right - left + 1) + region[2] * (region[1] - region[0] + 1)) / \
                    (right - left + region[1] - region[0] + 2)
                right = region[1]
            else:
                if right - left + 1 > min_cluster:
                    ref_merge_regions[tig].append((left, right))
                left = region[0]
                right = region[1]
                count = region[2]
        if left and (right - left + 1 > min_cluster):
            ref_merge_regions[tig].append((left, right))
        left = None
        right = None
        count = 0.
        for region in query_regions:
            if not left:
                left = region[0]
                right = region[1]
                count = region[2]
                continue
            if region[0] - right < overlap and abs(region[2] - count) < merge_diff + 1:
                count = (count * (right - left + 1) + region[2] * (region[1] - region[0] + 1)) / \
                    (right - left + region[1] - region[0] + 2)
                right = region[1]
            else:
                if right - left + 1 > min_cluster:
                    query_merge_regions[tig].append((left, right))
                left = region[0]
                right = region[1]
                count = region[2]
        if left:
            query_merge_regions[tig].append((left, right))
    # Write the merge file
    with open(ref_name + ".merge.region", "w") as refp, \
         open(query_name + ".merge.region", "w") as queryp:
        ref_lines = []
        query_lines = []
        for tig, regions in ref_merge_regions.iteritems():
            ref_tig = ref_name + "." + tig
            ref_line = ["%s\t%d\t%d\n" %((ref_tig,) + region) for region in regions]
            ref_lines.extend(ref_line)
        for tig, regions in query_merge_regions.iteritems():
            query_tig = query_name + "." + tig
            query_line = ["%s\t%d\t%d\n" %((query_tig,) + region) for region in regions]
            query_lines.extend(query_line)
        refp.write("".join(ref_lines))
        queryp.write("".join(query_lines))

    # Extract Repetitive Region and Flanking region for read count extraction

def tmp_extract_telomeric_region(size=50000, telo_seq=r"GGGGTCTGGGTGCTG", save_telo=2, safe_extend=5000):
    """
    Extract the telomeric region of 51-5n.
    """
    dirn = "/home/zhuwei/cglabarata/58x/genome/"
    os.chdir(dirn)
    if not os.path.isdir(dirn + "fos"):
        os.mkdir(dirn + "fos")
    fnames = ["51.fa", "52.fa", "53.fa", "54.fa"]
    for fname in fnames:
        names, seqs = simple_fasta_load(fname)
        for id_, name in enumerate(names):
            seq = seqs[id_]
            if len(seq) < size * 2:
                continue
            chr_pos = name.find("Chr")
            fos_left = fname[:2] + "." + name[chr_pos:] + ".Left"
            fos_right = fname[:2] + "." + name[chr_pos:] + ".Right"

            seq_left = seq[:size + safe_extend].upper()
            seq_right = seq[-size - safe_extend:].upper()
            seq_left = rc_seq(seq_left)
            left_telo_matches = re.finditer(telo_seq, seq_left)
            right_telo_matches = re.finditer(telo_seq, seq_right)
            pos = -2
            for count, match in enumerate(left_telo_matches, 0):
                if count < save_telo:
                    pos = match.end()
                else:
                    break
            seq_left = seq_left[:pos]
            seq_left = seq_left[-size:]
            pos = -2
            for count, match in enumerate(right_telo_matches, 0):
                if count < save_telo:
                    pos = match.end()
                else:
                    break
            seq_right = seq_right[:pos]
            seq_right = seq_right[-size:]
            fname_left = dirn + "fos/" + fos_left + ".fa"
            fname_right = dirn + "fos/" + fos_right + ".fa"
            simple_fasta_write(fname_left, [fos_left], [seq_left])
            simple_fasta_write(fname_right, [fos_right], [seq_right])

def tmp_tidy_clustal(dirn, fname, prefix="out" , pattern=r"Chr\w\.\w+"):
    """Tidy the alignment of Clustal Omega"""
    os.chdir(dirn)
    alignment = {}
    group = {}
    difference = {}
    tig_pattern = re.compile(pattern)

    with open(fname, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            name = entry[0]
            align = entry[1]
            pos = int(entry[2])
            alignment.setdefault(name, [])
            alignment[name].append((align, pos))
            match = tig_pattern.search(name)
            if not match:
                continue
            tig_name = match.group(0)
            strain_name = name[: match.start()]
            group.setdefault(tig_name, {})
            group[tig_name].setdefault(strain_name, [])
            group[tig_name][strain_name].append((align, pos))

    for tig_name, tig_align in group.iteritems():
        difference.setdefault(tig_name, {})
        for strain in tig_align.keys():
            difference[tig_name].setdefault(strain, [])

    for tig_name, tig_align in group.iteritems():
        strain_names = tig_align.keys()
        tmp_strain = strain_names[0]
        len_align = len(tig_align[tmp_strain])
        for id_ in xrange(len_align):
            cur_align = ""
            for strain in strain_names:
                align = tig_align[strain][id_][0]
                if not cur_align:
                    cur_align = align
                else:
                    if align != cur_align:
                        for strain_diff in strain_names:
                            difference[tig_name][strain_diff].append(tig_align[strain_diff][id_])
                        break
    # Save grouped result
    with open(prefix + ".group.align", "w") as filep:
        lines = []
        for tig_name, tig_align in group.iteritems():
            strain_names = tig_align.keys()
            tmp_strain = strain_names[0]
            len_align = len(tig_align[tmp_strain])
            for id_ in xrange(len_align):
                line = ["%s\t%s\t%d\n" % (strain + tig_name, tig_align[strain][id_][0],
                        tig_align[strain][id_][1]) for strain in strain_names]
                line.append("\n")
                lines.extend(line)
        filep.write("".join(lines))

    # Save difference alignment

    with open(prefix + ".diff.align", "w") as filep:
        lines = []
        for tig_name, tig_align in difference.iteritems():
            strain_names = tig_align.keys()
            tmp_strain = strain_names[0]
            len_align = len(tig_align[tmp_strain])
            for id_ in xrange(len_align):
                line = ["%s\t%s\t%d\n" % (strain + tig_name, tig_align[strain][id_][0],
                        tig_align[strain][id_][1]) for strain in strain_names]
                line.append("\n")
                lines.extend(line)
        filep.write("".join(lines))

def filter_gff_region(dirn, fname, prefix="out", region_dict=None):
    """Filter the gff file by region."""
    pass

def extract_mRNA_gff(dirn, mrna, gff, prefix="out", flank=1000):
    """Extract the genomic mRNA sequence by name."""
    os.chdir(dirn)
    # rna_list = $Name: [contig, start, end, direction, sequence, sequcne + flanking]
    rna_list = {}
    entries, names, seqs = simple_gff3_load(gff, return_fasta=True)
    gene_entries = extract_gff3_region(entries=entries, column=2, tag="gene")
    rna_entries = extract_gff3_region(entries=entries, column=2, tag="mRNA")

    with open(mrna, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            rna_list.setdefault(line, None)

    rna_names = []
    rna_seqs = []
    rna_seqs_flanking = []
    for tig, rna_entry in rna_entries.iteritems():
        for rna in rna_entry:
            rna_name = rna[3]
            if rna_name in rna_list:
                start = rna[0]
                end = rna[1]
                direction = rna[2]
                name_id = names.index(tig)
                if direction == 1:
                    seq = seqs[name_id][start-1: end]
                    seq_flank = seqs[name_id][start - flank -1: end + flank]
                    rna_names.append(rna_name)
                else:
                    seq = rc_seq(seqs[name_id][start-1: end])
                    seq_flank = rc_seq(seqs[name_id][start - flank -1: end + flank])
                    rna_names.append(rna_name + "_rc")
                rna_seqs.append(seq)
                rna_seqs_flanking.append(seq_flank)
                rna_list[rna_name] = (tig, start, end, direction)
    for tig, gene_entry in gene_entries.iteritems():
        for gene in gene_entry:
            gene_name = gene[3]
            if gene_name in rna_list and not rna_list[gene_name] :
                start = gene[0]
                end = gene[1]
                direction = gene[2]
                name_id = names.index(tig)
                if direction == 1:
                    seq = seqs[name_id][start-1: end]
                    seq_flank = seqs[name_id][start - flank -1: end + flank]
                    rna_names.append(gene_name)
                else:
                    seq = rc_seq(seqs[name_id][start-1: end])
                    seq_flank = rc_seq(seqs[name_id][start - flank -1: end + flank])
                    rna_names.append(gene_name + "_rc")
                rna_seqs.append(seq)
                rna_seqs_flanking.append(seq_flank)
                rna_list[gene_name] = (tig, start, end, direction)
    # write result
    with open(prefix + ".list", "w") as filep:
        lines = ["Gene\tTig_Name\tStart\tEnd\n"]
        for rna, rna_info in rna_list.iteritems():
            if not rna_info:
                print rna
                continue
            if rna_info[3] == 1:
                lines.append("%s\t%s\t%d\t%d\n" % (rna, rna_info[0], rna_info[1], rna_info[2]))
            else:
                lines.append("%s\t%s\t%d\t%d\n" % (rna + "_rc", rna_info[0], rna_info[2], rna_info[1]))
        filep.write("".join(lines))

    simple_fasta_write(prefix + ".rna.fa", names=rna_names, seqs=rna_seqs)
    simple_fasta_write(prefix + ".rna.flank.fa", names=rna_names, seqs=rna_seqs_flanking)

def extract_single_fasta_by_name(dirn, fasta, name_file, prefix=""):
    """Extract the multifasta file into single fasta file from the list."""
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fasta)
    names = [name.split()[0] for name in names]
    name_filter = []
    with open(name_file, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            name_filter.append(line)

    fasta = zip(names, seqs)
    filter_fasta = [fas for fas in fasta if fas[0] in name_filter]
    print len([nam for nam in name_filter if nam not in names])
    for fas in filter_fasta:
        if not prefix:
            simple_fasta_write(fas[0] + ".fa", [fas[0]], [fas[1]])
        else:
            simple_fasta_write(prefix + "." + fas[0] + ".fa", [prefix + "." + fas[0]], [fas[1]])

def tmp_extract_cds_from_mrna_pos(dirn, protein, rna_pos, rna, prefix="out"):
    """Extract the cds sequence before the given mrna positions.
    If the region exceeds the first exon, only the first exon is returned.
    """
    os.chdir(dirn)
    prot_names, prot_seqs = simple_fasta_load(protein)
    rna_names, rna_seqs = simple_fasta_load(rna)

    rna_poslist = []
    with open(rna_pos, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            entry = line.split()
            rna_poslist.append((entry[0], int(entry[1])))

    prot_filter_seqs = []
    prot_filter_names = []
    for pos in rna_poslist:
        rna_name = pos[0]
        prot_name = pos[0].replace("_rc", "")
        rna_id = rna_names.index(rna_name)
        prot_id = prot_names.index(prot_name)
        rna_seq = rna_seqs[rna_id][: pos[1]/3*3]
        rna_trans = translate_exon(rna_seq)

        prot_seq = prot_seqs[prot_id]

        seq = common_start(seq_a=prot_seq, seq_b=rna_trans)
        prot_filter_seqs.append(seq)
        prot_filter_names.append(prot_name)

    simple_fasta_write(fname=prefix+".filter.fa", names=prot_filter_names, seqs=prot_filter_seqs)


def common_start(seq_a, seq_b):
    """Return the longest common shared sequence from beginning."""
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()
    def _iter():
        for nt_a, nt_b in zip(seq_a, seq_b):
            if nt_a == nt_b:
                yield nt_a
            else:
                return

    return "".join(_iter())


def tmp_dna_repeat_mummer(dirn, suffix="fa"):
    """Run mummer to identify the DNA repeat region.
    All the single fasta file in the dirn is under inverstigation.
    The repeat-match and nucmer -maxmatch -nosimplify are used to locate DNA repeats."""
    os.chdir(dirn)
    for fname in os.listdir(dirn):
        if not fname.endswith(suffix):
            continue
        pos = fname.index(suffix)
        name = fname[: pos-1]
        os.system("nucmer -maxmatch -nosimplify -p %s %s %s" %(name, fname, fname))
        os.system("mummerplot -medium -png -color -p %s %s.delta" %(name, name))
        os.system("repeat-match %s > %s.repeat" %(fname, name))

def tmp_filter_fasta_by_name(fasta, namelist, dirn, prefix="out"):
    """Filter the multifasta file by namelist."""
    os.chdir(dirn)
    names, seqs = simple_fasta_load(fasta)
    names = [name.split()[0] for name in names]
    fasta_ = zip(names, seqs)

    name_filter = []
    with open(namelist, "r") as filep:
        for line in filep:
            line = line.strip()
            if not line:
                continue
            name_filter.append(line)
    fasta_filter = [fas for fas in fasta_ if fas[0] in name_filter]
    names, seqs = zip(*fasta_filter)
    simple_fasta_write(prefix + ".filter.fa", names=names, seqs=seqs)

def extract_telomeric_region(fasta, dirn, prefix="out", size=50000, telo_seq=r"GGGGTCTGGGTGCTG", save_telo=2, safe_extend=5000):
    """
    Extract the telomeric region of 51-5n.
    """
    os.chdir(dirn)
    if not os.path.isdir(dirn + "fos"):
        os.mkdir(dirn + "fos")
    names, seqs = simple_fasta_load(fasta)
    for id_, name in enumerate(names):
        seq = seqs[id_]
        if len(seq) < size * 2:
            continue
        chr_pos = name.find("Chr")
        fos_left = prefix + "." + name[chr_pos:] + ".Left"
        fos_right = prefix + "." + name[chr_pos:] + ".Right"

        seq_left = seq[:size + safe_extend].upper()
        seq_right = seq[-size - safe_extend:].upper()
        seq_left = rc_seq(seq_left)
        left_telo_matches = re.finditer(telo_seq, seq_left)
        right_telo_matches = re.finditer(telo_seq, seq_right)
        pos = -2
        for count, match in enumerate(left_telo_matches, 0):
            if count < save_telo:
                pos = match.end()
            else:
                break
        seq_left = seq_left[:pos]
        seq_left = seq_left[-size:]
        pos = -2
        for count, match in enumerate(right_telo_matches, 0):
            if count < save_telo:
                pos = match.end()
            else:
                break
        seq_right = seq_right[:pos]
        seq_right = seq_right[-size:]
        fname_left = dirn + "fos/" + fos_left + ".fa"
        fname_right = dirn + "fos/" + fos_right + ".fa"
        simple_fasta_write(fname_left, [fos_left], [seq_left])
        simple_fasta_write(fname_right, [fos_right], [seq_right])





if __name__ == '__main__':
    # adjust_direction(fname="/home/zhuwei/cglabarata/comp/chrcomp/cbs-2-annotated.fa")
    # filter_cbs_fos()
    # extract_fos()
    # compare_strain("BG2", "CBS",
    #                dirn="/home/zhuwei/cglabarata/comp/kmer/linkage_map/",
    #                ref_file="bg2-fosmid.fa", query_file="cbs-fosmid.fa", pattern=r"_(\w)_(\w)")
    # one2one_nucmer()
    # maker_gff_extract("/home/zhuwei/cglabarata/anno/cbs-fosmid/cbs-fosmid-fos-augustus/cbs-fosmid.all.gff",
    #                   "/home/zhuwei/cglabarata/anno/cbs-fosmid/",
    #                   "cbs-fos")
    # kmer_fingerprint()
    # compare2fos()

    # exract_cds()
    # sort_write_fasta()
    # get_junction(dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2-comp/",
    #              fn="bg2.cbs.comp.txt",
    #              query_fn="bg2.canu.fa",
    #              prefix="bg2-cbs-flank")
    # pairwise_dnadiff()
    # extract_gpi()
    # compare_snp_cluster_5x_58x()
    # protein_self_align(dirn="/home/zhuwei/cglabarata/anno/canu-fosmid/gpi",
    #                    fname="gpi.fa", window=10,
    #                    matrix="BLOSUM62", prefix="gpi")
    # protein_self_dotmatcher(dirn="/home/zhuwei/cglabarata/anno/canu-fosmid/gpi/",
                            # fname="gpi.fa", prefix="gpi.dotmatcher")
    # load_blastp_result(fname="gpi.s288c.txt",
    #                    dirn="/home/zhuwei/cglabarata/anno/canu-fosmid/gpi/gpi.blast.protein/",
    #                    prefix="gpi.yeast")
    # extract_fasta_namelist(dirn="/home/zhuwei/cglabarata/anno/canu-fosmid/gpi/gpi.blast.protein",
    #                        fasta="cbs.cds.v2.fa",
    #                        namelist="gpi.v3.list",
    #                        prefix="gpi.v3")
    # split_fasta(maxcount=500, fname="cbs.cds.v2.fa", dirn="/home/zhuwei/cglabarata/anno/canu-fosmid/gpi/gpi.blast.protein", prefix="cbs")
    # auto_maker_annotation(query_prot="bg2.maker.fix.protein.fa",
    #                       ref_prot="cbs.cds.v2.fa",
    #                       query_pattern="",
    #                       ref_pattern=r"CAGL0(\w)(\d*)",
    #                       query_blast_file="canu.q2r.tab",
    #                       ref_blast_file="canu.r2q.tab",
    #                       query_gff_file="bg2.maker.fix.gff",
    #                       ref_gff_file="",
    #                       dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output/auto_anoo_v2/",
    #                       prefix="canu.v2",
    #                       cutoff=0.1)
    # maker_gff_extract(fname="maker.fix.gff",
    #                   dirn="/home/zhuwei/cglabarata/anno/cbs-canu/cbs-canu.maker.output/", prefix="out.fix")
    # fix_maker_blastx(
    #     maker_gff="bg2.canu.all.gff",
    #     maker_protein="bg2.canu.all.maker.proteins.fasta",
    #     dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output",
    #     prefix="maker")
    # fix_maker_blastx(
    #     maker_gff="cbs-canu.all.gff",
    #     maker_protein="cbs-canu.all.maker.proteins.fasta",
    #     dirn="/home/zhuwei/cglabarata/anno/cbs-canu/cbs-canu.maker.output",
    #     prefix="maker")

    # bridge_map(map1="bg2.all.id.map", map2="./auto_anno/canu.reci.pair.txt",
    #            dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output/",
    #            prefix="paired")
    # find_share_snp_5x_58x()
    # adjust_direction(fname="/home/zhuwei/cglabarata/anno/canu_ref_cbs/genome/cbs-canu.fa")
    # match_strain(ref="cbs-ref.fasta",
    #              query="cbs-canu-adj.fa",
    #              dirn="/home/zhuwei/cglabarata/anno/canu_ref_cbs/genome/")
    # batch_assign()
    # comp_cbs_canu()
    # filter_cbs_snp_fosmid()
    # filter_cbs_canu_fosmid(fname="cbs.canu.v2.fa", dirn="/home/zhuwei/cglabarata/cbs/canu.v2/fos/")
    # comp_fos_in_pair()
    # filter_cbs_canu_snps()
    # analyze_mpileup(
    #     fname="cbs.ref.nfsnp.mpileup",
    #     dirn="/home/zhuwei/cglabarata/cbs/canu.v2/mpileup",
    #     prefix="cbs.ref.nfsnp", filter_bq=0, filter_mq=0)
    # mummer_pair()
    # tmp_extract_telomeric_region()
    # compare_strain("51", "54",
    #                dirn="/home/zhuwei/cglabarata/58x/genome/fos/pairwise/",
    #                ref_file="51.fos.fa", query_file="54.fos.fa", pattern=r"Chr(\w).(\w)", prefix="51.54.v2")
    # fix_maker_blastx(
    #     maker_gff="cbs.canu.v2.all.gff",
    #     maker_protein="cbs.canu.v2.all.maker.proteins.fasta",
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/",
    #     prefix="canu.v2.fix")
    # auto_maker_annotation(query_prot="canu.v2.fix.protein.fa",
    #                       ref_prot="cbs.cds.v2.fa",
    #                       query_pattern="",
    #                       ref_pattern=r"CAGL0(\w)(\d*)",
    #                       query_blast_file="canu.v2.q2r.tab",
    #                       ref_blast_file="canu.v2.r2q.tab",
    #                       query_gff_file="canu.v2.fix.gff",
    #                       ref_gff_file="",
    #                       dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/",
    #                       prefix="canu.v2",
    #                       cutoff=0.1)
    # extract_mRNA_gff(
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation",
    #     mrna="canu.gpi.txt",
    #     gff="canu.v2.rename.v2.gff", prefix="out", flank=1000)
    # extract_single_fasta_by_name(
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/gpi/",
    #     fasta="gpi.rna.fa",
    #     name_file="canu.gpi.dna.txt",
    #     prefix="cbs.canu.mrna")
    # tmp_dna_repeat_mummer(
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/gpi/canu.mrna/",
    #     suffix="fa")
    # mummer_pair(
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/gpi/canu.mrna/",
    #     min_repeat=3, merge_diff=200, min_length=20,
    #     min_cluster=20, overlap=50, ref_name="cbs.canu.mrna.",
    #     query_name="cbs.canu.mrna.", ref_tag="cbs.canu.mrna.", query_tag="cbs.canu.mrna.")
    # tmp_extract_cds_from_mrna_pos(
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/gpi/",
    #     protein="cbs.canu.gpi.fa",
    #     rna_pos="canu.gpi.mrna.repeat.start",
    #     rna="gpi.rna.fa",
    #     prefix="gpi.nterm")
    # fix_maker_blastx(
    #     maker_gff="51.fos.all.gff",
    #     maker_protein="51.fos.all.maker.proteins.fasta",
    #     dirn="/home/zhuwei/cglabarata/anno/5x/fos/51/51.fos/",
    #     prefix="51")
    # auto_maker_annotation(query_prot="51.fix.protein.fa",
    #                       ref_prot="cbs.cds.v2.fa",
    #                       query_pattern="",
    #                       ref_pattern=r"CAGL0(\w)(\d*)",
    #                       query_blast_file="",
    #                       ref_blast_file="",
    #                       query_gff_file="51.fix.gff",
    #                       ref_gff_file="",
    #                       dirn="/home/zhuwei/cglabarata/anno/5x/fos/51/51.fos/prot_annotation/",
    #                       prefix="51.fos",
    #                       cutoff=0.1)
    # auto_maker_annotation(query_prot="bg2.canu.maker.cds.v3.fasta",
    #                   ref_prot="cbs.cds.v2.fa",
    #                   query_pattern="",
    #                   ref_pattern=r"CAGL0(\w)(\d*)",
    #                   query_blast_file="",
    #                   ref_blast_file="",
    #                   query_gff_file="bg2.gene.gff",
    #                   ref_gff_file="",
    #                   dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output/prot_annotation/",
    #                   prefix="bg2",
    #                   cutoff=0.1)
    # tmp_filter_fasta_by_name(
    #     fasta="cbs.canu.gpi.fa",
    #     namelist="gpi.nonrepeat.list",
    #     dirn="/home/zhuwei/cglabarata/anno/cbs.canu.v2/cbs.canu.v2/protein_annotation/gpi/",
    #     prefix="gpi.nonrepeat")

    # tmp_filter_fasta_by_name(
    #     fasta="bg2.canu.maker.cds.v3.fasta",
    #     namelist="bg2.gpi.over5kb.v1.list",
    #     dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output/prot_annotation/gpi/",
    #     prefix="bg2.gpi.over5kb")

    # extract_mRNA_gff(
    #     dirn="/home/zhuwei/cglabarata/anno/bg2-canu/bg2.canu.maker.output/prot_annotation/gpi/",
    #     mrna="bg2.gpi.over5kb.v1.list",
    #     gff="bg2.gene.gff", prefix="bg2", flank=1000)
    # extract_telomeric_region(
    #     fasta="bg2.canu.fa",
    #     dirn="/home/zhuwei/cglabarata/bg2/genome/",
    #     prefix="bg2")
    align_contig_to_ref_fasta(
        ref="bg2_fosmid.fa",
        query="bg2.fosmid.fa",
        prefix="bg2.canu",
        dirn="/home/zhuwei/cglabarata/bg2/genome/fos/",
        query_strain="bg2.canu",
        ref_strain="bg2.fos",
        ref_cds_gff="", phred_threshold=40, homomer=6,
        remove_polish_tag="")