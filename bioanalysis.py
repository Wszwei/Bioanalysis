# This code aims to analyze the DNA seuqence by GC content window,
# rePetitive sequence(6mers), and frequent RE sites.


import os
import numpy
from fasta import Fasta
from genebank import GeneBank
from sequence_analysis import Seq_Analyzer
# from csv import
# Generate Bioanalyzer Window
# from biowindow import Window


def minichunk_for_Twist():
    dirn = "/Users/Zhuwei/Google_Drive/Project Data/ORDERS/Twist/"
    stat_name = "minichunk_stat.csv"
    GC_stat = []
    rep_stat = []
    stat_fn = open(dirn + stat_name, "w")
    fasta = Fasta("/Users/Zhuwei/Documents/Order_Twist_140607.fasta")
    for _n in range(fasta.size):
        _GC_ls = Seq_Analyzer(fasta.Seqs[_n]).GC_window_analyzer(100)
        print"pass GC\n"
        _mean = numpy.mean(_GC_ls)
        _min = min(_GC_ls)
        _max = max(_GC_ls)
        _std = numpy.std(_GC_ls)
        GC_stat.extend([fasta.Names[_n], _min, _max, _mean, _std])
        _rep_ls = Seq_Analyzer(fasta.Seqs[_n]).Repetitive_Seq_analyzer(3, 4, 2)
        print "pass rep\n"
        rep_stat.append(fasta.Names[_n])
        rep_stat.extend(_rep_ls)
        rep_stat.append(";;")
    # for fn in os.listdir(dirn):
    # 	if fn[-5:].upper() != "FASTA" and fn[-2:].upper() != "FA":
    # 		continue
    # 	fasta=Fasta(dirn+fn)
    # 	for _n in range(fasta.size):
    # 		_GC_ls=Seq_Analyzer(fasta.Seqs[_n]).GC_window_analyzer(100)
    # 		print"pass GC\n"
    # 		_mean=numpy.mean(_GC_ls)
    # 		_min=min(_GC_ls)
    # 		_max=max(_GC_ls)
    # 		_std=numpy.std(_GC_ls)
    # 		GC_stat.extend([fasta.Names[_n],_min,_max,_mean,_std])
    # 		_rep_ls=Seq_Analyzer(fasta.Seqs[_n]).Repetitive_Seq_analyzer(3,5,2)
    # 		print "pass rep\n"
    # 		rep_stat.append(fasta.Names[_n])
    # 		rep_stat.extend(_rep_ls)
    # 		rep_stat.append(";;")
    for _x in range(len(GC_stat)):
        stat_fn.write(str(GC_stat[_x]) + ";")
        if _x % 5 == 4:
            stat_fn.write("\n")
    for _x in rep_stat:
        if _x != ";;":
            stat_fn.write(str(_x) + ";")
        else:
            stat_fn.write(";\n")


def dig_finder():
    dirn = "./150126_BaG_Seq/"
    vector_seq = Fasta("./pZX-vector.fa").Seqs[0]
    stat_name = "chunk_dig_stat.csv"
    stat_dir = "/Users/xuz02/Google_Drive/workspace/Python/150126_BaG_Seq/"
    vector_seq = vector_seq.upper()
    vector_size = len(vector_seq)
    min_frag = 100
    max_dig = 5
    min_dig = 2

    enzyme_dict = {"EcoRI": "GAATTC",
                   "NotI": "GCGGCCGC",
                   "BamHI": "GGATCC",
                   "HindIII": "AAGCTT",
                   "KpnI": "GGTACC",
                   "SacI": "GAGCTC",
                   "SalI": "GTCGAC",
                   "SpeI": "ACTAGT",
                   "NheI": "GCTAGC",
                   "AgeI": "ACCGGT",
                   "BsaI": "GGTCTC",
                   "EcoRV": "GATATC",
                   "NcoI": "CCATGG",
                   "AgeI": "ACCGGT",
                   "PstI": "CTGCAG",
                   "XbaI": "TCTAGA",
                   }
    enzyme_dict_II = {"BsaI": "GGTCTC",
                      "BsmBI": "CGTCTC"}
    stat_fp = open(stat_dir + stat_name, "w")
    stat_fp.close()


def sticky_ends():
    enzyme_dict = {"CACCTG": "AarI",
                   "GCCGAG": "NmeAIII",
                   "CGTCTC": "BsmBI",
                   "GGTCTC": "BsaI"}
    enzyme_site_dict = {"AarI": "5/9",
                        "NmeAIII": "21/19",
                        "BsmBI": "1/5",
                        "BsaI": "1/5"}
    dirn = "/Users/Zhuwei/Google_Drive/Project Data/ORDERS/ChrIV_minichunks/"
    stat = "sticky_end_stat.csv"
    stat_fp = open(dirn + stat, "w")
    for fn in os.listdir(dirn):
        if fn[-5:].upper() != "FASTA":
            continue
        fasta = Fasta(dirn + fn)
        _seq_L = fasta.Seqs[0][:50]
        _seq_R = fasta.Seqs[0][-50:]
        if _seq_L[:6] in enzyme_dict:
            _enzyme = enzyme_dict[_seq_L[:6]]
            _site = enzyme_site_dict[_enzyme].split("/")
            _site.sort()
            _site_L = int(_site[0])
            _site_R = int(_site[1])
            _sticky_end = []
            _sticky_end.append(_seq_L[6 + _site_L:_site_R + 6])
            _sticky_end.append(_seq_R[-6 - _site_R:-6 - _site_L])
            stat_fp.write(fn[:-6] + "," + _enzyme + ",")
            for _end in _sticky_end:
                stat_fp.write(_end + ",")
            stat_fp.write("\n")
        else:
            print "ERRO in " + fn[:-5] + "\n"


def length_reporter():
    dirn = "/Users/xuz02/Google_Drive/Project Data/ORDERS/"
    stat = "length_stat_0819.csv"
    stat_fp = open(dirn + stat, "w")
    for fn in os.listdir(dirn):
        if fn != "20150819_TWIST_LabOrder.txt":
            continue
        fasta = Fasta(dirn + fn)
        for n in range(len(fasta.Seqs)):
            _len = len(fasta.Seqs[n])
            stat_fp.write(fasta.Names[n] + "," + str(_len) + "," + "\n")
    stat_fp.close()


def GC_analyzer(vis_flag, Seq, supth, infth):
    # fn=""
    # fasta=Fasta(fn)
    size = 20
    if vis_flag:
        Seq_Analyzer(Seq).GC_window_analyzer_visual(size, supth, infth)
    _GC_ls, _tmp = Seq_Analyzer(Seq).GC_window_analyzer(size, supth, infth)
    print min(_GC_ls), max(_GC_ls), numpy.mean(_GC_ls), numpy.std(_GC_ls)


def GC_tmp():
    dirn = "/Users/xuz02/Downloads/"
    fn = "dra_mt.fa"
    # for fn in os.listdir(dirn):
    # 	if fn[-3:].upper() != "TXT" :
    # 		continue
    fasta = Fasta(dirn + fn)
    n = len(fasta.Seqs)
    for i in range(n):
        seq = fasta.Seqs[i]
        name = fasta.Names[i]
        GC_content, outlets = [], []
        GC_content, outlets = Seq_Analyzer(
            seq).GC_window_analyzer(20, 0.70, 0.40)
        if len(outlets) > 0:
            Seq_Analyzer(seq).GC_window_analyzer_visual(20, 0.70, 0.30, name)
    # GC_content, outlets=[],[]
    # GC_content,outlets=Seq_Analyzer(seq).GC_window_analyzer(100,1.00,0.20)
    # if len(outlets)>0:
    # print fn
    # 	Seq_Analyzer(seq).GC_window_analyzer_visual(100,1.0,0.20)


def batch_PCR_Primers_at_End():

    folder = "/workspace/Python/170116_SynIV/minichunks/"
    savefile = "/workspace/Python/170116_SynIV/primers.csv"

    L_res = 2
    R_res = 2
    Tm_ls = range(52,59)
    minLen = 20
    maxLen = 50
    primer_ls = []
    output = open(savefile, "w+")

    _len = len(Tm_ls)

    for fn in os.listdir(folder):
        if "fasta" in fn:
            fa = Fasta(folder+fn)
            length = len(fa.Seqs)
            for n in range(length):
                seq = fa.Seqs[n]
                name = fa.Names[n]
                F_primer, R_primer, Tm_bin = Seq_Analyzer(seq).\
                    Find_Primer_at_Ends(L_res, R_res, Tm_ls, minLen, maxLen)
                primer_ls.append([name, F_primer, R_primer, Tm_bin])

    for primer_sub_ls in primer_ls:
        for n in range(_len):
            if primer_sub_ls[3][n]:
                output.write(str(primer_sub_ls[0][:-1]) + ",")
                for i in [1, 2]:
                    for j in range(4):
                        output.write(str(primer_sub_ls[i][n][j]) + ",")
                output.write(str(primer_sub_ls[2][n][4]) + ",")
                output.write("\n")


def Chara_analyzer():
    fasta = Fasta(
        "/Users/Zhuwei/Google_Drive/workspace/Python/141008_JIngchuan/\
        S288C_reference_genome_R64-1-1_20110203\
        /S288C_reference_sequence_R64-1-1_20110203.fsa")
    output = "/Users/Zhuwei/Google_Drive/workspace/Python/141008_Jingchuan/\
chara_stat_4.csv"
    fp = open(output, 'w')
    N = len(fasta.Names)
    word = {	"C": ["C"],
             "G": ["G"],
             "A": ["A"],
             "T": ["T"],
             "CpG": ["CG"],
             "CpT": ["CT"],
             "ApG": ["AG"],
             "CpC": ["CC"],
             "GpG": ["GG"],
             "CpA": ["CA"],
             "TpG": ['TG'],
             "CCGG": ["CCGG"]}
    chara_sites = []
    for i in range(N):
        seq = Seq_Analyzer(fasta.Seqs[i])
        name = fasta.Names[i]
        for idx in word:
            word_lst = word[idx]
            count = 0
            for site in word_lst:
                site_lst = seq.Degen_searcher(site)
                count += len(site_lst)
                chara_sites.append([name, site, len(site_lst), site_lst])
            fp.write(name + "," + idx + "," + str(count) + "\n")
        # for lst in chara_sites:
        # 	for item in lst:
        # 		fp.write(str(item)+",")
        # 	fp.write("\n")


def minichunk_csPCR_primers_two_pairs(dirn, minichunk, chunk, Tm,
                                      ref="", Lim=150, pcr=300):
    """ Function to make junction primers for csPCR.

    NAME OF THE MINICHUNKS = CHUNKNAME.01-0X
    Two pairs of csPCR is generated
    dirn: working dir
    minichunk: name of the multi fasta file with all the minichunks
    ref: reference genome
    Tm : Expected Tm
    """
    # upper limit of the length of the PCR product at each site
    half_PCR_length = pcr
    m_fasta = Fasta(dirn+minichunk)
    ref_fasta = Fasta(dirn+ref)
    c_fasta = Fasta(dirn+chunk)

    output_fn = "csPCR_primers1.csv"
    output = open(dirn + output_fn, "w")

    first_minichunk = []
    second_minichunk = []
    third_minichunk = []

    num_minichunks = len(m_fasta.Seqs)
    for i in range(num_minichunks):
        print m_fasta.Names[i][-1]
        if m_fasta.Names[i][-1] == "1":
            first_minichunk.append([m_fasta.Seqs[i], m_fasta.Names[i]])
        elif m_fasta.Names[i][-1] == "2":
            second_minichunk.append([m_fasta.Seqs[i], m_fasta.Names[i]])
        elif m_fasta.Names[i][-1] == "3":
            third_minichunk.append([m_fasta.Seqs[i], m_fasta.Names[i]])
        else:
            continue
    if (len(first_minichunk) != len(second_minichunk) or
       len(first_minichunk) != len(third_minichunk)):
        print("Not enough minichunks to make 2 pairs of csPCR primers!\n")
        return

    # forward primers at minichunk #1
    # return arguments for Seq_Analyzer.Find_Forward_Primer = [seq, ]
    forward_12 = []
    for minichunk in first_minichunk:
        seq = minichunk[0][- half_PCR_length:]
        primer = []
        primer_name = minichunk[1][:-3] + "_F_junc12"
        primer.append(primer_name)
        primer.extend(Seq_Analyzer(seq).
                      Find_Forward_Primer(
                      Tm, Lim, crossref=ref_fasta.Seqs[0]))
        if len(primer) == 1:
            print ("%s FAILED!\n" % (primer_name))
            primer.extend(["", -1, -1, -1, -1, -1])
        forward_12.append(primer)
    reverse_12 = []
    forward_23 = []
    for minichunk in second_minichunk:
        seqF = minichunk[0][- half_PCR_length:]
        seqR = minichunk[0][:half_PCR_length]
        primerF, primerR = [], []
        primer_name_F = minichunk[1][:-3] + "_F_junc23"
        primer_name_R = minichunk[1][:-3] + "_R_junc12"
        primerF.append(primer_name_F)
        primerR.append(primer_name_R)
        primerF.extend(Seq_Analyzer(seqF).
                       Find_Forward_Primer(
                       Tm, Lim, crossref=ref_fasta.Seqs[0]))
        primerR.extend(Seq_Analyzer(seqR).
                       Find_Reverse_Primer(
                       Tm, Lim, crossref=ref_fasta.Seqs[0]))
        if len(primerF) == 1:
            print ("%s FAILED!\n" % (primer_name_F))
            primerF.extend(["", -1, -1, -1, -1, -1])
        if len(primerR) == 1:
            print ("%s FAILED!\n" % (primer_name_R))
            primerR.extend(["", -1, -1, -1, -1, -1])
        reverse_12.append(primerR)
        forward_23.append(primerF)
    reverse_23 = []
    for minichunk in third_minichunk:
        seq = minichunk[0][:half_PCR_length]
        primer = []
        primer_name = minichunk[1][:-3] + "_R_junc23"
        primer.append(primer_name)
        primer.extend(Seq_Analyzer(seq).
                      Find_Reverse_Primer(
                      Tm, Lim, crossref=ref_fasta.Seqs[0]))
        if len(primer) == 1:
            print ("%s FAILED!\n" % (primer_name))
            primer.extend(["", -1, -1, -1, -1, -1])
        reverse_23.append(primer)
    # Verify the PCR products
    size = len(forward_12)
    if size != len(c_fasta.Seqs):
        print("# of reference chunks(%d) != # of primer request (%d)"
              % (len(c_fasta.Seqs), size))
        return
    refseq = ref_fasta.Seqs[0]

    output.write("Name, Sequence, 5' at minichunk, length, GC ratio\
 , flag_GC, flag_multi, Tm, 5' at chunk, 5' at genome, PCR at chunk\
 , PCR at genome\n")
    for n in range(size):
        pf12 = forward_12[n]
        pr12 = reverse_12[n]
        pf23 = forward_23[n]
        pr23 = reverse_23[n]
        chunkseq = c_fasta.Seqs[n]
        loc_f12 = chunkseq.find(pf12[1]) + 1
        loc_f23 = chunkseq.find(pf23[1]) + 1

        loc_r12 = chunkseq.find(Seq_Analyzer(pr12[1]).rcSeq()) + len(pr12[1])
        loc_r23 = chunkseq.find(Seq_Analyzer(pr23[1]).rcSeq()) + len(pr23[1])

        gen_f12 = refseq.find(pf12[1]) + 1
        gen_f23 = refseq.find(pf23[1]) + 1
        gen_r12 = refseq.find(Seq_Analyzer(pr12[1]).rcSeq()) + len(pr12[1])
        gen_r23 = refseq.find(Seq_Analyzer(pr23[1]).rcSeq()) + len(pr23[1])

        len_pcr_12 = loc_r12 - loc_f12
        len_pcr_23 = loc_r23 - loc_f23

        len_gen_12 = gen_r12 - gen_f12
        len_gen_23 = gen_r23 - gen_f23

        pf12.extend([loc_f12, gen_f12, len_pcr_12, len_gen_12])
        pf23.extend([loc_f23, gen_f23, len_pcr_23, len_gen_23])

        pr12.extend([loc_r12, gen_r12, len_pcr_12, len_gen_12])
        pr23.extend([loc_r23, gen_r23, len_pcr_23, len_gen_23])

        primers = [pf12, pr12, pf23, pr23]

        for p in primers:
            for t in p:
                output.write(str(t)+" ,")
            output.write("\n")
    output.close()


def tmp():
    import random
    dirn = "/Users/xuz02/Google_Drive/workspace/Python/150223_Jef/"
    fa = Fasta(dirn + "synIV_noTer.fa")
    output1 = dirn + "synIV_100seqments.fasta"
    op1 = open(output1, "w")
    name = "yeast_chr04_3_66.seg"

    seq = fa.Seqs[0]

    m = len(seq) / 1800

    index = range(m)
    random.shuffle(index)
    index = index[:100]
    index.sort()
    for i in index:
        sub_seq = seq[i * 1800:i * 1800 + 1800]
        sub_name = ">" + name + "%03d" % (i + 1)
        op1.write(sub_name + "\n")
        op1.write(sub_seq + "\n")
    op1.close()


def linkerPCR_primers():
    """ Make Linker PCR Primers."""
    chunks = "yeast_chr01_chunks.FA"
    dirn = "/Users/xuz02/Google_Drive/workspace/Python/150218_Leslie/"
    bacbone = "pZX4_lin.fa"
    output = "liner_primers.csv"
    v_fasta = Fasta(dirn + bacbone)
    c_fasta = Fasta(dirn + chunks)
    o_fp = open(dirn + output, "w")
    left = 58
    right = 57
    reverse = "ggccggccccagcttttgttc"
    forward = "cggccggccctatagtgagtcg"
    o_fp.write("Name, Forward Primer, Reverse Primer\n")
    for n in range(len(c_fasta.Seqs)):
        f_primer = c_fasta.Seqs[n][-right:] + forward
        r_primer = Seq_Analyzer(c_fasta.Seqs[n][:left]).rcSeq() + reverse
        name = c_fasta.Names[n]
        fn = "pZX4_" + name[:19] + ".fasta"
        fp = open(dirn+fn, "w")
        fp.write(">%s\n" % name)
        fp.write(c_fasta.Seqs[n] + v_fasta.Seqs[0])
        fp.close()
        o_fp.write("%s, %s, %s\n" % (name[:19], f_primer, r_primer))
    o_fp.close()


def chr01_pcrtag_stat():
    """Stat the pcrtags over the minichunks."""

    import re
    dirn = "./160521_JL/"
    primers = "csPCR_primers.csv"
    miniC = "synI_mini_chunks.fasta"
    output = "chr01_pcrtag_stat.csv"
    miniC_f = Fasta(dirn+miniC)
    primers_fp = open(dirn+primers, "r")
    op = open(dirn + output, "w")

    primer_ls = []
    csv = primers_fp.read()
    primer_ls = re.split("\r|,", csv)
    size = len(primer_ls)/3
    stats = []
    for i in range(size):
        p_seq = primer_ls[3*i+1].lower()
        for j in range(len(miniC_f.Names)):
            if p_seq in miniC_f.Seqs[j]:
                _pos = miniC_f.Seqs[j].find(p_seq)
                stats.append([miniC_f.Names[j], primer_ls[3*i],
                              primer_ls[3*i+2], _pos])
    for item in stats:
        op.write("%s, %s, %s, %d\n" % (item[0], item[1], item[2], item[3]))
    op.close()
    miniC_f.close()


def forLeslie():
    output = "./151101_leslie/order_v2.fa"
    fop = open(output, "w")
    dirn = "./151101_leslie/Aro1_final/"
    for fn in os.listdir(dirn):
        if fn[-3:] == "ape":
            _gb = GeneBank(dirn+fn)
            _name = fn[:-4]
            _seq = _gb.seq
            fop.write(">%s\n%s\n" % (_name, _seq))
            fop.write("\n")
            _gb.close()
    fop.close()


def tilling_PCR_primers(dist=500, Tm=58, disp=200, length=200,
                        seq="", minLen=20, maxLen=30, name="", ref=""):
    """Search Tilling primers along the sequence.
    dist: Distance between Forward primers
    disp: Maximum displacement around the expected distance
    length: Length of the qpcr product
    seq: Sequence."""

    size = len(seq)
    seq_list = []
    f_primer_list = []
    r_primer_list = []
    # First segemnt
    seq_list.append(seq[: 2*disp + length])
    # 2nd - last segment
    if size > dist + 2*disp:
        n = size / dist
        # 2nd through 2nd from last
        for i in range(n-2):
            seq_list.append(seq[(i+1)*dist-disp: (i+1)*dist + disp + length])
        if size - dist*n > length + disp:
            seq_list.append(seq[n*dist-disp:])
    for n in range(len(seq_list)):
        seg = seq_list[n]
        f_primer, r_primer = search_tilling_primer(seg=seg, ref=ref,
                                                   Tm=Tm, minLen=minLen,
                                                   maxLen=maxLen, disp=disp,
                                                   name=name, seq=seq, n=n,
                                                   length=length)
        if f_primer == []:
            Tm_trial_ls = [Tm-1, Tm+1, Tm-2, Tm+2, Tm-3,
                           Tm+3, Tm-4, Tm+4, Tm-5, Tm+5]
            for t in Tm_trial_ls:
                f_primer, r_primer = search_tilling_primer(seg=seg, ref=ref,
                                                           Tm=t, minLen=minLen,
                                                           maxLen=maxLen,
                                                           disp=disp,
                                                           name=name, seq=seq,
                                                           n=n, length=length)
                if f_primer:
                    f_primer_list.append(f_primer)
                    r_primer_list.append(r_primer)
                    break
            if not f_primer:
                print("No primers Found for %s_%d\n" % (name, n))
                f_primer_name = "%s_%d_F" % (name, n)
                r_primer_name = "%s_%d_R" % (name, n)
                f_primer_list.append([f_primer_name, "", "", "", ""])
                r_primer_list.append([r_primer_name, "", "", "", ""])
                continue
        f_primer_list.append(f_primer)
        r_primer_list.append(r_primer)
    return f_primer_list, r_primer_list


def search_tilling_primer(seg, ref, Tm, minLen, maxLen,
                          disp, name, seq, n, length):
    f_primer, r_primer = Seq_Analyzer(seg).Find_Primer_at_Ends_Tm(L_res=2*disp,
                                                                  R_res=2*disp,
                                                                  Tm=Tm,
                                                                  minLen=minLen,
                                                                  maxLen=maxLen)
    score = 10.
    best_p = [-1, -1]
    for i in range(len(f_primer)):
        fp = f_primer[i]
        if off_target_binding(seq=fp[0], ref=ref, rc=False):
            continue
        pos_sc = 2*((fp[2] - disp)/disp)**2
        GC_sc = (fp[3]-.5)**2
        _sc = pos_sc + GC_sc
        if _sc < score:
            score = _sc
            best_p[0] = i
    if best_p[0] == -1:
            print "%s_%d: Foward Primer Not FOUND at Tm=%d!\n" % (name, n, Tm)
            return [], []
    score = 10.
    fp = f_primer[i]
    for i in range(len(r_primer)):
        rp = r_primer[i]
        if off_target_binding(seq=rp[0], ref=ref, rc=True):
            continue
        if rp[2] < fp[2] + 100:
            # minimal product > 100 bp
            continue
        pos_sc = ((fp[2] - rp[2]-length)/length)**2
        GC_sc = (rp[3]-.5)**2
        _sc = pos_sc + GC_sc
        if _sc < score:
            score = _sc
            best_p[1] = i
    if best_p[1] == -1:
            print "%s_%d: Reverse Primer FOUND at Tm=%d!\n" % (name, n, Tm)
            return [], []
    else:
        fp = f_primer[best_p[0]]
        rp = r_primer[best_p[1]]
        f_primer_seq = fp[0]
        r_primer_seq = rp[0]
        f_primer_name = "%s_%d_F" % (name, n)
        print f_primer_name
        r_primer_name = "%s_%d_R" % (name, n)
        f_primer_pos = seq.find(f_primer_seq) + 1
        r_primer_pos = seq.find(Seq_Analyzer(r_primer_seq).rcSeq()) + 1
        f_primer_Tm = str(fp[1])
        r_primer_Tm = str(rp[1])
        f_primer_GC = str(fp[3])
        r_primer_GC = str(rp[3])
        f_primer = [f_primer_name, f_primer_seq,
                    f_primer_Tm, f_primer_pos, f_primer_GC]
        r_primer = [r_primer_name, r_primer_seq, r_primer_Tm,
                    r_primer_pos, r_primer_GC]
    return f_primer, r_primer


def off_target_binding(seq, ref, rc=False):
    seq_r = Seq_Analyzer(seq).rcSeq()
    if rc:
        if seq in ref:
            return True
        pos = ref.find(seq_r)
        ref = ref[:pos] + ref[pos+len(seq):]
        if seq_r in ref:
            return True
    else:
        if seq_r in ref:
            return True
        pos = ref.find(seq)
        ref = ref[:pos] + ref[pos+len(seq):]
        if seq in ref:
            return True
    if easy_blast(seq=seq, ref=ref):
        return True
    if easy_blast(seq=seq_r, ref=ref):
        return True
    return False


def SW_binding(seq, ref, threshold=.80):
    """Smith-Waterman algorithm."""
    mismatch = -1
    indel = -2
    match = 2
    len_seq = len(seq)
    len_ref = len(ref)
    if len_seq == 0:
        return 0
    sc_matrix = [0]*((len_seq+1)*(len_ref+1))
    # tra_matrix = [[0, 0]]((len_seq+1)*(len_ref+1))
    identity_matrix = [0]*((len_seq+1)*(len_ref+1))
    for i in xrange(len_ref):
        for j in range(len_seq):
            if ref[i] == seq[j]:
                ma = sc_matrix[i*len_seq+j] + match
                ide = True
            else:
                ma = sc_matrix[i*len_seq+j] + mismatch
                ide = False
            if sc_matrix[i*len_seq+j+1] > sc_matrix[(i+1)*len_seq+j]:
                ind = sc_matrix[i*len_seq+j+1] + indel
                ins = False
            else:
                ind = sc_matrix[(i+1)*len_seq+j] + indel
                ins = True
            if ma > ind:
                sc_matrix[(i+1)*len_seq+j+1] = ma
                if ide:
                    identity_matrix[(i+1)*len_seq+j+1] =\
                        identity_matrix[i*len_seq+j]+1
                else:
                    identity_matrix[(i+1)*len_seq+j+1] =\
                        identity_matrix[i*len_seq+j]
            else:
                sc_matrix[(i+1)*len_seq+j+1] = ind
                if ins:
                    identity_matrix[(i+1)*len_seq+j+1] =\
                        identity_matrix[(i+1)*len_seq+j]
                else:
                    identity_matrix[(i+1)*len_seq+j+1] =\
                        identity_matrix[i*len_seq+j+1]
    return identity_matrix[-1]

# def SW_binding_seq(seq, ref, threshold=.80):
#     """Smith-Waterman algorithm."""
#     """Alignment result is returned."""
#     mismatch = -1
#     indel = -2
#     match = 2
#     len_seq = len(seq)
#     len_ref = len(ref)
#     if len_seq == 0:
#         return 0
#     sc_matrix = [0]*((len_seq+1)*(len_ref+1))
#     # tra_matrix = [[0, 0]]((len_seq+1)*(len_ref+1))
#     identity_matrix = [0]*((len_seq+1)*(len_ref+1))
#     for i in xrange(len_ref):
#         for j in range(len_seq):
#             if ref[i] == seq[j]:
#                 ma = sc_matrix[i*len_seq+j] + match
#                 ide = True
#             else:
#                 ma = sc_matrix[i*len_seq+j] + mismatch
#                 ide = False
#             if sc_matrix[i*len_seq+j+1] > sc_matrix[(i+1)*len_seq+j]:
#                 ind = sc_matrix[i*len_seq+j+1] + indel
#                 ins = False
#             else:
#                 ind = sc_matrix[(i+1)*len_seq+j] + indel
#                 ins = True
#             if ma > ind:
#                 sc_matrix[(i+1)*len_seq+j+1] = ma
#                 if ide:
#                     identity_matrix[(i+1)*len_seq+j+1] =\
#                         identity_matrix[i*len_seq+j]+1
#                 else:
#                     identity_matrix[(i+1)*len_seq+j+1] =\
#                         identity_matrix[i*len_seq+j]
#             else:
#                 sc_matrix[(i+1)*len_seq+j+1] = ind
#                 if ins:
#                     identity_matrix[(i+1)*len_seq+j+1] =\
#                         identity_matrix[(i+1)*len_seq+j]
#                 else:
#                     identity_matrix[(i+1)*len_seq+j+1] =\
#                         identity_matrix[i*len_seq+j+1]
#     return identity_matrix[-1]


def easy_blast(seq, ref, threshold=.80):
    """Simplified Blast.
       Seed length = 8 bp
    """
    len_seed = 8
    len_seq = len(seq)
    for i in range(len_seq-len_seed+1):
        seed = seq[i:i+len_seed]
        pos = -1
        while seed in ref:
            pos = ref.find(seed)
            if i > 0:
                seq_l = seq[:i:-1]
                ref_l = ref[pos-2*i:pos:-1]
                ide_l = SW_binding(ref=ref_l, seq=seq_l)
            else:
                ide_l = 0
            if i < len_seq-len_seed:
                seq_r = seq[i+len_seed::]
                ref_r = ref[pos+len_seed:pos+len_seed+2*len(seq_r)]
                ide_r = SW_binding(ref=ref_r, seq=seq_r)
            else:
                ide_r = 0
            ide = len_seed + ide_l + ide_r
            if ide > len_seq * threshold:
                return True
            ref = ref[pos+len_seed:]
    return False


def DL_search_primers():
    dirn = "./150314_DL_PIG/"
    fn1 = "pign_3kbFLK.fa"
    fn2 = "piga_3kbFLK.fa"
    fn3 = "pigl_3kbFLK.fa"
    fn4 = "pigk_3kbFLK.fa"
    ref = Fasta(dirn+fn1).Seqs[0] + Fasta(dirn+fn2).Seqs[0]\
        + Fasta(dirn+fn3).Seqs[0] + Fasta(dirn+fn4).Seqs[0]

    name = "PigL"
    output = "pigL_primers_v3.csv"
    fp = open(dirn+output, "w")
    fa = Fasta(dirn+fn3)
    f_primer, r_primer = tilling_PCR_primers(seq=fa.Seqs[0],
                                             name=name, ref=ref)
    for p in f_primer:
        print p
        fp.write("%s,%s,%s,%s,%s\n" % (p[0], p[1], p[2], p[3], p[4]))
    for p in r_primer:
        fp.write("%s,%s,%s,%s,%s\n" % (p[0], p[1], p[2], p[3], p[4]))
    fp.close()


def load_enzymes():
    """Load list of REs."""
    # re_ls = [["Name of RE", "RE site"]]
    # re_ls_II = [["Name of RE", "RE site"]]
    # re_ls_II for REs with non_reciprocal binding site.
    re_ls = [["EcoRI", "GAATTC"],
             ["NotI", "GCGGCCGC"],
             ["BamHI", "GGATCC"],
             ["HindIII", "AAGCTT"],
             ["KpnI", "GGTACC"],
             ["SacI", "GAGCTC"],
             ["SalI", "GTCGAC"],
             ["SpeI", "ACTAGT"],
             ["NheI", "GCTAGC"],
             ["AgeI", "ACCGGT"],
             ["BsaI", "GGTCTC"],
             ["EcoRV", "GATATC"],
             ["NcoI", "CCATGG"],
             ["AgeI", "ACCGGT"],
             ["PstI", "CTGCAG"],
             ["XbaI", "TCTAGA"]]

    re_ls_II = [["BsaI", "GGTCTC"],
                ["BsmBI", "CGTCTC"]]

    return re_ls, re_ls_II


def analyze_re_site(seq, re_ls, re_ls_II):
    """Export list of RE restriction sites.

    REQUIRE load_enzymes to provide RE list

    Currently only report the 1st binding site of re.
    Note that TypeIIs will have RE site out of the binding site.
    Returns list of (1) Name of REs
                    (2) Sites in acscending seq
    """

    re_name_ls = []
    re_site_ls = []

    for re in re_ls:
        re_name_ls.append(re[0])
        seq = seq+seq[:len(re[1])]  # For Circular plasmid
        seq = seq.upper()
        site_ls = [i for i in range(len(seq)) if seq.startswith(re[1], i)]
        re_site_ls.append(site_ls)

    for re in re_ls_II:
        re_name_ls.append(re[0])
        re_site_1 = re[1]
        re_site_2 = Seq_Analyzer(re_site_1).rcSeq()
        site_ls_1 = [i for i in range(len(seq))
                     if seq.startswith(re_site_1, i)]
        site_ls_2 = [i for i in range(len(seq))
                     if seq.startswith(re_site_2, i)]
        site_ls = merge_ascend_ls(site_ls_1, site_ls_2)
        re_site_ls.append(site_ls)
    return re_name_ls, re_site_ls


def merge_ascend_ls(ls1, ls2):
    """Merge two ascending lists."""
    merge_ls = []
    i = 0
    j = 0
    len_1 = len(ls1)
    len_2 = len(ls2)
    while i < len_1:
        if j == len_2:
            merge_ls.extend(ls1[i:])
            break
        if ls2[j] < ls1[i]:
            merge_ls.append(ls2[j])
            j = j+1
        elif ls2[j] == ls1[i]:
            merge_ls.append(ls2[j])
            j = j+1
            i = i+1
        else:
            merge_ls.append(ls1[i])
            i = i+1
    if j < len_2:
        merge_ls.extend(ls2[j:])
    return merge_ls


def get_band_size(re_site_ls, plasmid_size):
    """Return size of RE bands in ascending seq.
       For digestion of circular plasmids.
    """
    if not re_site_ls:
        return []
    re_site_ls.append(re_site_ls[0] + plasmid_size)
    band_size_ls = [re_site_ls[i+1] - re_site_ls[i]
                    for i in range(len(re_site_ls)-1)]
    band_size_ls.sort()
    return band_size_ls


def filter_re(re_site_ls, plasmid_size):
    """Filter the REs from the analyze_re_site()."""

    # Filter Rules
    # Minimum site of dig band:
    min_dig_size = 100
    # Minimum distance between distinguishable  bands:
    min_dig_dist = 50
    # Minimum numbers of dig bands:
    min_dig_num = 2

    f_re_site_ls = []

    for re_site in re_site_ls:
        band_size_ls = get_band_size(re_site, plasmid_size)
        f_band_size_ls = [band for band in band_size_ls if band > min_dig_size]
        f_band_size_ls.append(plasmid_size+min_dig_dist)
        f_band_size_ls = [f_band_size_ls[i]
                          for i in range(len(f_band_size_ls)-1)
                          if (f_band_size_ls[i+1] - f_band_size_ls[i])
                          > min_dig_dist]
        if len(f_band_size_ls) > min_dig_num:
            f_re_site_ls.append(f_band_size_ls)
        else:
            f_re_site_ls.append([])
    return f_re_site_ls


def select_restriction_enyzme(seq_ls):
    """Select RE for plasmid digestion test.
    """

    # Number fo top REs
    top_num = 3
    re_ls, re_ls_II = load_enzymes()
    re_size = len(re_ls) + len(re_ls_II)
    seq_ls_size = len(seq_ls)
    dig_band_ls = []
    re_ls_for_seqs = []  # a_i_j=Ture if seq[i] could be verified by re[j]]]
    re_ls_for_seqs = [True] * ((len(re_ls)+len(re_ls_II))*seq_ls_size)
    num_verified = []

    for i in range(seq_ls_size):
        seq = seq_ls[i]
        plasmid_size = len(seq)
        re_name_ls, re_site_ls = analyze_re_site(seq, re_ls, re_ls_II)
        f_re_site_ls = filter_re(re_site_ls, plasmid_size)
        dig_band_ls.append(f_re_site_ls)
        for j in range(len(re_site_ls)):
            re_ls_for_seqs[j*seq_ls_size+i] = (len(f_re_site_ls[j]) > 0)
    re_name_ls, re_site_ls = analyze_re_site(seq_ls[0], re_ls, re_ls_II)

    # Find TOP RE (and its partener)
    for i in range(re_size):
        _num = re_ls_for_seqs[i*seq_ls_size: (i+1)*seq_ls_size].count(True)
        num_verified.append(_num)
    tops = []
    for i in range(len(num_verified)):
        if len(tops) < top_num:
                tops.append(i)
        else:
            nums = [num_verified[j] for j in tops]
            if num_verified[i] > max(nums):
                index_min = nums.index(min(nums))
                tops.pop(index_min)
                print ("%d poped\n" % index_min)
                tops.append(i)
            elif num_verified[i] == max(nums):
                if min(nums) == seq_ls_size:
                    top_num += 1
                    tops.append(i)
                else:
                    index_min = nums.index(min(nums))
                    tops.pop(index_min)
                    tops.append(i)

    # Find parteners for top RE
    parteners = []
    for i in range(len(tops)):
        num = num_verified[tops[i]]
        if num < seq_ls_size:
            n_verified = [j for j in range(seq_ls_size)
                          if not re_ls_for_seqs[tops[i]*seq_ls_size+j]]
            for m in range(len(re_site_ls)):
                _ver = True
                for n in n_verified:
                    _ver = (_ver and re_ls_for_seqs[m*seq_ls_size+n])
                if _ver:
                    if m in tops:
                        break
                    else:
                        parteners.append(m)
                        break
    enzymes = []
    enzymes.extend(tops)
    enzymes.extend(parteners)
    n_verified_ls = range(seq_ls_size)
    for re in enzymes:
        verified = [j for j in range(seq_ls_size)
                    if re_ls_for_seqs[re*seq_ls_size+j]]
        for v in verified:
            if v in n_verified_ls:
                n_verified_ls.remove(v)
    return enzymes, dig_band_ls, re_name_ls, n_verified_ls


def test_select_enzyme():
    dirn = "/Users/xuz02/Google_Drive/workspace/Python/test_dig/"
    seq_name_ls = []
    seq_ls = []
    for fn in os.listdir(dirn):
        if fn[-5:].upper() == "FASTA":
            fn = dirn+fn
            fa = Fasta(fn)
            seq_ls.append(fa.Seqs[0])
            seq_name_ls.append(fa.Names[0])
    enzymes, dig_band_ls, re_name_ls, n_verified_ls =\
        select_restriction_enyzme(seq_ls)
    # Write output csv file
    # Headers
    print enzymes
    for i in range(len(enzymes)):
        ouput = dirn + "re_" + re_name_ls[enzymes[i]] + ".csv"
        fop = open(ouput, "w")
        for n in range(len(seq_name_ls)):
            fop.write("%s, " % seq_name_ls[n])
            bands = dig_band_ls[n][enzymes[i]]
            for band in bands:
                fop.write("%d, " % band)
            fop.write("\n")
        fop.close()
    output = dirn+"non_verified.csv"
    fop = open(output, "w")
    for i in n_verified_ls:
        fop.write(seq_name_ls[i]+"\n")
    fop.close()

def annotate_variants():
    # Annotate the variants detected by Pacbio
    # Variants are labeled in LOWERCASE
    import re

    dirn = "/workspace/Python/170105_pacbio_var/"
    out1 = open(dirn + "150813_14_tst_variant.csv", "w")
    out2 = open(dirn + "150813_14_tst_verified.csv", "w")
    fn = dirn + "seq_plate_14.fasta"
    fa = Fasta(fn)
    pattern = r"[atcg]+"  #Lowercase that grouped together
    for i in range(len(fa.Names)):
        seq = fa.Seqs[i]
        name = fa.Names[i]
        if re.search(pattern, seq):
            out1.write("%s\n" %name)
            for m in re.finditer(pattern, seq):
                var = "%d to %d, %s\n" % (m.start()+1, m.end()+1, m.group(0))
                out1.write(var)
        else:
            out2.write("%s\n" %name)

def break_string_by_len(string, length):
    #Break down string in lines by fixed length
    return (string[0+i:length+i] for i in range(0, len(string), length))

def annotate_variants_v2():
    #Annotate Variants from consensus and ref file.
        # Annotate the variants detected by Pacbio
    # Variants are labeled in LOWERCASE
    import re
    from Bio import pairwise2

    line_length = 80  #Length for alignment shown in one line

    dirn = "/Workplace/Python/161111_Leslie/"
    out1 = open(dirn + "ZX_variant_1t_v2.txt", "w")
    out2 = open(dirn + "ZX_verified_1t.txt", "w")
    out3 = open(dirn + "ZX_low_read_count_1t_v2.txt", "w")
    con_fn = dirn + "012716-consensus.fasta"
    con_fa = Fasta(con_fn)  #Consensus file
    ref_fn = dirn + "ref_012716.fa"
    ref_fa = Fasta(ref_fn)

    pattern = r"[atcg]+"  #Lowercase that grouped together

    for i in range(len(con_fa.Names)):
        con_seq = con_fa.Seqs[i]
        name = con_fa.Names[i]
        if con_fa.Names[i][:-7] not in ref_fa.Names:
            continue
        j = ref_fa.Names.index(con_fa.Names[i][:-7])
        ref_seq = ref_fa.Seqs[j]
        _len = max([len(con_seq), len(ref_seq)])
        print _len
        alignment = pairwise2.align.globalms(con_seq.upper(), ref_seq.upper(), 2, -1, -.5, -.1)
        score = int(alignment[0][2])
        print score
        if score == _len*2:
            if re.search(pattern, con_seq):
                out3.write("%s\n" %name)
                out3.write("No variants found\n")
                out3.write("sequence with low read counts are annotated\n")
                for m in re.finditer(pattern, con_seq):
                    var = "%d to %d, %s\n" % (m.start()+1, m.end()+1, m.group(0))
                    out3.write(var)
            else:
                out2.write("%s\n" %name) 
        else:
            align1 = alignment[0][0]
            align2 = alignment[0][1]
            score = alignment[0][2]
            result = "Score = %g\n" %score
            res_1 = list(break_string_by_len(align1, line_length))
            res_2 = list(break_string_by_len(align2, line_length))
            for k in range(len(res_1)):
                if res_1[k] == res_2[k]:
                    continue
                cur = ("%d to %d:\n" % (k*line_length+1, (k+1)*line_length))
                cur = cur + res_1[k] + "\n"
                cur = cur + "|" * len(res_1[k]) + "\n"
                cur = cur + res_2[k] + "\n"
                result = result + cur
            out1.write("%s\n" %name)
            out1.write(result)
            out1.write("\n")
        print "%s: done\n" %name

def compare_similar():
    #Compare (multi)fasta file which are similar
    #Fragmentation before alignment
    import re
    from Bio import pairwise2

    line_length = 80  #Length for alignment shown in one line
    frag_length = 5000 #Fragmentation before align
    overlap_legnth = 500 #overlaps to join fragmented alignments

    dirn = "/Workplace/Python/161111_Leslie/"
    out1 = open(dirn + "Leslie_variant_1t_v2.txt", "w")
    out2 = open(dirn + "Leslie_verified_1t.txt", "w")
    out3 = open(dirn + "Leslie_low_read_count_1t_v2.txt", "w")
    con_fn = dirn + "Leslie_lib1.fasta"
    con_fa = Fasta(con_fn)  #Consensus file
    ref_fn = dirn + "111016_ref.fa"
    ref_fa = Fasta(ref_fn)

    def recstr(s1, s2, answer=''):
        #Combine the alignment result together for overlap test
        if not s1:
            return answer+s2
        if not s2:
            return answer+s1
        return recstr(s1[1:], s2[1:], answer+s1[0]+s2[0])

    for i in range(len(con_fa.Names)):
        con_seq = con_fa.Seqs[i].upper()
        name = con_fa.Names[i]
        j = ref_fa.Names.index(con_fa.Names[i][:-7])
        ref_seq = ref_fa.Seqs[j].upper()
        _len = max([len(con_seq), len(ref_seq)])
        print _len
        #Fragmentation
        con_frags = list(con_seq[0+i:frag_length+i+overlap_legnth] for i in range(0, len(con_seq), frag_length))
        ref_frags = list(ref_seq[0+i:frag_length+i+overlap_legnth] for i in range(0, len(ref_seq), frag_length))
        if len(con_frags) != len(ref_frags):
            print "Length difference is too LARGE to compare! %s\n" %name
            continue
        _is_same = True
        align_c = ""
        align_ref = ""
        for m in range(len(con_frags)):
            align_f = pairwise2.align.globalms(con_frags[m], ref_frags[m], 2, -1, -.5, -.1, one_alignment_only=True)
            _length = len(con_frags[m])
            _score = int(align_f[0][2])
            align_1 = align_f[0][0]
            align_2 = align_f[0][1]
            if not align_c: #First frag:
                align_c = align_c + align_1
                align_ref = align_ref + align_2
            else:
                overlap_1 = recstr(align_c[-overlap_legnth:], align_ref[-overlap_legnth:])
                overlap_2 = recstr(align_1[:overlap_legnth], align_2[:overlap_legnth])
                max_match = ""
                for n in range(overlap_legnth):
                    match =  ""
                    for l in range(overlap_legnth):
                        if (n + l < overlap_legnth and overlap_1[n+l] == overlap_2[l]):
                            match += overlap_2[l]
                        else:
                            if (len(match) > len(max_match)) : max_match = match
                            match = ""
                p_aligned = overlap_1.index(max_match)
                p_new = overlap_2.index(max_match)
                if p_aligned % 2 == 0:
                    p_aligned = p_aligned/2
                else: p_aligned = (p_aligned+1)/2
                if p_new % 2 == 0:
                    p_new = p_new/2
                else: p_new = (p_new+1)/2
                align_c = align_c[:-overlap_legnth + p_aligned] + align_1[p_new:]
                align_ref = align_ref[:-overlap_legnth + p_aligned] + align_2[p_new:]
                print align_c[-overlap_legnth + p_aligned:-overlap_legnth + p_aligned + len(max_match)/2]
                print align_1[p_new: p_new+len(max_match)/2]
                print overlap_1[p_aligned*2:p_aligned*2+len(max_match)]
                print overlap_2[p_new*2:p_new*2+len(max_match)]
                

        if align_c == align_ref:
            out2.write("%s \n" %name)
        else:
            res_1 = list(break_string_by_len(align_c, line_length))
            res_2 = list(break_string_by_len(align_ref, line_length))
            result = ""
            for k in range(len(res_1)):
                if res_1[k] == res_2[k]:
                    continue
                cur = ("%d to %d:\n" % (k*line_length+1, (k+1)*line_length))
                cur = cur + res_1[k] + "\n"
                cur = cur + "|" * len(res_1[k]) + "\n"
                cur = cur + res_2[k] + "\n"
                result = result + cur
            out1.write("%s\n" %name)
            out1.write(result)
            out1.write("\n")
        print "%s: done\n" %name

def count_loxP():
    import re
    dirn = "/Workplace/Python/161122_loxp/"
    synIV = Fasta(dirn+"synIV.fa") 
    loxP_seq = "ATAACTTCGTATAATGTACATTATACGAAGTTAT"
    pattern1 = r"[G]ATAACTTCGTATAATGTACATTATACGAAGTTAT[^G]"
    pattern2 = r"[^G]ATAACTTCGTATAATGTACATTATACGAAGTTAT[G]"
    pattern3 = r"GATAACTTCGTATAATGTACATTATACGAAGTTATG"
    ref_seq = synIV.Seqs[0].upper()
    loci_1 = [m.start() for m in re.finditer(pattern1, ref_seq)]
    loci_2 = [m.start() for m in re.finditer(pattern2, ref_seq)]
    loci_3 = [m.start() for m in re.finditer(pattern3, ref_seq)]
    loci = loci_1 + loci_2 + loci_3
    output = open(dirn + "loci_g.csv", "w")
    for l in loci:
        output.write(str(l)+"\n")

def compare_CDS():
    import Codon
    import re
    codon_table = Codon.c_table()
    dirn = "/Workplace/Python/PCRtagChange"
    wt = Fasta(dirn + "yeast_chr04_0_00_genes.fa")
    syn = Fasta(dirn + "yeast_chr04_3_66_genes.fa")
    # for cds_s in syn:

def chr04_pcrtag_stat():
    """Stat the pcrtags over the minichunks."""

    import re
    dirn = "/workspace/Python/161212_megachunk_csPCR/"
    primers = "PCRtags_syn.csv"
    mega = "synIV_mega.fa"
    output = "chr04_pcrtag_stat.csv"
    mega_f = Fasta(dirn+mega)
    primers_fp = open(dirn+primers, "r")
    op = open(dirn + output, "w")

    primer_ls = []
    csv = primers_fp.read()
    primer_ls = re.split("\r|,", csv)
    size = len(primer_ls)/2
    stats = []
    for i in range(size):
        if "synF" in primer_ls[2*i]:
            p_seq = primer_ls[2*i+1].upper()
        else:
            if "\n" not in primer_ls[2*i+1]:
                p_seq = Seq_Analyzer(primer_ls[2*i+1]).rcSeq().upper()
        for j in range(len(mega_f.Names)):
            sq = mega_f.Seqs[j].upper()
            if p_seq in sq:
                _pos = sq.find(p_seq)
                stats.append([mega_f.Names[j], primer_ls[2*i],
                              primer_ls[2*i], _pos])
        print "%s done!" %primer_ls[2*i]
    for item in stats:
        op.write("%s, %s, %s, %d\n" % (item[0], item[1], item[2], item[3]))
    op.close()
    mega_f.close()

def print_minichunk(chunk_ls):
    import minichunks as mc
    mc_dict = mc.minichunk_dict()
    out = open("out.txt", "w")

    for chunk in chunk_ls:
        for m in mc_dict[chunk]:
            ch = m[0].split()[0]
            out.write("%s \n" %ch)
    out.close()

def batch_PCR_Primers_for_minichunks():

    folder = "/workspace/Python/170116_SynIV/minichunks/"
    savefile = "/workspace/Python/170116_SynIV/primers.csv"

    L_res = 2
    R_res = 2
    Tm_ls = range(52,59)
    minLen = 20
    maxLen = 50
    primer_ls = []
    output = open(savefile, "w+")

    _len = len(Tm_ls)

    for fn in os.listdir(folder):
        if "fasta" in fn:
            fa = Fasta(folder+fn)
            length = len(fa.Seqs)
            for n in range(length):
                seq = fa.Seqs[n]
                name = fa.Names[n]
                F_primer, R_primer, Tm_bin = Seq_Analyzer(seq).\
                    Find_Primer_at_Ends(L_res, R_res, Tm_ls, minLen, maxLen)
                primer_ls.append([name, F_primer, R_primer, Tm_bin])

    for primer_sub_ls in primer_ls:
        for n in range(_len):
            if primer_sub_ls[3][n]:
                output.write(str(primer_sub_ls[0][:-1]) + ",")
                for i in [1, 2]:
                    for j in range(4):
                        output.write(str(primer_sub_ls[i][n][j]) + ",")
                output.write(str(primer_sub_ls[2][n][4]) + ",")
                output.write("\n")


def extract_terminial_seq(length=10000,prefix=""):
    dirn = "/home/zhuwei/cglabarata/comp/"
    fn = "52.fasta"
    outdirn = dirn + "52_ends/"
    fa = Fasta(dirn+fn)
    for n in range(len(fa.Names)):
        _name = fa.Names[n]
        if _name[:3].upper() != "CHR":
            continue
        _seq = fa.Seqs[n]
        _len = len(_seq)
        output_prefix = prefix+_name.split()[0]
        with open(outdirn+output_prefix+"_Left.fa", "w") as op:
            op.write(">%s_left\n" % output_prefix)
            op.write(_seq[:length])
            op.write("\n")
        with open(outdirn+output_prefix+"_Right.fa", "w") as op:
            op.write(">%s_right\n" % output_prefix)
            op.write(_seq[-length:])
            op.write("\n")

    


if __name__ == "__main__":
    extract_terminial_seq(prefix="52_")