from collections import deque
import copy

class Seq_Analyzer():

    def __init__(self, Seq):
        self.seq = Seq.upper()

    def rcSeq(self):
        dic = {"A": "T",
               "C": "G",
               "T": "A",
               "G": "C"}
        rcseq = ""
        for ch in self.seq:
            rcseq = dic[ch] + rcseq
        return rcseq

    def GC_window_analyzer(self, size, supth, infth):
        _GC_window = size
        _GC_contents = []
        _outbounds = []
        if _GC_window >= len(self.seq):
            return [(self.seq.count("C") + self.seq.count("G"))/len(self.seq)]
        _len = len(self.seq)
        for _i in range(_len - _GC_window + 1):
            _seq = self.seq[_i:_i + _GC_window]
            _GC_contents.append(
                (_seq.count("C") + _seq.count("G")) / float(_GC_window))
        _pos = 0
        for _ratio in _GC_contents:
            _pos += 1
            if _ratio > supth or _ratio < infth:
                _outbounds.append((_pos, _ratio))
        print len(self.seq)
        return _GC_contents, _outbounds

    def GC_window_analyzer_visual(self, size, supth, infth, name):
        """ size = size of the Window
            supth = upper GC limit
            infth = lower GC limit
            name = name of the seq
        """
        import matplotlib.pyplot as plt
        _GC_list, _outbounds = self.GC_window_analyzer(size, supth, infth)
        _len = len(_GC_list)
        # for _x in _outbounds:
        #   print _x
        x = range(_len)
        y1 = [supth] * _len
        y2 = [infth] * _len
        plt.figure(figsize=(16, 7))
        plt.plot(x, _GC_list, label="$GC content$", color="blue", linewidth=2)
        plt.plot(x, y1, color="red", linewidth=1)
        plt.plot(x, y2, color="red", linewidth=1)
        plt.xlabel("Position(bp)")
        plt.ylabel("GC content")
        plt.title("GC content:" + name)
        plt.ylim = (0, 1)
        plt.show()

    def Degen_generater(self):
        degen_nt = ["W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N"]
        degen_dict = {"W": ["A", "T"],
                      "S": ["C", "G"],
                      "M": ["A", "C"],
                      "K": ["G", "T"],
                      "R": ["A", "G"],
                      "Y": ["C", "T"],
                      "B": ["C", "G", "T"],
                      "D": ["A", "G", "T"],
                      "H": ["A", "C", "T"],
                      "V": ["A", "C", "G"],
                      "N": ["A", "T", "C", "G"]
                      }
        _degen_seq = deque([""])
        for nt in self.seq:
            if nt not in degen_nt:
                for _m in range(len(_degen_seq)):
                    _seq = (_degen_seq).popleft()
                    _degen_seq.append(_seq + nt)
                continue
            _degen_nt = degen_dict[nt]
            for _m in range(len(_degen_seq)):
                _seq = (_degen_seq).popleft()
                for _dnt in _degen_nt:
                    _degen_seq.append(_seq + _dnt)
        return _degen_seq

    def Degen_searcher(self, seq):
        Degen_seqs = Seq_Analyzer(seq).Degen_generater()
        site_list = []
        print Degen_seqs
        for _degen_seq in Degen_seqs:
            _pos = self.seq.find(_degen_seq)
            while _pos >= 0:
                site_list.append(_pos)
                _pos = self.seq.find(_degen_seq, _pos + 1)
        site_list.sort()
        return site_list

    def Repetitive_Seq_analyzer(self, mer, rep_offset, top_filter):
        _seq = ""
        for _m in range(mer):
            _seq += "N"
        _site_list = Seq_Analyzer(_seq).Degen_generater()
        _rep_site_list = []
        for _site in _site_list:
            _pos = self.seq.find(_site)
            _counter = 0
            _counter_list = []
            if _pos < 0:
                continue
            while _pos >= 0:
                _pos_next = self.seq.find(_site, _pos + mer)
                if _pos_next == _pos + mer:
                    _counter += 1
                    _pos = _pos_next
                else:
                    _counter_list.append(_counter + 1)
                    _counter = 0
                    _pos = _pos_next
            _max = max(_counter_list)
            if _max >= rep_offset:
                _rep_site_list.append([_site, _max])
        if len(_rep_site_list) <= top_filter:
            return _rep_site_list
        else:
            _tmp_list = []
            for _m in _rep_site_list:
                _tmp_list.append(_m[1])
            _counter_ls = _tmp_list[:top_filter]
            _min = min(_counter_ls)
            _max = max(_counter_ls)
            print "c\n"
            for _m in _tmp_list:
                if _m >= _max:
                    _counter_ls.remove(_min)
                    _counter_ls.append(_m)
                    _max = _m
                    _min = min(_counter_ls)
            _pos = 0
            _tmp = copy.copy(_rep_site_list)
            _rep_site_list = []
            for _m in _counter_ls:
                _pos = _tmp_list.index(_m, _pos)
                _rep_site_list.append(_tmp[_pos])
            _tmp = []
            return _rep_site_list

    def TmCalculator(self):
        """ Calculate Tm of a Sequence.

        primer: 50 nM ; Na+: 50 mM, pH = 7
        """
        _seq = self.seq
        _len = len(_seq)
        _GC = self.seq.count("G") + self.seq.count("C")

        if _len <= 13:
            _Tm = (_len-_GC) * 2 + _GC * 4
            return _Tm

        _Tm = 64.9 + 41 * (_GC - 16.4) / _len
        return _Tm

    def GC_ratio(self):
        """ Calculate GC ratio. """
        _len = len(self.seq)
        _GC = self.seq.count("G") + self.seq.count("C")
        return 1.0 * _GC / _len

    def Find_Forward_Primer(self, Tm, Lim, minLen=20, maxLen=30, crossref=""):
        """ Find the PCR primer from the beginning of a sequence.

        return:
            primer, h, t _GC, this_GC, this_multi
            primer, h postion, length , GC ratio,
            if 0.4<GC<0.6, if contains multi nt
        pos of 5' end: 1 - Lim
        default min length = 20
        default max length = 30
        crossref = cross reference sequence that the pcr primer also works
        rules:
            1. All primers end in G/C
            2. Avoid multi A/T/C/G (more than 4)
            2. GC ratio 0.4 - 0.6 preferred
            3. shorter primer preferred
        """
        isReffed = (len(crossref) != 0)
        crossref = crossref.upper()
        _length = len(self.seq)
        primer_ls = []
        multi_nt = ["AAAA", "TTTT", "CCCC", "GGGG"]
        # Flag to check the existance of primer with no multi nt
        _multi_flag = False
        # Flag to check the existance of primer with 0.4<GC<0.6
        _GC_flag = False
        for h in range(_length - Lim):
            for t in range(minLen, maxLen):
                this_GC = False
                this_multi = True
                if self.seq[h+t-1].upper() not in ["G", "C"]:
                    continue
                primer = self.seq[h:h+t]
                if isReffed:
                    if primer.upper() not in crossref:
                        continue
                _tm = int(round(Seq_Analyzer(primer).TmCalculator()))
                if _tm == Tm:
                    _GC = Seq_Analyzer(primer).GC_ratio()
                    if 0.4 <= _GC <= 0.6:
                        _GC_flag = True
                        this_GC = True
                    for nt in multi_nt:
                        if nt in primer:
                            this_multi = False
                    if this_multi:
                        _multi_flag = True
                    primer_ls.append([primer, h, t, _GC,
                                     this_GC, this_multi, Tm])
        primer_ls_filtered = []
        """ filter_condition:
        -1 : init
        0: all primers contain multi_nt and worse GC
        1: all primers with no multi_nt only
        2: all primers with no worse GC only
        3: all primers with both chara
        """
        filter_condition = -1
        for primer in primer_ls:
            if _multi_flag and _GC_flag:
                filter_condition = 3
                if primer[4] and primer[5]:
                    primer_ls_filtered.append(primer)
            elif (not _GC_flag) and (not _multi_flag):
                filter_condition = 0
                primer_ls_filtered.append(primer)
            elif not _GC_flag:
                filter_condition = 1
                if primer[5]:
                    primer_ls_filtered.append(primer)
            elif not _multi_flag:
                filter_condition = 2
                if primer[4]:
                    primer_ls_filtered.append(primer)
        size = len(primer_ls_filtered)
        if filter_condition == 2 or 3:
            if size == 0:
                primer_ls_filtered = primer_ls
                filter_condition = 0
                size = len(primer_ls_filtered)
        best_primer = []
        if filter_condition == 3:
            length = 100
            for primer in primer_ls_filtered:
                if len(primer[0]) < length:
                    best_primer = primer
                    length = len(primer[0])
        elif filter_condition == 1:
            count = 100
            for primer in primer_ls_filtered:
                c = 0
                for nt in multi_nt:
                    c += primer[0].count(nt)
                if c < count:
                    best_primer = primer
                    count = c
                elif c == count:
                    if len(primer[0]) < len(best_primer[0]):
                        best_primer = primer
        elif filter_condition == 2:
            gc = 0.5
            for primer in primer_ls_filtered:
                _GC = Seq_Analyzer(primer[0]).GC_ratio()
                delta = abs(_GC - 0.5)
                if delta < gc:
                    best_primer = primer
                    gc = delta
                elif delta == gc:
                    if len(primer[0]) < len(best_primer[0]):
                        best_primer = primer
        elif filter_condition == 0:
            count = 100
            gc = 0.5
            for primer in primer_ls_filtered:
                c = 0
                _GC = Seq_Analyzer(primer[0]).GC_ratio()
                delta = abs(_GC - 0.5)
                for nt in multi_nt:
                    c += primer[0].count(nt)
                if c < count:
                    best_primer = primer
                    count = c
                    gc = delta
                elif c == count:
                    if delta < gc:
                        best_primer = primer
                        gc = delta
                    elif delta == gc:
                        if len(primer[0]) < len(best_primer[0]):
                            best_primer = primer
        if not best_primer:
            # print("NO PRIMER for Tm = %s\n" % (Tm))
            if Tm >= 58:
                return Seq_Analyzer(self.seq).\
                    Find_Forward_Primer(Tm-1, Lim, minLen, maxLen, crossref)
            else:
                return []
        else:
            print("Success!\n")
        return best_primer

    def Find_Reverse_Primer(self, Tm, Lim, minLen=20, maxLen=30, crossref=""):
        """ Find the PCR primer from the end of a sequence.

        return:
            primer, h, t _GC, this_GC, this_multi
            primer, h postion, length , GC ratio,
            if 0.4<GC<0.6, if contains multi nt
        using Find_Forward_Primer to perform the function
        pos of 5' end from  the end: 1 - Lim
        default min length = 20
        default max length = 30
        crossref
        rules:
            1. All primers end in G/C
            2. Avoid multi A/T/C/G (more than 4)
            2. GC ratio 0.4 - 0.6 prefered
            3. shorter primer perfered
        """
        seq = Seq_Analyzer(self.seq).rcSeq()
        if len(crossref) != 0:
            crossref = Seq_Analyzer(crossref).rcSeq()
        return Seq_Analyzer(seq).Find_Forward_Primer(
            Tm, Lim, minLen, maxLen, crossref)

    def Find_Primer_at_Ends_Tm(self, L_res, R_res, Tm, minLen, maxLen):
        """Find all Primers around around 5' and 3' ends for specific Tm."""
        # Obtain all probable primers
        _length = len(self.seq)
        F_primer_ls = []
        R_primer_ls = []

        # Forward Primer
        for h in range(L_res):
            for t in range(minLen, maxLen):
                _seq = self.seq[h: h+t]
                _tm = int(round(Seq_Analyzer(_seq).TmCalculator()))
                if _tm == Tm:
                    _GC = Seq_Analyzer(_seq).GC_ratio()
                    F_primer_ls.append([_seq, _tm, h, _GC])

        for h in range(1, R_res):
            for t in range(minLen, maxLen):
                _seq = self.seq[-(h+t): -h]
                _tm = int(round(Seq_Analyzer(_seq).TmCalculator()))
                if _tm == Tm:
                    _seq = Seq_Analyzer(_seq).rcSeq()
                    _GC = Seq_Analyzer(_seq).GC_ratio()
                    R_primer_ls.append([_seq, _tm, _length - h, _GC, h])
        return F_primer_ls, R_primer_ls

    def Find_Primer_at_Ends(self, L_res, R_res, Tm_ls, minLen, maxLen):
        """ Find Primers around 5' and 3' ends.

        L_res, R_res = Int : search primers from 0-n from the ends
        Tm_ls = [Int] : List of accepted Tm
        minLen =  Int : minimum length of the primer
        maxLen = Int : maximum length of the primer
        All primers end in G/C
        return best_F_primer_ls, best_R_primer_ls, paired_Tm_bin
        best_R_primer_ls, best_R_primer_ls =[(Seq,Tm,pos,GC ratio)]
        pos = pos of 5'end
        paired_Tm_bin = [bool]
        """
        # First obtain all probable primers
        _length = len(self.seq)
        F_primer_ls = []
        R_primer_ls = []

        # Forward Primer
        for h in range(L_res):
            for t in range(minLen, maxLen):
                _seq = self.seq[h: h+t]
                _tm = int(round(Seq_Analyzer(_seq).TmCalculator()))
                if _tm in Tm_ls:
                    _GC = Seq_Analyzer(_seq).GC_ratio()
                    F_primer_ls.append([_seq, _tm, h, _GC])

        for h in range(1, R_res):
            for t in range(minLen, maxLen):
                _seq = self.seq[-(h+t): -h]
                _tm = int(round(Seq_Analyzer(_seq).TmCalculator()))
                if _tm in Tm_ls:
                    _seq = Seq_Analyzer(_seq).rcSeq()
                    _GC = Seq_Analyzer(_seq).GC_ratio()
                    R_primer_ls.append([_seq, _tm, _length - h, _GC, h])

        # remove Primers that don't end in G/C:
        F_primer_ls = filter(lambda x: x[0][-1] in ["G", "C"], F_primer_ls)
        R_primer_ls = filter(lambda x: x[0][-1] in ["G", "C"], R_primer_ls)

        # Select for best primer
        # Rule 1: better GC
        # Rule 2: longer product
        # Rule 3: shorter primer

        best_F_primer_ls = [("", -1, -1)] * len(Tm_ls)
        best_R_primer_ls = [("", -1, -1)] * len(Tm_ls)
        for primer in F_primer_ls:
            _p = Tm_ls.index(primer[1])
            currentbest = best_F_primer_ls[_p]
            if currentbest[0] == "":
                best_F_primer_ls[_p] = primer
            # GC range: 40% - 60%
            elif ((not (0.4 <= currentbest[3] <= 0.6))
                  and (0.4 <= primer[3] <= 0.6)):
                best_F_primer_ls[_p] = primer
            # longer product
            elif (primer[2] < currentbest[2]):
                best_F_primer_ls[_p] = primer
            # shorter primer
            elif (len(primer[0]) < len(currentbest[0])):
                best_F_primer_ls[_p] = primer
        for primer in R_primer_ls:
            _p = Tm_ls.index(primer[1])
            currentbest = best_R_primer_ls[_p]
            if currentbest[0] == "":
                best_R_primer_ls[_p] = primer
            # GC range: 40% - 60%
            elif ((not (0.4 <= currentbest[3] <= 0.6)) and
                  (0.4 <= primer[3] <= 0.6)):
                best_R_primer_ls[_p] = primer
            # longer product
            elif (primer[2] > currentbest[2]):
                best_R_primer_ls[_p] = primer
            # shorter primer
            elif (len(primer[0]) < len(currentbest[0])):
                best_R_primer_ls[_p] = primer

        paired_Tm_bin = [False] * len(Tm_ls)
        for n in range(len(Tm_ls)):
            if (best_F_primer_ls[n][0] != "" and best_R_primer_ls[n][0] != ""):
                paired_Tm_bin[n] = True

        return best_F_primer_ls, best_R_primer_ls, paired_Tm_bin

