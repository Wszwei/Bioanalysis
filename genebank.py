import os
import re
import datetime


class GenebBankFeature():

    def __init__(self, Type="", Pos=[-1, -1], isComp=False, Label=[]):
        self.Type = Type
        self.Pos = Pos
        self.isComp = isComp
        self.Label = Label


class GeneBank():

    def __init__(self, fn=""):
        self.seq = ""
        self.locus = ""
        self.length = 0
        self.type = ""
        self.circular = ""
        self.date = ""
        self.features = []
        self.fn = fn
        self.header = ""
        if os.path.isfile(fn):
            self.fp = open(fn, "ru")
            self.readfile()
            self.lengthcheck()
        else:
            self.fp = open(fn, "w")
            self._warn("New GBfile generated!")

    def close(self):
        if self.fp:
            self.fp.close()

    def readfile(self):
        self.fp.seek(0)
        self._content = self.fp.read()
        if self.header_reader() == -1:
            self._warn("Problem in reading header")
            return
        if self.feature_reader() == -1:
            self._warn("Problem in reading features!")
            return
        if self.seq_reader() == -1:
            self._warn("Problem in reading sequence!")
            return

    def header_reader(self):
        self.fp.seek(0)
        _line = self.nextline()
        if _line == -1:
            self._warn("Empty file!")
            return -1
        if _line[:5] != "LOCUS":
            self._warn("No LOCUS!")
            return -1
        _head = _line.split()
        if len(_head) < 2:
            self._warn("No Locus Info!")
        if "bp" not in _head:
            self.locus = self.merge_words(_head[1:])
        else:
            i = _head.index("bp")
            self.locus = self.merge_words(_head[1:i-1])
            self.length = int(_head[i-1])
            if len(_head) > i+1:
                self.type = _head[i+1]
            if len(_head) > i+2:
                self.circular = _head[i+2]
            if len(_head) > i+3:
                self.date = _head
        self.fp.seek(0)
        self.fp.readline()
        _i = self.fp.tell()
        _e = self._content.find("FEATURES")
        if _i < _e:
            self.header = self._content[_i: _e]
        return 1

    def feature_reader(self):
        _pos = self._content.find("FEATURES")
        if _pos == -1:
            return
        self.fp.seek(_pos)
        if not self.fp.readline():
            self._warn("No Features and Seq!")
            return -1
        _pos = self.fp.tell()
        if "ORIGIN" not in self._content:
            self._warn("No Seq Info!")
            return -1
        _end = self._content.find("ORIGIN")
        _feature = self._content[_pos-1: _end]
        _feature = "".join("%-4s" % i for i in _feature.split("\t"))
        _p_fheader = re.compile(r"^\s{5}\w+.*?")
        _feature_lines = _feature.split("\n")
        _n = len(_feature_lines)
        _i = 0
        _f = None
        _p_loc = re.compile(r"\d+\.\.\d+")
        _p_loc_r = re.compile(r"complement\(\d+\.\.\d+\)")
        _labels = []
        while _i < _n:
            _l = _feature_lines[_i]
            if _p_fheader.match(_l):
                _labels = []
                if _f:
                    _f.Label = _labels
                    self.features.append(_f)
                _t = _l.split()
                if len(_t) != 2:
                    self._warn("Problem in Feature line: %s" % _l)
                    return -1

                if _p_loc.match(_t[1]):
                    isComp = False
                    loc = [int(x) for x in _t[1].split("..")]
                elif _p_loc_r.match(_t[1]):
                    isComp = True
                    loc = [int(x) for x in _t[1][11:-1].split("..")]
                else:
                    self._warn("Problem in feature line %s" % _l)
                    return -1
                _f = GenebBankFeature(Type=_t[0], Pos=loc, isComp=isComp)
            else:
                _l = _l.strip()
                if _l:
                    if _l[0] == "/":
                        _labels.append(_l)
                    else:
                        if not _labels:
                            self._warn("Problem in feature line %s" % _l)
                            return -1
                        _labels[-1] += _l
            _i += 1
        _f.Label = _labels
        self.features.append(_f)

    def _warn(self, message):
        print "%s\t%s: %s\n"\
              % (datetime.datetime.now(), self.fn, message)
        # time.sleep(1)

    def seq_reader(self):
        if "ORIGIN" not in self._content or "//" not in self._content:
            self._warn("No Seq data in GB file!")
            return -1
        _s = self._content.find("ORIGIN")
        _e = self._content.find("//")
        if _e < _s:
            self._warn("Wrong format for Seq region!")
            return -1
        _content = self._content[_s+6: _e]
        _t = _content.split()
        self.seq = "".join([x for x in _t if not x.isdigit()])
        return 1

    def writefile(self, fn):
        if not os.path.isfile(fn):
            self._warn("New File Genereted")
        self.fp.seek(0)
        fp = open(fn, "w")
        top = "LOCUS    %s    %d bp %s     %s     %s\n"\
              % (self.locus, self.length, self.type,
                 self.circular, self.date)
        fp.write(top)
        fp.write(self.header)
        fp.write("FEATURES             Location/Qualifiers\n")
        for _f in self.features:
            fp.write("     %s          " % _f.Type)
            if _f.isComp:
                fp.write("complement(%d..%d)\n" % (_f.Pos[0], _f.Pos[1]))
            else:
                fp.write("%d..%d\n" % (_f.Pos[0], _f.Pos[1]))
            for _l in _f.Label:
                fp.write("                     %s\n" % _l)
        fp.write("ORIGIN\n")
        _len = len(self.seq)
        _line = _len/60
        _ori = ""
        if _line == 0:
            fp.write("       1 %s\n" % self.seq)
        else:
            for _l in range(_line):
                _ori += "       %d" % (_l*60+1)
                for i in range(6):
                    _ori += " "
                    _ori += self.seq[_l*60+i*10:_l*60+(i+1)*10]
                _ori += "\n"
            if _len % 60 != 0:
                _ori += "       %d " % (_line*60+1)
                _ori += self.seq[_line*60:]
                _ori += "\n"
        _ori += "//"
        fp.write(_ori)
        fp.close

    def lengthcheck(self):
        if self.length == -1:
            return
        else:
            if self.length != len(self.seq):
                self._warn("Incooridnate Length!")
                self.length = len(self.seq)
                self.writefile(self.fn)
            else:
                self._warn("Sequence Length verified!")

    def merge_words(self, wlist):
        s = wlist[0]
        for w in wlist[1:]:
            s = s + " " + w
        return s

    def nextline(self):
        _line = self.fp.readline()
        if not _line:
            return -1
        while "".join(_line.split()) == "":
            _line = self.fp.readline()
            if not _line:
                return -1
        return _line
