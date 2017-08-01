class Fasta(object):

    """Class to process fasta file."""

    def __init__(self, fn):
        """ fn = file name."""
        # if fn[-5:].upper() != "FASTA" and fn[-2:].upper() != "FA":
        #   print(fn+" is NOT a FASTA File!\n")
        #   return
        self.fp = open(fn, "ru")
        self.Names = []
        self.Seqs = []
        self.Fasta_reader()
        self.size = len(self.Names)
        self.close()
        print self.size

    def Fasta_reader(self):
        self.fp.seek(0)
        _line = self.fp.readline()
        _line_tmp = "".join(_line.split())
        _seq = ""
        while _line:
            if _line_tmp == "":
                _line = self.fp.readline()
                _line_tmp = "".join(_line.split())
                continue
            if _line_tmp[0] == ">":
                h = _line.find(">")
                self.Names.append(_line[h+1:-1])  # OMIT > \n
                print _line[h+1:]
                self.Seqs.append(_seq)
                _seq = ""
                _line = self.fp.readline()
                _line_tmp = "".join(_line.split())
                continue
            _seq += _line_tmp
            _line = self.fp.readline()
            _line_tmp = "".join(_line.split())
        self.Seqs.append(_seq)
        if "" in self.Seqs:
            self.Seqs.remove("")
        if len(self.Seqs) != len(self.Names):
            print "Error Reading Fasta File: \n\
Unequal number of Sequences and Names\n"

    def splitFasta(self, dirn):
        self.fp.seek(0)
        _fasta = self.fp.read()
        pos = -1

        for name in self.Names:
            pos = _fasta.find(name)
            _name = name.split()[0]
            pos_next = _fasta.find(">", pos + 1)
            print _name, pos, pos_next
            _fp = open(dirn + _name + ".fasta", "w")
            if pos_next == -1:
                _fp.write(_fasta[pos - 1:])
                _fp.close()
                continue
            _fp.write(_fasta[pos - 1:pos_next])
            _fp.close()

    def close(self):
        self.fp.close()
