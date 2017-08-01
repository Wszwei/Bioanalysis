# Analysis CSV files
# Number of the columns is controlled by the first line

import os
import datetime
import math


class csv(object):

    def __init__(self, fn="", isHeaded=False, sep=','):
        self.line = -1
        self.col = -1
        self.fn = fn
        self.isHeaded = isHeaded
        self.lined_content = []
        self.sep = sep
        self.fn = fn
        self.fp = None
        self.Header = []
        if fn != "":
            self.read()

    def read(self):
        if os.path.isfile(self.fn):
            with open(self.fn, "ru") as input:
                for line in input:
                    line = line.rstrip()
                    if line:
                        self.line += 1
                        line_col = line.split(self.sep)
                        if self.line == 0:
                            self.col = len(line_col)
                            if self.isHeaded:
                                self.Header = line_col
                            else:
                                self.lined_content.extend(line_col)
                        else:
                            self.lined_content.extend(line_col)
                            if len(line_col) != self.col:
                                self._warn("Error in line %d: Wrong col number"
                                           % line)
            self._warn("%s: Csv File loaded" % self.fn)
        else:
            self._warn("%s: File Not Found" % self.fn)

    def write(self):
        pass

    def _warn(self, message):
        print "%s\t%s\n"\
              % (datetime.datetime.now(), message)


class list_groups(object):
    # generate a list of Gourp Numbers
    # Max line and col = 65535
    # col_pattern: Int x: element in col a + t*x is Grouped for the same a<x
    # line_pattern: Int x: element in line a + t*x is Grouped for the same a<x
    # block_pattern: [Int x, Int y]: element in:
    #                line m: a*x <= m <(a+1)*x
    #                col n: b*y <= n < (b+1)*y is Grouped for the same a,b
    #                ONLY entact blocks are labeled
    # line and col: Number of lines and columns for the grouped list
    def __init__(self, line=0, col=0,
                 col_pattern=65535, line_pattern=65535,
                 block_pattern=[65535, 65535]):
        self.line = line
        self.col = col
        self.col_pattern = col_pattern
        self.line_pattern = line_pattern
        self.block_grouped = block_pattern
        self.col_grouped = [0] * line * col
        self.line_grouped = [0] * line * col
        self.block_grouped = [0] * line * col
        if col_pattern < col:
            self.self.group_column()
        if line_pattern < line:
            self.group_line()
        if block_pattern[0] < line or block_pattern[1] < col:
            self.group_block()

    def generate_group_list(self):
        # Return group labels
        # groups = [[line labeles], [col labels], [block labels], [All]]
        # [All]: block_line_col
        line_digit = int(math.log(self.line // self.line_pattern) /
                         math.log(10)) + 1
        col_digit = int(math.log(self.col // self.col_pattern) /
                        math.log(10)) + 1
        block_digit = (int(math.log(self.line // self.block_grouped[0]) /
                       math.log(10)) + 1) *\
                      (int(math.log(self.col // self.block_grouped[1]) /
                       math.log(10)) + 1)
        line_label = map(lambda x: str(x).zfill(line_digit), self.line_grouped)
        col_label = map(lambda x: str(x).zfill(col_digit), self.col_grouped)
        block_label = map(lambda x: str(x).zfill(block_digit),
                          self.block_grouped)
        all_label = map(lambda x, y, z: "Block_%s_LGroup_%s_CGroup_%s"
                        % (x, y, z), block_label, line_label, col_label)
        return [line_label, col_label, block_label, all_label]

    def group_column(self):
        _col_index = range(self.col_pattern)
        _number_entact_col_group = self.col // self.col_pattern
        _extra_col_per_line = self.col % self.col_pattern
        _col_label = (_col_index*_number_entact_col_group)
        _col_label = _col_label.extend(range(_extra_col_per_line))
        self.col_grouped = _col_label*self.line

    def group_line(self):
        _number_entact_line_group = self.line // self.line_pattern
        _extra_line_per_line = self.line % self.line_pattern

        def _line_block(_line):
            _block = []
            for i in range(_line):
                _block.extend([i]*self.col)
            return _block
        self.line_grouped = _line_block(self.line_pattern) * \
            _number_entact_line_group
        if _extra_line_per_line > 0:
            self.line_grouped.extend(_line_block(_extra_line_per_line))

    def group_block(self):
        # Only Entact blocks are labeled
        _line = self.block_pattern[0]
        _col = self.block_grouped[1]
        _height = self.line // _line
        _length = self.col // _col
        _odd_col = self.col % _col
        _odd_line = self.line % _line
        self.block_grouped = []
        for i in _height:
            _block = []
            for j in _length:
                _no = i*_length + j
                _block.extend([_no] * _col)
            if _odd_col > 0:
                _block.extend([-1] * _odd_col)
            _block = _block * _line
            self.block_grouped.extend(_block)
        if _odd_line > 0:
            self.block_grouped.extend([-1]*self.col*_odd_line)

    def _warn(self, message):
        print "%s\t%s\n"\
              % (datetime.datetime.now(), message)

