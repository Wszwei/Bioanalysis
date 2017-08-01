from PyQt4 import QtGui
from PyQt4 import QtCore
import os
from sequence_analysis import Seq_Analyzer
from fasta import Fasta
import sys


class Window(QtGui.QMainWindow):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.initUI()

    def initUI(self):
        _ver = "version 0.1 alpha"
        self.setGeometry(300, 300, 400, 400)
        self.setWindowTitle("Sequence Processor " + _ver)
        self.tabWidget = QtGui.QTabWidget(self)
        self.tab_seq_analyzer = QtGui.QFrame()
        self.tab_seq_processor = QtGui.QFrame()
        self.tabWidget.addTab(self.tab_seq_analyzer, "Sequence Analyzer")
        self.tabWidget.addTab(self.tab_seq_processor, "Sequences Processor")
        self.initSeq_Analyzer()
        self.initSeq_Processor()
        self.initMenu()
        self.setCentralWidget(self.tabWidget)
        self.dirn = "/home"
        self.currentSeq = ""
        self.show()

    def initSeq_Analyzer(self):
        self.lName = QtGui.QLabel("")
        self.initSeqEdit()
        self.SelectEdit = QtGui.QTextEdit(self)
        self.findSeqLine = QtGui.QLineEdit(self)
        self.findSeqLine.returnPressed.connect(self.findSeq)
        _vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
        _vsplitter.addWidget(self.SeqEdit)
        _vsplitter.addWidget(self.SelectEdit)
        _vsplitter.setCollapsible(0, False)
        _vsplitter.setCollapsible(1, False)
        _layout = QtGui.QGridLayout()
        _layout.addWidget(self.lName, 0, 0)
        _layout.addWidget(_vsplitter, 1, 0)
        _layout.addWidget(self.findSeqLine, 2, 0)
        self.SeqEdit.setFocus()
        self.tab_seq_analyzer.setLayout(_layout)
        _vsplitter.setSizes([300, 100])

    def initSeqEdit(self):
        class SeqEdit(QtGui.QTextEdit):
            def __init__(self, parent=None):
                super(SeqEdit, self).__init__(parent)
                self.cursor = self.textCursor()
                self.sel_format = QtGui.QTextCharFormat()
                self.sel_color = "red"
                self.sel_format.setBackground(
                    QtGui.QBrush(QtGui.QColor(self.sel_color)))

            def hlightRegion(self, pos_init_lst, length):
                # current_cursor_pos=self.cursor.position()
                for _pos in pos_init_lst:
                    self.cursor.setPosition(_pos)
                    self.cursor.movePosition(
                        QtGui.QTextCursor.NextCharacter, 1, length)
                    self.cursor.mergeCharFormat(self.sel_format)

            def clean_background(self, slot):
                self.cursor.setPosition(QtGui.QTextCursor.Start)
                self.cursor.movePosition(QtGui.QTextCursor.End, 1)
                bg_color_format = QtGui.QTextCharFormat()
                bg_color_format.setBackground(
                    QtGui.QBrush(QtGui.QColor("white")))
                self.cursor.mergeCharFormat(bg_color_format)
        self.SeqEdit = SeqEdit(self)

    def initSeq_Processor(self):
        self.fnlist = QtGui.QListWidget(self)
        self.fnlist.itemActivated.connect(self.openfasta)
        self.featurelist = QtGui.QListWidget(self)
        self.initFunclst()
        self.initSeq_Processor_Layout()

    def initFunclst(self):
        self.funclst = QtGui.QListWidget(self)
        func_dic = {"Restriction Enzyme Selector": 1000}
        for _func in func_dic:
            self.funclst.addItem(
                QtGui.QListWidgetItem(_func, self.funclst, func_dic[_func]))
        self.funclst.itemActivated.connect(self.perform_function)

    def initMenu(self):
        exitAction = self.addAction("Exit", "Ctrl+Q", "Exit", self.close)
        openAction = self.addAction(
            "Open", "Ctrl+O", "Open", self.showOpenDialog)
        openDirAction = self.addAction(
            "Open Directory", "Ctrl+D", "Open Directory of Sequence files",
            self.showOpenDir)
        cleanHighlights_Action = self.addAction(
            "Clean all the highlights", "Ctrl+Alt+C",
            "Clean all the Highlighted Sequence",
            self.SeqEdit.clean_background)
        # findAction=self.addAction("Find","Ctrl+F","Find",self.showFindDialog)
        self.menubar = self.menuBar()
        # self.menuBar.setNativeMenuBar(False)
        fileMenu = self.menubar.addMenu('&File')
        editMenu = self.menubar.addMenu('&Edit')
        fileMenu.addAction(exitAction)
        fileMenu.addAction(openAction)
        fileMenu.addAction(openDirAction)
        editMenu.addAction(cleanHighlights_Action)

    def initSeq_Processor_Layout(self):
        lftlst = QtGui.QLabel("Feature List", self)
        lfnlst = QtGui.QLabel("File List", self)
        lfunclst = QtGui.QLabel("Function List", self)

        _hsplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
        _vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)

        _frfn = QtGui.QFrame(self)
        _frft = QtGui.QFrame(self)
        _frfunc = QtGui.QFrame(self)
        _vlayout = QtGui.QVBoxLayout(self)
        _vlayout.addWidget(lfnlst)
        _vlayout.addWidget(self.fnlist)
        _frfn.setLayout(_vlayout)

        _vlayout2 = QtGui.QVBoxLayout(self)
        _vlayout2.addWidget(lftlst)
        _vlayout2.addWidget(self.featurelist)
        _frft.setLayout(_vlayout2)

        _vlayout3 = QtGui.QVBoxLayout(self)
        _vlayout3.addWidget(lfunclst)
        _vlayout3.addWidget(self.funclst)
        _frfunc.setLayout(_vlayout3)

        _vsplitter.addWidget(_frft)
        _vsplitter.addWidget(_frfunc)
        _vsplitter.setCollapsible(0, False)
        _vsplitter.setCollapsible(1, False)
        _vsplitter.setHandleWidth(0)

        _hsplitter.addWidget(_frfn)
        _hsplitter.addWidget(_vsplitter)
        _hsplitter.setCollapsible(0, False)
        _hsplitter.setCollapsible(1, False)
        _hsplitter.setHandleWidth(0)
        # _layout=QtGui.QGridLayout(self)
        # _layout.addWidget(lfnlst,0,0)
        # _layout.addWidget(self.fnlist,1,0)
        _layout = QtGui.QVBoxLayout(self)
        _layout.addWidget(_hsplitter)
        self.tab_seq_processor.setLayout(_layout)

    def perform_function(self, emitter):
        func_dic = {1000: self.re_selector}
        func_dic[emitter.type()]()

    def re_selector(self):
        window = QtGui.QMainWindow(self)
        window.setWindowTitle("Restriction Enzyme Selector")
        _layout = QtGui.QGridLayout(window)
        _lcoef_min = QtGui.QLabel("Minimum Effective Size of thr Band (bp):")
        _lcoef_diff = QtGui.QLabel("Minimun Difference in Size (bp):")
        lEdit_Size_Limit = QtGui.QLineEdit(window)
        lEdit_Diff_limit = QtGui.QLineEdit(window)
        _layout.addWidget(_lcoef_min, 0, 0)
        _layout.addWidget(_lcoef_diff, 1, 0)
        _layout.addWidget(lEdit_Size_Limit, 0, 1)
        _layout.addWidget(lEdit_Diff_limit, 1, 1)
        _frame = QtGui.QFrame(window)
        _frame.setLayout(_layout)
        window.setCentralWidget(_frame)
        window.show()

    def findSeq(self):
        Seq = str(self.findSeqLine.text())
        print Seq
        if not str(self.currentSeq):
            return
        self.seqFindLst = Seq_Analyzer(self.currentSeq).Degen_searcher(Seq)
        _tmp_seq_str = ""
        for _pos in self.seqFindLst:
            _tmp_seq_str += (str(_pos) + " ; ")
        self.SeqEdit.hlightRegion(self.seqFindLst, len(Seq))
        self.SelectEdit.setText(_tmp_seq_str)

    def openfasta(self, fn):
        if fn.__class__.__name__ != "str":
            fn = self.dirn + "/" + str(fn.text())
        fasta = Fasta(fn)
        data = fasta.Seqs[0]
        self.lName.setText(fasta.Names[0])
        self.currentSeq = data
        self.SeqEdit.setPlainText(data)
        self.dirn = os.path.dirname(fn)
        self.tabWidget.setCurrentWidget(self.tab_seq_analyzer)

    def showOpenDir(self):
        dirn = str(QtGui.QFileDialog.getExistingDirectory(
            self, 'SelectDirectory', '/home'))
        self.fnlist.clear()
        flist = os.listdir(dirn)
        self.dirn = dirn
        for fn in flist:
            if fn.upper()[-5:] == "FASTA" or fn.upper()[-2:] == "FA":
                self.fnlist.addItem(fn)
        self.tabWidget.setCurrentWidget(self.tab_seq_processor)

    def showGC_Analyzer(self):
        pass

    def showTags_Labeler(self):
        pass

    def showOpenDialog(self):
        fn = str(QtGui.QFileDialog.getOpenFileName(
            self, 'Open file', self.dirn, "Fasta Files (*.fasta *.fa)"))
        self.openfasta(fn)

    def addAction(self, Name, ShortCut, Tip, Function):
        _Action = QtGui.QAction("&" + Name, self)
        if ShortCut:
            _Action.setShortcut(ShortCut)
        if Tip:
            _Action.setStatusTip(Tip)
        _Action.triggered.connect(Function)
        return _Action


def showapp():
    app = QtGui.QApplication(sys.argv)
    Window()
    sys.exit(app.exec_())


