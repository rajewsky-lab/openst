import sys, os
import numpy as np

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from matplotlib import pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


def cpselect(_img1, _img2):
    """
    Tool for selection a individual number of control points in any two pictures

    :param img_path1: first image
    :param img_path2: second image

    :return: tuple with the (2d) points of image 1 and image 2
    """
    global img1
    global img2

    img1, img2 = _img1, _img2

    app = QApplication(sys.argv)
    cps = _MainWindow()
    cps.raise_()
    app.exec_()

    points1, points2 = [], []
    for cp in cps.wp.canvas.CPlist:
        points1.append(cp[0])
        points2.append(cp[1])

    del img1, img2

    points1 = np.array(points1)
    points2 = np.array(points2)

    return points1, points2


class _MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # os.path.join(__file__, relative)

        self.createWindow()
        self.createConn()

        self.showMaximized()
        self.show()
        self.wp.canvas.fig.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99)

        q1 = QTimer(self)
        q1.setSingleShot(False)
        q1.timeout.connect(self.updateGUI)
        q1.start(100)

    def createWindow(self):
        self.setWindowTitle('Control Point Selection Tool')
        self.setWindowIcon(QIcon(os.path.join(os.path.dirname(__file__), 'img/cpselect32.ico')))

        widget = QWidget(self)

        self.setCentralWidget(widget)

        hlay = QHBoxLayout(widget)
        vlay = QVBoxLayout()
        vlay2 = QVBoxLayout()
        vlay2.setSpacing(20)
        hlay_buttons = QHBoxLayout()

        hlay.addLayout(vlay)
        hlay.addLayout(vlay2)

        self.wp = _WidgetPlot(self)
        vlay.addWidget(self.wp)

        self.help = QTextEdit()
        self.help.setReadOnly(True)
        self.help.setMaximumWidth(400)
        self.help.setMinimumHeight(540)

        help_header = '<!DOCTYPE html><html lang="de" id="main"><head><meta charset="UTF-8"><title>Manual selection for open-ST alignment</title><style>td,th{font-size:14px;}p{font-size: 14px;}</style></head>'
        help_body = '<body><h1>Manual selection for open-ST alignment&emsp;</h1><h2>Navigation Toolbar</h2><img src="{}" alt="navbuttons"><br/><table cellspacing="20px"><tr><th valign="middle" height="20px">Tool</th><th valign="middle" height="20px">how to use</th></tr><tr><td><img src="{}" alt="homebutton"></td><td valign="middle">For all Images, reset to the original view.</td></tr><tr><td><img src="{}" alt="backwardforwardbutton"></td><td valign="middle">Go back to the last or forward to the next view.</td></tr><tr><td><img src="{}" alt="panzoombutton"></td><td valign="middle">Activate the pan/zoom tool. Pan with left mouse button, zoom with right</td></tr><tr><td><img src="{}" alt="backwardforwardbutton"></td><td valign="middle">Zoom with drawing a rectangle</td></tr></table><h2>Pick Mode</h2><p>Change into pick mode to pick up your control points. You have to pick the control points in both images before you can start to pick the next point.</p><p>If you have activated a navigation command (pan or zoom), you have to deactivate the command, before you can pick again</p><p>Press the red button below to start pick mode.</p><h2>Control Point list</h2><p>Below in the table, all your control points are listed. You can delete one ore more selected control points with the <b>delete</b> button.</p><h2>Return</h2><p>If you are finished, please press the <b>return</b> button below. You will come back to wherever you have been.</p></body></html>'
        help_html = help_header + help_body.format(os.path.join(os.path.dirname(__file__), 'img/navbuttons.PNG'),
                                                   os.path.join(os.path.dirname(__file__), 'img/homebutton.png'),
                                                   os.path.join(os.path.dirname(__file__), 'img/backforwardbutton.png'),
                                                   os.path.join(os.path.dirname(__file__), 'img/panzoombutton.png'),
                                                   os.path.join(os.path.dirname(__file__), 'img/zoomboxbutton.png'))
        self.help.insertHtml(help_html)
        self.cpTabelModel = QStandardItemModel(self)
        self.cpTable = QTableView(self)
        self.cpTable.setModel(self.cpTabelModel)
        self.cpTable.setMaximumWidth(400)

        self.delButton = QPushButton('Delete selected Control Point')
        self.delButton.setStyleSheet("font-size: 16px")

        self.pickButton = QPushButton("pick mode")
        self.pickButton.setFixedHeight(60)
        self.pickButton.setStyleSheet("color: red; font-size: 16px;")

        self.exitButton = QPushButton('Return')
        self.exitButton.setFixedHeight(60)
        self.exitButton.setStyleSheet("font-size: 16px;")

        vlay2.addWidget(self.help)
        vlay2.addWidget(self.cpTable)
        vlay2.addWidget(self.delButton)

        vlay2.addLayout(hlay_buttons)
        hlay_buttons.addWidget(self.pickButton)
        hlay_buttons.addWidget(self.exitButton)

        self.updateCPtable()
        self.statusBar().showMessage('Ready')

    def createConn(self):
        self.pickButton.clicked.connect(self.pickmodechange)
        self.exitButton.clicked.connect(self.menu_quit)
        self.delButton.clicked.connect(self.delCP)

    def menu_quit(self):
        self.close()

    def pickmodechange(self):
        if self.wp.canvas.toolbar.mode.value in ['', None]:
            if self.wp.canvas.pickmode == True:
                self.wp.canvas.pickMode_changed = True
                self.wp.canvas.pickmode = False
                self.statusBar().showMessage('Pick Mode deactivate.')
                self.wp.canvas.cursorGUI = 'arrow'
                self.wp.canvas.cursorChanged = True
            else:
                self.wp.canvas.pickMode_changed = True
                self.wp.canvas.pickmode = True
                #self.wp.canvas.toolbar._active = ''
                self.statusBar().showMessage('Pick Mode activate. Select Control Points.')

    def delCP(self):

        rows = self.cpTable.selectionModel().selectedRows()
        for row in rows:
            try:
                idp = int(row.data())
                for cp in self.wp.canvas.CPlist:
                    if cp.idp == idp:
                        index = self.wp.canvas.CPlist.index(cp)
                        self.wp.canvas.CPlist.pop(index)
            except:
                pass

        self.wp.canvas.updateCanvas()
        self.wp.canvas.cpChanged = True

    def updateGUI(self):

        if self.wp.canvas.toolbar.mode.value not in ['', None]:
            self.wp.canvas.pickmode = False
            self.wp.canvas.pickMode_changed = True

        if self.wp.canvas.pickMode_changed:
            if not self.wp.canvas.pickmode:
                self.pickButton.setStyleSheet("color: red; font-size: 20px;")
            elif self.wp.canvas.pickmode:
                self.pickButton.setStyleSheet("color: green; font-size: 20px;")
            self.wp.canvas.pickMode_changed = False

        if self.wp.canvas.cursorChanged:
            if self.wp.canvas.cursorGUI == 'cross':
                QApplication.setOverrideCursor(QCursor(Qt.CrossCursor))
            elif self.wp.canvas.cursorGUI == 'arrow':
                QApplication.restoreOverrideCursor()
            self.wp.canvas.cursorChanged = False

        if self.wp.canvas.cpChanged:
            self.updateCPtable()

    def updateCPtable(self):
        self.wp.canvas.cpChanged = False
        self.cpTable.clearSelection()
        self.cpTabelModel.clear()
        self.cpTabelModel.setHorizontalHeaderLabels(
            ['Point Number', 'x (Img 1)', 'y (Img 1)', 'x (Img 2)', 'y (Img 2)'])

        for cp in self.wp.canvas.CPlist:
            idp, x1, y1, x2, y2 = cp.coordText

            c1 = QStandardItem(idp)
            c2 = QStandardItem(x1)
            c3 = QStandardItem(y1)
            c4 = QStandardItem(x2)
            c5 = QStandardItem(y2)

            row = [c1, c2, c3, c4, c5]

            for c in row:
                c.setTextAlignment(Qt.AlignCenter)
                c.setFlags(Qt.ItemIsEditable)
                c.setFlags(Qt.ItemIsSelectable)

            self.cpTabelModel.appendRow(row)

        self.cpTable.resizeColumnsToContents()


class _WidgetPlot(QWidget):
    def __init__(self, *args, **kwargs):
        QWidget.__init__(self, *args, **kwargs)
        self.setLayout(QVBoxLayout())
        self.canvas = _PlotCanvas(self)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.layout().addWidget(self.toolbar)
        self.layout().addWidget(self.canvas)


class _PlotCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self.fig = Figure()
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self, QSizePolicy.Expanding, QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.plot()
        self.createConn()

        self.figureActive = False
        self.axesActive = None
        self.CPactive = None
        self.pickmode = False
        self.pickMode_changed = True
        self.cpChanged = False
        self.cursorGUI = 'arrow'
        self.cursorChanged = False
        self.CPlist = []
        self.lastIDP = 0

    def plot(self):
        gs0 = self.fig.add_gridspec(1, 2)

        self.ax11 = self.fig.add_subplot(gs0[0], xticks=[], yticks=[], title='Image 1: select Control Points')
        self.ax12 = self.fig.add_subplot(gs0[1], xticks=[], yticks=[], title='Image 2: select Control Points')

        self.ax11.imshow(img1)
        self.ax12.imshow(img2)

    def updateCanvas(self, event=None):
        ax11_xlim = self.ax11.get_xlim()
        ax11_xvis = ax11_xlim[1] - ax11_xlim[0]
        ax12_xlim = self.ax12.get_xlim()
        ax12_xvis = ax12_xlim[1] - ax12_xlim[0]

        while len(self.ax11.patches) > 0:
            [p.remove() for p in self.ax11.patches]
        while len(self.ax12.patches) > 0:
            [p.remove() for p in self.ax12.patches]
        while len(self.ax11.texts) > 0:
            [t.remove() for t in self.ax11.texts]
        while len(self.ax12.texts) > 0:
            [t.remove() for t in self.ax12.texts]

        ax11_units = ax11_xvis * 0.003
        ax12_units = ax12_xvis * 0.003

        for cp in self.CPlist:
            x1 = cp.img1x
            y1 = cp.img1y
            x2 = cp.img2x
            y2 = cp.img2y
            idp = str(cp.idp)

            if x1:
                symb1 = plt.Circle((x1, y1), ax11_units * 8, fill=False, color='red')
                symb2 = plt.Circle((x1, y1), ax11_units * 1, fill=True, color='red')
                self.ax11.text(x1 + ax11_units * 5, y1 + ax11_units * 5, idp)
                self.ax11.add_patch(symb1)
                self.ax11.add_patch(symb2)

            if x2:
                symb1 = plt.Circle((x2, y2), ax12_units * 8, fill=False, color='red')
                symb2 = plt.Circle((x2, y2), ax12_units * 1, fill=True, color='red')
                self.ax12.text(x2 + ax12_units * 5, y2 + ax12_units * 5, idp)
                self.ax12.add_patch(symb1)
                self.ax12.add_patch(symb2)

        self.fig.canvas.draw()

    def createConn(self):
        self.fig.canvas.mpl_connect('figure_enter_event', self.activeFigure)
        self.fig.canvas.mpl_connect('figure_leave_event', self.leftFigure)
        self.fig.canvas.mpl_connect('axes_enter_event', self.activeAxes)
        self.fig.canvas.mpl_connect('button_press_event', self.mouseClicked)
        self.ax11.callbacks.connect('xlim_changed', self.updateCanvas)
        self.ax12.callbacks.connect('xlim_changed', self.updateCanvas)

    def activeFigure(self, event):

        self.figureActive = True
        if self.pickmode and self.cursorGUI != 'cross':
            self.cursorGUI = 'cross'
            self.cursorChanged = True

    def leftFigure(self, event):

        self.figureActive = False
        if self.cursorGUI != 'arrow':
            self.cursorGUI = 'arrow'
            self.cursorChanged = True

    def activeAxes(self, event):
        self.axesActive = event.inaxes

    def mouseClicked(self, event):
        x = event.xdata
        y = event.ydata

        if self.toolbar.mode != '':
            self.pickmode = False

        if self.pickmode and (event.inaxes == self.ax11 or event.inaxes == self.ax12):

            if self.CPactive and not self.CPactive.status_complete:
                self.CPactive.appendCoord(x, y)
                self.cpChanged = True
            else:
                idp = self.lastIDP + 1
                cp = _ControlPoint(idp, x, y, self)
                self.CPlist.append(cp)
                self.cpChanged = True
                self.lastIDP += 1

            self.updateCanvas()


class _ControlPoint:
    def __init__(self, idp, x, y, other):
        self.img1x = None
        self.img1y = None
        self.img2x = None
        self.img2y = None
        self.status_complete = False
        self.idp = idp

        self.mn = other
        self.mn.CPactive = self

        self.appendCoord(x, y)

    def appendCoord(self, x, y):

        if self.mn.axesActive == self.mn.ax11 and self.img1x is None:
            self.img1x = x
            self.img1y = y
        elif self.mn.axesActive == self.mn.ax12 and self.img2x is None:
            self.img2x = x
            self.img2y = y

        else:
            raise Exception("Please, select the control point in the other image")

        if self.img1x and self.img2x:
            self.status_complete = True
            self.mn.cpActive = None

    @property
    def coord(self):
        return self.idp, self.img1x, self.img1y, self.img2x, self.img2y

    @property
    def coordText(self):
        if self.img1x and not self.img2x:
            return str(round(self.idp, 2)), str(round(self.img1x, 2)), str(round(self.img1y, 2)), '', ''
        elif not self.img1x and self.img2x:
            return str(round(self.idp, 2)), '', '', str(round(self.img2x, 2)), str(round(self.img2y, 2))
        else:
            return str(round(self.idp, 2)), str(round(self.img1x, 2)), str(round(self.img1y, 2)), str(
                round(self.img2x, 2)), str(round(self.img2y, 2))

    def __str__(self):
        return f"CP {self.idp}: {self.coord}"

    @property
    def getdict(self):
        dict = {
            'point_id': self.idp,
            'img1_x': self.img1x,
            'img1_y': self.img1y,
            'img2_x': self.img2x,
            'img2_y': self.img2y,
        }

        return dict
    

    @property
    def getarray(self):
        return [[self.img1x, self.img1y], [self.img2x, self.img2y]]