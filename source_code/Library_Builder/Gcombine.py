import os
import sys

import pandas as pd
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, QThread, Qt, QCoreApplication, QTimer
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QApplication, QDesktopWidget
from PyQt5.QtGui import QGuiApplication
from combine import Ui_MainWindow
from Library_builder import CombineRtMsp
from qt_material import apply_stylesheet
import images_rc


class WorkerThread(QThread):
    finished = pyqtSignal(int)

    def __init__(self, all_msp_path, all_rt_path, RI_path, RIalertmin, RIalertmax,
                 RI_threshold_value, ri_window_scale, RTmin, RTmax, RImin, RImax, check_RT, check_latin, out_path,
                 use_unknown, unknow_msp_path, unknow_rt_path, rt_window_unknown, similarity_score_threshold_unknown):
        super().__init__()
        self.all_msp_path = all_msp_path
        self.all_rt_path = all_rt_path
        self.RI_path = RI_path
        self.RIalertmin = RIalertmin
        self.RIalertmax = RIalertmax
        self.RI_threshold_value = RI_threshold_value
        self.ri_window_scale = ri_window_scale
        self.RTmin = RTmin
        self.RTmax = RTmax
        self.RImin = RImin
        self.RImax = RImax
        self.check_RT = check_RT
        self.check_latin = check_latin
        self.out_path = out_path
        self.use_unknown = use_unknown
        self.unknow_msp_path = unknow_msp_path
        self.unknow_rt_path = unknow_rt_path
        self.rt_window_unknown = rt_window_unknown
        self.similarity_score_threshold_unknown = similarity_score_threshold_unknown


    def run(self):
        mymainwindow = MyMainWindow()

        try:
            mymainwindow.Main(self.all_msp_path, self.all_rt_path, self.RI_path, self.RIalertmin, self.RIalertmax,
                              self.RI_threshold_value, self.ri_window_scale,
                              self.RTmin, self.RTmax, self.RImin, self.RImax, self.check_RT, self.check_latin,
                              self.out_path, self.use_unknown, self.unknow_msp_path, self.unknow_rt_path,
                              self.rt_window_unknown, self.similarity_score_threshold_unknown)
            self.finished.emit(0)
        except:
            self.finished.emit(1)




class MyMainWindow(QMainWindow, Ui_MainWindow, CombineRtMsp):

    def __init__(self, parent=None):
        super(MyMainWindow, self).__init__(parent)
        self.setupUi(self)
        self.progressBar.setRange(0, 100)
        self.progress_value = 0
        self.progressBar.hide()
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_progress_bar)
        self.frame_5.setHidden(True)
        self.listWidget.setHidden(True)
        self.label_22.setHidden(True)
        self.label_21.setHidden(True)
        self.pushButton_16.setHidden(True)
        self.pushButton_17.setHidden(True)
        self.frame_19.setHidden(True)
        self.listWidget.itemDoubleClicked.connect(self.delete_item)
        self.listWidget_2.setHidden(True)
        self.listWidget_2.itemDoubleClicked.connect(self.delete_item_2)
        self.listWidget_3.setHidden(True)
        self.listWidget_3.itemDoubleClicked.connect(self.delete_item_3)
        self.listWidget_4.setHidden(True)
        self.listWidget_4.itemDoubleClicked.connect(self.delete_item_4)
        self.msp_input = []
        self.rt_input = []
        self.RI_path = ''
        self.out_path = ''
        self.RTmax = 69.0
        self.RTmin = 0
        self.RIalertmin = 600
        self.RIalertmax = 2000
        self.RImax = 3000
        self.RImin = 0
        self.RI_threshold_value = 100
        self.check_latin = False
        self.check_RT = False
        self.use_unknown = False
        self.ri_window_scale = 5
        self.rt_window_unknown = 2
        self.similarity_score_threshold_unknown = 0.85
        self.unknow_msp_path = ''
        self.unknow_rt_path = ''
        self.rt_window_unknown = 2
        self.similarity_score_threshold_unknown = 0.85

    def open_file1(self):

        msp_file, _ = QFileDialog.getOpenFileNames(self, 'Open MSP File', './', 'MSP files (*.msp)')

        if msp_file:
            self.msp_input.extend(msp_file)
            self.listWidget.setHidden(False)
            for i in self.msp_input:
                self.listWidget.addItem(os.path.basename(i))

    def open_file2(self):

        rt_file, _ = QFileDialog.getOpenFileNames(self, 'Open RT File', './', 'All files (*.csv)')
        if rt_file:
            self.rt_input.extend(rt_file)
            self.listWidget_2.setHidden(False)
            for i in self.rt_input:
                self.listWidget_2.addItem(os.path.basename(i))


    def open_file3(self):

        self.RI_path, _ = QFileDialog.getOpenFileName(self, 'Open RI File', './', 'All files (*.csv)')
        self.listWidget_3.setHidden(False)
        self.listWidget_3.addItem(os.path.basename(self.RI_path))

    def is_use_unknown(self, bool):

        self.use_unknown = bool
        if self.use_unknown:
            self.label_22.setHidden(False)
            self.label_21.setHidden(False)
            self.pushButton_16.setHidden(False)
            self.pushButton_17.setHidden(False)
            self.frame_19.setHidden(False)
        else:
            self.label_22.setHidden(True)
            self.label_21.setHidden(True)
            self.pushButton_16.setHidden(True)
            self.pushButton_17.setHidden(True)
            self.frame_19.setHidden(True)

    def open_unknown_rt(self):
        self.unknow_rt_path, _ = QFileDialog.getOpenFileName(self, 'Open RT File', './', 'All files (*.csv)')

    def open_unknown_msp(self):
        self.unknow_msp_path, _ = QFileDialog.getOpenFileName(self, 'Open MSP File', './', 'MSP files (*.msp)')

    def rt_window_unknown(self, value):
        self.rt_window_unknown = value

    def similarity_score_threshold_unknown(self, value):
        self.similarity_score_threshold_unknown = value

    def check1(self, bool):

        self.check_RT = bool

    def check2(self, bool):

        self.check_latin = bool

    def chooseRTmin(self, value):

        self.RTmin = value

    def chooseRTmax(self, value):
        self.RTmax = value

    def chooseRIalertmin(self, value):
        self.RIalertmin = value

    def chooseRIalertmax(self, value):
        self.RIalertmax = value

    def chooseRIerror(self, value):
        self.RI_threshold_value = value

    def chooseRImin(self, value):
        self.RImin = value

    def chooseRImax(self, value):
        self.RImax = value

    def choose_ri_window(self, value):
        self.ri_window_scale = value

    def save(self):
        self.out_path = QFileDialog.getExistingDirectory(None, "Choose Output Path", './')
        self.listWidget_4.setHidden(False)
        self.listWidget_4.addItem(self.out_path)

    def delete_item(self, item):
        row = self.listWidget.row(item)
        if row != 0:
            self.listWidget.takeItem(row)
            del self.msp_input[row - 1]

    def delete_item_2(self, item):
        row = self.listWidget_2.row(item)
        if row != 0:
            self.listWidget_2.takeItem(row)
            del self.rt_input[row - 1]

    def delete_item_3(self, item):
        row = self.listWidget_3.row(item)
        if row != 0:
            self.listWidget_3.takeItem(row)
            self.RI_path = ''

    def delete_item_4(self, item):
        row = self.listWidget_4.row(item)
        if row != 0:
            self.listWidget_4.takeItem(row)
            self.out_path = ''

    def ridemo(self):

        ri_demo_example = pd.DataFrame(
            {'RI': [427, 517, 600, 700, 800],
             'RT (min)': [1.575, 1.584, 2.036, 3.319, 5.199]})

        try:
            ri_demo_example.to_csv('RI_demo.csv', index=False)
            QMessageBox.about(
                None,
                'Help',
                'An RI sample file has been stored in the same directory as the software')
        except:
            pass

    def rtdemo(self):

        rt_demo_example = pd.DataFrame(
            {'Name': ['Ethylene oxide', 'Acetonitrile', 'Ethanethiol', 'Ethanethiol', '2-Propen-1-ol'],
             'RT': [1.486, 1.648, 1.683, 1.8055875, 1.825]})

        try:
            rt_demo_example.to_excel('RT_demo.xlsx', index=False)
            QMessageBox.about(
                None,
                'Help',
                'An RT sample file has been stored in the same directory as the software')
        except:
            pass

    def tooltip1(self):
        QMessageBox.about(
            None,
            'Help',
            'Import a list of compounds and their RT        ')

    def tooltip2(self):
        QMessageBox.about(
            None,
            'Help',
            'Import RI calibration data     ')

    def tooltip3(self):
        QMessageBox.about(
            None,
            'Help',
            'Enter the allowed minimum and maximum RT       ')

    def tooltip4(self):
        QMessageBox.about(
            None,
            'Help',
            'Enter the allowed minimum and maximum RI. default: 0 to 3000       ')

    def tooltip5(self):
        QMessageBox.about(
            None,
            'Help',
            'Enter the RI range, in which compounds with user-defined differences between the measured RI and the library RI will be marked. default: 600 to 2000       ')

    def tooltip6(self):
        QMessageBox.about(
            None,
            'Help',
            'Compound will be marked if the difference between the measured RI and the database RI exceeds this threshold. default: 100     ')

    def tooltip7(self):
        QMessageBox.about(
            None,
            'Help',
            'The larger this value, the larger the RI alerting threshold when the RI is larger. Setting it to 0 disables this feature. default: 5      ')

    def tooltip8(self):
        QMessageBox.about(
            None,
            'Help',
            'e. g. replace “.beta.” with “beta”       ')

    def tooltip9(self):
        QMessageBox.about(
            None,
            'Help',
            'Unknowns are compared for similarity with compounds within the user-defined RT window when selecting qualitative ions      ')

    def tooltip10(self):
        QMessageBox.about(
            None,
            'Help',
            'Only the unknowns with spectral similarity scores lower than the threshold are incorporated into the library. default: 0.85')

    def run(self):

        if self.out_path == '':
            QMessageBox.critical(
                None,
                'Error',
                'Please select the output path！')
        else:
            self.frame_5.setHidden(False)
            self.progressBar.show()
            self.progress_value = 0
            self.timer.start(100)

            self.pushButton_4.setEnabled(False)

            self.worker_thread = WorkerThread(self.msp_input, self.rt_input, self.RI_path, self.RIalertmin, self.RIalertmax,
                                              self.RI_threshold_value, self.ri_window_scale, self.RTmin,
                                              self.RTmax, self.RImin, self.RImax, self.check_RT, self.check_latin,
                                              self.out_path,
                                              self.use_unknown, self.unknow_msp_path, self.unknow_rt_path,
                                              self.rt_window_unknown, self.similarity_score_threshold_unknown)
            self.worker_thread.finished.connect(self.hide_progress_bar)

            self.worker_thread.start()

    def update_progress_bar(self):
        self.progress_value += 1
        if self.progress_value == 101:
            self.progress_value = 0
        else:
            self.progressBar.setValue(self.progress_value)


    def hide_progress_bar(self, run_stat):
        self.frame_5.setHidden(True)
        self.progressBar.hide()
        self.timer.stop()
        self.pushButton_4.setEnabled(True)
        self.pushButton_5.setEnabled(True)
        if run_stat == 0:
            QMessageBox.information(self, 'Success', 'Run successfully!')
        elif run_stat == 1:
            QMessageBox.critical(self, 'Error', 'Run Failed!')

    def center(self):

        screen = QDesktopWidget().availableGeometry()
        size = self.geometry()

        x = (screen.width() - size.width()) // 2
        y = (screen.height() - size.height()) // 2

        if screen.width() < size.width() or screen.height() < size.height():
            x = 0
            y = 0

        self.move(x, y)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)
    myWin = MyMainWindow()
    myWin.center()
    apply_stylesheet(app, theme='light_cyan_500.xml', invert_secondary=True)
    myWin.show()
    sys.exit(app.exec_())
