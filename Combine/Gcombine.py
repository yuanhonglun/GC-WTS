import os
import sys

import pandas as pd
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, QThread, Qt, QCoreApplication
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QApplication, QDesktopWidget
from PyQt5.QtGui import QGuiApplication
from combine import Ui_MainWindow
from combine_rt_msp_final import CombineRtMsp
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

        # try:
        #     mymainwindow.Main(self.all_msp_path, self.all_rt_path, self.RI_path, self.RIalertmin, self.RIalertmax,
        #                       self.RI_threshold_value, self.ri_window_scale,
        #                       self.RTmin, self.RTmax, self.RImin, self.RImax, self.check_RT, self.check_latin,
        #                       self.out_path,self.use_unknown, self.unknow_msp_path, self.unknow_rt_path, self.rt_window_unknown, self.similarity_score_threshold_unknown)
        #     self.finished.emit(0)
        # except:
        #     self.finished.emit(1)

        mymainwindow.Main(self.all_msp_path, self.all_rt_path, self.RI_path, self.RIalertmin, self.RIalertmax,
                          self.RI_threshold_value, self.ri_window_scale,
                          self.RTmin, self.RTmax, self.RImin, self.RImax, self.check_RT, self.check_latin,
                          self.out_path, self.use_unknown, self.unknow_msp_path, self.unknow_rt_path,
                          self.rt_window_unknown, self.similarity_score_threshold_unknown)
        self.finished.emit(0)


class MyMainWindow(QMainWindow, Ui_MainWindow, CombineRtMsp):  # 继承 QMainWindow类和 Ui_MainWindow界面类

    def __init__(self, parent=None):
        super(MyMainWindow, self).__init__(parent)  # 初始化父类
        self.setupUi(self)  # 继承 Ui_MainWindow 界面类
        self.progressBar.setRange(0, 0)
        self.progressBar.hide()
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
        self.RTmax = 68.8
        self.RTmin = 0
        self.RIalertmin = 600
        self.RIalertmax = 2000
        self.RImax = 3000
        self.RImin = 0
        self.RI_threshold_value = 40
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
        '''
            打开全部msp所在的文件目录
        '''
        # MyMainWindow.all_msp_path = QFileDialog.getExistingDirectory(None, "选择msp文件路径", './')
        msp_file, _ = QFileDialog.getOpenFileNames(self, 'Open MSP File', './', 'MSP files (*.msp)')
        if msp_file:
            self.msp_input.extend(msp_file)
            self.listWidget.setHidden(False)
            self.listWidget.addItem(os.path.basename(msp_file[0]))

    def open_file2(self):
        '''
            打开全部RT所在的文件目录
        '''
        # MyMainWindow.all_rt_path = QFileDialog.getExistingDirectory(None, "选择msp文件路径", './')
        rt_file, _ = QFileDialog.getOpenFileNames(self, 'Open RT File', './', 'All files (*.xlsx)')
        if rt_file:
            self.rt_input.extend(rt_file)
            self.listWidget_2.setHidden(False)
            self.listWidget_2.addItem(os.path.basename(rt_file[0]))

    def open_file3(self):
        '''
            输入RI文件
        '''
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
        self.unknow_rt_path, _ = QFileDialog.getOpenFileName(self, 'Open RT File', './', 'All files (*.xlsx)')

    def open_unknown_msp(self):
        self.unknow_msp_path, _ = QFileDialog.getOpenFileName(self, 'Open MSP File', './', 'MSP files (*.msp)')

    def rt_window_unknown(self, value):
        self.rt_window_unknown = value

    def similarity_score_threshold_unknown(self, value):
        self.similarity_score_threshold_unknown = value

    def check1(self, bool):
        '''
            是否导出理论RT
        '''
        self.check_RT = bool

    def check2(self, bool):
        '''
            是否替换拉丁文i
        '''
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
            'The larger this value, the larger the RI error threshold when the RI is larger. Setting it to 0 disables this feature. default: 5      ')

    def tooltip8(self):
        QMessageBox.about(
            None,
            'Help',
            'e. g. “.beta.” to “beta”       ')

    def tooltip9(self):
        QMessageBox.about(
            None,
            'Help',
            'Target compound are compared for similarity with substances within the user-defined RT window when selecting qualitative ions      ')

    def tooltip10(self):
        QMessageBox.about(
            None,
            'Help',
            'two spectra are considered distinguishable when similarity score below the threshold. default: 0.85')

    def run(self):
        self.frame_5.setHidden(False)
        self.progressBar.show()

        self.pushButton_4.setEnabled(False)

        self.worker_thread = WorkerThread(self.msp_input, self.rt_input, self.RI_path, self.RIalertmin, self.RIalertmax,
                                          self.RI_threshold_value, self.ri_window_scale, self.RTmin,
                                          self.RTmax, self.RImin, self.RImax, self.check_RT, self.check_latin,
                                          self.out_path,
                                          self.use_unknown, self.unknow_msp_path, self.unknow_rt_path,
                                          self.rt_window_unknown, self.similarity_score_threshold_unknown)
        self.worker_thread.finished.connect(self.hide_progress_bar)

        self.worker_thread.start()

    def hide_progress_bar(self, run_stat):
        self.frame_5.setHidden(True)
        self.progressBar.hide()
        self.pushButton_4.setEnabled(True)
        self.pushButton_5.setEnabled(True)
        if run_stat == 0:
            QMessageBox.information(self, 'Success', 'Run successfully!')
        elif run_stat == 1:
            QMessageBox.critical(self, 'Error', 'Run Failed!')

    def center(self):
        # 获取屏幕的分辨率信息
        screen = QDesktopWidget().availableGeometry()

        # 获取窗口的大小
        size = self.geometry()
        print(screen)
        # 计算窗口的左上角位置
        x = (screen.width() - size.width()) // 2
        y = (screen.height() - size.height()) // 2

        # 如果分辨率太小，窗口显示在左上角
        if screen.width() < size.width() or screen.height() < size.height():
            x = 0
            y = 0

        # 移动窗口到指定位置
        self.move(x, y)


if __name__ == '__main__':
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication(sys.argv)  # 在 QApplication 方法中使用，创建应用程序对象
    myWin = MyMainWindow()  # 实例化 MyMainWindow 类，创建主窗口
    myWin.center()
    apply_stylesheet(app, theme='light_cyan_500.xml', invert_secondary=True)
    myWin.show()  # 在桌面显示控件 myWin
    sys.exit(app.exec_())  # 结束进程，退出程序
