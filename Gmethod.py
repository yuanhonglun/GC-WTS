import os
import sys
from PyQt5.QtCore import pyqtSignal, QThread
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QWidget, QTextEdit, QPushButton, QVBoxLayout, \
    QApplication
from method import Ui_MainWindow
from get_collect_method_final import GetMethod
from qt_material import apply_stylesheet
import pandas as pd

class WorkerThread(QThread):
    finished = pyqtSignal(int)

    def __init__(self, msp_path, rt_data_path,set_name_list, name_list_path,
                      mz_min, mz_max, outpath, rt_window,
                      min_ion_intensity_percent, min_ion_num,  prefer_mz_threshold, similarity_threshold, fr_factor,
                      retention_time_max, solvent_delay, sim_sig_max,
                      min_dwell_time,min_point_per_s,min_point_per_s_limit, convert_to_ag_method):
        super().__init__()
        self.msp_path = msp_path
        self.rt_data_path = rt_data_path
        self.set_name_list = set_name_list
        self.name_list_path = name_list_path
        self.mz_min = mz_min
        self.mz_max = mz_max
        self.outpath = outpath
        self.rt_window = rt_window
        self.min_ion_intensity_percent = min_ion_intensity_percent
        self.min_ion_num = min_ion_num
        self.prefer_mz_threshold = prefer_mz_threshold
        self.similarity_threshold = similarity_threshold
        self.fr_factor = fr_factor
        self.retention_time_max = retention_time_max
        self.solvent_delay = solvent_delay
        self.sim_sig_max = sim_sig_max
        self.min_dwell_time = min_dwell_time
        self.min_point_per_s = min_point_per_s
        self.min_point_per_s_limit = min_point_per_s_limit
        self.convert_to_ag_method = convert_to_ag_method

        self.mymainwindow = MyMainWindow()



    def run(self):

        try:
            self.mymainwindow.Main(self.msp_path, self.rt_data_path, self.set_name_list, self.name_list_path,
                              self.mz_min, self.mz_max, self.outpath, self.rt_window,
                              self.min_ion_intensity_percent, self.min_ion_num, self.prefer_mz_threshold, self.similarity_threshold, self.fr_factor,
                              self.retention_time_max, self.solvent_delay, self.sim_sig_max,
                              self.min_dwell_time, self.min_point_per_s, self.min_point_per_s_limit, self.convert_to_ag_method)
            self.finished.emit(0)
        except:
            self.finished.emit(1)



class MyMainWindow(QMainWindow, Ui_MainWindow, GetMethod):  # 继承 QMainWindow类和 Ui_MainWindow界面类


    def __init__(self, parent=None):
        super(MyMainWindow, self).__init__(parent)  # 初始化父类

        self.setupUi(self)  # 继承 Ui_MainWindow 界面类
        self.progressBar.setRange(0, 0)
        self.progressBar.hide()
        self.listWidget.setHidden(True)
        self.listWidget.itemDoubleClicked.connect(self.delete_item)
        self.listWidget_2.setHidden(True)
        self.listWidget_2.itemDoubleClicked.connect(self.delete_item_2)
        self.pushButton_3.setHidden(True)
        self.listWidget_3.setHidden(True)
        self.listWidget_3.itemDoubleClicked.connect(self.delete_item_3)
        self.listWidget_4.setHidden(True)
        self.listWidget_4.itemDoubleClicked.connect(self.delete_item_4)



        self.msp_path = ''
        self.rt_data_path = ''
        self.set_name_list = False
        self.name_list_path = ''
        self.mz_min = 35
        self.mz_max = 400
        self.outpath = ''
        self.rt_window = 0.5
        self.min_ion_intensity_percent = 7
        self.min_ion_num = 2
        self.retention_time_max = 68.8
        self.sim_sig_max = 99
        self.min_dwell_time = 10
        self.min_point_per_s = 2
        self.min_point_per_s_limit = 0.5
        self.prefer_mz_threshold = 60
        self.similarity_threshold = 0.85
        self.fr_factor = 2
        self.solvent_delay = 0
        self.convert_to_ag_method = False



    def open1(self):
        '''
            导入msp文件
        '''

        self.listWidget.takeItem(1)
        self.msp_path, _ = QFileDialog.getOpenFileNames(self, 'Open MSP File', './', 'MSP files (*.msp)')
        self.msp_path = self.msp_path[0]
        self.listWidget.setHidden(False)
        self.listWidget.addItem(os.path.basename(self.msp_path))


    def open2(self):
        '''
            输入RT文件
        '''
        self.listWidget_2.takeItem(1)
        self.rt_data_path, _ = QFileDialog.getOpenFileNames(self, 'Open RT File', './', 'All files (*.xlsx)')
        self.rt_data_path = self.rt_data_path[0]
        self.listWidget_2.setHidden(False)
        self.listWidget_2.addItem(os.path.basename(self.rt_data_path))

    def open3(self):
        '''
            输入物质名列表文件
        '''
        self.listWidget_3.takeItem(1)
        self.name_list_path, _ = QFileDialog.getOpenFileNames(self, 'Open Metabolite File', './', 'All files (*.txt)')
        self.name_list_path = self.name_list_path[0]
        self.listWidget_3.setHidden(False)
        self.listWidget_3.addItem(os.path.basename(self.name_list_path))

    def save(self):
        self.listWidget_4.takeItem(1)
        self.outpath = QFileDialog.getExistingDirectory(None, "Choose Output Path", './')
        self.listWidget_4.setHidden(False)
        self.listWidget_4.addItem(self.outpath)


    def ismateriallist(self, bool):
        '''
            是否导入物质名列表
        '''
        self.set_name_list = bool
        if self.set_name_list:
            self.pushButton_3.setHidden(False)
        else:
            self.pushButton_3.setHidden(True)
            self.listWidget_3.setHidden(True)


    def isanjielun(self, bool):
        '''
            是否输入为安捷伦方法
        '''
        self.convert_to_ag_method = bool


    def mzmin(self, value):
        '''
            m/z范围0-500， 默认35
        :param value:
        :return:
        '''
        self.mz_min = value


    def mzmax(self, value):
        '''
            m/z范围0-500， 默认400
        :param value:
        :return:
        '''
        self.mz_max = value


    def maxrt(self, value):
        '''
            最大保留时间，大于0， 默认空
        :param value:
        :return:
        '''
        self.retention_time_max = value

    def rtwindow(self, value):
        '''
            RT窗口范围，大于0小于rt_max， 默认0.5
        :param value:
        :return:
        '''
        self.rt_window = value

    def ionmax(self, value):
        '''
            离子响应阈值，0-100，默认7
        :param value:
        :return:
        '''
        self.min_ion_intensity_percent = value

    def ionnum(self, value):
        '''
            最少离子数，大于2，默认2
        :param value:
        :return:
        '''
        self.min_ion_num = value

    def simnum(self, value):
        '''
            SIM最多段数，大于0，默认99
        :param value:
        :return:
        '''
        self.sim_sig_max = value


    def simtime(self, value):
        '''
            SIM最小驻留时间，1-1000，默认10
        :param value:
        :return:
        '''
        self.min_dwell_time = value

    def simmincount(self, value):
        '''
            SIM最少打点次数，大于sim_ionmincount，默认2， 过大易报错， 不建议超过5
        :param value:
        :return:
        '''
        self.min_point_per_s = value

    def simionmincount(self, value):
        '''
            SIM离子密集处最少打点次数，小于sim_mincount，默认0.5, 过大易报错，不建议超过2
        :param value:
        :return:
        '''
        self.min_point_per_s_limit = value



    def get_prefer_mz_threshold(self, value):
        self.prefer_mz_threshold = value


    def get_similarity_threshold(self, value):
        self.similarity_threshold = value


    def get_fr_factor(self, value):
        self.fr_factor = value

    def get_solvent_delay(self, value):
        self.solvent_delay = value

    def delete_item(self, item):
        row = self.listWidget.row(item)
        if row != 0:
            self.listWidget.takeItem(row)
            self.msp_path = None

    def delete_item_2(self, item):
        row = self.listWidget_2.row(item)
        if row != 0:
            self.listWidget_2.takeItem(row)
            self.rt_data_path = None

    def delete_item_3(self, item):
        row = self.listWidget_3.row(item)
        if row != 0:
            self.listWidget_3.takeItem(row)
            self.name_list_path = None

    def delete_item_4(self, item):
        row = self.listWidget_4.row(item)
        if row != 0:
            self.listWidget_4.takeItem(row)
            self.outpath = None

    def rtdemo(self):
        rt_demo_example = pd.DataFrame(
            {'Name': ['Ethylene oxide', 'Acetonitrile', 'Ethanethiol', 'Ethanethiol', '2-Propen-1-ol'],
             'RT': [1.486, 1.648, 1.683, 1.8055875, 1.825]})

        try:
            rt_demo_example.to_excel('RT_demo.xlsx', index=False)
            QMessageBox.about(
                None,
                'Help',
                '一个RT示例文件存储在软件同目录下')
        except:
            pass

    def namelistdemo(self):
        namelist_demo_example = pd.DataFrame(
            {'Name': ['3-Heptanone', 'Acetonitrile', 'Ethanethiol', 'Ethanethiol', '2-Propen-1-ol']})

        try:
            namelist_demo_example.to_csv('Namelist_demo.txt', sep="\t", index=False)
            QMessageBox.about(
                None,
                'Help',
                '一个namelist示例文件存储在软件同目录下')
        except:
            pass

    def tooltip1(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip2(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip3(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip4(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip5(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip6(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip7(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip8(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip9(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip10(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip11(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip12(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def tooltip13(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')


    def tooltip14(self):
        QMessageBox.about(
            None,
            'Help',
            'ssssss')

    def run(self):

        self.progressBar.show()

        self.pushButton_4.setEnabled(False)
        print("msp_path", self.msp_path)
        print("rt_data_path", self.rt_data_path)
        print("set_name_list", self.set_name_list)
        print("name_list_path", self.name_list_path)
        print("mz_min", self.mz_min)
        print("mz_max", self.mz_max)
        print("outpath", self.outpath)
        print("rt_window", self.rt_window)
        print("min_ion_intensity_percent", self.min_ion_intensity_percent)
        print("min_ion_num", self.min_ion_num)
        print("retention_time_max", self.retention_time_max)
        print("sim_sig_max", self.sim_sig_max)
        print("min_dwell_time", self.min_dwell_time)
        print("min_point_per_s_limit", self.min_point_per_s_limit)
        print("convert_to_ag_method", self.convert_to_ag_method)

        self.worker_thread = WorkerThread(self.msp_path, self.rt_data_path, self.set_name_list, self.name_list_path,
                          self.mz_min, self.mz_max, self.outpath, self.rt_window,
                          self.min_ion_intensity_percent, self.min_ion_num, self.prefer_mz_threshold, self.similarity_threshold, self.fr_factor,
                          self.retention_time_max, self.solvent_delay, self.sim_sig_max,
                          self.min_dwell_time, self.min_point_per_s, self.min_point_per_s_limit, self.convert_to_ag_method)
        self.worker_thread.finished.connect(self.hide_progress_bar)

        # self.log_window = ChildWindow()
        # self.log_window.show()
        #

        self.worker_thread.start()




    def hide_progress_bar(self, run_stat):
        self.progressBar.hide()
        self.pushButton_4.setEnabled(True)

        if run_stat == 0:
            QMessageBox.information(self, 'Success', 'Run successfully!')
        elif run_stat == 1:
            QMessageBox.critical(self, 'Error', 'Run Failed!')


# class ChildWindow(QWidget):
#
#     def __init__(self):
#         super().__init__()
#         self.init_ui()
#
#     def init_ui(self):
#         self.setWindowTitle('Log Window')
#         self.resize(400, 300)
#
#         self.log_text_edit = QTextEdit(self)
#         self.log_text_edit.setReadOnly(True)
#
#         clear_button = QPushButton('Clear', self)
#         clear_button.clicked.connect(self.log_text_edit.clear)
#
#         layout = QVBoxLayout(self)
#         layout.addWidget(self.log_text_edit)
#         layout.addWidget(clear_button)
#
#     def add_log(self, log):
#         self.log_text_edit.append(log)
#






if __name__ == '__main__':
    app = QApplication(sys.argv)  # 在 QApplication 方法中使用，创建应用程序对象
    myWin = MyMainWindow()  # 实例化 MyMainWindow 类，创建主窗口
    apply_stylesheet(app, theme='light_cyan_500.xml', invert_secondary=True)
    myWin.show()  # 在桌面显示控件 myWin
    sys.exit(app.exec_())  # 结束进程，退出程序
