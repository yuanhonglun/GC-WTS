# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'step_3_RT.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form_RT(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(857, 654)
        self.gridLayout = QtWidgets.QGridLayout(Form)
        self.gridLayout.setObjectName("gridLayout")
        self.frame_3 = QtWidgets.QFrame(Form)
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame_3)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(self.frame_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(22)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label.setFont(font)
        self.label.setStyleSheet("font: 22pt;")
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        self.gridLayout.addWidget(self.frame_3, 0, 0, 1, 2)
        self.frame_2 = QtWidgets.QFrame(Form)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.frame_2)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.frame_6 = QtWidgets.QFrame(self.frame_2)
        self.frame_6.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_6.setObjectName("frame_6")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.frame_6)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.label_6 = QtWidgets.QLabel(self.frame_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy)
        self.label_6.setMinimumSize(QtCore.QSize(0, 0))
        self.label_6.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_6.setFont(font)
        self.label_6.setStyleSheet("font: 18pt;")
        self.label_6.setObjectName("label_6")
        self.gridLayout_6.addWidget(self.label_6, 0, 0, 1, 1)
        self.radioButton = QtWidgets.QRadioButton(self.frame_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton.sizePolicy().hasHeightForWidth())
        self.radioButton.setSizePolicy(sizePolicy)
        self.radioButton.setStyleSheet("font: 18pt;")
        self.radioButton.setText("")
        self.radioButton.setChecked(False)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout_6.addWidget(self.radioButton, 0, 1, 1, 1, QtCore.Qt.AlignHCenter)
        self.gridLayout_3.addWidget(self.frame_6, 1, 0, 1, 1)
        self.widget_3 = QtWidgets.QWidget(self.frame_2)
        self.widget_3.setObjectName("widget_3")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.widget_3)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_12 = QtWidgets.QLabel(self.widget_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(20)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_12.setFont(font)
        self.label_12.setStyleSheet("font: 20pt;")
        self.label_12.setObjectName("label_12")
        self.horizontalLayout_5.addWidget(self.label_12)
        self.gridLayout_3.addWidget(self.widget_3, 0, 0, 1, 1)
        self.frame_5 = QtWidgets.QFrame(self.frame_2)
        self.frame_5.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_5.setObjectName("frame_5")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.frame_5)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_7 = QtWidgets.QLabel(self.frame_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_7.sizePolicy().hasHeightForWidth())
        self.label_7.setSizePolicy(sizePolicy)
        self.label_7.setMinimumSize(QtCore.QSize(0, 30))
        self.label_7.setMaximumSize(QtCore.QSize(16777215, 30))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_7.setFont(font)
        self.label_7.setStyleSheet("font: 18pt;")
        self.label_7.setObjectName("label_7")
        self.gridLayout_4.addWidget(self.label_7, 0, 0, 1, 1)
        self.doubleSpinBox_12 = QtWidgets.QDoubleSpinBox(self.frame_5)
        self.doubleSpinBox_12.setMinimumSize(QtCore.QSize(150, 30))
        self.doubleSpinBox_12.setMaximumSize(QtCore.QSize(150, 30))
        self.doubleSpinBox_12.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_12.setMaximum(9999.0)
        self.doubleSpinBox_12.setSingleStep(0.1)
        self.doubleSpinBox_12.setProperty("value", 0.3)
        self.doubleSpinBox_12.setObjectName("doubleSpinBox_12")
        self.gridLayout_4.addWidget(self.doubleSpinBox_12, 0, 1, 1, 1)
        self.pushButton_6 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_6.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_6.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_6.setText("")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("C:/Users/86724/Desktop/R-C.jpg"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.pushButton_6.setIcon(icon)
        self.pushButton_6.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_6.setObjectName("pushButton_6")
        self.gridLayout_4.addWidget(self.pushButton_6, 0, 2, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.frame_5)
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_8.setFont(font)
        self.label_8.setStyleSheet("font: 18pt;")
        self.label_8.setWordWrap(True)
        self.label_8.setObjectName("label_8")
        self.gridLayout_4.addWidget(self.label_8, 1, 0, 1, 1)
        self.doubleSpinBox_7 = QtWidgets.QDoubleSpinBox(self.frame_5)
        self.doubleSpinBox_7.setMinimumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_7.setMaximumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_7.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_7.setDecimals(2)
        self.doubleSpinBox_7.setSingleStep(0.01)
        self.doubleSpinBox_7.setProperty("value", 0.05)
        self.doubleSpinBox_7.setObjectName("doubleSpinBox_7")
        self.gridLayout_4.addWidget(self.doubleSpinBox_7, 1, 1, 1, 1)
        self.pushButton_7 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_7.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_7.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_7.setText("")
        self.pushButton_7.setIcon(icon)
        self.pushButton_7.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_7.setObjectName("pushButton_7")
        self.gridLayout_4.addWidget(self.pushButton_7, 1, 2, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.frame_5)
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_10.setFont(font)
        self.label_10.setStyleSheet("font: 18pt;")
        self.label_10.setWordWrap(True)
        self.label_10.setObjectName("label_10")
        self.gridLayout_4.addWidget(self.label_10, 2, 0, 1, 1)
        self.doubleSpinBox_8 = QtWidgets.QDoubleSpinBox(self.frame_5)
        self.doubleSpinBox_8.setMinimumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_8.setMaximumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_8.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_8.setDecimals(2)
        self.doubleSpinBox_8.setMaximum(1.0)
        self.doubleSpinBox_8.setSingleStep(0.01)
        self.doubleSpinBox_8.setProperty("value", 0.1)
        self.doubleSpinBox_8.setObjectName("doubleSpinBox_8")
        self.gridLayout_4.addWidget(self.doubleSpinBox_8, 2, 1, 1, 1)
        self.pushButton_8 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_8.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_8.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_8.setText("")
        self.pushButton_8.setIcon(icon)
        self.pushButton_8.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_8.setObjectName("pushButton_8")
        self.gridLayout_4.addWidget(self.pushButton_8, 2, 2, 1, 1)
        self.label_11 = QtWidgets.QLabel(self.frame_5)
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_11.setFont(font)
        self.label_11.setStyleSheet("font: 18pt;")
        self.label_11.setWordWrap(True)
        self.label_11.setObjectName("label_11")
        self.gridLayout_4.addWidget(self.label_11, 3, 0, 1, 1)
        self.doubleSpinBox_9 = QtWidgets.QDoubleSpinBox(self.frame_5)
        self.doubleSpinBox_9.setMinimumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_9.setMaximumSize(QtCore.QSize(120, 30))
        self.doubleSpinBox_9.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_9.setMaximum(1.0)
        self.doubleSpinBox_9.setSingleStep(0.01)
        self.doubleSpinBox_9.setProperty("value", 0.05)
        self.doubleSpinBox_9.setObjectName("doubleSpinBox_9")
        self.gridLayout_4.addWidget(self.doubleSpinBox_9, 3, 1, 1, 1)
        self.pushButton_9 = QtWidgets.QPushButton(self.frame_5)
        self.pushButton_9.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_9.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_9.setText("")
        self.pushButton_9.setIcon(icon)
        self.pushButton_9.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_9.setObjectName("pushButton_9")
        self.gridLayout_4.addWidget(self.pushButton_9, 3, 2, 1, 1)
        self.gridLayout_3.addWidget(self.frame_5, 2, 0, 1, 1)
        self.gridLayout_3.setRowStretch(0, 1)
        self.gridLayout_3.setRowStretch(1, 1)
        self.gridLayout_3.setRowStretch(2, 10)
        self.gridLayout.addWidget(self.frame_2, 1, 1, 1, 1)
        self.frame = QtWidgets.QFrame(Form)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.frame)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.widget_2 = QtWidgets.QWidget(self.frame)
        self.widget_2.setObjectName("widget_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.widget_2)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label_9 = QtWidgets.QLabel(self.widget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(20)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_9.setFont(font)
        self.label_9.setStyleSheet("font: 20pt;")
        self.label_9.setObjectName("label_9")
        self.horizontalLayout.addWidget(self.label_9)
        self.gridLayout_5.addWidget(self.widget_2, 0, 0, 1, 1)
        self.frame_4 = QtWidgets.QFrame(self.frame)
        self.frame_4.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_4.setObjectName("frame_4")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.frame_4)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label_3 = QtWidgets.QLabel(self.frame_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setMinimumSize(QtCore.QSize(0, 0))
        self.label_3.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_3.setFont(font)
        self.label_3.setStyleSheet("font: 18pt;")
        self.label_3.setObjectName("label_3")
        self.gridLayout_2.addWidget(self.label_3, 0, 0, 1, 1)
        self.doubleSpinBox_5 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_5.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_5.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_5.setStyleSheet("font: 16pt;")
        self.doubleSpinBox_5.setDecimals(2)
        self.doubleSpinBox_5.setMaximum(9999999999999.0)
        self.doubleSpinBox_5.setSingleStep(0.1)
        self.doubleSpinBox_5.setProperty("value", 1.5)
        self.doubleSpinBox_5.setObjectName("doubleSpinBox_5")
        self.gridLayout_2.addWidget(self.doubleSpinBox_5, 0, 1, 1, 1)
        self.pushButton_10 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_10.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_10.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_10.setText("")
        self.pushButton_10.setIcon(icon)
        self.pushButton_10.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_10.setObjectName("pushButton_10")
        self.gridLayout_2.addWidget(self.pushButton_10, 0, 2, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.frame_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setMinimumSize(QtCore.QSize(0, 30))
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 30))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_2.setFont(font)
        self.label_2.setStyleSheet("font: 18pt;")
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 1, 0, 1, 1)
        self.doubleSpinBox_3 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_3.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_3.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_3.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_3.setDecimals(2)
        self.doubleSpinBox_3.setMaximum(99999.0)
        self.doubleSpinBox_3.setSingleStep(0.1)
        self.doubleSpinBox_3.setProperty("value", 0.7)
        self.doubleSpinBox_3.setObjectName("doubleSpinBox_3")
        self.gridLayout_2.addWidget(self.doubleSpinBox_3, 1, 1, 1, 1)
        self.pushButton_11 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_11.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_11.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_11.setText("")
        self.pushButton_11.setIcon(icon)
        self.pushButton_11.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_11.setObjectName("pushButton_11")
        self.gridLayout_2.addWidget(self.pushButton_11, 1, 2, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.frame_4)
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_4.setFont(font)
        self.label_4.setStyleSheet("font: 18pt;")
        self.label_4.setWordWrap(True)
        self.label_4.setObjectName("label_4")
        self.gridLayout_2.addWidget(self.label_4, 2, 0, 1, 1)
        self.doubleSpinBox_4 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_4.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_4.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_4.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_4.setDecimals(2)
        self.doubleSpinBox_4.setSingleStep(0.1)
        self.doubleSpinBox_4.setProperty("value", 0.3)
        self.doubleSpinBox_4.setObjectName("doubleSpinBox_4")
        self.gridLayout_2.addWidget(self.doubleSpinBox_4, 2, 1, 1, 1)
        self.pushButton_12 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_12.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_12.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_12.setText("")
        self.pushButton_12.setIcon(icon)
        self.pushButton_12.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_12.setObjectName("pushButton_12")
        self.gridLayout_2.addWidget(self.pushButton_12, 2, 2, 1, 1)
        self.label_5 = QtWidgets.QLabel(self.frame_4)
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_5.setFont(font)
        self.label_5.setStyleSheet("font: 18pt;")
        self.label_5.setWordWrap(True)
        self.label_5.setObjectName("label_5")
        self.gridLayout_2.addWidget(self.label_5, 3, 0, 1, 1)
        self.spinBox_2 = QtWidgets.QSpinBox(self.frame_4)
        self.spinBox_2.setMinimumSize(QtCore.QSize(100, 30))
        self.spinBox_2.setMaximumSize(QtCore.QSize(100, 30))
        self.spinBox_2.setStyleSheet("font: 18pt;")
        self.spinBox_2.setMaximum(99999)
        self.spinBox_2.setProperty("value", 1)
        self.spinBox_2.setObjectName("spinBox_2")
        self.gridLayout_2.addWidget(self.spinBox_2, 3, 1, 1, 1)
        self.pushButton_13 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_13.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_13.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_13.setText("")
        self.pushButton_13.setIcon(icon)
        self.pushButton_13.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_13.setObjectName("pushButton_13")
        self.gridLayout_2.addWidget(self.pushButton_13, 3, 2, 1, 1)
        self.label_16 = QtWidgets.QLabel(self.frame_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setMinimumSize(QtCore.QSize(0, 0))
        self.label_16.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_16.setFont(font)
        self.label_16.setStyleSheet("font: 18pt;")
        self.label_16.setWordWrap(True)
        self.label_16.setObjectName("label_16")
        self.gridLayout_2.addWidget(self.label_16, 4, 0, 1, 1)
        self.doubleSpinBox_6 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_6.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_6.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_6.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_6.setDecimals(2)
        self.doubleSpinBox_6.setSingleStep(0.1)
        self.doubleSpinBox_6.setProperty("value", 0.2)
        self.doubleSpinBox_6.setObjectName("doubleSpinBox_6")
        self.gridLayout_2.addWidget(self.doubleSpinBox_6, 4, 1, 1, 1)
        self.pushButton_14 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_14.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_14.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_14.setText("")
        self.pushButton_14.setIcon(icon)
        self.pushButton_14.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_14.setObjectName("pushButton_14")
        self.gridLayout_2.addWidget(self.pushButton_14, 4, 2, 1, 1)
        self.label_17 = QtWidgets.QLabel(self.frame_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_17.sizePolicy().hasHeightForWidth())
        self.label_17.setSizePolicy(sizePolicy)
        self.label_17.setMinimumSize(QtCore.QSize(0, 30))
        self.label_17.setMaximumSize(QtCore.QSize(16777215, 30))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_17.setFont(font)
        self.label_17.setStyleSheet("font: 18pt;")
        self.label_17.setWordWrap(False)
        self.label_17.setObjectName("label_17")
        self.gridLayout_2.addWidget(self.label_17, 5, 0, 1, 1)
        self.doubleSpinBox_10 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_10.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_10.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_10.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_10.setSuffix("")
        self.doubleSpinBox_10.setDecimals(2)
        self.doubleSpinBox_10.setMaximum(99.0)
        self.doubleSpinBox_10.setSingleStep(0.1)
        self.doubleSpinBox_10.setProperty("value", 0.8)
        self.doubleSpinBox_10.setObjectName("doubleSpinBox_10")
        self.gridLayout_2.addWidget(self.doubleSpinBox_10, 5, 1, 1, 1)
        self.pushButton_15 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_15.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_15.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_15.setText("")
        self.pushButton_15.setIcon(icon)
        self.pushButton_15.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_15.setObjectName("pushButton_15")
        self.gridLayout_2.addWidget(self.pushButton_15, 5, 2, 1, 1)
        self.label_18 = QtWidgets.QLabel(self.frame_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy)
        self.label_18.setMinimumSize(QtCore.QSize(0, 30))
        self.label_18.setMaximumSize(QtCore.QSize(16777215, 30))
        font = QtGui.QFont()
        font.setPointSize(18)
        font.setBold(False)
        font.setItalic(False)
        font.setWeight(50)
        self.label_18.setFont(font)
        self.label_18.setStyleSheet("font: 18pt;")
        self.label_18.setWordWrap(False)
        self.label_18.setObjectName("label_18")
        self.gridLayout_2.addWidget(self.label_18, 6, 0, 1, 1)
        self.doubleSpinBox_11 = QtWidgets.QDoubleSpinBox(self.frame_4)
        self.doubleSpinBox_11.setMinimumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_11.setMaximumSize(QtCore.QSize(100, 30))
        self.doubleSpinBox_11.setStyleSheet("font: 18pt;")
        self.doubleSpinBox_11.setDecimals(2)
        self.doubleSpinBox_11.setMaximum(1.0)
        self.doubleSpinBox_11.setSingleStep(0.1)
        self.doubleSpinBox_11.setProperty("value", 0.4)
        self.doubleSpinBox_11.setObjectName("doubleSpinBox_11")
        self.gridLayout_2.addWidget(self.doubleSpinBox_11, 6, 1, 1, 1)
        self.pushButton_16 = QtWidgets.QPushButton(self.frame_4)
        self.pushButton_16.setMinimumSize(QtCore.QSize(20, 20))
        self.pushButton_16.setMaximumSize(QtCore.QSize(20, 20))
        self.pushButton_16.setText("")
        self.pushButton_16.setIcon(icon)
        self.pushButton_16.setIconSize(QtCore.QSize(20, 20))
        self.pushButton_16.setObjectName("pushButton_16")
        self.gridLayout_2.addWidget(self.pushButton_16, 6, 2, 1, 1)
        self.gridLayout_5.addWidget(self.frame_4, 1, 0, 1, 1)
        self.gridLayout_5.setRowStretch(0, 1)
        self.gridLayout_5.setRowStretch(1, 10)
        self.gridLayout.addWidget(self.frame, 1, 0, 1, 1)

        self.retranslateUi(Form)
        self.pushButton_10.clicked.connect(Form.slot1) # type: ignore
        self.pushButton_11.clicked.connect(Form.slot2) # type: ignore
        self.pushButton_12.clicked.connect(Form.slot3) # type: ignore
        self.pushButton_13.clicked.connect(Form.slot4) # type: ignore
        self.pushButton_14.clicked.connect(Form.slot5) # type: ignore
        self.pushButton_15.clicked.connect(Form.slot6) # type: ignore
        self.pushButton_16.clicked.connect(Form.slot7) # type: ignore
        self.pushButton_6.clicked.connect(Form.slot8) # type: ignore
        self.pushButton_7.clicked.connect(Form.slot9) # type: ignore
        self.pushButton_8.clicked.connect(Form.slot10) # type: ignore
        self.pushButton_9.clicked.connect(Form.slot11) # type: ignore
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label.setText(_translate("Form", "Step3:选择RT"))
        self.label_6.setText(_translate("Form", "RT参与打分："))
        self.label_12.setText(_translate("Form", "保留信息参数"))
        self.label_7.setText(_translate("Form", "罚分窗口："))
        self.doubleSpinBox_12.setPrefix(_translate("Form", "±"))
        self.doubleSpinBox_12.setSuffix(_translate("Form", "min"))
        self.label_8.setText(_translate("Form", "罚分系数："))
        self.label_10.setText(_translate("Form", "最大罚分："))
        self.label_11.setText(_translate("Form", "保留信息缺失罚分："))
        self.label_9.setText(_translate("Form", "定性参数"))
        self.label_3.setText(_translate("Form", "搜库范围："))
        self.label_2.setText(_translate("Form", "峰组Match系数："))
        self.label_4.setText(_translate("Form", "峰组R Match系数："))
        self.label_5.setText(_translate("Form", "参与定性group最少离子数："))
        self.label_16.setText(_translate("Form", "峰组打分权重："))
        self.label_17.setText(_translate("Form", "直接搜索打分权重："))
        self.label_18.setText(_translate("Form", "相似性阈值："))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form_RT()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())
