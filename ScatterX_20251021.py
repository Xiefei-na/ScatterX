# -*- coding: utf-8 -*-
"""
Created on Sun Oct 19 23:48:40 2025

@author: admin
"""
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                            QPushButton, QListWidget, QTextEdit, QFileDialog, QGridLayout,
                            QLineEdit, QSplitter, QLabel, QFrame, QSizePolicy, QScrollArea, 
                            QComboBox, QMessageBox, QMenu, QAction, QSpinBox, QDoubleSpinBox,
                            QSplashScreen, QTabWidget, QGroupBox)   # ← 新增



import sys
import os
import re
import numpy as np
import csv
from PyQt5.QtCore import Qt, QPoint
from PyQt5.QtGui import (QPixmap, QImage, QPainter, QColor, QFont, QPen, QCursor, 
                        QGuiApplication, QIcon)
from PIL import Image
import matplotlib.cm as cm
from scipy.stats import pearsonr, describe
from scipy.signal import find_peaks
import math
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import datetime
import time

# ===================== 替换 1：__init__ =====================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        # 默认字体大小
        self.current_font_size = 8
        self.setWindowTitle("ScatterX - 散射数据分析系统")
        self._tab2_widgets_initialized = False
        self.curve_data = None
        self.apply_params_btn = None  # 显式置空，避免误调用
        self.resize(1200, 700)
        self.axes_tab2 = [None] * 6
        # 设置窗体图标
        self.set_window_icon()
        self.region_lines = [None] * 6   # 0-5 对应 4 张图 6 条竖线
        # -------------------------- 核心数据存储（原有+新增）--------------------------
        # 原有图像与分析数据
        self.current_image_array = None
        self.current_image_mode = None  # 'color' 或 'grayscale'
        self.image_matrix = None  # 原始图像数据矩阵（会被mask修改）
        self.original_image_matrix = None  # 原始图像备份（用于取消mask）
        self.horizontal_center = None  # 水平中心（列，x坐标，0-based）
        self.vertical_center = None    # 竖直中心（行，y坐标，0-based）
        self.polar_expand_matrix = None
        self.diff_matrix = None
        self.current_file_path = None  # 当前加载的文件路径
        self.background_file_path = None  # 本底图像文件路径

        # 本底相关数据
        self.background_image_array = None
        self.background_image_mode = None
        self.background_matrix = None
        self.background_coeff = 1.0

        # 能量与距离参数
        self.energy = None  # 能量，单位eV
        self.detector_distance = None  # 样品到探测器距离，单位mm

        # 曲线裁剪参数
        self.crop_front_nan = 0  # 前n个点设为NaN
        self.crop_back_nan = 0   # 后m个点设为NaN

        # q变换参数（支持小数能量值）
        self.q_energy = 1000.0  # 用于q变换的能量（默认1000.0 eV）
        self.q_distance = 200   # 用于q变换的探测器距离（默认200 mm）

        # Mask功能相关变量
        self.mask_mode = None  # None: 正常模式, 'circle': 圆形mask, 'square': 方形mask
        self.mask_start_point = None  # 拖拽起始点（图像坐标系）
        self.mask_current_point = None  # 拖拽当前点（图像坐标系）
        self.mask_active = False  # 是否正在进行mask绘制
        self.masked_regions = []  # 已创建的mask区域记录 [(type, params), ...]
        self.mask_threshold = 0  # 阈值Mask的阈值
        # 对比度手动调整标记
        self.contrast_manual_set = False

        # -------------------------- 新增：数据分析选项卡数据 --------------------------
        self.analyze_curve_list = []  # 已加载的分析曲线列表 [{"file": "", "q": [], "intensity": [], "label": ""}, ...]
        self.analyze_figure = Figure(figsize=(8, 5), dpi=100)  # 数据分析曲线画布
        self.analyze_canvas = FigureCanvas(self.analyze_figure)

        # 峰值检测参数
        self.peak_height = 0.1  # 峰值最小高度（相对最大值的比例）
        self.peak_distance = 5  # 峰值间最小水平距离（像素）

        # -------------------------- 初始化选项卡布局 --------------------------
        self._ensure_tab2_widgets()  # ✅ 先初始化控件

        self.init_tab_widget()

        # 应用样式
        self.apply_styles()

        

        

        self._tab2_widgets_initialized = False



            
            
    def _ensure_tab2_widgets(self):
        """确保 tab2 的所有控件已实例化"""
        if self._tab2_widgets_initialized:
            return
    
        # 初始化Guinier范围控件
        self.guinier_qmin = QDoubleSpinBox()
        self.guinier_qmin.setMinimum(0)
        self.guinier_qmin.setMaximum(1000)
        self.guinier_qmin.setDecimals(3)
        self.guinier_qmin.setSingleStep(0.01)
        self.guinier_qmin.setToolTip("Guinier分析的q²最小值（单位：nm⁻²）")
    
        self.guinier_qmax = QDoubleSpinBox()
        self.guinier_qmax.setMinimum(0)
        self.guinier_qmax.setMaximum(1000)
        self.guinier_qmax.setDecimals(3)
        self.guinier_qmax.setSingleStep(0.01)
        self.guinier_qmax.setToolTip("Guinier分析的q²最大值（单位：nm⁻²）")
    
        self.guinier_result = QLabel("斜率: --, 截距: --")
        self.guinier_result.setStyleSheet("color: #2c3e50; font-weight: bold;")
    
        # 初始化Debye范围控件
        self.debye_qmin = QDoubleSpinBox()
        self.debye_qmin.setMinimum(0)
        self.debye_qmin.setMaximum(1000)
        self.debye_qmin.setDecimals(3)
        self.debye_qmin.setSingleStep(0.01)
        self.debye_qmin.setToolTip("Debye分析的q²最小值（单位：nm⁻²）")
    
        self.debye_qmax = QDoubleSpinBox()
        self.debye_qmax.setMinimum(0)
        self.debye_qmax.setMaximum(1000)
        self.debye_qmax.setDecimals(3)
        self.debye_qmax.setSingleStep(0.01)
        self.debye_qmax.setToolTip("Debye分析的q²最大值（单位：nm⁻²）")
    
        self.debye_result = QLabel("斜率: --, 截距: --")
        self.debye_result.setStyleSheet("color: #2c3e50; font-weight: bold;")
    
        # 初始化Porod范围控件
        self.porod_qmin = QDoubleSpinBox()
        self.porod_qmin.setMinimum(0)
        self.porod_qmin.setMaximum(1000)
        self.porod_qmin.setDecimals(3)
        self.porod_qmin.setSingleStep(0.01)
        self.porod_qmin.setToolTip("Porod分析的q²最小值（单位：nm⁻²）")
    
        self.porod_qmax = QDoubleSpinBox()
        self.porod_qmax.setMinimum(0)
        self.porod_qmax.setMaximum(1000)
        self.porod_qmax.setDecimals(3)
        self.porod_qmax.setSingleStep(0.01)
        self.porod_qmax.setToolTip("Porod分析的q²最大值（单位：nm⁻²）")
    
        self.porod_result = QLabel("斜率: --, 截距: --")
        self.porod_result.setStyleSheet("color: #2c3e50; font-weight: bold;")
    
        # 初始化按钮控件
        self.apply_params_btn = QPushButton("应用所有范围并拟合")
        self.guinier_fit_btn = QPushButton("Guinier拟合")
        self.debye_fit_btn = QPushButton("Debye拟合")
        self.porod_fit_btn = QPushButton("Porod拟合")
    
        # 初始化曲线列表控件
        self.curve_list_tab2 = QListWidget()
        self.curve_list_tab2.setStyleSheet("border: 1px solid #d0d0d0; border-radius: 4px;")
        self.curve_list_tab2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    
        # 初始化加载曲线按钮
        self.load_curve_btn_tab2 = QPushButton("导入一维曲线文件夹")
        self.load_curve_btn_tab2.setStyleSheet("background-color: #27ae60;")
    
        # 新增：物理量计算结果标签（含误差标签）
        self.guinier_rg = QLabel("Rg: --")
        self.guinier_rg.setStyleSheet("color: #e74c3c; font-weight: bold;")
        self.guinier_rg_error = QLabel("误差: --")
        self.guinier_rg_error.setStyleSheet("color: #3498db; font-weight: bold;")
        
        self.debye_rc = QLabel("Rc: --")
        self.debye_rc.setStyleSheet("color: #e74c3c; font-weight: bold;")
        self.debye_rc_error = QLabel("误差: --")
        self.debye_rc_error.setStyleSheet("color: #3498db; font-weight: bold;")
        
        self.porod_deviation = QLabel("Porod偏离: --")
        self.porod_deviation.setStyleSheet("color: #e74c3c; font-weight: bold;")
        self.porod_deviation_error = QLabel("误差: --")
        self.porod_deviation_error.setStyleSheet("color: #3498db; font-weight: bold;")
    
        self._tab2_widgets_initialized = True
            
        
    def set_window_icon(self):
        """设置软件窗体图标（使用ScatterX_logo.png）"""
        icon_path = "ScatterX_logo.png"  # 图标文件路径（与脚本同目录）
        try:
            icon_pixmap = QPixmap(icon_path)
            if not icon_pixmap.isNull():
                self.setWindowIcon(QIcon(icon_pixmap))
        except Exception as e:
            pass  # 图标加载失败不影响主程序运行

    def init_tab_widget(self):
        """初始化选项卡控件：选项卡1（原界面）、选项卡2（数据分析）"""
        # 1. 创建主选项卡容器
        self.tab_widget = QTabWidget()
        self.setCentralWidget(self.tab_widget)
        
        # 2. 创建选项卡1：图像采集与极坐标分析（原界面）
        self.tab1 = QWidget()
        self.create_tab1_layout()  # 加载原界面布局
        self.tab_widget.addTab(self.tab1, "图像处理")
        
        # 3. 创建选项卡2：数据分析（新增）
        self.tab2 = QWidget()
        self.create_tab2_layout()  # 加载数据分析布局
        self.tab_widget.addTab(self.tab2, "曲线分析")
        
        # 4. 选项卡切换事件（可选：用于刷新界面）
        self.tab_widget.currentChanged.connect(self.on_tab_switch)

    
 

    def create_tab1_layout(self):
        """选项卡1布局：原有的三栏式界面（文件选择+图像显示+参数设置）"""
        main_layout = QHBoxLayout(self.tab1)
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(20, 20, 20, 20)
        
        # 初始化三栏布局（复用原有逻辑）
        self.create_first_column()
        self.create_second_column()
        self.create_third_column()
        
        main_layout.addWidget(self.first_column, 1)    # 第一栏权重1
        main_layout.addWidget(self.second_column, 2)   # 第二栏权重2
        main_layout.addWidget(self.third_column, 1)    # 第三栏权重2
    # ----------------------------  新增：自动 q 范围  ----------------------------

    


    # 1) 自动范围函数：把原来的 q-range 换成 q²-range
    def auto_set_q_range_controls(self, q):
        """根据 q 实际范围自动设置 q²-range"""
        q_min, q_max = float(np.nanmin(q)), float(np.nanmax(q))
        if not (np.isfinite(q_min) and np.isfinite(q_max)) or q_max <= q_min:
            self.log_print("警告：无法获取有效 q 范围")
            return
        q2_min, q2_max = q_min**2, q_max**2
    
        # Guinier：低 q² 区（前 1/3）
        guinier_q2max = q2_min + (q2_max - q2_min)/3.0
        self.guinier_qmin.setValue(q2_min)
        self.guinier_qmax.setValue(guinier_q2max)
    
        # Porod：高 q² 区（后 1/2）
        porod_q2min = q2_min + (q2_max - q2_min)/2.0
        self.porod_qmin.setValue(porod_q2min)
        self.porod_qmax.setValue(q2_max)
    
        # Debye：中间 q² 区（1/3~2/3）
        debye_q2min = q2_min + (q2_max - q2_min)/3.0
        debye_q2max = q2_min + 2*(q2_max - q2_min)/3.0
        self.debye_qmin.setValue(debye_q2min)
        self.debye_qmax.setValue(debye_q2max)
    
        self.log_print(f"已自动设置 q² 范围：{q2_min:.4f} – {q2_max:.4f} nm⁻²")
        
        
        
    # ===================== 替换 2：create_tab2_layout =====================
    def create_tab2_layout(self):
        """选项卡2布局：曲线分析功能区（分三栏）"""
        from PyQt5.QtWidgets import (QVBoxLayout, QHBoxLayout, QLabel,
                                     QPushButton, QFrame, QSizePolicy,
                                     QScrollArea, QWidget, QGridLayout,
                                     QGroupBox, QSpinBox, QDoubleSpinBox)
    
        main_layout = QHBoxLayout(self.tab2)
        main_layout.setSpacing(15)
        main_layout.setContentsMargins(20, 20, 20, 20)
    
        # -------------------------- 第一栏：数据加载区 --------------------------
        self.first_column_tab2 = QFrame()
        self.first_column_tab2.setFrameStyle(QFrame.StyledPanel)
        first_layout = QVBoxLayout(self.first_column_tab2)
        first_layout.setContentsMargins(10, 10, 10, 10)
        first_layout.setSpacing(10)
    
        # 绑定已初始化的加载按钮
        first_layout.addWidget(self.load_curve_btn_tab2)
        self.load_curve_btn_tab2.clicked.connect(self.load_analyze_curve_folder)
    
        # 绑定已初始化的曲线列表
        self.curve_list_tab2.itemDoubleClicked.connect(self.on_curve_double_clicked)
        first_layout.addWidget(self.curve_list_tab2)
        main_layout.addWidget(self.first_column_tab2, 1)
    
        # -------------------------- 第二栏：图表显示区（上下结构3:1） --------------------------
        self.second_column_tab2 = QFrame()
        self.second_column_tab2.setFrameStyle(QFrame.StyledPanel)
        second_layout = QVBoxLayout(self.second_column_tab2)
        second_layout.setContentsMargins(10, 10, 10, 10)
        second_layout.setSpacing(10)
    
        # 上半部分：6张图表（占3/4）
        upper_widget = QWidget()
        upper_layout = QGridLayout(upper_widget)
    
        plot_titles = ["原始 q~I", "LogLog 曲线", "Lorentz 曲线 (q~q²·I)",
                      "Guinier 曲线 (q²~ln I)", "Debye 曲线 (q²~I^(-½))", "Porod 曲线 (q²~ln[q⁴·I])"]
        self.figures_tab2 = []
        self.canvases_tab2 = []
        self.axes_tab2 = [None] * 6
        # 存储各区域的选中范围和拟合线
        self.region_masks = {
            'guinier': None,
            'debye': None,
            'porod': None
        }
        self.region_lines = {
            'guinier': [None, None],  # 左右边界线
            'debye': [None, None],
            'porod': [None, None]
        }
        self.fit_lines = {
            'guinier': None,
            'debye': None,
            'porod': None
        }
        self.fit_results = {
            'guinier': {'slope': 0, 'intercept': 0, 'r_squared': 0, 'error': 0},
            'debye': {'slope': 0, 'intercept': 0, 'r_squared': 0, 'error': 0},
            'porod': {'slope': 0, 'intercept': 0, 'r_squared': 0, 'error': 0}
        }
    
        for i in range(6):
            fig = Figure(figsize=(4, 3), dpi=80)
            canvas = FigureCanvas(fig)
            ax = fig.add_subplot(111)
            ax.set_title(plot_titles[i], fontsize=self.current_font_size)
            ax.tick_params(labelsize=self.current_font_size - 1)
            self.figures_tab2.append(fig)
            self.canvases_tab2.append(canvas)
            self.axes_tab2[i] = ax
            row, col = divmod(i, 3)
            upper_layout.addWidget(canvas, row, col)
    
        # 下半部分：范围选择和拟合结果（占1/4）
        lower_widget = QWidget()
        lower_layout = QVBoxLayout(lower_widget)
    
        # 创建分组框
        group_box = QGroupBox("分析范围设置 (q² 单位)")
        group_layout = QGridLayout(group_box)
        group_box.setFont(QFont("Arial", self.current_font_size))
    
        # -------------------------- Guinier范围与结果布局（含误差） --------------------------
        group_layout.addWidget(QLabel("Guinier 范围:"), 0, 0)
        group_layout.addWidget(self.guinier_qmin, 0, 1)
        self.guinier_qmin.valueChanged.connect(self.update_guinier_range)
        group_layout.addWidget(QLabel("至"), 0, 2)
        group_layout.addWidget(self.guinier_qmax, 0, 3)
        self.guinier_qmax.valueChanged.connect(self.update_guinier_range)
        group_layout.addWidget(self.guinier_result, 0, 4)
        self.guinier_fit_btn.clicked.connect(lambda: self._fit_range('guinier', 3))
        group_layout.addWidget(self.guinier_fit_btn, 0, 5)
        group_layout.addWidget(self.guinier_rg, 0, 6)
        group_layout.addWidget(self.guinier_rg_error, 0, 7)
    
        # -------------------------- Debye范围与结果布局（含误差） --------------------------
        group_layout.addWidget(QLabel("Debye 范围:"), 1, 0)
        group_layout.addWidget(self.debye_qmin, 1, 1)
        self.debye_qmin.valueChanged.connect(self.update_debye_range)
        group_layout.addWidget(QLabel("至"), 1, 2)
        group_layout.addWidget(self.debye_qmax, 1, 3)
        self.debye_qmax.valueChanged.connect(self.update_debye_range)
        group_layout.addWidget(self.debye_result, 1, 4)
        self.debye_fit_btn.clicked.connect(lambda: self._fit_range('debye', 4))
        group_layout.addWidget(self.debye_fit_btn, 1, 5)
        group_layout.addWidget(self.debye_rc, 1, 6)
        group_layout.addWidget(self.debye_rc_error, 1, 7)
    
        # -------------------------- Porod范围与结果布局（含误差） --------------------------
        group_layout.addWidget(QLabel("Porod 范围:"), 2, 0)
        group_layout.addWidget(self.porod_qmin, 2, 1)
        self.porod_qmin.valueChanged.connect(self.update_porod_range)
        group_layout.addWidget(QLabel("至"), 2, 2)
        group_layout.addWidget(self.porod_qmax, 2, 3)
        self.porod_qmax.valueChanged.connect(self.update_porod_range)
        group_layout.addWidget(self.porod_result, 2, 4)
        self.porod_fit_btn.clicked.connect(lambda: self._fit_range('porod', 5))
        group_layout.addWidget(self.porod_fit_btn, 2, 5)
        group_layout.addWidget(self.porod_deviation, 2, 6)
        group_layout.addWidget(self.porod_deviation_error, 2, 7)
    
        # 应用按钮与保存按钮布局（适配误差列）
        self.apply_params_btn.clicked.connect(self.apply_all_ranges)
        group_layout.addWidget(self.apply_params_btn, 3, 0, 1, 6)
        self.save_data_btn = QPushButton("保存分析结果")
        self.save_data_btn.clicked.connect(self.save_analysis_results)
        group_layout.addWidget(self.save_data_btn, 3, 6, 1, 2)
    
        lower_layout.addWidget(group_box)
        lower_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
    
        # 按3:1比例添加上下部分
        second_layout.addWidget(upper_widget, 3)
        second_layout.addWidget(lower_widget, 1)
    
        main_layout.addWidget(self.second_column_tab2, 3)



    def save_analysis_results(self):
        """保存当前分析结果到.saxs文件（含误差）"""
        # 获取当前选中的曲线
        idx = self.curve_list_tab2.currentRow()
        if idx < 0 or idx >= len(self.analyze_curve_list):
            QMessageBox.warning(self, "保存失败", "请先选择要保存的曲线")
            return
        
        curve = self.analyze_curve_list[idx]
        original_file = curve["file"]
        base_name = os.path.splitext(original_file)[0]
        save_path = f"{base_name}.saxs"
        
        try:
            with open(save_path, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                
                # 写入基本信息
                writer.writerow(["文件名:", os.path.basename(original_file)])
                writer.writerow(["数据点数量:", len(curve["q"])])
                writer.writerow(["保存时间:", datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")])
                writer.writerow([])
                
                # 写入拟合结果（新增误差列）
                writer.writerow(["分析类型", "q²最小值", "q²最大值", "斜率", "截距", "物理量结果", "误差结果"])
                
                # Guinier结果（含误差）
                writer.writerow([
                    "Guinier",
                    f"{self.guinier_qmin.value():.6f}",
                    f"{self.guinier_qmax.value():.6f}",
                    f"{self.fit_results['guinier']['slope']:.6f}",
                    f"{self.fit_results['guinier']['intercept']:.6f}",
                    self.guinier_rg.text(),
                    self.guinier_rg_error.text()
                ])
                
                # Debye结果（含误差）
                writer.writerow([
                    "Debye",
                    f"{self.debye_qmin.value():.6f}",
                    f"{self.debye_qmax.value():.6f}",
                    f"{self.fit_results['debye']['slope']:.6f}",
                    f"{self.fit_results['debye']['intercept']:.6f}",
                    self.debye_rc.text(),
                    self.debye_rc_error.text()
                ])
                
                # Porod结果（含误差）
                writer.writerow([
                    "Porod",
                    f"{self.porod_qmin.value():.6f}",
                    f"{self.porod_qmax.value():.6f}",
                    f"{self.fit_results['porod']['slope']:.6f}",
                    f"{self.fit_results['porod']['intercept']:.6f}",
                    self.porod_deviation.text(),
                    self.porod_deviation_error.text()
                ])
                
                writer.writerow([])
                writer.writerow(["q值", "强度值"])
                
                # 写入原始数据
                for q_val, i_val in zip(curve["q"], curve["intensity"]):
                    writer.writerow([q_val, i_val])
            
            QMessageBox.information(self, "保存成功", f"分析结果已保存至:\n{save_path}")
            self.log_print(f"分析结果保存成功: {save_path}")
            
        except Exception as e:
            QMessageBox.warning(self, "保存失败", f"保存文件时出错:\n{str(e)}")
            self.log_print(f"分析结果保存失败: {str(e)}")     
            
            
            
    def update_guinier_range(self):
        """更新Guinier范围显示"""
        self._update_range('guinier', 3)  # 3是Guinier图的索引
    
    def update_debye_range(self):
        """更新Debye范围显示"""
        self._update_range('debye', 4)  # 4是Debye图的索引
    
    def update_porod_range(self):
        """更新Porod范围显示"""
        self._update_range('porod', 5)  # 5是Porod图的索引
    
    def _update_range(self, region, ax_index):
        """通用的范围更新函数：清除旧标记并绘制新的选中区域标记"""
        if ax_index >= len(self.axes_tab2) or self.axes_tab2[ax_index] is None:
            return
    
        # 获取当前选中的曲线
        idx = self.curve_list_tab2.currentRow()
        if idx < 0 or idx >= len(self.analyze_curve_list):
            return
    
        curve = self.analyze_curve_list[idx]
        q = curve["q"].astype(float)
        I = curve["intensity"].astype(float)
    
        # 获取当前范围
        qmin_widget = getattr(self, f"{region}_qmin")
        qmax_widget = getattr(self, f"{region}_qmax")
        q2_min = qmin_widget.value()
        q2_max = qmax_widget.value()
    
        # 清除之前的标记（核心修复部分）
        ax = self.axes_tab2[ax_index]
    
        # 清除区域掩码
        if self.region_masks[region] is not None:
            try:
                self.region_masks[region].remove()
            except Exception as e:
                self.log_print(f"清除{region}区域掩码失败: {str(e)}")
            finally:
                self.region_masks[region] = None  # 重置引用
    
        # 清除边界线
        for i in range(len(self.region_lines[region])):
            line = self.region_lines[region][i]
            if line is not None:
                try:
                    if hasattr(line, 'remove'):
                        line.remove()
                except Exception as e:
                    self.log_print(f"清除{region}边界线失败: {str(e)}")
                finally:
                    self.region_lines[region][i] = None  # 重置引用
    
        # 清除拟合线（处理列表类型）
        if self.fit_lines[region] is not None:
            try:
                if isinstance(self.fit_lines[region], list):
                    for fit_line in self.fit_lines[region]:
                        fit_line.remove()
                else:
                    self.fit_lines[region].remove()
            except Exception as e:
                self.log_print(f"清除{region}拟合线失败: {str(e)}")
            finally:
                self.fit_lines[region] = None  # 重置引用
    
        # 计算有效数据范围
        q2 = q**2
        valid_idx = (q2 >= q2_min) & (q2 <= q2_max)
    
        # 在图上标记选中区域
        if region == 'guinier':
            y_data = np.log(I)
            valid_data_idx = valid_idx & (I > 0)  # 确保log有效
            x_display = q2[valid_data_idx]
            y_display = y_data[valid_data_idx]
        elif region == 'debye':
            y_data = 1/np.sqrt(I)
            valid_data_idx = valid_idx & (I > 0)  # 确保sqrt有效
            x_display = q2[valid_data_idx]
            y_display = y_data[valid_data_idx]
        elif region == 'porod':
            y_data = np.log(q**4 * I)
            valid_data_idx = valid_idx & (q != 0) & (I > 0)  # 确保log有效
            x_display = q2[valid_data_idx]
            y_display = y_data[valid_data_idx]
    
        # 绘制选中区域
        if len(x_display) > 0:
            self.region_masks[region] = ax.scatter(x_display, y_display, color='red', s=10, alpha=0.8)
            # 绘制边界线并保存引用
            self.region_lines[region][0] = ax.axvline(x=q2_min, color='red', linestyle='--', alpha=0.5)
            self.region_lines[region][1] = ax.axvline(x=q2_max, color='red', linestyle='--', alpha=0.5)
    
        self.canvases_tab2[ax_index].draw()
    
    def apply_all_ranges(self):
        """应用所有范围并执行线性拟合"""
        self._fit_range('guinier', 3)
        self._fit_range('debye', 4)
        self._fit_range('porod', 5)
    
    def _fit_range(self, region, ax_index):
        """对指定区域执行线性拟合并计算物理量与误差"""
        if ax_index >= len(self.axes_tab2) or self.axes_tab2[ax_index] is None:
            return
    
        # 获取当前选中的曲线
        idx = self.curve_list_tab2.currentRow()
        if idx < 0 or idx >= len(self.analyze_curve_list):
            return
    
        curve = self.analyze_curve_list[idx]
        q = curve["q"].astype(float)
        I = curve["intensity"].astype(float)
    
        # 获取当前范围
        qmin_widget = getattr(self, f"{region}_qmin")
        qmax_widget = getattr(self, f"{region}_qmax")
        q2_min = qmin_widget.value()
        q2_max = qmax_widget.value()
    
        if q2_min >= q2_max:
            getattr(self, f"{region}_result").setText("错误: 起始值需小于结束值")
            return
    
        # 计算有效数据范围
        q2 = q**2
        valid_idx = (q2 >= q2_min) & (q2 <= q2_max)
    
        # 准备拟合数据
        ax = self.axes_tab2[ax_index]
        try:
            if region == 'guinier':
                y_data = np.log(I)
                valid_data_idx = valid_idx & (I > 0)  # 确保log有效
                x_fit = q2[valid_data_idx]
                y_fit = y_data[valid_data_idx]
            elif region == 'debye':
                y_data = 1/np.sqrt(I)
                valid_data_idx = valid_idx & (I > 0)  # 确保sqrt有效
                x_fit = q2[valid_data_idx]
                y_fit = y_data[valid_data_idx]
            elif region == 'porod':
                y_data = np.log(q**4 * I)
                valid_data_idx = valid_idx & (q != 0) & (I > 0)  # 确保log有效
                x_fit = q2[valid_data_idx]
                y_fit = y_data[valid_data_idx]
    
            if len(x_fit) < 2:
                getattr(self, f"{region}_result").setText("数据点不足，无法拟合")
                return
    
            # 清除旧拟合线
            if self.fit_lines[region] is not None:
                try:
                    if isinstance(self.fit_lines[region], list):
                        for line in self.fit_lines[region]:
                            line.remove()
                    else:
                        self.fit_lines[region].remove()
                except Exception as e:
                    self.log_print(f"清除{region}旧拟合线失败: {str(e)}")
                finally:
                    self.fit_lines[region] = None
    
            # 执行线性拟合 y = kx + b
            coefficients = np.polyfit(x_fit, y_fit, 1)
            slope = coefficients[0]
            intercept = coefficients[1]
    
            # 新增：计算拟合误差（R²决定系数）
            y_fit_pred = slope * x_fit + intercept
            tss = np.sum((y_fit - np.mean(y_fit)) ** 2)
            rss = np.sum((y_fit - y_fit_pred) ** 2)
            r_squared = 1 - (rss / tss) if tss != 0 else 0.0
            fit_error = 1 - r_squared  # 拟合误差（0~1范围）
    
            # 保存拟合结果（含R²和误差）
            self.fit_results[region]['slope'] = slope
            self.fit_results[region]['intercept'] = intercept
            self.fit_results[region]['r_squared'] = r_squared
            self.fit_results[region]['error'] = fit_error
    
            # 绘制新拟合线
            x_line = np.linspace(np.min(x_fit), np.max(x_fit), 100)
            y_line = slope * x_line + intercept
            self.fit_lines[region] = ax.plot(x_line, y_line, 'r-', linewidth=2)
    
            # 更新结果显示
            result_text = f"斜率: {slope:.6f}, 截距: {intercept:.6f}"
            getattr(self, f"{region}_result").setText(result_text)
    
            # ---------------- 计算物理量与误差并更新标签 ----------------
            if region == 'guinier':
                # Guinier公式: Rg = sqrt(-3k)，误差 = Rg × fit_error
                if slope < 0:
                    rg = np.sqrt(-3 * slope)
                    rg_error = rg * fit_error
                    self.guinier_rg.setText(f"Rg: {rg:.4f} nm")
                    self.guinier_rg_error.setText(f"误差: {rg_error:.6f} nm")
                else:
                    self.guinier_rg.setText("Rg: 无效(斜率非负)")
                    self.guinier_rg_error.setText("误差: 无效")
    
            elif region == 'debye':
                # Debye公式: Rc = 3k/b，误差 = Rc × fit_error
                if intercept != 0:
                    rc = 3 * slope / intercept
                    rc_error = rc * fit_error
                    self.debye_rc.setText(f"Rc: {rc:.4f} nm")
                    self.debye_rc_error.setText(f"误差: {rc_error:.6f} nm")
                else:
                    self.debye_rc.setText("Rc: 无效(截距为0)")
                    self.debye_rc_error.setText("误差: 无效")
    
            elif region == 'porod':
                # Porod公式: 偏离值 = |k|，误差 = 偏离值 × fit_error
                porod_deviation = abs(slope)
                porod_error = porod_deviation * fit_error
                self.porod_deviation.setText(f"偏离: {porod_deviation:.6f}")
                self.porod_deviation_error.setText(f"误差: {porod_error:.8f}")
    
            self.canvases_tab2[ax_index].draw()
    
        except Exception as e:
            error_msg = f"拟合失败: {str(e)[:50]}..."
            getattr(self, f"{region}_result").setText(error_msg)
            # 异常时清空物理量与误差标签
            if region == 'guinier':
                self.guinier_rg.setText("Rg: 无效")
                self.guinier_rg_error.setText("误差: 无效")
            elif region == 'debye':
                self.debye_rc.setText("Rc: 无效")
                self.debye_rc_error.setText("误差: 无效")
            elif region == 'porod':
                self.porod_deviation.setText("Porod偏离: 无效")
                self.porod_deviation_error.setText("误差: 无效")
            self.log_print(f"{region}拟合异常: {str(e)}")
       
    def load_analyze_curve(self):
        """加载一维曲线数据（txt/dat格式）"""
        file_paths, _ = QFileDialog.getOpenFileNames(
            self, "选择一维曲线数据", "", "文本文件 (*.txt *.dat);;所有文件 (*)"
        )
        if not file_paths:
            return
        
        for file_path in file_paths:
            try:
                # 读取曲线数据（q值+强度值）
                data = np.loadtxt(file_path, skiprows=0, delimiter=None)
                if data.shape[0] < 2:
                    raise ValueError("数据行数不足（至少2行有效数据）")
                
                # 自动检测列格式（支持q-I或仅I）
                if data.shape[1] >= 2:
                    q = data[:, 0]
                    intensity = data[:, 1]
                else:
                    intensity = data[:, 0]
                    q = np.arange(len(intensity))
                
                file_name = os.path.basename(file_path)
                label = f"曲线_{len(self.analyze_curve_list)+1}_{os.path.splitext(file_name)[0]}"
                
                # 添加到曲线列表
                self.analyze_curve_list.append({
                    "file": file_path,
                    "q": q,
                    "intensity": intensity,
                    "label": label
                })
                
                # 更新列表显示
                self.curve_list_tab2.addItem(f"{label} - {file_name} ({len(q)}点)")
                self.log_print(f"加载曲线数据成功：{file_path}（{len(q)}个数据点）")
            
            except Exception as e:
                QMessageBox.warning(self, "加载失败", f"曲线{file_path}加载错误：{str(e)}")
                self.log_print(f"曲线加载失败：{file_path} - {str(e)}")

    def on_curve_double_clicked(self, item):
        """双击曲线文件 -> 6 张坐标变换图同步刷新，并在第一幅图标题显示文件名"""
        if not self.analyze_curve_list:
            return
        idx = self.curve_list_tab2.currentRow()
        if idx < 0 or idx >= len(self.analyze_curve_list):
            return
    
        curve = self.analyze_curve_list[idx]
        q = curve["q"].astype(float)
        I = curve["intensity"].astype(float)
        
        # 获取文件名（不含路径和扩展名）
        file_name = os.path.splitext(os.path.basename(curve["file"]))[0]
        
        # 刷新第一幅图（原始 q~I）并更新标题
        self.axes_tab2[0].clear()
        self.axes_tab2[0].plot(q, I)
        self.axes_tab2[0].set_title(f"Original q~I - {file_name}")  # 英文标题
        self.axes_tab2[0].set_xlabel("q")  # 英文标签
        self.axes_tab2[0].set_ylabel("Intensity")  # 英文标签
        self.axes_tab2[0].grid(True)
        
        # 第二幅图：LogLog 曲线
        self.axes_tab2[1].clear()
        self.axes_tab2[1].loglog(q, I)
        self.axes_tab2[1].set_title("LogLog Curve")  # 英文标题
        self.axes_tab2[1].set_xlabel("q (log)")  # 英文标签
        self.axes_tab2[1].set_ylabel("Intensity (log)")  # 英文标签
        self.axes_tab2[1].grid(True)
        
        # 第三幅图：Lorentz 曲线 (q~q²·I)
        self.axes_tab2[2].clear()
        self.axes_tab2[2].plot(q, q**2 * I)
        self.axes_tab2[2].set_title("Lorentz Curve (q vs q^2*I)")  # 英文标题，用^表示平方
        self.axes_tab2[2].set_xlabel("q")  # 英文标签
        self.axes_tab2[2].set_ylabel("q^2*I")  # 英文标签，用^表示平方
        self.axes_tab2[2].grid(True)
        
        # 第四幅图：Guinier 曲线 (q²~ln I)
        self.axes_tab2[3].clear()
        valid_idx = I > 0  # 避免log(0)或负数
        self.axes_tab2[3].plot(q[valid_idx]**2, np.log(I[valid_idx]))
        self.axes_tab2[3].set_title("Guinier Curve (q^2 vs ln I)")  # 英文标题，用^表示平方
        self.axes_tab2[3].set_xlabel("q^2")  # 英文标签，用^表示平方
        self.axes_tab2[3].set_ylabel("ln(I)")  # 英文标签
        self.axes_tab2[3].grid(True)
        
        # 第五幅图：Debye 曲线 (q²~I^(-½))
        self.axes_tab2[4].clear()
        valid_idx = I > 0  # 避免开方负数
        self.axes_tab2[4].plot(q[valid_idx]**2, 1/np.sqrt(I[valid_idx]))
        self.axes_tab2[4].set_title("Debye Curve (q^2 vs I^(-1/2))")  # 英文标题，用^表示指数
        self.axes_tab2[4].set_xlabel("q^2")  # 英文标签，用^表示平方
        self.axes_tab2[4].set_ylabel("I^(-1/2)")  # 英文标签，用^表示指数
        self.axes_tab2[4].grid(True)
        
        # 第六幅图：Porod 曲线 (q²~ln[q⁴·I])
        self.axes_tab2[5].clear()
        valid_idx = (q != 0) & (I > 0)  # 避免log(0)或负数
        self.axes_tab2[5].plot(q[valid_idx]**2, np.log(q[valid_idx]**4 * I[valid_idx]))
        self.axes_tab2[5].set_title("Porod Curve (q^2 vs ln(q^4*I))")  # 英文标题，用^表示指数
        self.axes_tab2[5].set_xlabel("q^2")  # 英文标签，用^表示平方
        self.axes_tab2[5].set_ylabel("ln(q^4*I)")  # 英文标签，用^表示指数
        self.axes_tab2[5].grid(True)
        
        # 更新所有画布
        for canvas in self.canvases_tab2:
            canvas.draw()
    
        



 

    def apply_tab2_font(self, size):
        """一次性把选项卡2所有文字改成同一字号"""
        font = QFont("Arial", size)
        # 列表
        self.curve_list_tab2.setFont(font)
        # 所有 QLabel / QSpinBox / QDoubleSpinBox
        for w in self.third_column_tab2.findChildren((QLabel, QSpinBox, QDoubleSpinBox)):
            w.setFont(font)
        # 按钮
        self.apply_params_btn.setFont(font)
        # 4 张图默认标题字体
        for ax in self.axes_tab2:
            ax.title.set_fontsize(size)
            ax.xaxis.label.set_fontsize(size)
            ax.yaxis.label.set_fontsize(size)
        for canv in self.canvases_tab2:
            canv.draw()
        
      


    def on_tab_switch(self, index):
        """选项卡切换事件"""
        pass  # 可以在这里添加切换时的刷新逻辑

    # -------------------------- 原有：选项卡1（原界面）其他方法（复用，未修改） --------------------------
    def create_first_column(self):
        """第一栏：文件选择区域（复用原有逻辑）"""
        self.first_column = QFrame()
        self.first_column.setFrameStyle(QFrame.StyledPanel)
        layout = QVBoxLayout(self.first_column)
        
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(10)
        
        self.select_folder_btn = QPushButton("选择二维图像文件夹")
        self.select_folder_btn.clicked.connect(self.select_folder)
        layout.addWidget(self.select_folder_btn)
        
        self.file_list = QListWidget()
        self.file_list.itemDoubleClicked.connect(self.on_file_double_clicked)
        self.file_list.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.file_list.setFont(QFont("Arial", self.current_font_size))
        layout.addWidget(self.file_list)
        
        # 取消所有Mask按钮
        self.clear_all_mask_btn = QPushButton("取消所有Mask")
        self.clear_all_mask_btn.clicked.connect(self.clear_all_masks)
        self.clear_all_mask_btn.setEnabled(False)
        layout.addWidget(self.clear_all_mask_btn)
        
        self.first_column.setMinimumWidth(240)
        self.first_column.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)

    def create_second_column(self):
        """第二栏：图像显示与曲线分析（复用原有逻辑）"""
        self.second_column = QFrame()
        self.second_column.setFrameStyle(QFrame.StyledPanel)
        layout = QVBoxLayout(self.second_column)
        
        # 对比度调节 + 颜色图选择区域
        contrast_widget = QWidget()
        contrast_layout = QHBoxLayout(contrast_widget)
        contrast_layout.addWidget(QLabel("对比度下限:"))
        self.lower_combo = QComboBox()
        self.lower_combo.addItems(["100", "500", "1000", "3000", "5000", "10000"])
        self.lower_combo.currentTextChanged.connect(self.on_contrast_manual_change)
        self.lower_combo.currentTextChanged.connect(self.apply_contrast)
        contrast_layout.addWidget(self.lower_combo)
        
        contrast_layout.addWidget(QLabel("对比度上限:"))
        self.upper_combo = QComboBox()
        self.upper_combo.addItems(["400", "1000", "3000", "5000", "10000", 
                                 "20000", "50000", "100000", "200000", "300000"])
        self.upper_combo.currentTextChanged.connect(self.on_contrast_manual_change)
        self.upper_combo.currentTextChanged.connect(self.apply_contrast)
        contrast_layout.addWidget(self.upper_combo)
        
        # 颜色图选择（支持viridis/plasma/jet）
        contrast_layout.addSpacing(15)
        contrast_layout.addWidget(QLabel("颜色图:"))
        self.cmap_combo = QComboBox()
        self.cmap_combo.addItems(["viridis", "plasma", "jet"])
        self.cmap_combo.currentTextChanged.connect(self.apply_contrast)
        self.cmap_combo.currentTextChanged.connect(self.update_background_cmap)
        self.cmap_combo.currentTextChanged.connect(self.update_diff_cmap)
        contrast_layout.addWidget(self.cmap_combo)
        
        layout.addWidget(contrast_widget)
        
        # Mask控制按钮区域
        mask_control_widget = QWidget()
        mask_control_layout = QHBoxLayout(mask_control_widget)
        mask_control_layout.setContentsMargins(0, 5, 0, 5)
        
        # 圆形Mask按钮
        self.circle_mask_btn = QPushButton("圆形Mask模式")
        self.circle_mask_btn.clicked.connect(lambda: self.set_mask_mode('circle'))
        self.circle_mask_btn.setStyleSheet("background-color: #6c757d;")
        mask_control_layout.addWidget(self.circle_mask_btn)
        
        # 方形Mask按钮
        self.square_mask_btn = QPushButton("方形Mask模式")
        self.square_mask_btn.clicked.connect(lambda: self.set_mask_mode('square'))
        self.square_mask_btn.setStyleSheet("background-color: #6c757d;")
        mask_control_layout.addWidget(self.square_mask_btn)
        
        # 阈值Mask相关控件 - 新增
        mask_control_layout.addSpacing(10)
        self.threshold_mask_label = QLabel("阈值:")
        mask_control_layout.addWidget(self.threshold_mask_label)
        self.threshold_mask_input = QSpinBox()
        self.threshold_mask_input.setMinimum(0)
        self.threshold_mask_input.setMaximum(100000)
        self.threshold_mask_input.setValue(0)
        self.threshold_mask_input.valueChanged.connect(self.set_mask_threshold)
        mask_control_layout.addWidget(self.threshold_mask_input)
        self.threshold_mask_btn = QPushButton("阈值Mask")
        self.threshold_mask_btn.clicked.connect(lambda: self.set_mask_mode('threshold'))
        self.threshold_mask_btn.setStyleSheet("background-color: #6c757d;")
        mask_control_layout.addWidget(self.threshold_mask_btn)
        
        # 退出Mask模式按钮
        self.exit_mask_btn = QPushButton("退出Mask模式")
        self.exit_mask_btn.clicked.connect(lambda: self.set_mask_mode(None))
        self.exit_mask_btn.setEnabled(False)
        mask_control_layout.addWidget(self.exit_mask_btn)
        
        layout.addWidget(mask_control_widget)
        
        # 原始图像与本底图像显示（水平平分）
        top_split_widget = QWidget()
        top_split_layout = QHBoxLayout(top_split_widget)
        
        # 左侧：原始图像
        self.left_image_widget = QWidget()
        self.left_image_layout = QVBoxLayout(self.left_image_widget)
        self.left_image_label = QLabel("原始图像（Mask操作区）")
        self.left_image_label.setAlignment(Qt.AlignCenter)
        self.left_image_layout.addWidget(self.left_image_label)
        
        self.left_image_display = QLabel()
        self.left_image_display.setAlignment(Qt.AlignCenter)
        self.left_image_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.left_image_display.mousePressEvent = self.on_image_mouse_press
        self.left_image_display.mouseMoveEvent = self.on_image_mouse_move
        self.left_image_display.mouseReleaseEvent = self.on_image_mouse_release
        self.left_scroll_area = QScrollArea()
        self.left_scroll_area.setWidgetResizable(True)
        self.left_scroll_area.setWidget(self.left_image_display)
        self.left_image_layout.addWidget(self.left_scroll_area)
        top_split_layout.addWidget(self.left_image_widget)
        
        # 右侧：本底图像
        self.right_image_widget = QWidget()
        self.right_image_layout = QVBoxLayout(self.right_image_widget)
        self.right_image_label = QLabel("探测器本底图像 (右键选择)")
        self.right_image_label.setAlignment(Qt.AlignCenter)
        self.right_image_layout.addWidget(self.right_image_label)
        
        self.right_image_display = QLabel()
        self.right_image_display.setAlignment(Qt.AlignCenter)
        self.right_image_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.right_image_display.setContextMenuPolicy(Qt.CustomContextMenu)
        self.right_image_display.customContextMenuRequested.connect(self.show_background_menu)
        self.right_scroll_area = QScrollArea()
        self.right_scroll_area.setWidgetResizable(True)
        self.right_scroll_area.setWidget(self.right_image_display)
        self.right_image_layout.addWidget(self.right_scroll_area)
        top_split_layout.addWidget(self.right_image_widget)
        
        top_split_layout.setStretch(0, 1)
        top_split_layout.setStretch(1, 1)
        layout.addWidget(top_split_widget)
        
        # 上下分区分割器
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(top_split_widget)
        
        # 曲线分析区
        self.curve_widget = QWidget()
        self.curve_layout = QVBoxLayout(self.curve_widget)
        
        # 曲线标题
        self.curve_title = QLabel("散射曲线（每列对应一个径向距离的角度平均 | 横坐标：散射矢量q）")
        self.curve_title.setAlignment(Qt.AlignCenter)
        self.curve_layout.addWidget(self.curve_title)
        
        # 合并参数调节区和保存按钮为同一行
        self.control_parameters_widget = QWidget()
        self.control_parameters_layout = QHBoxLayout(self.control_parameters_widget)
        self.control_parameters_layout.setContentsMargins(0, 5, 0, 10)
        
        # q变换参数调节区（能量E和距离L）
        self.q_param_widget = QWidget()
        self.q_param_layout = QHBoxLayout(self.q_param_widget)
        self.q_param_layout.setContentsMargins(0, 0, 0, 0)
        
        # 能量E调节（支持0.1eV步长）
        self.energy_label = QLabel("能量 E (eV):")
        self.energy_spin = QDoubleSpinBox()
        self.energy_spin.setMinimum(100.0)
        self.energy_spin.setMaximum(10000.0)
        self.energy_spin.setSingleStep(0.1)
        self.energy_spin.setValue(self.q_energy)
        self.energy_spin.valueChanged.connect(self.on_q_param_changed)
        self.energy_spin.setToolTip("用于散射矢量q计算的X射线能量（步长0.1eV）")
        
        # 距离L调节（支持1mm步长）
        self.distance_label = QLabel("探测器距离 L (mm):")
        self.distance_spin = QSpinBox()
        self.distance_spin.setMinimum(10)
        self.distance_spin.setMaximum(1000)
        self.distance_spin.setSingleStep(1)
        self.distance_spin.setValue(self.q_distance)
        self.distance_spin.valueChanged.connect(self.on_q_param_changed)
        self.distance_spin.setToolTip("样品到探测器的物理距离（步长1mm）")
        
        self.q_param_layout.addWidget(self.energy_label)
        self.energy_spin.setFixedWidth(120)
        self.q_param_layout.addWidget(self.energy_spin)
        self.q_param_layout.addSpacing(15)
        self.q_param_layout.addWidget(self.distance_label)
        self.distance_spin.setFixedWidth(100)
        self.q_param_layout.addWidget(self.distance_spin)
        
        # 裁剪控件
        self.crop_control_widget = QWidget()
        self.crop_control_layout = QHBoxLayout(self.crop_control_widget)
        self.crop_control_layout.setContentsMargins(0, 0, 0, 0)
        
        self.crop_front_label = QLabel("前n个点设为NaN:")
        self.crop_front_spin = QSpinBox()
        self.crop_front_spin.setMinimum(0)
        self.crop_front_spin.setMaximum(10000)
        self.crop_front_spin.setValue(self.crop_front_nan)
        self.crop_front_spin.valueChanged.connect(self.on_crop_param_changed)
        
        self.crop_back_label = QLabel("后m个点设为NaN:")
        self.crop_back_spin = QSpinBox()
        self.crop_back_spin.setMinimum(0)
        self.crop_back_spin.setMaximum(10000)
        self.crop_back_spin.setValue(self.crop_back_nan)
        self.crop_back_spin.valueChanged.connect(self.on_crop_param_changed)
        
        self.crop_control_layout.addWidget(self.crop_front_label)
        self.crop_front_spin.setFixedWidth(80)
        self.crop_control_layout.addWidget(self.crop_front_spin)
        self.crop_control_layout.addSpacing(15)
        self.crop_control_layout.addWidget(self.crop_back_label)
        self.crop_back_spin.setFixedWidth(80)
        self.crop_control_layout.addWidget(self.crop_back_spin)
        
        # 保存数据按钮
        self.save_data_btn = QPushButton("保存数据")
        self.save_data_btn.clicked.connect(self.save_data)
        self.save_data_btn.setToolTip("保存当前分析数据（差值矩阵、曲线数据等）")
        
        # 将所有控件添加到同一行布局
        self.control_parameters_layout.addWidget(self.q_param_widget)
        self.control_parameters_layout.addSpacing(20)
        self.control_parameters_layout.addWidget(self.crop_control_widget)
        self.control_parameters_layout.addSpacing(20)
        self.control_parameters_layout.addWidget(self.save_data_btn)
        self.control_parameters_layout.addStretch(1)
        
        self.curve_layout.addWidget(self.control_parameters_widget)
        
        # Matplotlib曲线画布
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.curve_layout.addWidget(self.canvas)
        
        # 日志输出区
        self.log_text_edit = QTextEdit()
        self.log_text_edit.setPlaceholderText("操作日志显示区域")
        self.log_text_edit.setMaximumHeight(100)
        self.curve_layout.addWidget(self.log_text_edit)
        
        splitter.addWidget(self.curve_widget)
        splitter.setSizes([1, 1])
        layout.addWidget(splitter)

    def create_third_column(self):
        """第三栏：参数设置与差值结果显示（上下同宽）"""
        self.third_column = QFrame()
        self.third_column.setFrameStyle(QFrame.StyledPanel)
        # 关键：让第三栏自身可以横向扩展
        self.third_column.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
    
        main_layout = QVBoxLayout(self.third_column)
    
        # -------------------------- 1. 参数区 --------------------------
        top_widget = QWidget()          # 不再设置固定宽度
        top_layout = QGridLayout(top_widget)
        top_layout.setContentsMargins(0, 0, 0, 0)
        top_layout.setVerticalSpacing(8)
    
        # 字体大小
        top_layout.addWidget(QLabel("界面字体大小:"), 0, 0)
        self.font_size_spin = QSpinBox()
        self.font_size_spin.setRange(4, 14)
        self.font_size_spin.setValue(self.current_font_size)
        self.font_size_spin.setToolTip("调整界面所有元素的字体大小（4-14pt）")
        top_layout.addWidget(self.font_size_spin, 0, 1)
        self.apply_font_btn = QPushButton("应用字体大小")
        self.apply_font_btn.setFixedWidth(240)
        self.apply_font_btn.clicked.connect(self.change_font_size)
        top_layout.addWidget(self.apply_font_btn, 0, 2, 1, 2)
    
        # 能量 + 距离
        top_layout.addWidget(QLabel("能量 (eV)"), 1, 0)
        self.energy_edit = QLineEdit()
        self.energy_edit.setPlaceholderText("能量值")
        top_layout.addWidget(self.energy_edit, 1, 1)
        top_layout.addWidget(QLabel("样品到探测器距离 (mm)"), 1, 2)
        self.distance_edit = QLineEdit()
        self.distance_edit.setPlaceholderText("距离值")
        top_layout.addWidget(self.distance_edit, 1, 3)
    
        # 行号 + 行展宽
        top_layout.addWidget(QLabel("行号（y坐标）"), 2, 0)
        self.row_edit = QLineEdit()
        self.row_edit.setPlaceholderText("输入行号")
        top_layout.addWidget(self.row_edit, 2, 1)
        top_layout.addWidget(QLabel("行展宽"), 2, 2)
        self.row_width_edit = QLineEdit()
        self.row_width_edit.setPlaceholderText("行展宽")
        self.row_width_edit.setText("50")
        top_layout.addWidget(self.row_width_edit, 2, 3)
    
        # 列号 + 列展宽
        top_layout.addWidget(QLabel("列号（x坐标）"), 3, 0)
        self.col_edit = QLineEdit()
        self.col_edit.setPlaceholderText("输入列号")
        top_layout.addWidget(self.col_edit, 3, 1)
        top_layout.addWidget(QLabel("列展宽"), 3, 2)
        self.col_width_edit = QLineEdit()
        self.col_width_edit.setPlaceholderText("列展宽")
        self.col_width_edit.setText("50")
        top_layout.addWidget(self.col_width_edit, 3, 3)
    
        # 水平/竖直中心按钮
        self.calc_horizontal_btn = QPushButton("计算水平中心（列）")
        self.calc_horizontal_btn.clicked.connect(self.get_horizencenter)
        top_layout.addWidget(self.calc_horizontal_btn, 4, 0, 1, 2)
        self.calc_vertical_btn = QPushButton("计算竖直中心（行）")
        self.calc_vertical_btn.clicked.connect(self.get_verticalcenter)
        top_layout.addWidget(self.calc_vertical_btn, 4, 2, 1, 2)
    
        # 中心结果显示
        self.horizontal_result_label = QLineEdit()
        self.horizontal_result_label.setPlaceholderText("水平中心(1-based)")
        self.horizontal_result_label.editingFinished.connect(self.on_horizontal_center_edit)
        top_layout.addWidget(self.horizontal_result_label, 5, 0, 1, 2)
        self.vertical_result_label = QLineEdit()
        self.vertical_result_label.setPlaceholderText("竖直中心(1-based)")
        self.vertical_result_label.editingFinished.connect(self.on_vertical_center_edit)
        top_layout.addWidget(self.vertical_result_label, 5, 2, 1, 2)
    
        # 清除中心线
        self.clear_center_btn = QPushButton("清除中心行列线")
        self.clear_center_btn.clicked.connect(self.clear_center_lines)
        self.clear_center_btn.setEnabled(False)
        top_layout.addWidget(self.clear_center_btn, 6, 0, 1, 4)
    
        # 探测器像素 + Binning
        top_layout.addWidget(QLabel("探测器像素"), 7, 0)
        self.detector_pixel_edit = QLineEdit()
        self.detector_pixel_edit.setPlaceholderText("输入像素尺寸(μm)")
        self.detector_pixel_edit.setText("15")
        self.detector_pixel_edit.textChanged.connect(self.calculate_equivalent_pixel)
        top_layout.addWidget(self.detector_pixel_edit, 7, 1)
        top_layout.addWidget(QLabel("Binning模式"), 7, 2)
        self.binning_combo = QComboBox()
        self.binning_combo.addItems(["1*1", "2*2", "4*4"])
        self.binning_combo.currentTextChanged.connect(self.calculate_equivalent_pixel)
        top_layout.addWidget(self.binning_combo, 7, 3)
    
        # 本底系数 + 等效像素值
        top_layout.addWidget(QLabel("本底系数"), 8, 0)
        self.background_coeff_edit = QLineEdit()
        self.background_coeff_edit.setPlaceholderText("本底缩放系数")
        self.background_coeff_edit.setText("1.0")
        self.background_coeff_edit.textChanged.connect(self.update_background_coeff)
        top_layout.addWidget(self.background_coeff_edit, 8, 1)
        top_layout.addWidget(QLabel("等效像素值"), 8, 2)
        self.equivalent_pixel_edit = QLineEdit()
        self.equivalent_pixel_edit.setPlaceholderText("等效像素尺寸(μm)")
        self.equivalent_pixel_edit.setReadOnly(True)
        top_layout.addWidget(self.equivalent_pixel_edit, 8, 3)
    
        # 积分角度范围
        top_layout.addWidget(QLabel("积分起始角度(°)"), 9, 0)
        self.angle_start_edit = QLineEdit()
        self.angle_start_edit.setPlaceholderText("0-360")
        self.angle_start_edit.setText("0")
        top_layout.addWidget(self.angle_start_edit, 9, 1)
        top_layout.addWidget(QLabel("积分终止角度(°)"), 9, 2)
        self.angle_end_edit = QLineEdit()
        self.angle_end_edit.setPlaceholderText("0-360")
        self.angle_end_edit.setText("360")
        top_layout.addWidget(self.angle_end_edit, 9, 3)
    
        # 极坐标展开按钮
        self.polar_expand_btn = QPushButton("极坐标展开与差值计算")
        self.polar_expand_btn.clicked.connect(self.polar_coordinate_expand_with_diff)
        top_layout.addWidget(self.polar_expand_btn, 10, 0, 1, 4)
    
        main_layout.addWidget(top_widget)
    
        # -------------------------- 2. 差值结果显示 --------------------------
        self.diff_title_label = QLabel("极坐标展开差值结果（数据 - 系数×本底）")
        self.diff_title_label.setAlignment(Qt.AlignCenter)
        main_layout.addWidget(self.diff_title_label)
    
        diff_container = QWidget()      # 不再设置固定宽度
        diff_container_layout = QVBoxLayout(diff_container)
        diff_container_layout.setContentsMargins(0, 0, 0, 0)
    
        self.diff_scroll_area = QScrollArea()
        self.diff_scroll_area.setWidgetResizable(True)
        self.diff_image_label = QLabel()
        self.diff_image_label.setAlignment(Qt.AlignCenter)
        self.diff_image_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.diff_scroll_area.setWidget(self.diff_image_label)
        diff_container_layout.addWidget(self.diff_scroll_area)
        main_layout.addWidget(diff_container)
    
        # 比例分配
        main_layout.setStretch(0, 1)   # 参数区
        main_layout.setStretch(1, 0)   # 标题
        main_layout.setStretch(2, 5)   # 差值图

    # -------------------------- 原有：选项卡1（原界面）其他方法（复用，未修改） --------------------------
    def display_image(self, file_path):
        try:
            pil_image = Image.open(file_path)
            ext = os.path.splitext(file_path)[1].lower()
            is_grayscale = pil_image.mode in ['L', 'I', 'I;16', 'I;16L', 'I;16B'] or \
                          (ext in ['.tif', '.tiff'] and len(pil_image.getbands()) == 1)
            
            img_array = np.array(pil_image)
            self.current_image_array = img_array
            self.original_image_matrix = img_array.astype(np.float32).copy()
            self.image_matrix = self.original_image_matrix.copy()
            self.horizontal_center = self.vertical_center = None
            self.horizontal_result_label.clear()
            self.vertical_result_label.clear()
            self.diff_image_label.clear()
            self.diff_matrix = self.polar_expand_matrix = None
            
            self.crop_front_nan = self.crop_back_nan = 0
            self.crop_front_spin.setValue(0)
            self.crop_back_spin.setValue(0)
            self.crop_front_spin.setMaximum(10000)
            self.crop_back_spin.setMaximum(10000)
            
            self.masked_regions.clear()
            self.clear_all_mask_btn.setEnabled(False)
            if self.mask_mode is not None:
                self.set_mask_mode(None)
            
            self.clear_center_btn.setEnabled(False)
            
            height, width = img_array.shape[:2]
            self.row_edit.setText(str(int(height / 2)))
            self.col_edit.setText(str(int(width / 2)))
            
            self.clear_curve()
            
            if self.background_image_array is not None:
                orig_h, orig_w = img_array.shape[:2]
                bg_h, bg_w = self.background_image_array.shape[:2]
                if orig_h != bg_h or orig_w != bg_w:
                    QMessageBox.warning(self, "尺寸不匹配", 
                                      f"原始与本底图像尺寸不同！\n原始: {orig_w}×{orig_h} (宽×高)\n本底: {bg_w}×{bg_h} (宽×高)")
            
            if is_grayscale:
                self.current_image_mode = 'grayscale'
                max_val = np.max(img_array)
                
                if not self.contrast_manual_set:
                    self.set_closest_combo_value(self.lower_combo, max_val / 5)
                    self.set_closest_combo_value(self.upper_combo, max_val / 2)
                
                self.apply_contrast()
            else:
                self.current_image_mode = 'color'
                if pil_image.mode != 'RGB':
                    pil_image = pil_image.convert('RGB')
                    img_array = np.array(pil_image)
                q_image = QImage(img_array.data, pil_image.width, pil_image.height, 
                               3 * pil_image.width, QImage.Format_RGB888)
                self.display_left_pixmap_with_center(QPixmap.fromImage(q_image))
            
            self.log_print(f"加载图像：{file_path}，尺寸：{img_array.shape[1]}×{img_array.shape[0]} (宽×高)，模式：{'灰度' if is_grayscale else '彩色'}")
            self.log_print(f"对比度设置：下限={self.lower_combo.currentText()}, 上限={self.upper_combo.currentText()} ({'使用当前值' if self.contrast_manual_set else '自动设置'})")
        except Exception as e:
            self.log_print(f"加载图像失败: {str(e)}")

    def change_font_size(self):
        new_size = self.font_size_spin.value()
        if new_size == self.current_font_size:
            return
            
        self.current_font_size = new_size
        
        self.file_list.setFont(QFont("Arial", new_size))
        self.curve_list_tab2.setFont(QFont("Arial", new_size))
        
        self.apply_styles()
        
        if hasattr(self, 'curve_data') and self.curve_data is not None:
            self.plot_col_average_curve()
        else:
            self.clear_curve()
            
        self.update()
        
        self.log_print(f"界面字体大小已调整为: {new_size}pt")

    def clear_center_lines(self):
        if self.horizontal_center is None and self.vertical_center is None:
            QMessageBox.information(self, "提示", "没有可清除的中心行列线")
            return
            
        self.horizontal_center = None
        self.vertical_center = None
        
        self.horizontal_result_label.clear()
        self.vertical_result_label.clear()
        
        self.clear_center_btn.setEnabled(False)
        
        self.redraw_images_with_center_lines()
        
        self.clear_analysis_results()
        
        self.log_print("已清除光斑中心行列线")

    def on_horizontal_center_edit(self):
        if self.image_matrix is None:
            QMessageBox.warning(self, "警告", "请先加载图像")
            self.horizontal_result_label.clear()
            return
        try:
            input_val = self.horizontal_result_label.text().strip()
            if not input_val:
                self.horizontal_center = None
                self.clear_center_lines()
                return
            center_col_1based = int(input_val)
            total_cols = self.image_matrix.shape[1]
            if center_col_1based < 1 or center_col_1based > total_cols:
                raise ValueError(f"必须在1-{total_cols}之间（当前图像共{total_cols}列）")
            self.horizontal_center = center_col_1based - 1
            self.log_print(f"手动设置水平中心：{center_col_1based} (1-based) → 列索引{self.horizontal_center} (0-based)")
            self.clear_center_btn.setEnabled(True)
            self.redraw_images_with_center_lines()
            self.clear_analysis_results()
        except ValueError as e:
            QMessageBox.warning(self, "输入错误", f"水平中心输入无效：{str(e)}")
            if self.horizontal_center is not None:
                self.horizontal_result_label.setText(f"{self.horizontal_center + 1}")
            else:
                self.horizontal_result_label.clear()

    def on_vertical_center_edit(self):
        if self.image_matrix is None:
            QMessageBox.warning(self, "警告", "请先加载图像")
            self.vertical_result_label.clear()
            return
        try:
            input_val = self.vertical_result_label.text().strip()
            if not input_val:
                self.vertical_center = None
                self.clear_center_lines()
                return
            center_row_1based = int(input_val)
            total_rows = self.image_matrix.shape[0]
            if center_row_1based < 1 or center_row_1based > total_rows:
                raise ValueError(f"必须在1-{total_rows}之间（当前图像共{total_rows}行）")
            self.vertical_center = center_row_1based - 1
            self.log_print(f"手动设置竖直中心：{center_row_1based} (1-based) → 行索引{self.vertical_center} (0-based)")
            self.clear_center_btn.setEnabled(True)
            self.redraw_images_with_center_lines()
            self.clear_analysis_results()
        except ValueError as e:
            QMessageBox.warning(self, "输入错误", f"竖直中心输入无效：{str(e)}")
            if self.vertical_center is not None:
                self.vertical_result_label.setText(f"{self.vertical_center + 1}")
            else:
                self.vertical_result_label.clear()

    
    
    # 添加阈值设置方法
    def set_mask_threshold(self, value):
        self.mask_threshold = value
        # 如果当前是阈值mask模式，立即应用
        if self.mask_mode == 'threshold':
            self.apply_threshold_mask()
    
    
    # 修改set_mask_mode方法以支持阈值mask模式
    def set_mask_mode(self, mode):
        self.mask_mode = mode
        self.mask_active = False
        self.mask_start_point = None
        self.mask_current_point = None
        
        # 更新按钮样式
        self.circle_mask_btn.setStyleSheet("background-color: #6c757d;" if mode != 'circle' else "background-color: #28a745;")
        self.square_mask_btn.setStyleSheet("background-color: #6c757d;" if mode != 'square' else "background-color: #28a745;")
        self.threshold_mask_btn.setStyleSheet("background-color: #6c757d;" if mode != 'threshold' else "background-color: #28a745;")  # 新增
        self.exit_mask_btn.setEnabled(mode is not None)
        self.clear_all_mask_btn.setEnabled(len(self.masked_regions) > 0)
        
        # 如果是阈值mask模式，立即应用
        if mode == 'threshold':
            self.apply_threshold_mask()
    
    
    # 添加阈值Mask应用方法
    def apply_threshold_mask(self):
        if self.original_image_matrix is None:
            return
            
        # 创建阈值mask
        mask = self.original_image_matrix < self.mask_threshold
        self.image_matrix = self.original_image_matrix.copy()
        self.image_matrix[mask] = np.nan
        
        # 记录mask区域
        self.masked_regions.append(('threshold', {'value': self.mask_threshold}))
        
        # 更新图像显示
        self.update_image_display()
        self.clear_all_mask_btn.setEnabled(True)
    
    
    # 修改clear_all_masks方法以支持清除阈值mask
    def clear_all_masks(self):
        if self.original_image_matrix is not None:
            self.image_matrix = self.original_image_matrix.copy()
        self.masked_regions = []
        self.update_image_display()
        self.clear_all_mask_btn.setEnabled(False)
    
    
    # 添加自动检测表头的函数
    def detect_header_lines(self, file_path):
        """检测文件中表头的行数"""
        header_lines = 0
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:  # 跳过空行
                    header_lines += 1
                    continue
                # 尝试将行按空白分割并转换为数值
                parts = re.split(r'\s+', line)
                try:
                    # 尝试转换为数值
                    [float(part) for part in parts]
                    # 如果成功转换，说明不是表头
                    break
                except ValueError:
                    # 转换失败，认为是表头
                    header_lines += 1
        return header_lines
    
    
    
    def load_analyze_curve_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "选择一维曲线文件夹")
        if not folder:
            return
    
        import glob
        files = glob.glob(os.path.join(folder, "*.txt")) + \
                glob.glob(os.path.join(folder, "*.dat"))
        if not files:
            QMessageBox.information(self, "提示", "文件夹内未找到 txt/dat 文件")
            return
    
        loaded = 0
        for fp in files:
            try:
                # 自动检测表头行数
                skip_rows = self.detect_header_lines(fp)
                
                data = np.loadtxt(fp, skiprows=skip_rows)
                if data.shape[0] < 2:
                    continue          # 跳过无效文件
                if data.shape[1] >= 2:
                    q, I = data[:, 0], data[:, 1]
                else:
                    I = data[:, 0]
                    q = np.arange(len(I))
    
                label = f"曲线_{len(self.analyze_curve_list)+1}_{os.path.splitext(os.path.basename(fp))[0]}"
                self.analyze_curve_list.append({"file": fp, "q": q, "intensity": I, "label": label})
                self.curve_list_tab2.addItem(f"{label} - {os.path.basename(fp)} ({len(q)}点, 跳过{skip_rows}行表头)")
                loaded += 1
            except Exception as e:
                self.log_print(f"跳过文件 {fp} ：{e}")
        self.log_print(f"文件夹扫描完成，共加载 {loaded} 条曲线")



    def on_contrast_manual_change(self):
        self.contrast_manual_set = True
        self.log_print("对比度已手动调整，后续加载图像将保持当前设置")

    def update_background_cmap(self):
        if self.background_file_path is not None:
            try:
                self.load_background_image(self.background_file_path)
                self.log_print(f"本底图像颜色图已更新为：{self.cmap_combo.currentText()}")
            except Exception as e:
                self.log_print(f"更新本底图像颜色图失败: {str(e)}")

    def update_diff_cmap(self):
        if self.diff_matrix is not None:
            try:
                start_angle = float(self.angle_start_edit.text().strip())
                end_angle = float(self.angle_end_edit.text().strip())
                angle_range = np.arange(start_angle, end_angle + 1, 1)
                img_h, img_w = self.image_matrix.shape
                max_radius = int(round(math.sqrt(img_w**2 + img_h**2)))
                self.display_diff_result(angle_range, max_radius)
                self.log_print(f"差值矩阵颜色图已更新为：{self.cmap_combo.currentText()}")
            except Exception as e:
                self.log_print(f"更新差值矩阵颜色图失败: {str(e)}")

    def set_mask_mode(self, mode):
        self.mask_mode = mode
        self.mask_active = False
        
        if mode == 'circle':
            self.circle_mask_btn.setStyleSheet("background-color: #dc3545;")
            self.square_mask_btn.setStyleSheet("background-color: #6c757d;")
            self.exit_mask_btn.setEnabled(True)
            self.left_image_display.setCursor(QCursor(Qt.CrossCursor))
            QMessageBox.information(self, "Mask模式", "已进入圆形Mask模式：\n在图像上拖拽鼠标绘制圆形区域（松开后设为NaN）")
            self.log_print("已进入圆形Mask模式：在图像上拖拽创建圆形区域（松开后设为NaN）")
        elif mode == 'square':
            self.square_mask_btn.setStyleSheet("background-color: #dc3545;")
            self.circle_mask_btn.setStyleSheet("background-color: #6c757d;")
            self.exit_mask_btn.setEnabled(True)
            self.left_image_display.setCursor(QCursor(Qt.CrossCursor))
            QMessageBox.information(self, "Mask模式", "已进入方形Mask模式：\n在图像上拖拽鼠标绘制方形区域（松开后设为NaN）")
            self.log_print("已进入方形Mask模式：在图像上拖拽创建方形区域（松开后设为NaN）")
        else:
            self.circle_mask_btn.setStyleSheet("background-color: #6c757d;")
            self.square_mask_btn.setStyleSheet("background-color: #6c757d;")
            self.exit_mask_btn.setEnabled(False)
            self.left_image_display.setCursor(QCursor(Qt.ArrowCursor))
            if self.mask_mode is not None:
                self.log_print("已退出Mask模式")
        
        self.set_other_buttons_state(mode is None)

    def set_other_buttons_state(self, enabled):
        buttons = [
            self.select_folder_btn, self.calc_horizontal_btn, self.calc_vertical_btn,
            self.polar_expand_btn, self.save_data_btn, self.clear_center_btn,
            self.apply_font_btn, self.load_curve_btn_tab2
        ]
        for btn in buttons:
            btn.setEnabled(enabled)
        
        widgets = [
            self.lower_combo, self.upper_combo, self.cmap_combo,
            self.energy_spin, self.distance_spin,
            self.crop_front_spin, self.crop_back_spin, self.energy_edit, self.distance_edit,
            self.row_edit, self.row_width_edit, self.col_edit, self.col_width_edit,
            self.detector_pixel_edit, self.binning_combo, self.background_coeff_edit,
            self.angle_start_edit, self.angle_end_edit,
            self.font_size_spin, self.horizontal_result_label, self.vertical_result_label
        ]
        for widget in widgets:
            widget.setEnabled(enabled)

    def on_image_mouse_press(self, event):
        if self.mask_mode is None:
            return
            
        if self.current_image_array is None:
            QMessageBox.warning(self, "无图像", "请先加载原始图像再使用Mask功能")
            self.log_print("错误：未加载图像，无法绘制Mask")
            return
            
        if event.button() == Qt.LeftButton:
            widget_pos = event.pos()
            self.mask_start_point = self.widget_to_image_coords(widget_pos)
            self.mask_current_point = self.mask_start_point
            self.mask_active = True

    def on_image_mouse_move(self, event):
        if not self.mask_active or self.mask_mode is None or self.current_image_array is None:
            return
            
        self.mask_current_point = self.widget_to_image_coords(event.pos())
        self.redraw_image_with_mask_preview()

    def on_image_mouse_release(self, event):
        if not self.mask_active or self.mask_mode is None or self.current_image_array is None:
            return
            
        if event.button() == Qt.LeftButton:
            self.mask_active = False
            end_point = self.widget_to_image_coords(event.pos())
            
            start_x, start_y = self.mask_start_point
            end_x, end_y = end_point
            
            img_height, img_width = self.image_matrix.shape[:2]
            
            start_x = max(0, min(start_x, img_width - 1))
            start_y = max(0, min(start_y, img_height - 1))
            end_x = max(0, min(end_x, img_width - 1))
            end_y = max(0, min(end_y, img_height - 1))
            
            if self.mask_mode == 'circle':
                self.apply_circle_mask(start_x, start_y, end_x, end_y)
            elif self.mask_mode == 'square':
                self.apply_square_mask(start_x, start_y, end_x, end_y)
            
            self.redraw_images_with_center_lines()
            self.clear_all_mask_btn.setEnabled(len(self.masked_regions) > 0)

    def widget_to_image_coords(self, widget_pos):
        """将控件坐标转换为原始图像坐标"""
        pixmap = self.left_image_display.pixmap()
        if pixmap.isNull() or self.current_image_array is None:
            return (0, 0)
        
        orig_h, orig_w = self.current_image_array.shape[:2]
        scaled_w, scaled_h = pixmap.width(), pixmap.height()
        widget_w, widget_h = self.left_image_display.width(), self.left_image_display.height()
        
        scale_x = orig_w / scaled_w if scaled_w > 0 else 1.0
        scale_y = orig_h / scaled_h if scaled_h > 0 else 1.0
        
        pixmap_offset_x = (widget_w - scaled_w) // 2
        pixmap_offset_y = (widget_h - scaled_h) // 2
        scaled_x = widget_pos.x() - pixmap_offset_x
        scaled_y = widget_pos.y() - pixmap_offset_y
        scaled_x = max(0, min(scaled_x, scaled_w - 1))
        scaled_y = max(0, min(scaled_y, scaled_h - 1))
        
        scroll_h = self.left_scroll_area.horizontalScrollBar()
        scroll_v = self.left_scroll_area.verticalScrollBar()
        scroll_x_orig = scroll_h.value() * scale_x
        scroll_y_orig = scroll_v.value() * scale_y
        
        img_x = int(scaled_x * scale_x + scroll_x_orig)
        img_y = int(scaled_y * scale_y + scroll_y_orig)
        
        img_x = max(0, min(img_x, orig_w - 1))
        img_y = max(0, min(img_y, orig_h - 1))
        
        return (img_x, img_y)

    def redraw_image_with_mask_preview(self):
        """重新绘制图像并叠加当前Mask预览"""
        if self.current_image_array is None or not self.mask_active:
            return
            
        if self.current_image_mode == 'grayscale':
            lower = int(self.lower_combo.currentText())
            upper = int(self.upper_combo.currentText())
            img = self.current_image_array.astype(np.float32)
            img = np.clip(img, lower, upper)
            img = (img - lower) / (upper - lower) * 255 if upper > lower else np.zeros_like(img)
            img = np.clip(img, 0, 255).astype(np.uint8)
            colored_array = (cm.get_cmap(self.cmap_combo.currentText())(img) * 255).astype(np.uint8)
            q_image = QImage(colored_array.data, colored_array.shape[1], colored_array.shape[0], 
                           4 * colored_array.shape[1], QImage.Format_RGBA8888)
        else:
            pil_image = Image.fromarray(self.current_image_array)
            if pil_image.mode != 'RGB':
                pil_image = pil_image.convert('RGB')
                img_array = np.array(pil_image)
            q_image = QImage(self.current_image_array.data, pil_image.width, pil_image.height, 
                           3 * pil_image.width, QImage.Format_RGB888)
        
        pixmap = QPixmap.fromImage(q_image)
        painter = QPainter(pixmap)
        
        start_x, start_y = self.mask_start_point
        end_x, end_y = self.mask_current_point
        
        pen = QPen(QColor(255, 0, 0, 180))
        pen.setWidth(2)
        pen.setStyle(Qt.DashLine)
        painter.setPen(pen)
        
        if self.mask_mode == 'circle':
            center_x = (start_x + end_x) / 2
            center_y = (start_y + end_y) / 2
            radius = math.hypot(end_x - center_x, end_y - center_y)
            painter.drawEllipse(QPoint(int(center_x), int(center_y)), int(radius), int(radius))
        elif self.mask_mode == 'square':
            x1 = min(start_x, end_x)
            y1 = min(start_y, end_y)
            x2 = max(start_x, end_x)
            y2 = max(start_y, end_y)
            painter.drawRect(x1, y1, x2 - x1, y2 - y1)
        
        self.draw_center_lines(painter, pixmap.width(), pixmap.height())
        
        painter.end()
        
        avail_size = self.left_scroll_area.size()
        scaled = pixmap.scaled(avail_size.width()-20, avail_size.height()-20, 
                             Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.left_image_display.setPixmap(scaled)
        self.left_image_display.update()

    def apply_circle_mask(self, start_x, start_y, end_x, end_y):
        """应用圆形Mask"""
        center_x = (start_x + end_x) / 2
        center_y = (start_y + end_y) / 2
        radius = math.hypot(end_x - center_x, end_y - center_y)
        
        if radius < 2:
            QMessageBox.warning(self, "区域过小", "忽略过小的圆形Mask区域（半径<2像素），请拖拽更大范围")
            self.log_print("忽略过小的圆形Mask区域（半径<2像素）")
            return
        
        img_height, img_width = self.image_matrix.shape[:2]
        y, x = np.ogrid[0:img_height, 0:img_width]
        distance = np.sqrt((x - center_x)** 2 + (y - center_y)**2)
        mask = distance <= radius
        self.image_matrix[mask] = np.nan
        
        self.masked_regions.append(('circle', {
            'center_x': center_x,
            'center_y': center_y,
            'radius': radius,
            'count': np.sum(mask)
        }))
        
        self.log_print(f"已应用圆形Mask：圆心({int(center_x)},{int(center_y)})，半径{int(radius)}像素，覆盖{np.sum(mask)}个像素")
        self.clear_analysis_results()

    def apply_square_mask(self, start_x, start_y, end_x, end_y):
        """应用方形Mask"""
        x_min = int(min(start_x, end_x))
        x_max = int(max(start_x, end_x))
        y_min = int(min(start_y, end_y))
        y_max = int(max(start_y, end_y))
        
        if (x_max - x_min) < 2 or (y_max - y_min) < 2:
            QMessageBox.warning(self, "区域过小", "忽略过小的方形Mask区域（<2x2像素），请拖拽更大范围")
            self.log_print("忽略过小的方形Mask区域（<2x2像素）")
            return
        
        x_min = max(0, x_min)
        x_max = min(self.image_matrix.shape[1] - 1, x_max)
        y_min = max(0, y_min)
        y_max = min(self.image_matrix.shape[0] - 1, y_max)
        
        self.image_matrix[y_min:y_max+1, x_min:x_max+1] = np.nan
        
        self.masked_regions.append(('square', {
            'x_min': x_min,
            'x_max': x_max,
            'y_min': y_min,
            'y_max': y_max,
            'width': x_max - x_min + 1,
            'height': y_max - y_min + 1,
            'count': (x_max - x_min + 1) * (y_max - y_min + 1)
        }))
        
        self.log_print(f"已应用方形Mask：区域({x_min},{y_min})-({x_max},{y_max})，尺寸{self.masked_regions[-1]['width']}×{self.masked_regions[-1]['height']}像素，覆盖{self.masked_regions[-1]['count']}个像素")
        self.clear_analysis_results()

    def clear_all_masks(self):
        """清除所有Mask"""
        if self.original_image_matrix is not None and self.image_matrix is not None:
            self.image_matrix = self.original_image_matrix.copy()
            mask_count = len(self.masked_regions)
            self.masked_regions.clear()
            self.apply_contrast()
            self.redraw_images_with_center_lines()
            self.clear_all_mask_btn.setEnabled(False)
            self.log_print(f"已清除所有({mask_count}个)Mask区域，恢复原始图像数据")
            self.clear_analysis_results()
        else:
            self.log_print("没有可清除的Mask区域")

    def clear_analysis_results(self):
        """清空现有分析结果"""
        if self.polar_expand_matrix is not None or self.diff_matrix is not None:
            self.polar_expand_matrix = None
            self.diff_matrix = None
            self.diff_image_label.clear()
            self.diff_title_label.setText("极坐标展开差值结果（数据 - 系数×本底）")
            self.clear_curve()
            self.log_print("已清空现有分析结果，请重新进行极坐标展开计算")

    def draw_center_lines(self, painter, pixmap_width, pixmap_height):
        """在指定painter上绘制中心线"""
        if self.image_matrix is None:
            return
            
        pen = QPen(QColor(255, 0, 0))
        pen.setWidth(2)
        painter.setPen(pen)
        
        if self.vertical_center is not None:
            y = self.vertical_center * pixmap_height / self.image_matrix.shape[0]
            painter.drawLine(0, int(y), pixmap_width, int(y))
        
        if self.horizontal_center is not None:
            x = self.horizontal_center * pixmap_width / self.image_matrix.shape[1]
            painter.drawLine(int(x), 0, int(x), pixmap_height)

    def extract_parameters_from_filename(self, filename):
        """从文件名中提取能量和距离参数"""
        energy = None
        distance = None
        
        base_name = os.path.splitext(filename)[0]
        
        energy_patterns = [
            r'(\d+\.?\d*)eV',
            r'E=?(\d+\.?\d*)',
            r'energy=?(\d+\.?\d*)',
            r'(\d+\.?\d*)\s*电子伏特'
        ]
        
        for pattern in energy_patterns:
            match = re.search(pattern, base_name, re.IGNORECASE)
            if match:
                try:
                    energy = float(match.group(1))
                    break
                except ValueError:
                    continue
        
        distance_patterns = [
            r'(\d+\.?\d*)mm',
            r'D=?(\d+\.?\d*)',
            r'dist=?(\d+\.?\d*)',
            r'距离=?(\d+\.?\d*)',
            r'(\d+\.?\d*)\s*毫米'
        ]
        
        for pattern in distance_patterns:
            match = re.search(pattern, base_name, re.IGNORECASE)
            if match:
                try:
                    distance = float(match.group(1))
                    break
                except ValueError:
                    continue
        
        return energy, distance

    def show_background_menu(self, position):
        """本底图像右键菜单"""
        menu = QMenu()
        select_action = QAction("选择探测器本底图像", self)
        select_action.triggered.connect(self.select_background_image)
        menu.addAction(select_action)
        
        if self.background_image_array is not None:
            clear_action = QAction("清除本底图像", self)
            clear_action.triggered.connect(self.clear_background_image)
            menu.addAction(clear_action)
            
        menu.exec_(self.right_image_display.mapToGlobal(position))

    def select_background_image(self):
        """选择本底图像"""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "选择探测器本底图像", "", "TIFF图像 (*.tif *.tiff);;所有图像 (*.png *.jpg *.jpeg *.bmp)"
        )
        if file_path:
            try:
                self.background_file_path = file_path  # 保存本底文件路径用于颜色图更新
                self.load_background_image(file_path)
            except Exception as e:
                QMessageBox.warning(self, "加载失败", f"无法加载本底图像: {str(e)}")
                self.log_print(f"本底图像加载失败: {str(e)}")

    def load_background_image(self, file_path):
        """加载并验证本底图像（使用当前选择的颜色图）"""
        pil_image = Image.open(file_path)
        img_array = np.array(pil_image)
        self.background_image_array = img_array
        self.background_matrix = img_array.astype(np.float32)
        
        is_grayscale = pil_image.mode in ['L', 'I', 'I;16', 'I;16L', 'I;16B'] or \
                      (os.path.splitext(file_path)[1].lower() in ['.tif', '.tiff'] and len(pil_image.getbands()) == 1)
        self.background_image_mode = 'grayscale' if is_grayscale else 'color'
        
        if is_grayscale:
            max_val, min_val = np.max(img_array), np.min(img_array)
            normalized = (img_array - min_val) / (max_val - min_val) * 255 if max_val > min_val else np.zeros_like(img_array)
            normalized = normalized.astype(np.uint8)
            # 使用当前选择的颜色图
            colored_array = (cm.get_cmap(self.cmap_combo.currentText())(normalized) * 255).astype(np.uint8)
            q_image = QImage(colored_array.data, colored_array.shape[1], colored_array.shape[0], 
                           4 * colored_array.shape[1], QImage.Format_RGBA8888)
        else:
            if pil_image.mode != 'RGB':
                pil_image = pil_image.convert('RGB')
                img_array = np.array(pil_image)
            q_image = QImage(img_array.data, pil_image.width, pil_image.height, 
                           3 * pil_image.width, QImage.Format_RGB888)
        
        self.display_background_pixmap_with_center(QPixmap.fromImage(q_image))
        
        if self.image_matrix is not None:
            orig_h, orig_w = self.image_matrix.shape[:2]
            bg_h, bg_w = img_array.shape[:2]
            if orig_h != bg_h or orig_w != bg_w:
                QMessageBox.warning(self, "尺寸不匹配", 
                                  f"本底图像尺寸与原始图像不同！\n原始: {orig_w}×{orig_h} (宽×高)\n本底: {bg_w}×{bg_h} (宽×高)")
                self.log_print(f"本底与原始图像尺寸不匹配（原始：{orig_w}×{orig_h}，本底：{bg_w}×{bg_h}）")
            else:
                self.log_print(f"本底图像加载成功（尺寸：{bg_w}×{bg_h}，颜色图：{self.cmap_combo.currentText()}，与原始图像匹配）")
        else:
            self.log_print(f"本底图像加载成功（尺寸：{img_array.shape[1]}×{img_array.shape[0]}，颜色图：{self.cmap_combo.currentText()}）")

    def clear_background_image(self):
        """清除本底图像"""
        self.background_image_array = self.background_matrix = None
        self.background_file_path = None  # 清除文件路径记录
        self.right_image_display.clear()
        self.log_print("本底图像已清除")

    def update_background_coeff(self):
        """更新本底系数"""
        try:
            coeff = float(self.background_coeff_edit.text().strip()) if self.background_coeff_edit.text().strip() else 1.0
            if coeff < 0:
                raise ValueError("系数不能为负数")
            self.background_coeff = coeff
            self.log_print(f"本底系数更新为: {coeff:.4f}")
        except ValueError as e:
            QMessageBox.warning(self, "输入错误", f"本底系数设置错误：{str(e)}")
            self.background_coeff_edit.setText(f"{self.background_coeff:.4f}")

    def calculate_equivalent_pixel(self):
        """计算等效像素尺寸"""
        try:
            pixel = float(self.detector_pixel_edit.text().strip()) if self.detector_pixel_edit.text().strip() else 0
            if pixel <= 0:
                raise ValueError("像素值必须为正数")
            binning = int(self.binning_combo.currentText().split("*")[0])
            eq_pixel = pixel * binning
            self.equivalent_pixel_edit.setText(f"{eq_pixel:.2f}")
            self.log_print(f"等效像素计算：{pixel}μm × {binning} = {eq_pixel:.2f}μm")
        except ValueError as e:
            self.equivalent_pixel_edit.setText("输入无效")
            QMessageBox.warning(self, "输入错误", f"探测器像素设置错误：{str(e)}")

    def select_folder(self):
        """选择图像文件夹"""
        folder_path = QFileDialog.getExistingDirectory(self, "选择文件夹")
        if folder_path:
            self.file_list.clear()
            self.current_folder = folder_path
            try:
                img_exts = ['.jpg', '.jpeg', '.png', '.bmp', '.gif', '.tif', '.tiff']
                for f in os.listdir(folder_path):
                    if os.path.isfile(os.path.join(folder_path, f)) and os.path.splitext(f)[1].lower() in img_exts:
                        self.file_list.addItem(f)
                self.log_print(f"加载文件夹：{folder_path}，找到{len(self.file_list)}个图像文件")
            except Exception as e:
                self.file_list.addItem(f"错误: {str(e)}")
                self.log_print(f"加载文件夹失败：{str(e)}")

    def on_file_double_clicked(self, item):
        """双击加载图像"""
        if hasattr(self, 'current_folder'):
            file_name = item.text()
            self.current_file_path = os.path.join(self.current_folder, file_name)
            
            energy, distance = self.extract_parameters_from_filename(file_name)
            self.energy = energy
            self.detector_distance = distance
            
            if energy is not None:
                self.energy_edit.setText(f"{energy:.2f}")
                self.energy_spin.setValue(energy)
                self.q_energy = energy
                self.log_print(f"从文件名提取能量: {energy:.2f} eV（同步更新q变换参数）")
            else:
                self.energy_edit.clear()
                self.log_print("未从文件名提取到能量信息（使用默认E=1000.0 eV）")
                
            if distance is not None:
                self.distance_edit.setText(f"{distance:.2f}")
                self.distance_spin.setValue(int(round(distance)))
                self.q_distance = int(round(distance))
                self.log_print(f"从文件名提取样品到探测器距离: {distance:.2f} mm（同步更新q变换参数）")
            else:
                self.distance_edit.clear()
                self.log_print("未从文件名提取到距离信息（使用默认L=200 mm）")
            
            self.display_image(self.current_file_path)

    def set_closest_combo_value(self, combo, value):
        """设置ComboBox最接近值"""
        items = [int(combo.itemText(i)) for i in range(combo.count())]
        combo.setCurrentText(str(min(items, key=lambda x: abs(x - value))))

    def apply_contrast(self):
        """应用对比度（使用当前选择的颜色图）"""
        if self.current_image_mode != 'grayscale' or self.current_image_array is None:
            return
        try:
            lower = int(self.lower_combo.currentText())
            upper = int(self.upper_combo.currentText())
            if lower >= upper:
                self.log_print("警告: 对比度下限必须小于上限")
                return
            
            img = self.current_image_array.astype(np.float32)
            img = np.clip(img, lower, upper)
            img = (img - lower) / (upper - lower) * 255
            img = np.clip(img, 0, 255).astype(np.uint8)
            # 使用当前选择的颜色图
            colored = (cm.get_cmap(self.cmap_combo.currentText())(img / 255.0) * 255).astype(np.uint8)
            q_image = QImage(colored.data, colored.shape[1], colored.shape[0], 
                           4 * colored.shape[1], QImage.Format_RGBA8888)
            self.display_left_pixmap_with_center(QPixmap.fromImage(q_image))
            self.log_print(f"应用对比度：下限{lower}，上限{upper}，颜色图：{self.cmap_combo.currentText()}")
        except Exception as e:
            self.log_print(f"对比度调整失败: {str(e)}")
    
    def get_horizencenter(self):
        """计算水平中心（列，x坐标）"""
        try:
            m = int(self.row_edit.text().strip())
            mc = int(self.row_width_edit.text().strip())
            if m <= 0 or mc <= 0:
                raise ValueError("行数和展宽必须为正整数")
            if self.image_matrix is None:
                raise ValueError("请先加载图像")
            
            base_row = m - 1
            total_rows = self.image_matrix.shape[0]
            if base_row < 0 or base_row >= total_rows:
                raise IndexError(f"行数超出范围（有效1-{total_rows}）")
            
            start_row = max(0, base_row - mc)
            end_row = min(total_rows - 1, base_row + mc)
            avg_row = np.mean(self.image_matrix[start_row:end_row+1, :], axis=0)
            
            if len(avg_row) < 4:
                raise ValueError("曲线长度不足，无法分析")
            
            corr_coeffs = []
            splits = []
            for split in range(1, len(avg_row)-1):
                p1 = avg_row[:split][-min(len(avg_row[:split]), len(avg_row[split+1:])):]
                p2 = avg_row[split+1:][:min(len(avg_row[:split]), len(avg_row[split+1:])):]
                if len(p1) >= 2 and len(p2) >= 2:
                    corr, _ = pearsonr(p1[::-1], p2)
                    corr_coeffs.append(corr)
                    splits.append(split)
            
            if not corr_coeffs:
                raise RuntimeError("无法计算有效相关系数")
            
            exclude = int(len(corr_coeffs) * 0.1)
            if len(corr_coeffs) > 2 * exclude:
                mid_corr = corr_coeffs[exclude:-exclude]
                mid_splits = splits[exclude:-exclude]
                center_col = mid_splits[np.argmax(mid_corr)]
            else:
                center_col = splits[np.argmax(corr_coeffs)]
            
            self.horizontal_center = center_col
            self.horizontal_result_label.setText(f"{center_col + 1}")
            self.log_print(f"水平中心计算完成：行{m}，展宽{mc}，中心列{center_col + 1} (1-based)")
            
            self.clear_center_btn.setEnabled(True)
            self.redraw_images_with_center_lines()
        except (ValueError, IndexError) as e:
            QMessageBox.warning(self, "输入错误", str(e))
            self.log_print(f"水平中心计算错误: {str(e)}")
        except Exception as e:
            QMessageBox.warning(self, "错误", str(e))
            self.log_print(f"水平中心计算错误: {str(e)}")
    
    def get_verticalcenter(self):
        """计算竖直中心（行，y坐标）"""
        try:
            n = int(self.col_edit.text().strip())
            nc = int(self.col_width_edit.text().strip())
            if n <= 0 or nc <= 0:
                raise ValueError("列数和展宽必须为正整数")
            if self.image_matrix is None:
                raise ValueError("请先加载图像")
            
            base_col = n - 1
            total_cols = self.image_matrix.shape[1]
            if base_col < 0 or base_col >= total_cols:
                raise IndexError(f"列数超出范围（有效1-{total_cols}）")
            
            start_col = max(0, base_col - nc)
            end_col = min(total_cols - 1, base_col + nc)
            avg_col = np.mean(self.image_matrix[:, start_col:end_col+1], axis=1)
            
            if len(avg_col) < 4:
                raise ValueError("曲线长度不足，无法分析")
            
            corr_coeffs = []
            splits = []
            for split in range(1, len(avg_col)-1):
                p1 = avg_col[:split][-min(len(avg_col[:split]), len(avg_col[split+1:])):]
                p2 = avg_col[split+1:][:min(len(avg_col[:split]), len(avg_col[split+1:])):]
                if len(p1) >= 2 and len(p2) >= 2:
                    corr, _ = pearsonr(p1[::-1], p2)
                    corr_coeffs.append(corr)
                    splits.append(split)
            
            if not corr_coeffs:
                raise RuntimeError("无法计算有效相关系数")
            
            exclude = int(len(corr_coeffs) * 0.1)
            if len(corr_coeffs) > 2 * exclude:
                mid_corr = corr_coeffs[exclude:-exclude]
                mid_splits = splits[exclude:-exclude]
                center_row = mid_splits[np.argmax(mid_corr)]
            else:
                center_row = splits[np.argmax(corr_coeffs)]
            
            self.vertical_center = center_row
            self.vertical_result_label.setText(f"{center_row + 1}")
            self.log_print(f"竖直中心计算完成：列{n}，展宽{nc}，中心行{center_row + 1} (1-based)")
            
            self.clear_center_btn.setEnabled(True)
            self.redraw_images_with_center_lines()
        except (ValueError, IndexError) as e:
            QMessageBox.warning(self, "输入错误", str(e))
            self.log_print(f"竖直中心计算错误: {str(e)}")
        except Exception as e:
            QMessageBox.warning(self, "错误", str(e))
            self.log_print(f"竖直中心计算错误: {str(e)}")
    
    def redraw_images_with_center_lines(self):
        """重新绘制图像以显示中心线"""
        if self.left_image_display.pixmap():
            self.display_left_pixmap_with_center(self.left_image_display.pixmap())
        
        if self.right_image_display.pixmap():
            self.display_background_pixmap_with_center(self.right_image_display.pixmap())
    
    def display_left_pixmap_with_center(self, pixmap):
        """显示原始图像并叠加中心线"""
        if pixmap.isNull():
            return
            
        new_pixmap = QPixmap(pixmap)
        painter = QPainter(new_pixmap)
        
        self.draw_center_lines(painter, new_pixmap.width(), new_pixmap.height())
        
        painter.end()
        
        avail_size = self.left_scroll_area.size()
        scaled = new_pixmap.scaled(avail_size.width()-20, avail_size.height()-20, 
                                 Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.left_image_display.setPixmap(scaled)
    
    def display_background_pixmap_with_center(self, pixmap):
        """显示本底图像并叠加中心线"""
        if pixmap.isNull() or self.image_matrix is None:
            return
            
        new_pixmap = QPixmap(pixmap)
        painter = QPainter(new_pixmap)
        
        self.draw_center_lines(painter, new_pixmap.width(), new_pixmap.height())
        
        painter.end()
        
        avail_size = self.right_scroll_area.size()
        scaled = new_pixmap.scaled(avail_size.width()-20, avail_size.height()-20, 
                                 Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.right_image_display.setPixmap(scaled)
    
    def polar_coordinate_expand_with_diff(self):
        """极坐标展开与差值计算"""
        try:
            if self.image_matrix is None:
                raise ValueError("请先加载TIFF图像")
            if self.horizontal_center is None or self.vertical_center is None:
                raise ValueError("请先计算水平和竖直中心")
            
            start_angle = float(self.angle_start_edit.text().strip())
            end_angle = float(self.angle_end_edit.text().strip())
            if not (0 <= start_angle < 360 and 0 < end_angle <= 360 and start_angle < end_angle):
                raise ValueError("角度范围必须为0-360°且起始<终止")
            
            center_x = self.horizontal_center
            center_y = self.vertical_center
            img_height, img_width = self.image_matrix.shape
            angle_range = np.arange(start_angle, end_angle + 1, 1)
            angle_count = len(angle_range)
            diagonal = math.sqrt(img_width**2 + img_height**2)
            fixed_len = int(round(diagonal))
            self.log_print(f"极坐标展开：中心({center_x+1},{center_y+1}) (1-based列,1-based行)，角度{start_angle}-{end_angle}°（{angle_count}个角度），采样长度{fixed_len}")
            
            if fixed_len <= 0:
                raise RuntimeError("采样长度无效")
            
            angles_rad = np.radians(angle_range)
            cos_thetas = np.cos(angles_rad)
            sin_thetas = -np.sin(angles_rad)
            
            radii = np.arange(fixed_len + 1)
            
            x_coords = center_x + np.outer(cos_thetas, radii)
            y_coords = center_y + np.outer(sin_thetas, radii)
            
            x_indices = np.round(x_coords).astype(int)
            y_indices = np.round(y_coords).astype(int)
            
            valid_mask = (x_indices >= 0) & (x_indices < img_width) & (y_indices >= 0) & (y_indices < img_height)
            
            self.polar_expand_matrix = np.zeros((angle_count, fixed_len + 1), dtype=np.float32)
            self.polar_expand_matrix[valid_mask] = self.image_matrix[
                y_indices[valid_mask], x_indices[valid_mask]
            ]
            
            if self.background_matrix is not None:
                bg_expand_matrix = np.zeros((angle_count, fixed_len + 1), dtype=np.float32)
                bg_expand_matrix[valid_mask] = self.background_matrix[
                    y_indices[valid_mask], x_indices[valid_mask]
                ] * self.background_coeff
                self.log_print(f"本底图像展开完成，应用系数: {self.background_coeff}")
            else:
                bg_expand_matrix = np.zeros_like(self.polar_expand_matrix)
                self.log_print("未加载本底，本底按全0处理")
            
            self.diff_matrix = self.polar_expand_matrix - bg_expand_matrix
            self.diff_matrix[self.diff_matrix < 0] = np.nan
            
            valid_values = self.diff_matrix[~np.isnan(self.diff_matrix)]
            if len(valid_values) > 0:
                self.log_print(f"差值矩阵：尺寸{self.diff_matrix.shape}，有效取值范围{np.min(valid_values):.2f}~{np.max(valid_values):.2f}，"
                              f"负值已设为NaN（共{np.sum(np.isnan(self.diff_matrix))}个NaN值）")
            else:
                self.log_print(f"差值矩阵：尺寸{self.diff_matrix.shape}，所有值均为NaN")
            
            self.display_diff_result(angle_range, fixed_len)
            self.plot_col_average_curve()
        except ValueError as e:
            QMessageBox.warning(self, "参数错误", str(e))
            self.log_print(f"极坐标展开错误: {str(e)}")
        except Exception as e:
            QMessageBox.warning(self, "错误", str(e))
            self.log_print(f"极坐标展开错误: {str(e)}")
    
    def display_diff_result(self, angle_range, max_radius):
        """显示差值矩阵图像（使用当前选择的颜色图）"""
        if self.diff_matrix is None:
            return
        try:
            mat = self.diff_matrix.copy()
            mask = np.isnan(mat)
            mat[mask] = 0
            
            min_val, max_val = np.min(mat), np.max(mat)
            if max_val > min_val:
                normalized = (mat - min_val) / (max_val - min_val) * 255
            else:
                normalized = np.zeros_like(mat)
                
            normalized = normalized.astype(np.uint8)
            
            # 使用当前选择的颜色图，NaN区域保持红色
            colored = (cm.get_cmap(self.cmap_combo.currentText())(normalized / 255.0) * 255).astype(np.uint8)
            colored[mask] = [255, 0, 0, 255]
            
            q_image = QImage(colored.data, colored.shape[1], colored.shape[0], 
                           4 * colored.shape[1], QImage.Format_RGBA8888)
            
            pixmap = QPixmap.fromImage(q_image)
            avail_size = self.diff_scroll_area.size()
            scaled = pixmap.scaled(avail_size.width()-20, avail_size.height()-20, 
                                 Qt.IgnoreAspectRatio, Qt.SmoothTransformation)
            self.diff_image_label.setPixmap(scaled)
            
            # 更新标题，显示当前颜色图
            self.diff_title_label.setText(
                f"极坐标展开差值结果（数据 - {self.background_coeff:.4f}×本底，负值已设为NaN），角度{angle_range[0]}-{angle_range[-1]}°，颜色图：{self.cmap_combo.currentText()}"
            )
        except Exception as e:
            self.log_print(f"显示差值结果失败: {str(e)}")
    
    def on_q_param_changed(self):
        """能量E或距离L改变时更新曲线"""
        self.q_energy = self.energy_spin.value()
        self.q_distance = self.distance_spin.value()
        
        if self.diff_matrix is not None:
            self.plot_col_average_curve()
            self.log_print(f"q变换参数更新：能量E={self.q_energy:.1f} eV，探测器距离L={self.q_distance} mm")
    
    def on_crop_param_changed(self):
        """裁剪参数变更时更新曲线"""
        self.crop_front_nan = self.crop_front_spin.value()
        self.crop_back_nan = self.crop_back_spin.value()
        
        if self.diff_matrix is not None:
            self.plot_col_average_curve()
            self.log_print(f"曲线裁剪更新：前{self.crop_front_nan}个点设为NaN，后{self.crop_back_nan}个点设为NaN")
    
    def plot_col_average_curve(self):
        """绘制列平均值曲线"""
        if self.diff_matrix is None:
            self.clear_curve()
            self.log_print("无差值矩阵，无法绘制曲线")
            return
        
        try:
            col_avg = np.nanmean(self.diff_matrix, axis=0)
            total_points = len(col_avg)
            radii_x = np.arange(total_points)
            
            cropped_avg = col_avg.copy().astype(np.float64)
            
            if self.crop_front_nan > 0 and self.crop_front_nan < total_points:
                cropped_avg[:self.crop_front_nan] = np.nan
            elif self.crop_front_nan >= total_points:
                cropped_avg[:] = np.nan
                self.log_print(f"警告：前{self.crop_front_nan}个点≥总点数{total_points}，曲线全部设为NaN")
            
            if self.crop_back_nan > 0 and self.crop_back_nan < total_points:
                back_start_idx = total_points - self.crop_back_nan
                cropped_avg[back_start_idx:] = np.nan
            elif self.crop_back_nan >= total_points:
                cropped_avg[:] = np.nan
                self.log_print(f"警告：后{self.crop_back_nan}个点≥总点数{total_points}，曲线全部设为NaN")
            
            if self.crop_front_nan + self.crop_back_nan > total_points:
                cropped_avg[:] = np.nan
                self.log_print(f"警告：前{self.crop_front_nan} + 后{self.crop_back_nan} > 总点数{total_points}，无有效数据点")
            
            self.crop_front_spin.setMaximum(total_points)
            self.crop_back_spin.setMaximum(total_points)
            
            try:
                equivalent_pixel_μm = float(self.equivalent_pixel_edit.text().strip())
                equivalent_pixel_mm = equivalent_pixel_μm / 1000.0
                if equivalent_pixel_mm <= 0:
                    raise ValueError("等效像素尺寸必须为正数")
            except ValueError as e:
                self.log_print(f"等效像素尺寸无效：{str(e)}，使用默认值15μm")
                equivalent_pixel_mm = 0.015
                equivalent_pixel_μm = 15.0
            
            L = self.q_distance
            E = self.q_energy
            arctan_term = np.arctan((radii_x * equivalent_pixel_mm) / L)
            sin_term = np.sin(0.5 * arctan_term)
            denominator = 1.239 / (E / 1000)
            q = (4 * np.pi * sin_term) / denominator
            
            self.figure.clear()
            ax = self.figure.add_subplot(111)
            ax.plot(q, cropped_avg, color='#2E86AB', linewidth=2, marker='o', markersize=1)
            
            ax.set_xlabel('Scattering vector q (nm^-1)', fontsize=self.current_font_size)
            ax.set_ylabel('强度', fontsize=self.current_font_size)
            ax.set_yscale('log')
            ax.set_title(
                f'散射曲线（E={E:.1f}eV, L={L}mm, 等效像素={equivalent_pixel_μm:.2f}μm | '
                f'前{self.crop_front_nan}后{self.crop_back_nan}个点设为NaN）',
                fontsize=self.current_font_size, fontweight='bold'
            )
            
            ax.grid(True, alpha=0.3, linestyle='--')
            self.figure.tight_layout()
            self.canvas.draw()
            
            self.curve_data = np.column_stack((q, cropped_avg))
            
            valid_points = np.count_nonzero(~np.isnan(cropped_avg))
            self.log_print(
                f"曲线绘制完成：总点数{total_points}，有效点数{valid_points}，"
                f"q范围{np.nanmin(q):.4f}~{np.nanmax(q):.4f} nm^-1，"
                f"使用等效像素尺寸：{equivalent_pixel_μm:.2f}μm"
            )
        except Exception as e:
            self.clear_curve()
            self.log_print(f"绘制曲线失败: {str(e)}")
    
    def clear_curve(self):
        """清空曲线"""
        try:
            equivalent_pixel_μm = float(self.equivalent_pixel_edit.text().strip())
        except:
            equivalent_pixel_μm = 15.0
            
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.set_xlabel('Scattering vector q (nm^-1)', fontsize=self.current_font_size)
        ax.set_ylabel('强度', fontsize=self.current_font_size)
        ax.set_yscale('log')
        ax.set_title(
            f'散射曲线（待计算 | E={self.q_energy:.1f}eV, L={self.q_distance}mm, 等效像素={equivalent_pixel_μm:.2f}μm | '
            f'前{self.crop_front_nan}后{self.crop_back_nan}个点设为NaN）',
            fontsize=self.current_font_size, fontweight='bold'
        )
        ax.grid(True, alpha=0.3, linestyle='--')
        self.canvas.draw()
        
        self.curve_data = None
    
    def save_data(self):
        """保存分析数据"""
        if self.current_file_path is None:
            QMessageBox.information(self, "提示", "请先加载图像并进行分析")
            return
            
        if self.diff_matrix is None and self.polar_expand_matrix is None and not hasattr(self, 'curve_data'):
            QMessageBox.information(self, "提示", "没有可保存的分析数据，请先进行极坐标展开")
            return
            
        base_name = os.path.splitext(os.path.basename(self.current_file_path))[0]
        default_dir = os.path.dirname(self.current_file_path)
        default_filename = f"{base_name}_analysis"
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "保存分析数据", os.path.join(default_dir, default_filename), 
            "文本文件 (*.txt);;CSV文件 (*.csv);;所有文件 (*)"
        )
        
        if not file_path:
            return
            
        try:
            save_dir = os.path.dirname(file_path)
            if not os.path.exists(save_dir):
                os.makedirs(save_dir)
            
            if hasattr(self, 'curve_data') and self.curve_data is not None:
                curve_file = f"{os.path.splitext(file_path)[0]}_curve.txt"
                with open(curve_file, 'w', encoding='utf-8') as f:
                    f.write(f"# Analysis Data - {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    f.write(f"# Original file: {self.current_file_path}\n")
                    
                    if self.horizontal_center is not None and self.vertical_center is not None:
                        f.write(f"# Beam center coordinates:\n")
                        f.write(f"# Horizontal center (column/x): {self.horizontal_center + 1} pixels (1-based)\n")
                        f.write(f"# Vertical center (row/y): {self.vertical_center + 1} pixels (1-based)\n")
                    else:
                        f.write(f"# Beam center coordinates not calculated\n")
                    
                    f.write(f"# Energy: {self.q_energy} eV\n")
                    f.write(f"# Detector distance: {self.q_distance} mm\n")
                    f.write(f"# Equivalent pixel size: {self.equivalent_pixel_edit.text()} μm\n")
                    f.write(f"# Background coefficient: {self.background_coeff}\n")
                    f.write(f"# Crop parameters: front={self.crop_front_nan}, back={self.crop_back_nan}\n")
                    
                    f.write("\n# Data:\n")
                    f.write("q(nm^-1)\tIntensity\n")
                    
                    for q_val, intensity in self.curve_data:
                        if not np.isnan(intensity):
                            f.write(f"{q_val:.6f}\t{intensity:.6f}\n")
            
                self.log_print(f"散射曲线数据已保存至: {curve_file}")
            
            if self.diff_matrix is not None:
                diff_file = f"{os.path.splitext(file_path)[0]}_diff_matrix.csv"
                np.savetxt(diff_file, self.diff_matrix, delimiter=',', fmt='%.6f')
                self.log_print(f"差值矩阵已保存至: {diff_file}")
            
            if self.polar_expand_matrix is not None:
                polar_file = f"{os.path.splitext(file_path)[0]}_polar_matrix.csv"
                np.savetxt(polar_file, self.polar_expand_matrix, delimiter=',', fmt='%.6f')
                self.log_print(f"极坐标展开矩阵已保存至: {polar_file}")
            
            param_file = f"{os.path.splitext(file_path)[0]}_parameters.txt"
            with open(param_file, 'w', encoding='utf-8') as f:
                f.write(f"Analysis parameters - {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Original file: {self.current_file_path}\n\n")
                f.write(f"Beam center coordinates:\n")
                f.write(f"Horizontal center (column/x): {self.horizontal_center + 1 if self.horizontal_center is not None else 'Not set'} pixels (1-based)\n")
                f.write(f"Vertical center (row/y): {self.vertical_center + 1 if self.vertical_center is not None else 'Not set'} pixels (1-based)\n\n")
                f.write(f"Energy: {self.q_energy} eV\n")
                f.write(f"Detector distance: {self.q_distance} mm\n")
                f.write(f"Equivalent pixel size: {self.equivalent_pixel_edit.text()} μm\n")
                f.write(f"Background coefficient: {self.background_coeff}\n")
                f.write(f"Angle range: {self.angle_start_edit.text()}° - {self.angle_end_edit.text()}°\n")
                f.write(f"Curve cropping: front {self.crop_front_nan} points, back {self.crop_back_nan} points\n")
            
            self.log_print(f"分析参数已保存至: {param_file}")
            QMessageBox.information(self, "成功", "数据保存成功")
            
        except Exception as e:
            QMessageBox.warning(self, "保存失败", f"无法保存数据: {str(e)}")
            self.log_print(f"数据保存失败: {str(e)}")
    
    def log_print(self, text):
        """日志输出"""
        print(text)
        self.log_text_edit.append(text)
        self.log_text_edit.verticalScrollBar().setValue(self.log_text_edit.verticalScrollBar().maximum())
    
    def resizeEvent(self, event):
        """窗口 resize 时更新图像显示"""
        super().resizeEvent(event)
        if hasattr(self.left_image_display, 'pixmap') and self.left_image_display.pixmap():
            self.display_left_pixmap_with_center(self.left_image_display.pixmap())
        if hasattr(self.right_image_display, 'pixmap') and self.right_image_display.pixmap():
            self.display_background_pixmap_with_center(self.right_image_display.pixmap())
        if self.diff_matrix is not None:
            try:
                start_angle = float(self.angle_start_edit.text().strip())
                end_angle = float(self.angle_end_edit.text().strip())
                angle_range = np.arange(start_angle, end_angle + 1, 1)
                img_h, img_w = self.image_matrix.shape
                max_radius = int(round(math.sqrt(img_w**2 + img_h**2)))
                self.display_diff_result(angle_range, max_radius)
            except:
                pass
        self.canvas.resize(self.canvas.size())
    
    def apply_styles(self):
        """界面样式（使用当前字体大小）"""
        style = f"""
            QMainWindow {{ background-color: #f5f5f5; }}
            QFrame {{ background-color: white; border-radius: 8px; border: 1px solid #d0d0d0; }}
            QPushButton {{ 
                background-color: #4a86e8; color: white; border: none; 
                padding: 8px 16px; border-radius: 4px; font-weight: bold;
                font-size: {self.current_font_size}pt;
            }}
            QPushButton:hover {{ background-color: #3a76d8; }}
            QPushButton:pressed {{ background-color: #2a66c8; }}
            QListWidget {{ 
                border: 1px solid #d0d0d0; border-radius: 4px; 
                background-color: white; alternate-background-color: #f8f8f8;
            }}
            QTextEdit, QLineEdit, QSpinBox, QDoubleSpinBox {{ 
                border: 1px solid #d0d0d0; border-radius: 4px; 
                padding: 5px; background-color: white;
                font-size: {self.current_font_size}pt;
            }}
            QLineEdit, QSpinBox, QDoubleSpinBox {{ padding: 8px; }}
            QLineEdit:read-only {{ background-color: #f5f5f5; color: #666666; }}
            QScrollArea {{ 
                border: 1px solid #d0d0d0; border-radius: 4px; 
                background-color: #f8f8f8;
            }}
            QComboBox {{ 
                border: 1px solid #d0d0d0; border-radius: 4px; 
                padding: 5px; background-color: white;
                font-size: {self.current_font_size}pt;
            }}
            QLabel {{ padding: 2px; font-size: {self.current_font_size}pt; }}
        """
        self.setStyleSheet(style)
if __name__ == "__main__":
# 解决Matplotlib中文显示
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    app = QApplication(sys.argv)
    
    # 启动logo显示
    splash = None  # 启动画面实例
    try:
        splash_pixmap = QPixmap("ScatterX_logo.png")
        if not splash_pixmap.isNull():
            max_width = 600
            max_height = 400
            scaled_pixmap = splash_pixmap.scaled(
                max_width,         
                max_height,        
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
            splash = QSplashScreen(scaled_pixmap, Qt.WindowStaysOnTopHint)
            splash.show()
            time.sleep(2)
            app.processEvents()
        else:
            print("启动logo加载失败：图片为空（可能路径错误或图片损坏）")
    except Exception as e:
        print(f"启动logo显示异常：{str(e)}")
    
    # 初始化主窗口
    window = MainWindow()
    window.showMaximized()
    
    # 关闭启动画面
    if splash is not None:
        splash.finish(window)
    
    sys.exit(app.exec_())        