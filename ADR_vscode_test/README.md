# ADR_Fortran_my - 一维对流扩散方程的增幅色散关系分析

## 项目简介

本项目是一个用于分析数值格式在求解一维对流扩散方程（ADR）时的增幅色散关系（Amplification Dispersion Relation）的Fortran程序。通过计算不同波数下的增幅因子和色散关系，可以评估数值格式的精度、稳定性和耗散/色散特性。

## 代码功能概述

该程序主要用于分析WCNS（Weighted Compact Nonlinear Scheme）系列格式在求解对流项时的性能表现。程序支持多种WCNS格式，包括：

- WCNS-MR (Mapped WENO with Renormalization)
- WCNS-JS (Jiang-Shu WENO)
- WCNS-Z (WENO-Z)

通过傅里叶分析方法，程序计算不同精度阶数（3阶、5阶、7阶、9阶、11阶）的数值格式在不同网格分辨率下的增幅因子和相位误差，从而评估格式的性能。

## 代码文件夹结构

```Fortran
.
├── Code_BC/                    # 边界条件处理模块
│   ├── ImposeBC.f90           # 传统SD边界处理方法
│   └── ImposeBC_CBM.f90       # 特征边界方法（当前未使用）
├── Code_Mesh/                 # 网格处理模块
│   └── 2_CalGridMetrics.f90   # 网格度量计算（当前未使用）
├── Code_Scheme/               # 数值格式核心模块
│   ├── Coe_Interface.f90      # 接口系数计算
│   ├── Limiter_PP_ConsVar.f90 # 保守变量限制器
│   ├── Limiter_PP_Flux.f90    # 通量限制器
│   ├── RHS.f90                # 右端项计算
│   ├── Time_RK3.f90           # 三阶TVD Runge-Kutta时间推进
│   ├── ghostNode_Diff.f90     # 差分格式中的虚拟节点处理
│   ├── ghostNode_Int.f90      # 插值格式中的虚拟节点处理
│   └── parameters.mod         # 参数模块
├── Main.f90                   # 主程序
└── README.md                  # 本说明文件

```

## 核心模块功能详解

### 1. 主程序 (Main.f90)

主程序实现了对一维对流扩散方程的傅里叶分析：

- 设置计算参数和网格
- 初始化不同波数的初始条件
- 使用TVD-RK3时间推进求解方程
- 通过傅里叶变换计算增幅因子和相位误差
- 输出结果到文件用于后续分析

### 2. 数值格式模块 (Code_Scheme/)

#### 插值方法 (ghostNode_Int.f90)

实现了WCNS-MR系列格式的不同阶数实现：

- 3阶、5阶、7阶、9阶、11阶WCNS-MR格式
- 使用重构-映射方法提高格式的精度和稳定性

#### 差分方法 (ghostNode_Diff.f90)

实现了不同精度的差分格式：

- E4, E6, E8, E10等偶数阶中心差分格式

#### 时间推进 (Time_RK3.f90)

实现了三阶TVD Runge-Kutta时间推进方法：

- 三步Runge-Kutta格式，保证数值稳定性

#### 右端项计算 (RHS.f90)

计算方程的右端项（空间导数）：

- 调用插值和差分模块处理虚拟节点
- 实现完整的空间离散化过程

### 3. 边界条件模块 (Code_BC/)

虽然包含边界条件处理文件，但在当前的傅里叶分析程序中未使用。这些模块主要用于实际CFD计算中的边界处理：

- [ImposeBC.f90](file:///d:/whf00/Documents/GitHub/CFD_inWin/ADR_Fortran_my/Code_BC/ImposeBC.f90)：传统SD边界处理方法
- [ImposeBC_CBM.f90](file:///d:/whf00/Documents/GitHub/CFD_inWin/ADR_Fortran_my/Code_BC/ImposeBC_CBM.f90)：特征边界方法（当前未使用）

## 编译命令

gfortran -O2 Main.f90 Code_Scheme/Time_RK3.f90 Code_Scheme/RHS.f90 Code_Scheme/ghostNode_Int.f90 Code_Scheme/ghostNode_Diff.f90 Code_Scheme/Coe_Interface.f90 Code_Scheme/Limiter_PP_ConsVar.f90 Code_Scheme/Limiter_PP_Flux.f90 -o Main.exe

gfortran -O2 Main.f90 Code_Scheme/Coe_Interface.f90 Code_Scheme/Time_RK3.f90 Code_Scheme/RHS.f90 Code_Scheme/ghostNode_Int.f90 Code_Scheme/ghostNode_Diff.f90 Code_Scheme/Limiter_PP_ConsVar.f90 Code_Scheme/Limiter_PP_Flux.f90 -o Main.exe

gfortran -O2 Main.f90
Code_Scheme/Coe_Interface.f90
Code_Scheme/Time_RK3.f90
Code_Scheme/RHS.f90
Code_Scheme/ghostNode_Int.f90
Code_Scheme/ghostNode_Diff.f90
Code_Scheme/Limiter_PP_ConsVar.f90
Code_Scheme/Limiter_PP_Flux.f90 -o Main.exe

./Main.exe
