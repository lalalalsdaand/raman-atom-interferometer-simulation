# Raman Atom Interferometer Simulation

<p align="center">
  <img src="assets/logo.png" alt="A-KNOWS logo" width="420">
</p>

## Ownership

- Author: 赵博士
- Affiliation: A-KNOWS 实验室，清华大学
- Contact: yingpeng-zhao@mail.tsinghua.edu.cn

这是一个用于 **Raman 原子干涉仿真** 的 MATLAB 程序，包含单维原子干涉和多维原子干涉两类模型。程序面向冷原子干涉实验中的 Raman 光脉冲过程，可用于模拟原子初始化、自由飞行、微波态制备、Raman 跃迁、AC Stark shift、Zeeman 修正、速度选择、Rabi 振荡、Raman 谱线和干涉条纹等过程。

当前整理版的主代码为：

- `double_raman_test0507.m`：主仿真脚本，串联完整实验流程。
- `MW_transition.m`：微波跃迁演化函数，用于原子内态制备和谱线扫描。

## Features

- 单维 Raman 原子干涉仿真
- 多维/多通道 Raman 原子干涉仿真
- Rb87 D2 线相关偶极跃迁矩阵构造
- AC Stark shift 计算与光强比扫描
- Zeeman 哈密顿量构造
- 原子自由下落、速度选择和态选择
- Raman 谱、微波谱、Rabi 振荡和干涉条纹绘图

## Requirements

- MATLAB
- 代码主要依赖 MATLAB 基础矩阵运算、`ode45` 和绘图函数
- 若使用 `parfor`，需要 Parallel Computing Toolbox；也可以将 `parfor` 改成普通 `for` 运行

## Quick Start

在 MATLAB 中进入项目目录：

```matlab
cd("path/to/Raman-Atom-Interferometer-Simulation")
addpath(genpath(pwd))
```

运行主脚本：

```matlab
double_raman_test0507
```

如果只是快速测试流程，建议先把 `double_raman_test0507.m` 中的原子数调小：

```matlab
N = 10;   % 原代码中为 1E3
```

## Main Workflow

`double_raman_test0507.m` 中的典型仿真流程包括：

1. 读取默认参数并生成初始原子云。
2. 进行自由下落演化。
3. 通过 `MW_transition.m` 模拟微波跃迁和态制备。
4. 构造 Raman 光束参数、偏振、AC Stark shift 和 Zeeman 哈密顿量。
5. 分别调用单维和多维 Raman 演化函数。
6. 输出 Raman 谱、Rabi 振荡、速度选择结果和干涉条纹。

## Code Structure

### Main scripts

- `double_raman_test0507.m`：主仿真脚本。
- `MW_transition.m`：微波跃迁演化函数。

### Atom state and time evolution

- `default_aii_config.m`：默认实验和原子参数。
- `generate_atom.m`：生成初始原子位置、速度和密度矩阵。
- `state_eval.m`、`state_eval_2.m`、`state_eval_3.m`：封装 `ode45` 的状态演化函数。
- `free_fall.m`：自由下落演化。
- `blow_state.m`：按内态筛选原子。
- `init_blocks_state.m`：构造多动量块初态。
- `statistic_plot.m`：原子云统计绘图。
- `noisegen.m`：给仿真信号加入噪声。

### Raman and microwave transitions

- `microwave_delta_rabi.m`：微波失谐和 Rabi 耦合矩阵。
- `Raman_transition_new.m`：单维 Raman 演化函数。
- `Raman_transition_5blocks_40.m`：多维/多通道 Raman 演化函数。
- `build_ac_eff_coprop.m`：构造 Raman 有效 AC 和耦合项。
- `init_prePair_terms.m`：预计算 Raman 光束组合项。
- `scan_Iratio_cancel_ac.m`：扫描光强比并寻找差分 AC Stark shift 抵消点。
- `make_polarization.m`：根据光束方向和偏振角生成偏振矢量。

### Atomic physics helpers

- `dipole_T_FF_auto.m`：构造偶极耦合张量。
- `rb87_D2_dq_steck.m`：Rb87 D2 线跃迁矩阵元。
- `get_zeeman_hamiltonian.m`：构造 Zeeman 哈密顿量。
- `get_F_operators.m`、`get_spin_operators.m`：角动量与电子自旋算符。
- `Wigner3j.m`、`Wigner6j.m`、`Wigner9j.m`：Wigner 符号计算。
- `if3jc.m`、`ArranA.m`、`sper.m`：Wigner 计算辅助函数。

### Dimension conversion

- `rho8_to_rho56.m`、`rho56_to_rho8.m`：8 维内态与 56 维多块密度矩阵转换。
- `convert_y1606_to_y70.m`：40x40 多块演化历史压缩。
- `convert_Y3142_to_Y70.m`：56x56 多块演化历史压缩。

## Notes

- 原始代码中部分中文注释存在编码显示异常，但主要 MATLAB 调用关系已经整理保留。
- `double_raman_test0507.m` 是实验探索型脚本，后半部分包含多段可选绘图和扫描代码。分段运行时，请保证前序变量已经生成。
- 整理版删除了 `.asv` 自动备份、旧版函数、早期测试脚本、独立检查脚本和临时 `.mat` 数据文件。
- 整理版修正了 `state_eval_3.m` 首行函数名与文件名不一致的问题，计算逻辑未改。

## Repository Status

This repository is an organized research-code snapshot for Raman atom interferometer simulation. It is intended for experiment-oriented numerical exploration rather than a packaged MATLAB toolbox.
