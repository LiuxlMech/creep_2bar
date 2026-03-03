# creep_2bar

两杆模型（A 常温弹性杆 + B 高温蠕变杆）热应力-蠕变耦合计算示例，支持截面积比 `k=A_A/A_B`，并支持从外部文件加载蠕变模型（方程 + 参数）。

## 1. 问题与假设

- 两根杆 A/B 材料相同，长度同为 `L=100 mm`，但截面积可不同。定义 `k=A_A/A_B`（例如 `k=2`）。
- A 杆保持室温（默认 20°C），仅弹性；B 杆升温到高温（默认 800°C），发生弹性 + 热膨胀 + 蠕变。
- 并联连接，左端固定右端相连：
1. 两杆总应变相同；
2. 接头处无外载，故内力平衡 `σ_A A_A + σ_B A_B = 0`。
- 小变形、单轴应力、B 杆温度恒定（持温阶段）。

## 2. 控制方程

### 2.1 应变协调与力平衡

- A 杆：
  \[
  \varepsilon_A = \frac{\sigma_A}{E_A}
  \]
- B 杆：
  \[
  \varepsilon_B = \frac{\sigma_B}{E_B} + \varepsilon_{th,B} + \varepsilon_{creep,B}
  \]
- 条件：\(\varepsilon_A=\varepsilon_B\)，且 \(\sigma_A=-\sigma_B/k\)。

代入后可得：
\[
\sigma_B = -\frac{\varepsilon_{th,B}+\varepsilon_{creep,B}}{1/(kE_A)+1/E_B}
\]
记
\[
C = \left(\frac{1}{kE_A}+\frac{1}{E_B}\right)^{-1}, \quad \sigma_B = -C\,(\varepsilon_{th,B}+\varepsilon_{creep,B}), \quad \sigma_A=-\sigma_B/k
\]

### 2.2 热应变（按“平均热膨胀系数”解释）

本题给定的是从 20°C 到各温度点的平均热膨胀系数，不是瞬时系数。  
程序做法：
- 对每个表格温度点 `T_i` 先构造累计热应变  
  `eps(20->T_i) = alpha_mean(20->T_i) * (T_i - 20)`
- 再对累计热应变进行分段线性插值，最终得到
  `eps_th_B = eps(20->T_B) - eps(20->T_A)`

### 2.3 蠕变模型（外部文件）

主程序不再内置 `--K --n` 这类蠕变参数，而是通过 `--creep-model` 指定一个 Python 文件。  
该文件中定义方程与参数，主程序按统一接口调用。

## 3. 蠕变模型文件接口

`--creep-model` 指向一个 `.py` 文件，文件至少包含：

- `MODEL_NAME`：模型名称字符串
- `PARAMS`：参数字典
- `creep_rate(...)`：返回 `d(eps_creep)/dt`，单位 `1/h`

函数签名约定：

```python
def creep_rate(*, sigma_mpa: float, eps_creep: float, t_h: float, temp_c: float, params: dict) -> float:
    ...
```

仓库内置模型示例：

- `creep_models/norton.py`
- `creep_models/time_hardened.py`
- `creep_models/garofalo.py`
- `creep_models/custom_template.py`（自定义模板）

默认 `creep_models/norton.py` 已改为温度相关 `K(T), n(T)`，采用 Inconel 625 的公开蠕变数据（650/900/980°C）进行拟合插值：

- 参考：Ozgowicz et al., *Materials* 2023, 16, 7724（Table 3 与文中 `n=10.8` 说明）  
  https://www.mdpi.com/1996-1944/16/24/7724
- `K(T)` 在对数空间线性插值，`n(T)` 线性插值
- 文献数据应力范围约 20-400 MPa，模型中对 `|sigma|` 做同范围钳制，减少超范围外推失真

## 4. 运行方法

### 4.1 默认示例（Norton, 500h）

```bash
python3 two_bar_creep.py --hours 500 --k 2 --creep-model creep_models/norton.py
```

### 4.2 使用其它模型（例如 Garofalo）

```bash
python3 two_bar_creep.py --hours 500 --k 2 --creep-model creep_models/garofalo.py
```

### 4.3 自定义参数

```bash
python3 two_bar_creep.py --T-a 20 --T-b 800 --hours 500 --k 2 --dt 0.01 --creep-model creep_models/time_hardened.py
```

参数说明：

- `--T-a`：A 杆温度（°C）
- `--T-b`：B 杆温度（°C）
- `--hours`：总时长（h）
- `--k`：截面积比 `A_A/A_B`
- `--dt`：积分步长（h）
- `--creep-model`：蠕变模型文件路径

## 5. 输出解释

程序会输出：

- `E_A, E_B`：按温度插值后的杨氏模量（MPa）
- `k`：截面积比
- `creep_model`、`creep_model_file`、`creep_model_params`
- `eps_th_B`：B 杆自由热应变
- `sigma_B(0)`：初始热约束应力（通常为压应力）
- `sigma_B(t_end)`：蠕变松弛后的 B 杆应力
- `sigma_A(t_end)` 可由 `-sigma_B/k` 得到
- `eps_creep_B(t_end)`：目标时刻蠕变应变

## 6. 工程使用建议

1. 用试验数据标定各模型参数后再用于预测。
2. 若温度随时间变化，可改造为每步更新 `E(T)` 与热应变，并将 `temp_c` 传入模型方程。
3. 若进入显著塑性区，建议使用弹塑-蠕变统一模型。
4. 本脚本用于快速一维估算，适合方案比选与参数敏感性分析。
