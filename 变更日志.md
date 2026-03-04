# creep_2bar 变更日志

更新日期：2026-03-03

## 一、已完成改动

1. 截面积比参数化  
- 模型由等截面扩展为可变截面积比 `k = A_A / A_B`。  
- 命令行支持 `--k`。  
- 力平衡改为 `sigma_A * A_A + sigma_B * A_B = 0`，并据此更新应力关系。

2. 蠕变模型插件化（主程序与模型解耦）  
- 主程序不再通过 `--K --n` 直接输入 Norton 参数。  
- 改为 `--creep-model <path>` 指定模型文件。  
- 模型文件统一接口：  
  - `MODEL_NAME`  
  - `PARAMS`  
  - `creep_rate(*, sigma_mpa, eps_creep, t_h, temp_c, params) -> float`

3. 新增模型文件  
- `creep_models/norton.py`  
- `creep_models/time_hardened.py`  
- `creep_models/garofalo.py`  
- `creep_models/custom_template.py`

4. 文档合并  
- 原说明书内容合并到 `README.md`。  
- `说明书.md` 保留为跳转说明（提示查看 README）。

5. A/B 温度变量化  
- 已支持命令行变量：`--T-a`、`--T-b`。

6. Norton 模型温度相关化（Inconel 625）  
- `creep_models/norton.py` 改为 `K(T), n(T)` 随温度变化。  
- 采用文献数据点（650/900/980°C）拟合插值：  
  - `K(T)`：对数空间线性插值  
  - `n(T)`：线性插值  
- 新增应力有效区间钳制（默认 `20~400 MPa`），减少超范围外推导致的数值失真。

## 二、数据来源与说明

1. 文献来源  
- Ozgowicz et al., Materials 2023, 16, 7724  
- 链接：https://www.mdpi.com/1996-1944/16/24/7724  
- 公开镜像：https://pmc.ncbi.nlm.nih.gov/articles/PMC10673153/

2. 参数说明（重要）  
- 650°C 点：`n=10.8` 取自文中说明，`K` 由表格应力-稳态蠕变速率点反算。  
- 900°C/980°C 点：`n` 与 `K` 由各两组 `(sigma, eps_dot)` 数据反算。  
- 以上反算属于工程化重建参数，不是论文直接给出的完整 `K,n(T)` 参数表。

## 三、注意事项

1. 模型适用范围  
- 当前温度相关 Norton 参数来自 Inconel 625 数据，不应直接用于其它合金。  
- 若换材料，建议建立新的温度-参数表（或新模型文件）。

2. 应力范围限制  
- 文献数据应力范围约为 `20~400 MPa`。  
- 当计算应力超出该范围，`norton.py` 会钳制到该范围内计算速率。  
- 这能提高稳定性，但也意味着超范围结果是“保守近似”，非严格外推。

3. 数值稳定性  
- RK4 积分仍可能在激进参数下产生不稳定。  
- 若出现非有限值/异常，优先减小 `--dt`（如 `0.01` 或 `0.001`）。

4. 自定义模型接口  
- 自定义文件必须实现 `creep_rate(...)` 并返回有限浮点数。  
- 返回 `NaN/Inf` 或抛异常会被主程序识别并报错。

## 四、当前推荐运行方式

```bash
python two_bar_creep.py --T-a 20 --T-b 800 --hours 500 --k 2 --dt 0.01 --creep-model creep_models/norton.py
```

切换模型示例：

```bash
python two_bar_creep.py --hours 500 --k 2 --creep-model creep_models/time_hardened.py
python two_bar_creep.py --hours 500 --k 2 --creep-model creep_models/garofalo.py
```
