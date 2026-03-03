# creep_2bar

两杆模型（A 常温弹性杆 + B 高温蠕变杆）热应力-蠕变耦合计算示例，支持截面积比 k=A_A/A_B。

## 文件

- `two_bar_creep.py`：计算脚本（温度相关 `E(T)`、按“20°C到T平均热膨胀系数”重建热应变 + Norton 蠕变 + RK4 积分）
- `说明书.md`：详细公式、参数和运行说明

## 快速开始

```bash
python3 two_bar_creep.py --hours 500 --k 2
```
