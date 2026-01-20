# 如何运行和获得结果 / How to Run and Get Results

**快速参考指南 / Quick Reference Guide**

---

## 🚀 最简单的运行方式 / Simplest Way to Run

### 中文步骤：

```bash
# 1. 进入模型目录
cd models/three_species_crossfeeding

# 2. 安装依赖（只需要运行一次）
pip3 install numpy scipy matplotlib seaborn

# 3. 运行快速入门脚本
python3 quick_start.py
```

**运行时间**: 约15秒
**生成结果**: 2张图片在 `figures/` 文件夹

---

### English Steps:

```bash
# 1. Navigate to model directory
cd models/three_species_crossfeeding

# 2. Install dependencies (only once)
pip3 install numpy scipy matplotlib seaborn

# 3. Run quick start script
python3 quick_start.py
```

**Runtime**: ~15 seconds
**Output**: 2 figures in `figures/` folder

---

## 📊 查看结果 / View Results

### 生成的图片文件 / Generated Figure Files

运行完成后，检查以下文件夹：
After running, check this folder:

```bash
ls figures/
```

你会看到 / You will see:

```
quick_start_timeseries.png         # 时间序列图 / Time series
quick_start_omega_comparison.png   # ω参数对比 / ω comparison
```

### 如何打开图片 / How to Open Figures

**方法1 - 在文件管理器中**:
```bash
# Linux
xdg-open figures/quick_start_timeseries.png

# macOS
open figures/quick_start_timeseries.png

# Windows
start figures/quick_start_timeseries.png
```

**方法2 - 直接浏览**:
在文件浏览器中打开 `figures/` 文件夹，双击图片查看

---

## 🎓 三种运行级别 / Three Levels of Usage

### Level 1: 快速入门 (初学者) / Quick Start (Beginners)

**命令** / **Command**:
```bash
python3 quick_start.py
```

**适合** / **For**:
- ✅ 第一次使用
- ✅ 快速验证模型工作正常
- ✅ 了解基本功能

**获得结果** / **Results**:
- 控制台输出：平衡点、稳定性分析
- 2张图片：时间序列和ω影响

**示例输出** / **Example Output**:
```
最终种群密度 (Final populations):
  N_S = 99.69    ← S-specialist 存活
  N_M = 0.00     ← M-specialist 灭绝
  N_G = 5.04     ← Generalist 低密度共存

找到 3 个平衡点
平衡点 3: S-G coexistence → 稳定 (Stable) ✓
```

---

### Level 2: 交互式示例 (进阶) / Interactive Examples (Advanced)

**命令** / **Command**:
```bash
python3 交互式示例.py
```

**适合** / **For**:
- ✅ 系统化探索参数影响
- ✅ 理解模型行为
- ✅ 测试不同场景

**获得结果** / **Results**:
- 5个详细示例的分析
- 4张对比图片
- 参数敏感性结果

**示例包括** / **Examples Include**:
1. **基础模拟** - 不同初始条件
2. **ω参数扫描** - 路径权重影响
3. **合作强度** - M-specialist生存条件
4. **初始条件** - 系统鲁棒性
5. **共存搜索** - 三物种共存参数

**示例输出** / **Sample Output**:
```
ω值 | N_S最终 | N_M最终 | N_G最终
----------------------------------------
0.1 |   99.91 |    0.00 |    0.00
0.3 |   99.95 |    0.00 |    0.00
0.5 |   99.74 |    0.00 |    4.99
0.7 |   93.98 |    0.00 |   31.54  ← ω越大，G越强
0.9 |   79.59 |    0.00 |   61.66
```

---

### Level 3: Jupyter笔记本 (专业) / Jupyter Notebook (Professional)

**命令** / **Command**:
```bash
jupyter notebook notebooks/three_species_phase_analysis.ipynb
```

**适合** / **For**:
- ✅ 完整的相平面分析
- ✅ 发表级图片
- ✅ 深度参数探索
- ✅ 学术研究

**获得结果** / **Results**:
- 7张发表质量的图片
- 完整的数学分析
- 生物学解释
- 分岔图、相平面图、3D可视化

**包含内容** / **Includes**:
1. 模型初始化和参数探索
2. 平衡点寻找和稳定性分析
3. 时间序列动力学
4. 2D相平面图（S-M, S-G, M-G）
5. 3D相空间可视化
6. ω分岔分析
7. 合作-竞争敏感性分析
8. 生态学洞见总结

---

## 📖 如何理解结果 / How to Interpret Results

### 种群密度解释 / Population Density Interpretation

| 最终值 | 含义 |
|--------|------|
| > 50 | 物种占主导地位 |
| 5-50 | 物种共存但非主导 |
| 1-5 | 物种低密度存活 |
| < 1 | 物种灭绝 |

### 平衡点类型 / Equilibrium Types

| 类型 | 稳定性 | 含义 |
|------|--------|------|
| 灭绝 Extinction | 通常不稳定 | 所有物种灭绝 |
| S单独 S-only | 取决于参数 | 只有S存活 |
| S-G共存 S-G coexistence | **常见且稳定** | S主导，G低密度 |
| 三物种共存 Three-species | 罕见 | 需要特殊参数 |

### ω参数的影响 / Effect of ω Parameter

```
ω = 0    → Generalist表现像M-specialist（代谢物路径）
         → 竞争力弱，易被排除

ω = 0.5  → Generalist平衡使用两条路径
         → 可以与S共存

ω = 1    → Generalist表现像S-specialist（底物路径）
         → 竞争力强，种群大
```

---

## 🔧 常见问题排查 / Troubleshooting

### 问题1：找不到模块 / Module not found

```bash
# 错误信息 / Error:
ModuleNotFoundError: No module named 'numpy'

# 解决方法 / Solution:
pip3 install numpy scipy matplotlib seaborn
```

### 问题2：无法打开Jupyter / Cannot open Jupyter

```bash
# 安装Jupyter / Install Jupyter:
pip3 install jupyter

# 启动 / Launch:
jupyter notebook
```

### 问题3：M-specialist总是灭绝 / M-specialist always dies

这是**正常现象**！因为：
This is **expected**! Because:

- M是强制性互养生物（基础生长率 = -1）
- 需要 σ_MS > 某个阈值才能生存
- 默认参数下无法满足条件

**如何让M存活** / **How to make M survive**:
```python
# 修改参数 / Modify parameters:
model.params['sigma_MS'] = 0.9  # 增大合作强度
model.params['K_S'] = 150.0     # 增加S的携带容量
```

### 问题4：中文显示乱码 / Chinese characters display incorrectly

**不影响结果！** / **Does not affect results!**

图片中的中文可能显示为方块，但：
- 数据完全正确
- 英文标签清晰
- 可以忽略这些警告

---

## 📁 文件结构参考 / File Structure Reference

```
three_species_crossfeeding/
├── src/
│   ├── three_species_model.py         ← 核心模型 / Core model
│   └── phase_plane_analysis.py        ← 分析工具 / Analysis tools
│
├── notebooks/
│   └── three_species_phase_analysis.ipynb  ← 完整分析 / Full analysis
│
├── figures/                            ← 生成的图片 / Generated figures
│   ├── quick_start_*.png
│   ├── example*.png
│   └── (运行notebook后生成更多 / More after running notebook)
│
├── quick_start.py                      ← ⭐ 从这里开始 / START HERE
├── 交互式示例.py                       ← 交互示例 / Interactive examples
├── example_analysis.py                 ← 完整示例 / Full example
│
├── 使用说明_中文.md                    ← 详细中文文档 / Detailed Chinese docs
├── README.md                           ← 英文文档 / English documentation
└── requirements.txt                    ← 依赖清单 / Dependencies
```

---

## 🎯 推荐学习路径 / Recommended Learning Path

### 第1天 / Day 1: 基础入门
```bash
python3 quick_start.py
```
- 理解模型输出
- 查看生成的图片
- 阅读控制台结果

### 第2天 / Day 2: 深入探索
```bash
python3 交互式示例.py
```
- 理解参数影响
- 对比不同场景
- 修改参数重新运行

### 第3天 / Day 3: 专业分析
```bash
jupyter notebook notebooks/three_species_phase_analysis.ipynb
```
- 运行所有单元格
- 理解相平面分析
- 生成发表级图片

### 第4天 / Day 4: 自定义研究
- 修改参数探索新场景
- 设计自己的实验
- 比较模型预测和实验数据

---

## 💡 快速示例代码 / Quick Example Code

### Python交互式使用 / Python Interactive Usage

```python
# 启动Python / Start Python
python3

# 运行以下代码 / Run this code:
import sys
sys.path.append('src')

from three_species_model import ThreeSpeciesModel
import numpy as np

# 创建模型 / Create model
model = ThreeSpeciesModel()

# 设置初始条件 / Set initial conditions
N0 = np.array([50.0, 50.0, 50.0])

# 模拟 / Simulate
sol = model.simulate(N0, (0, 100))

# 查看结果 / View results
print(f"Final: S={sol['N_S'][-1]:.1f}, M={sol['N_M'][-1]:.1f}, G={sol['N_G'][-1]:.1f}")

# 寻找平衡点 / Find equilibria
equilibria = model.find_equilibria()
print(f"Found {len(equilibria)} equilibria")

# 稳定性分析 / Stability analysis
for eq in equilibria:
    stability = model.stability_analysis(eq)
    print(f"Equilibrium: {eq}")
    print(f"Stable: {stability['stable']}\n")
```

---

## ✅ 检查清单 / Checklist

运行前确认 / Before running:
- [ ] 已进入正确目录 `cd models/three_species_crossfeeding`
- [ ] 已安装依赖 `pip3 install -r requirements.txt`
- [ ] Python版本 ≥ 3.8 (`python3 --version`)

运行后确认 / After running:
- [ ] 控制台有输出结果
- [ ] `figures/` 文件夹有图片
- [ ] 图片可以正常打开
- [ ] 理解了输出的含义

---

## 📞 获取帮助 / Get Help

1. **查看详细文档** / **Read detailed docs**:
   - 中文：`使用说明_中文.md`
   - English: `README.md`

2. **查看代码示例** / **Check code examples**:
   - `quick_start.py` - 带注释的简单示例
   - `交互式示例.py` - 5个详细场景

3. **检查函数文档** / **Check function docs**:
   ```python
   from three_species_model import ThreeSpeciesModel
   help(ThreeSpeciesModel)
   help(ThreeSpeciesModel.simulate)
   ```

---

## 🎉 成功标志 / Success Indicators

运行成功的标志 / Signs of successful run:

✅ 控制台输出类似：
```
✓ 模拟成功！
✓ Simulation successful!
最终种群密度: N_S = 99.69, N_M = 0.00, N_G = 5.04
✓ 图片已保存到: figures/quick_start_timeseries.png
```

✅ 生成的图片显示：
- 清晰的曲线
- 坐标轴标签
- 图例说明

✅ 可以理解：
- 哪些物种存活
- 为什么某些物种灭绝
- 参数如何影响结果

---

**祝您使用愉快！**
**Happy modeling!**

🔬 三物种交叉喂养模型团队
   Three-Species Cross-Feeding Model Team
