# 三物种系统解析分析 - 简化参数形式

## 参数简化

定义**净相互作用参数**（Net Interaction Parameters）：

$$
\begin{align}
a &:= (1-\omega)\sigma_{SG} - \omega \cdot \alpha_{SG} \quad &\text{(G → S net effect)} \\
c &:= (1-\omega)\sigma_{GS} - \omega \cdot \alpha_{GS} \quad &\text{(S → G net effect)} \\
b &:= \omega \cdot \sigma_{MG} - (1-\omega) \cdot \alpha_{MG} \quad &\text{(G → M net effect)} \\
e &:= \omega \cdot \sigma_{GM} - (1-\omega) \cdot \alpha_{GM} \quad &\text{(M → G net effect)} \\
d &:= 2\omega - 1 \quad &\text{(generalist pathway balance)}
\end{align}
$$

### 参数意义

- **a, c**: S-G之间的相互作用
  - 当 $\omega$ 小（G偏向代谢物途径）：合作项 $\sigma$ 占主导
  - 当 $\omega$ 大（G偏向底物途径）：竞争项 $\alpha$ 占主导

- **b, e**: M-G之间的相互作用
  - $b$: G对M的影响（$\omega$ 大时合作，小时竞争）
  - $e$: M对G的影响

- **d**: Generalist的"身份认同"
  - $d < 0$ ($\omega < 0.5$): G更像M（代谢物专家）
  - $d = 0$ ($\omega = 0.5$): 真正的通才
  - $d > 0$ ($\omega > 0.5$): G更像S（底物专家）

---

## 系统方程（简化形式）

### 原始形式

$$
\begin{align}
\dot{s} &= r_S \cdot s \left(1 + \sigma_{SM} \cdot m + a \cdot g - s\right) \\
\dot{m} &= r_M \cdot m \left(-1 + \sigma_{MS} \cdot s + b \cdot g - m\right) \\
\dot{g} &= r_G \cdot g \left(d + c \cdot s + e \cdot m - g\right)
\end{align}
$$

### 向量场形式

定义生长率函数：
$$
\begin{align}
f_S(s,m,g) &= 1 + \sigma_{SM} \cdot m + a \cdot g - s \\
f_M(s,m,g) &= -1 + \sigma_{MS} \cdot s + b \cdot g - m \\
f_G(s,m,g) &= d + c \cdot s + e \cdot m - g
\end{align}
$$

则：
$$\dot{\mathbf{N}} = \begin{pmatrix} \dot{s} \\ \dot{m} \\ \dot{g} \end{pmatrix} = \begin{pmatrix} r_S \cdot s \cdot f_S \\ r_M \cdot m \cdot f_M \\ r_G \cdot g \cdot f_G \end{pmatrix}$$

---

## 平衡点推导

### E1: 灭绝 (0, 0, 0)
总是存在。

---

### E2: S-only (1, 0, 0)
总是存在。

**G入侵条件**：检查 $f_G(1, 0, 0)$ 的符号
$$f_G(1, 0, 0) = d + c \cdot 1 + e \cdot 0 - 0 = d + c$$

$$\boxed{\text{G can invade} \Leftrightarrow d + c > 0}$$

即：
$$(2\omega - 1) + [(1-\omega)\sigma_{GS} - \omega\alpha_{GS}] > 0$$

---

### E3: S-M coexistence (s*, m*, 0)

从S-M子系统（之前已推导）：
$$
\begin{align}
s^* &= \frac{1 - \sigma_{SM}}{1 - \sigma_{MS}\cdot\sigma_{SM}} \\
m^* &= \frac{\sigma_{MS} - 1}{1 - \sigma_{MS}\cdot\sigma_{SM}}
\end{align}
$$

**G入侵条件**：检查 $f_G(s^*, m^*, 0)$

$$f_G(s^*, m^*, 0) = d + c \cdot s^* + e \cdot m^*$$

$$\boxed{\text{G can invade S-M} \Leftrightarrow d + c \cdot s^* + e \cdot m^* > 0}$$

代入 $s^*, m^*$：
$$d + c \cdot \frac{1-\sigma_{SM}}{1-\sigma_{MS}\cdot\sigma_{SM}} + e \cdot \frac{\sigma_{MS}-1}{1-\sigma_{MS}\cdot\sigma_{SM}} > 0$$

$$\frac{d(1-\sigma_{MS}\cdot\sigma_{SM}) + c(1-\sigma_{SM}) + e(\sigma_{MS}-1)}{1-\sigma_{MS}\cdot\sigma_{SM}} > 0$$

分子：
$$\boxed{d(1-\sigma_{MS}\cdot\sigma_{SM}) + c(1-\sigma_{SM}) + e(\sigma_{MS}-1)}$$

这就是**G入侵S-M平衡点的条件**！

---

### E4: S-G coexistence (s*, 0, g*)

平衡条件：
$$
\begin{align}
f_S(s, 0, g) &= 1 + a \cdot g - s = 0 \quad &\Rightarrow s = 1 + a \cdot g \\
f_G(s, 0, g) &= d + c \cdot s - g = 0 \quad &\Rightarrow g = d + c \cdot s
\end{align}
$$

代入求解：
$$g = d + c(1 + a \cdot g) = d + c + ac \cdot g$$

$$g(1 - ac) = d + c$$

$$\boxed{g^* = \frac{d + c}{1 - ac}}$$

$$\boxed{s^* = 1 + a \cdot g^* = \frac{1 - ac + a(d+c)}{1 - ac} = \frac{1 + ad}{1 - ac}}$$

**存在条件**：
- $g^* > 0$: 要求 $(d+c)$ 和 $(1-ac)$ 同号
- $s^* > 0$: 要求 $(1+ad)$ 和 $(1-ac)$ 同号

**M入侵条件**：
$$f_M(s^*, 0, g^*) = -1 + \sigma_{MS} \cdot s^* + b \cdot g^*$$

$$\boxed{\text{M can invade} \Leftrightarrow \sigma_{MS} \cdot s^* + b \cdot g^* > 1}$$

---

### E5: M-G coexistence (0, m*, g*)

M无法独立于S生存（$f_M$ 基础项为-1），所以此平衡点**不存在**。

---

### E6: 三物种共存 (s*, m*, g*)

平衡条件：
$$
\begin{align}
1 + \sigma_{SM} \cdot m + a \cdot g - s &= 0 \quad &\text{...(1)} \\
-1 + \sigma_{MS} \cdot s + b \cdot g - m &= 0 \quad &\text{...(2)} \\
d + c \cdot s + e \cdot m - g &= 0 \quad &\text{...(3)}
\end{align}
$$

这是一个**3×3线性系统**（在s, m, g中）。

**矩阵形式**：
$$\begin{pmatrix}
-1 & \sigma_{SM} & a \\
\sigma_{MS} & -1 & b \\
c & e & -1
\end{pmatrix} \begin{pmatrix} s \\ m \\ g \end{pmatrix} = \begin{pmatrix} -1 \\ 1 \\ -d \end{pmatrix}$$

**解析解**：使用Cramer法则或符号计算可得（表达式很长）。

关键观察：
- 系数矩阵的行列式 $\Delta$ 包含所有相互作用参数
- $s^*, m^*, g^*$ 都是 $a, b, c, d, e, \sigma_{SM}, \sigma_{MS}$ 的函数
- 特别地，都是 $\omega$ 的函数（通过 $a, b, c, d, e$）

---

## Jacobian矩阵

$$
J = \begin{pmatrix}
\frac{\partial \dot{s}}{\partial s} & \frac{\partial \dot{s}}{\partial m} & \frac{\partial \dot{s}}{\partial g} \\
\frac{\partial \dot{m}}{\partial s} & \frac{\partial \dot{m}}{\partial m} & \frac{\partial \dot{m}}{\partial g} \\
\frac{\partial \dot{g}}{\partial s} & \frac{\partial \dot{g}}{\partial m} & \frac{\partial \dot{g}}{\partial g}
\end{pmatrix}
$$

**在三物种平衡点 $(s^*, m^*, g^*)$ 处**：

利用平衡条件 $f_S = f_M = f_G = 0$，可简化：

$$
J(s^*, m^*, g^*) = \begin{pmatrix}
-r_S \cdot s^* & r_S \cdot s^* \cdot \sigma_{SM} & r_S \cdot s^* \cdot a \\
r_M \cdot m^* \cdot \sigma_{MS} & -r_M \cdot m^* & r_M \cdot m^* \cdot b \\
r_G \cdot g^* \cdot c & r_G \cdot g^* \cdot e & -r_G \cdot g^*
\end{pmatrix}
$$

**迹（Trace）**：
$$\mathrm{Tr}(J) = -(r_S \cdot s^* + r_M \cdot m^* + r_G \cdot g^*) < 0$$

总是负的！（假设所有种群密度 > 0）

**行列式（Determinant）**：
$$\det(J) = -r_S \cdot r_M \cdot r_G \cdot s^* \cdot m^* \cdot g^* \cdot [1 - \sigma_{SM}\cdot\sigma_{MS} + \text{(G相关项)}]$$

（完整表达式非常复杂，但可以符号计算）

---

## 分岔分析（ω作为分岔参数）

### 关键观察

所有净相互作用参数都是 $\omega$ 的函数：
$$
\begin{align}
a(\omega) &= (1-\omega)\sigma_{SG} - \omega\alpha_{SG} = \sigma_{SG} - \omega(\sigma_{SG} + \alpha_{SG}) \\
c(\omega) &= (1-\omega)\sigma_{GS} - \omega\alpha_{GS} = \sigma_{GS} - \omega(\sigma_{GS} + \alpha_{GS}) \\
b(\omega) &= \omega\sigma_{MG} - (1-\omega)\alpha_{MG} = -\alpha_{MG} + \omega(\sigma_{MG} + \alpha_{MG}) \\
e(\omega) &= \omega\sigma_{GM} - (1-\omega)\alpha_{GM} = -\alpha_{GM} + \omega(\sigma_{GM} + \alpha_{GM}) \\
d(\omega) &= 2\omega - 1
\end{align}
$$

都是 $\omega$ 的**线性函数**！

### 分岔点的位置

**第一个临界点 $\omega_{crit1}$**：G开始能入侵S-M平衡点

$$d + c \cdot s^*_{SM} + e \cdot m^*_{SM} = 0$$

这是关于 $\omega$ 的线性方程（因为 $d, c, e$ 都是 $\omega$ 的线性函数）：

$$(2\omega - 1) + [\sigma_{GS} - \omega(\sigma_{GS}+\alpha_{GS})] \cdot s^*_{SM} + [-\alpha_{GM} + \omega(\sigma_{GM}+\alpha_{GM})] \cdot m^*_{SM} = 0$$

整理后可解出：
$$\boxed{\omega_{crit1} = \frac{1 - \sigma_{GS} \cdot s^*_{SM} + \alpha_{GM} \cdot m^*_{SM}}{2 - (\sigma_{GS}+\alpha_{GS}) \cdot s^*_{SM} + (\sigma_{GM}+\alpha_{GM}) \cdot m^*_{SM}}}$$

**第二个临界点 $\omega_{crit2}$**：M开始无法在三物种系统中生存

在三物种平衡点，$m^* \to 0$ 的条件（需要完整的三物种解析解）。

---

## 分岔类型：Transcritical

### 在 $\omega = \omega_{crit1}$

**平衡点交换**：
- **Before** ($\omega < \omega_{crit1}$):
  - E3 (S-M) 稳定
  - E6 (三物种) 不存在或 $g^* < 0$

- **At** ($\omega = \omega_{crit1}$):
  - G的入侵率 = 0
  - E3和E6"相遇"

- **After** ($\omega > \omega_{crit1}$):
  - E3 (S-M) 不稳定（G可入侵）
  - E6 (三物种) 稳定且 $g^* > 0$

这是**经典的跨临界分岔**！

---

## 生态学解释

### 参数d的核心作用

$$d = 2\omega - 1$$

- **$\omega = 0$**: G完全依赖代谢物途径 → $d = -1$
  - G的基础生长率 = -1（像M一样，无法独立生存）
  - G与M竞争同一生态位

- **$\omega = 0.5$**: G真正的通才 → $d = 0$
  - G的基础生长率 = 0（刚好自给自足）
  - 平衡两种途径

- **$\omega = 1$**: G完全依赖底物途径 → $d = 1$
  - G的基础生长率 = 1（像S一样，可独立生存）
  - G与S竞争同一生态位

### 共存窗口

三物种共存需要 $\omega$ 在中间范围：
$$\boxed{\omega_{crit1} < \omega < \omega_{crit2}}$$

**生态学直觉**：
- $\omega$ 太小：G太像M，竞争排斥M
- $\omega$ 太大：G太像S，可能排斥S或被M排斥
- **中间值**：G有自己的生态位，促进共存

---

## 关键预测

### 实验可验证的预测

1. **共存窗口**：
   - 改变G的代谢途径比例（通过基因工程或环境条件）
   - 预测：存在最优 $\omega$ 范围使三物种共存

2. **分岔点可测量**：
   - 在 $\omega_{crit1}$：G从无到有
   - 在 $\omega_{crit2}$：M从有到无

3. **参数依赖性**：
   - $\omega_{crit1}$ 依赖于 S-M 互惠强度 $(\sigma_{SM}, \sigma_{MS})$
   - 更强的S-M互惠 → 更难被G入侵 → $\omega_{crit1}$ 更大

---

## 总结

使用简化参数 $(a, b, c, d, e)$ 的优势：

1. ✅ **数学简洁**：方程形式更对称
2. ✅ **生态学直观**：净相互作用意义清晰
3. ✅ **分岔分析简化**：都是 $\omega$ 的线性函数
4. ✅ **符合Jeff Gore风格**：强调机制而非复杂公式

**下一步**：
- 数值验证解析预测
- 绘制完整的分岔图
- 参数空间的二维切片
- 与实验数据对比
