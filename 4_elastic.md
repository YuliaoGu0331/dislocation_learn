# 4. 位错的弹性性质（Elastic Properties of Dislocations）

> **本章概要**：位错使原子偏离完整晶格位置，并在其周围产生长程内应力场。除位错核心附近外，这些应力与应变通常足够小，可用小变形线弹性理论描述。本章在各向同性近似下，讨论直位错的应力场、位错能、位错受力、位错间相互作用、攀移的化学力以及自由表面引起的镜像力。

---

## 4.1 基本认识与适用范围

位错周围的晶格畸变会产生**内应力（internal stress）**。在外载荷作用下，晶体中某点的总应力为内应力与外应力之和：

$$
\boldsymbol{\sigma}^{\mathrm{tot}}
=\boldsymbol{\sigma}^{\mathrm{int}}
+\boldsymbol{\sigma}^{\mathrm{ext}}.
$$

在线弹性叠加近似下，外载荷不会改变给定位错构型本身的固有弹性场，而是与该内应力场线性叠加；若外载荷驱动位错改变位置或形状，内应力场才会随构型重新分布。

本章采用以下近似：

- 材料为均匀、各向同性的线弹性体；
- 变形与应变均为小量；
- 位错线近似为无限长直线；
- 连续介质弹性理论仅用于位错核心半径 $r_0$ 之外，即 $r\ge r_0$ 的区域。

> [!question] 为什么位错产生塑性变形，却能用弹性理论分析？
> 位错的不可逆运动产生宏观塑性变形，但在任一给定构型下，位错周围除核心区外的局部晶格畸变仍可视为小的弹性畸变。因此，可以用弹性理论求位错的瞬时应力场；该应力场又会反过来驱动位错滑移或攀移。换言之，**塑性过程由位错运动承担，而位错之间的长程相互作用主要由弹性场描述**。

### 4.1.1 主要符号

| 符号 | 含义 |
| --- | --- |
| $\mathbf{u}$ | 位移矢量（displacement vector） |
| $e_{ij}$ | 小应变张量分量（small-strain tensor component） |
| $\sigma_{ij}$ | 应力张量分量（stress tensor component） |
| $G$ | 剪切模量（shear modulus） |
| $E$ | 杨氏模量（Young's modulus） |
| $K$ | 体积模量（bulk modulus） |
| $\lambda$ | Lamé 第一常数（first Lamé parameter） |
| $\nu$ | 泊松比（Poisson's ratio） |
| $\mathbf b$、$b=|\mathbf b|$ | 伯格斯矢量（Burgers vector）及其模 |
| $E_{\mathrm{el}}$ | 位错单位长度的弹性应变能（除非另有说明） |
| $r_0$ | 位错核心的有效截断半径 |
| $R$ | 弹性场的外截断半径；在线张力问题中也表示曲率半径，需结合上下文判断 |

---

## 4.2 小变形各向同性线弹性理论

### 4.2.1 位移与应变

设位移场为 $\mathbf u=(u_x,u_y,u_z)$。法向应变为：

$$
e_{xx}=\frac{\partial u_x}{\partial x},\qquad
e_{yy}=\frac{\partial u_y}{\partial y},\qquad
e_{zz}=\frac{\partial u_z}{\partial z}.
\tag{4.2}
$$

剪切应变为：

$$
\begin{aligned}
e_{yz}=e_{zy}
&=\frac12\left(\frac{\partial u_y}{\partial z}+\frac{\partial u_z}{\partial y}\right),\\
e_{zx}=e_{xz}
&=\frac12\left(\frac{\partial u_z}{\partial x}+\frac{\partial u_x}{\partial z}\right),\\
e_{xy}=e_{yx}
&=\frac12\left(\frac{\partial u_x}{\partial y}+\frac{\partial u_y}{\partial x}\right).
\end{aligned}
\tag{4.3}
$$

体积应变（dilatation）为：

$$
\Delta=\frac{\Delta V}{V}=e_{xx}+e_{yy}+e_{zz}.
\tag{4.4}
$$

### 4.2.2 应力、静水压力与本构关系

在无体偶力的经典连续介质中，应力张量对称：

$$
\sigma_{yz}=\sigma_{zy},\qquad
\sigma_{zx}=\sigma_{xz},\qquad
\sigma_{xy}=\sigma_{yx}.
\tag{4.5}
$$

作用于体积元的有效静水压力（压应力取正）定义为：

$$
p=-\frac13\left(\sigma_{xx}+\sigma_{yy}+\sigma_{zz}\right).
\tag{4.6}
$$

各向同性线弹性体的 Hooke 定律为：

$$
\begin{aligned}
\sigma_{xx}&=2Ge_{xx}+\lambda(e_{xx}+e_{yy}+e_{zz}),\\
\sigma_{yy}&=2Ge_{yy}+\lambda(e_{xx}+e_{yy}+e_{zz}),\\
\sigma_{zz}&=2Ge_{zz}+\lambda(e_{xx}+e_{yy}+e_{zz}),\\
\sigma_{xy}&=2Ge_{xy},\qquad
\sigma_{yz}=2Ge_{yz},\qquad
\sigma_{zx}=2Ge_{zx}.
\end{aligned}
\tag{4.7}
$$

常用弹性常数之间的关系为：

$$
E=2G(1+\nu),\qquad
\nu=\frac{\lambda}{2(\lambda+G)},\qquad
K=\frac{E}{3(1-2\nu)}.
\tag{4.8}
$$

### 4.2.3 弹性应变能

体积元 $\mathrm dV$ 中储存的弹性应变能为：

$$
\mathrm dE_{\mathrm{el}}
=\frac12\left(\sum_{i=x,y,z}\sum_{j=x,y,z}\sigma_{ij}e_{ij}\right)\mathrm dV
=\frac12\sigma_{ij}e_{ij}\,\mathrm dV,
\tag{4.9}
$$

其中最后一个等号采用 Einstein 求和约定。

---

## 4.3 直位错的弹性场

统一取位错线沿 $z$ 轴。对刃型位错，另取 $\mathbf b=b\mathbf e_x$，滑移面为 $y=0$；并记

$$
r=\sqrt{x^2+y^2}.
$$

### 4.3.1 螺型位错（Screw Dislocation）

考虑伯格斯矢量为 $\mathbf b=b\mathbf e_z$ 的螺型位错。由于问题绕 $z$ 轴呈轴对称，采用圆柱坐标 $(r,\theta,z)$ 最方便。其位移场可写为：

$$
u_z=\frac{b\theta}{2\pi},\qquad u_r=u_\theta=0.
$$

非零应变与应力分量为：

$$
e_{\theta z}=e_{z\theta}=\frac{b}{4\pi r},\qquad
\sigma_{\theta z}=\sigma_{z\theta}=\frac{Gb}{2\pi r}.
\tag{4.15}
$$

$\sigma_{\theta z}$ 与 $\sigma_{z\theta}$ 是相互对应的剪切应力分量。由上式可见：

- 螺型位错只产生剪切应力，不产生法向应力或静水压力；
- 应力大小仅与到位错线的距离 $r$ 有关，同一圆柱面上大小相同；
- 应力按 $1/r$ 衰减，因此是长程应力场。

这组剪切应力会在垂直于位错线的端面上产生扭矩；因此，对有限长圆柱中的螺型位错，若边界允许，晶体还会出现附加的整体扭转以满足力矩平衡。

> [!note] 为什么在核心外采用圆柱区域分析？
> 螺型位错具有轴对称性，圆柱坐标能直接体现其 $1/r$ 应力场。公式在 $r\to0$ 时发散，但真实晶体中的应力不可能无限大，这意味着线弹性理论在位错核心区失效，需由原子尺度的非线性理论描述。

完整晶体的理论临界剪切应力可估为：

$$
\tau_{\mathrm{th}}=\frac{bG}{2\pi a},
\tag{1.5}
$$

其中 $a$ 为原子尺度的晶格周期。将 $\sigma_{\theta z}=Gb/(2\pi r)$ 与 $\tau_{\mathrm{th}}$ 比较，可得线弹性理论开始失效的尺度 $r\sim a\sim b$。实际常取位错核心半径

$$
r_0\sim b\text{--}4b.
$$

### 4.3.2 刃型位错（Edge Dislocation）

刃型位错是一个**平面应变（plane strain）**问题，即 $u_z=0$ 且 $e_{zz}=0$，但由于泊松效应一般有 $\sigma_{zz}\ne0$。其应力场为：

$$
\begin{aligned}
\sigma_{xx}&=-D\,y\frac{3x^2+y^2}{(x^2+y^2)^2},\\
\sigma_{yy}&=D\,y\frac{x^2-y^2}{(x^2+y^2)^2},\\
\sigma_{xy}=\sigma_{yx}&=D\,x\frac{x^2-y^2}{(x^2+y^2)^2},\\
\sigma_{zz}&=\nu(\sigma_{xx}+\sigma_{yy}),\\
\sigma_{xz}=\sigma_{zx}&=\sigma_{yz}=\sigma_{zy}=0,
\end{aligned}
\qquad
D=\frac{Gb}{2\pi(1-\nu)}.
\tag{4.16}
$$

与螺型位错不同，刃型位错同时产生拉压应力与剪切应力。对应的有效静水压力为：

$$
p=\frac23(1+\nu)D\frac{y}{x^2+y^2}.
\tag{4.17}
$$

因此刃型位错滑移面上、下两侧分别形成拉伸区和压缩区，其静水应力可与点缺陷或溶质原子发生相互作用。上述连续介质解仍仅适用于 $r\ge r_0$。

### 4.3.3 混合型位错（Mixed Dislocation）

对于混合型位错，可将伯格斯矢量分解为平行和垂直于位错线方向 $\boldsymbol\xi$ 的分量：

$$
\mathbf b=\mathbf b_{\mathrm{s}}+\mathbf b_{\mathrm{e}},\qquad
\mathbf b_{\mathrm{s}}\parallel\boldsymbol\xi,
\qquad
\mathbf b_{\mathrm{e}}\perp\boldsymbol\xi.
$$

其中 $\mathbf b_{\mathrm{s}}$ 和 $\mathbf b_{\mathrm{e}}$ 分别对应螺型与刃型分量。由于小变形线弹性方程是线性的，两部分弹性场可分别计算后线性叠加。

---

## 4.4 位错的弹性应变能

### 4.4.1 螺型与刃型位错

对单位长度螺型位错，在 $r_0$ 到外截断半径 $R$ 之间积分可得：

$$
E_{\mathrm{el}}^{\mathrm{screw}}
=\frac{Gb^2}{4\pi}\int_{r_0}^{R}\frac{\mathrm dr}{r}
=\frac{Gb^2}{4\pi}\ln\left(\frac{R}{r_0}\right).
\tag{4.20}
$$

更复杂的几何中，可把位错的形成过程看作剪切面逐步错开一个伯格斯矢量，并通过外力做功计算能量。螺型与刃型位错的单位长度弹性能分别为：

$$
\begin{aligned}
E_{\mathrm{el}}^{\mathrm{screw}}
&=\frac{Gb^2}{4\pi}\ln\left(\frac{R}{r_0}\right),\\
E_{\mathrm{el}}^{\mathrm{edge}}
&=\frac{Gb^2}{4\pi(1-\nu)}\ln\left(\frac{R}{r_0}\right).
\end{aligned}
\tag{4.22}
$$

这里的 $R$ 不能真正取无穷大。对有限晶体，$R$ 可取到自由表面或几何边界的特征距离；对含有大量位错的晶体，其他位错会屏蔽或抵消长程应力场，$R$ 常取位错平均间距量级，有时近似为平均间距的一半。

位错核心能（core energy）的估计较粗略，通常小于核心外的弹性应变能，但并非总可忽略。核心能随位错位置变化会形成晶格周期势垒，从而贡献位错运动的晶格阻力。

### 4.4.2 混合型位错

设 $\theta$ 为 $\mathbf b$ 与位错线方向 $\boldsymbol\xi$ 的夹角，则螺型分量为 $b\cos\theta$，刃型分量为 $b\sin\theta$。混合型位错的单位长度弹性能为：

$$
\begin{aligned}
E_{\mathrm{el}}^{\mathrm{mixed}}
&=\left[
\frac{Gb^2\sin^2\theta}{4\pi(1-\nu)}
+\frac{Gb^2\cos^2\theta}{4\pi}
\right]\ln\left(\frac{R}{r_0}\right)\\
&=\frac{Gb^2(1-\nu\cos^2\theta)}{4\pi(1-\nu)}
\ln\left(\frac{R}{r_0}\right).
\end{aligned}
\tag{4.23}
$$

由于位错能对位错类型及 $R/r_0$ 的依赖较弱，常近似写为：

$$
E_{\mathrm{el}}\approx\alpha Gb^2,
\qquad \alpha\approx0.5\text{--}1.0.
\tag{4.24}
$$

该近似给出了判断位错反应是否有利的简单能量判据。例如反应 $\mathbf b_1+\mathbf b_2\to\mathbf b_3$ 在其他条件相近时，若

$$
b_3^2<b_1^2+b_2^2,
$$

则反应后能量降低，因而在能量上更有利；实际能否发生还取决于反应势垒、几何条件和外加应力。

---

## 4.5 位错所受的力

### 4.5.1 虚功与 Peach–Koehler 力

可通过**虚功原理（principle of virtual work）**定义作用于位错的等效力：位错线元 $\mathrm dl$ 发生虚位移 $\delta\mathbf r$ 时，应力场对它所做的功为

$$
\delta W
=\mathbf f_{\mathrm{PK}}\cdot\delta\mathbf r\,\mathrm dl,
$$

其中 $\mathbf f_{\mathrm{PK}}$ 是**单位长度位错所受的力**。若位错线元在滑移面内移动 $\mathrm ds$，扫过面积 $\mathrm dA=\mathrm ds\,\mathrm dl$，则简单滑移几何下有：

$$
F=\frac{\mathrm dW}{\mathrm ds\,\mathrm dl}
=\frac{\mathrm dW}{\mathrm dA}
=\tau b.
\tag{4.27}
$$

其中 $\tau$ 是滑移面上沿 $\mathbf b$ 方向的分切应力。对任意位错线方向和应力状态，上式推广为 **Peach–Koehler 力**：

$$
\boxed{
\mathbf f_{\mathrm{PK}}
=\left(\boldsymbol\sigma\cdot\mathbf b\right)\times\boldsymbol\xi
}
$$

其分量形式为：

$$
f_i=\varepsilon_{ijk}\sigma_{jl}b_l\xi_k,
$$

其中 $\varepsilon_{ijk}$ 为 Levi–Civita 符号。计算时，$\boldsymbol\sigma$ 应取位错线位置处由**外载荷、其他位错、点缺陷和边界镜像场**等产生的局部应力；不能直接把该位错自身发散的直线弹性应力代入。弯曲位错的自作用通常通过线张力或经过核心正则化的自力另行处理。

### 4.5.2 公式中各项的物理意义

| 项 | 数学含义 | 物理意义 |
| --- | --- | --- |
| $\mathbf f_{\mathrm{PK}}$ | 单位长度力矢量，单位为 $\mathrm{N\,m^{-1}}$ | 应力场驱动位错线元改变位置的机械驱动力；长度为 $L$ 且受力近似均匀的线段所受总力约为 $\mathbf F=\mathbf f_{\mathrm{PK}}L$ |
| $\boldsymbol\sigma$ | 二阶应力张量，单位为 $\mathrm{Pa}$ | 描述位错所在位置的局部应力状态；既可以来自外载荷，也可以来自其他缺陷或边界 |
| $\mathbf b$ | 伯格斯矢量，单位为 $\mathrm m$ | 表征位错造成的晶格闭合差；其方向选取与应力耦合的分量，其大小决定耦合强弱，所以简单情形下有 $f\propto b$ |
| $\boldsymbol\xi$ | 位错线的单位切向矢量 | 规定局部线方向，并通过叉乘把应力—伯格斯矢量耦合转换成垂直于位错线的驱动力 |
| $\boldsymbol\sigma\cdot\mathbf b$ | 分量为 $q_i=\sigma_{ij}b_j$ 的矢量，单位为 $\mathrm{N\,m^{-1}}$ | 可理解为“局部应力对晶格闭合差的耦合”；它决定应力场中哪些分量能对该伯格斯矢量做功 |
| $\times\boldsymbol\xi$ | 与线方向作叉乘 | 决定力的方向，并保证 $\mathbf f_{\mathrm{PK}}\cdot\boldsymbol\xi=0$，即 Peach–Koehler 力没有沿位错线切向的分量 |

量纲上：

$$
[\mathbf f_{\mathrm{PK}}]
=[\boldsymbol\sigma][\mathbf b]
=\mathrm{Pa\cdot m}
=\mathrm{N\,m^{-1}},
$$

因此它是**每单位位错线长度的力**，而不是作用于整条位错的总力。

> [!note] 位错线方向与伯格斯矢量的符号约定
> $\boldsymbol\xi$ 的线向与 $\mathbf b$ 的符号通过 Burgers 回路约定相联系。只反转其中一个矢量会改变位错符号，并使 Peach–Koehler 力反向；同时作变换 $(\mathbf b,\boldsymbol\xi)\to(-\mathbf b,-\boldsymbol\xi)$ 则仍表示同一物理位错，计算所得的力不变。

### 4.5.3 滑移分量、攀移分量与典型情形

对非纯螺型位错，滑移面由 $\mathbf b$ 与 $\boldsymbol\xi$ 张成。定义滑移面单位法向和面内滑移方向：

$$
\mathbf n
=\frac{\boldsymbol\xi\times\mathbf b}
{|\boldsymbol\xi\times\mathbf b|},
\qquad
\mathbf g=\mathbf n\times\boldsymbol\xi.
$$

Peach–Koehler 力可分解为：

$$
\mathbf f_{\mathrm{PK}}
=\underbrace{(\mathbf f_{\mathrm{PK}}\cdot\mathbf g)\mathbf g}_{\mathbf f_{\mathrm g}\text{：滑移力}}
+\underbrace{(\mathbf f_{\mathrm{PK}}\cdot\mathbf n)\mathbf n}_{\mathbf f_{\mathrm c}\text{：攀移力}}.
$$

- **滑移力（glide force）** $\mathbf f_{\mathrm g}$ 位于滑移面内并垂直于位错线，驱动保守滑移；简单情形下其大小为 $f_{\mathrm g}=\tau b$。
- **攀移力（climb force）** $\mathbf f_{\mathrm c}$ 垂直于滑移面，驱动刃型分量发生非保守攀移；实际攀移还需要空位或间隙原子的扩散。

对本章采用的刃型位错坐标 $\mathbf b=b\mathbf e_x$、$\boldsymbol\xi=\mathbf e_z$，有 $\mathbf n=\mathbf e_y$、$\mathbf g=\mathbf e_x$，因此：

$$
\mathbf f_{\mathrm{PK}}
=b\sigma_{xy}\mathbf e_x
-b\sigma_{xx}\mathbf e_y.
$$

这表明 $\sigma_{xy}$ 提供滑移力，而 $\sigma_{xx}$ 提供攀移力，正好对应式 (4.35)。例如静水应力 $\boldsymbol\sigma=-p\mathbf I$ 对纯螺型位错不给出 Peach–Koehler 力，但会对刃型位错产生攀移驱动力。

对纯螺型位错 $\mathbf b\parallel\boldsymbol\xi$，由 $\mathbf b$ 与 $\boldsymbol\xi$ 无法唯一确定滑移面，需要结合晶体学选择具体滑移面。若 $\mathbf b=b\mathbf e_z$、$\boldsymbol\xi=\mathbf e_z$，则：

$$
\mathbf f_{\mathrm{PK}}
=b\sigma_{yz}\mathbf e_x
-b\sigma_{xz}\mathbf e_y.
$$

两个分量分别可驱动螺型位错在包含 $z$ 轴的不同滑移面上运动；对纯螺型位错通常不定义唯一的攀移方向。

> [!important] 力不等于实际运动
> Peach–Koehler 公式给出的是热力学/力学**驱动力**。位错是否运动以及运动多快，还取决于 Peierls–Nabarro 晶格阻力、溶质或析出物钉扎、声子/电子阻尼、温度以及攀移所需的点缺陷扩散。因此，$\mathbf f_{\mathrm{PK}}\ne0$ 并不必然意味着位错立即运动。

### 4.5.4 位错线张力与弯曲位错

位错线张力（line tension）定义为位错线每增加单位长度所需增加的能量，本质上与单位长度位错能同量级：

$$
T\approx\alpha Gb^2.
\tag{4.28}
$$

线张力倾向于缩短并拉直弯曲位错；外加滑移力则使位错弓出。设弯曲位错局部曲率半径为 $R$，几何关系为 $\mathrm d\theta=\mathrm dl/R$。平衡时：

$$
T\,\mathrm d\theta=\tau_0b\,\mathrm dl,
\qquad
\tau_0=\frac{T}{bR}.
\tag{4.29}
$$

采用 $T\approx\alpha Gb^2$ 可得：

$$
\tau_0\approx\frac{\alpha Gb}{R}.
\tag{4.30}
$$

该近似把不同取向位错的单位长度能量视为相同，因此弓出的位错线近似为圆弧；只有在 $\nu=0$ 时这一点才严格成立。考虑位错能对线方向的依赖后，真实线张力为：

$$
T(\theta)=E_{\mathrm{el}}(\theta)
+\frac{\mathrm d^2E_{\mathrm{el}}(\theta)}{\mathrm d\theta^2}.
\tag{4.31}
$$

在均匀应力下，局部曲率半径仍满足式 (4.29)，但整体形状更接近椭圆，长轴平行于伯格斯矢量，轴比约为 $1/(1-\nu)$。多数估算中，式 (4.30) 已足够准确。

---

## 4.6 位错之间的相互作用力

位错 1 的应力场作用于位错 2，可通过式 (4.27) 或 Peach–Koehler 公式求得相互作用力。也可先计算位错 2 在位错 1 应力场中移动所引起的能量变化，再由力等于势能的负梯度得到结果。

### 4.6.1 两条平行刃型位错

考虑两条平行刃型位错，伯格斯矢量均沿 $x$ 方向，位错 2 相对位错 1 的坐标为 $(x,y)$。单位长度位错所受力为：

$$
F_x=\sigma_{xy}b,
\qquad
F_y=-\sigma_{xx}b,
\tag{4.35}
$$

即

$$
\begin{aligned}
F_x&=\frac{Gb^2}{2\pi(1-\nu)}
\frac{x(x^2-y^2)}{(x^2+y^2)^2},\\
F_y&=\frac{Gb^2}{2\pi(1-\nu)}
\frac{y(3x^2+y^2)}{(x^2+y^2)^2}.
\end{aligned}
\tag{4.36}
$$

$F_x$ 是滑移方向的力，$F_y$ 是垂直于滑移面的攀移方向力。对 $y>0$ 且两位错同号的情况，$F_x$ 的规律为：

| $x$ 的范围 | $F_x$ 的符号 | 滑移相互作用 |
| --- | ---: | --- |
| $-\infty<x<-y$ | $-$ | 排斥 |
| $-y<x<0$ | $+$ | 吸引 |
| $0<x<y$ | $-$ | 吸引 |
| $y<x<\infty$ | $+$ | 排斥 |

$F_x=0$ 的位置为 $x=0,\pm y$ 以及 $|x|\to\infty$。其中，同号刃型位错的稳定滑移平衡位置为 $x=0$；异号位错则在 $x=\pm y$ 处具有稳定滑移平衡。因而同号刃型位错倾向于沿垂直方向排列；位于相邻平行滑移面上的异号位错则可形成**位错偶极子（dislocation dipole）**。若异号位错处在同一滑移面上，它们会相互吸引并可能湮灭。

### 4.6.2 两条平行螺型位错

对两条同号、平行的螺型位错，有：

$$
F_r=\sigma_{z\theta}b,
\qquad
F_\theta=-\sigma_{zr}b,
\tag{4.37}
$$

这里采用右手圆柱坐标 $(\mathbf e_r,\mathbf e_\theta,\mathbf e_z)$ 和 $\boldsymbol\xi=\mathbf e_z$；若改变线向或角向定义，分量符号也会相应改变。

代入式 (4.15) 得：

$$
F_r=\frac{Gb^2}{2\pi r},
\qquad
F_\theta=0.
\tag{4.38}
$$

因此同号螺型位错相互排斥，异号螺型位错相互吸引。对于相互平行且分别为纯刃型和纯螺型的两条位错，其 Peach–Koehler 相互作用力为零。

以上结论也说明：无论应力来自外载荷还是其他位错，均可用同一个 Peach–Koehler 公式计算位错受力。

---

## 4.7 攀移力与空位化学势

式 (4.35) 中的 $F_y$ 扮演机械攀移力的角色。攀移力既可来自外加应力，也可来自其他缺陷产生的内部应力。外部攀移力在高温蠕变中十分重要；内部攀移力则可驱动某些保守攀移过程。

位错攀移需要吸收或发射点缺陷，因此除机械力外，还必须考虑空位浓度变化产生的化学力。若单位长度位错受到攀移力 $F$，其平衡空位浓度为：

$$
\begin{aligned}
c
&=\exp\left[-\frac{E_f^{\mathrm v}+F\Omega/b}{kT}\right]\\
&=c_0\exp\left(-\frac{F\Omega}{bkT}\right),
\end{aligned}
\qquad
c_0=\exp\left(-\frac{E_f^{\mathrm v}}{kT}\right),
\tag{4.39}
$$

其中 $E_f^{\mathrm v}$ 为空位形成能，$\Omega$ 为原子体积，$k$ 为 Boltzmann 常数，$T$ 为绝对温度。反过来，空位浓度偏离热平衡值 $c_0$ 时，对位错产生的单位长度化学力为：

$$
f_{\mathrm{chem}}=\frac{bkT}{\Omega}\ln\left(\frac{c}{c_0}\right).
\tag{4.40}
$$

力的正负取决于所采用的攀移方向和空位吸收符号约定。即使是适度的空位过饱和，也可能产生显著大于常见外加机械攀移力的化学驱动力。

攀移速率最终取决于：

- 机械攀移力与化学力的合力；
- 割阶（jog）的浓度与迁移率；
- 空位向位错扩散的速率。

---

## 4.8 自由表面的镜像力（Image Force）

位错通常被自由表面吸引，因为靠近自由表面时其弹性场被削弱，系统弹性能降低；相反，高刚度覆盖层或刚性界面通常排斥位错。

### 4.8.1 镜像位错法

为满足自由表面上的零牵引边界条件，可在表面另一侧引入虚构的**镜像位错（image dislocation）**。考虑自由表面 $x=0$，真实螺型位错位于 $(d,0)$，在 $(-d,0)$ 放置一条符号相反的镜像螺型位错。定义

$$
x_-=x-d,
\qquad
x_+=x+d,
\qquad
A=\frac{Gb}{2\pi}.
$$

真实位错与镜像位错叠加后的应力场为：

$$
\sigma_{zx}
=-\frac{Ay}{x_-^2+y^2}
+\frac{Ay}{x_+^2+y^2},
\tag{4.41}
$$

$$
\sigma_{zy}
=\frac{Ax_-}{x_-^2+y^2}
-\frac{Ax_+}{x_+^2+y^2}.
\tag{4.42}
$$

镜像位错在真实螺型位错处产生的单位长度作用力为：

$$
F_x=-\frac{Gb^2}{4\pi d}.
\tag{4.43}
$$

负号表示力指向自由表面。对距自由表面为 $d$ 的刃型位错，相应镜像力为：

$$
F_x=-\frac{Gb^2}{4\pi(1-\nu)d}.
\tag{4.45}
$$

对螺型位错，一个异号镜像位错即可满足平面自由表面的边界条件；对刃型位错，除镜像位错外通常还需加入额外的表面修正应力场，才能完整满足边界条件。复杂表面、有限尺寸或多位错体系通常需要数值方法求解。

---

## 4.9 本章小结

- 位错核心外的弹性场可由小变形线弹性理论描述；核心内应采用原子尺度理论。
- 螺型位错产生纯剪切应力，刃型位错同时产生拉压与剪切应力；二者应力均按 $1/r$ 衰减。
- 位错单位长度弹性能近似与 $Gb^2$ 成正比，并对外、内截断半径呈对数依赖。
- 位错在应力场中的一般受力由 Peach–Koehler 公式 $\mathbf f=(\boldsymbol\sigma\cdot\mathbf b)\times\boldsymbol\xi$ 给出；该力是垂直于位错线的单位长度驱动力，可进一步分解为滑移力与攀移力。
- 位错线张力倾向于缩短位错线，外加滑移力可使位错弓出；二者平衡给出 $\tau\approx\alpha Gb/R$。
- 同号平行位错通常排斥，异号平行位错通常吸引；刃型位错因应力场具有方向性，可形成稳定的垂直排列或位错偶极子。
- 攀移同时受机械力与点缺陷化学力控制；自由表面则通过镜像力吸引位错。
