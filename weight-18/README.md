# 权重 18 · $(9;7,6)$ 的三阶小高度比较

取 $Y=4y=4\pi\Im\tau$、$\Delta=\Delta_{18}$，并记 $\Lambda(j)=\Lambda(\Delta,j)$。DGW 的 $(s;s_1,s_2)=(9;7,6)$ 对应 DKS 的偶解 $(s,m,k)=(9,6,7)$。本文从 FKL 的 $(r_1,r_2,d,k)=(13,11,-4,18)$ 算出

$$
F-\check F=\frac{143}{891813888000}\,
 \frac{\Lambda(10)\Lambda(21)}{\langle\Delta,\Delta\rangle}\,H_\Delta^+ .
\tag{1}
$$

证明沿用[权重 12 的小高度方法](../weight-12/s-cocycle-proof.md)。本例需比较的是 $Y^{11}[\check F]_n$ 的三次项。特解取自 DKS 附件 $ \texttt{cFFp[6,7,9]}$，对照值取自 DGW 猜想公式。[DKS 特解附件](https://arxiv.org/src/2109.05017v3/anc/checkFplus.txt) · [FKL 原文](https://arxiv.org/html/2512.21413v1#S1.SS2) · [DGW 一般系数](https://arxiv.org/html/2501.03996v2#S4.SS2)

## 1. FKL 归一化与四阶留数

将 FKL 式 (33) 的正则化结果乘 $\Gamma(21)/\pi$，记为 $W_n$。设 $U\sim\mathrm{Beta}(10,8)$ 及 $\mu_j=\mathbb E(1/2-U)^j$。前五个矩依次为

$$
(\mu_0,\mu_1,\mu_2,\mu_3,\mu_4)=
\left(1,-\frac1{18},\frac{11}{684},-\frac1{456},\frac{13}{19152}\right).
$$

正正频率中 $p+r=n$，取 $p$ 为 $\sigma_{-13}$ 的下标。式 (32) 的整个核是

$$
P_{++}(p,r)=-2(r/n)^{11}P_9^{(11,-13)}(2p/n-1)
-\frac{1140}{n^3}\mathbb E(p-nU)^3.
\tag{2}
$$

混合频率仍写 $x=h+n/2$。$\varepsilon=+1$ 取 $\sigma_{-13}$ 带正频 $h+n$，$\varepsilon=-1$ 取 $\sigma_{-11}$ 带正频。记对应除数权重为 $a_\varepsilon(h)$，$Z_\varepsilon(v)=\sum_{h\ge1}a_\varepsilon(h)x^{-v}$ 的极点留数为

$$
A_n=\operatorname{Res}_{v=1} Z_\varepsilon(v)
=\frac{\zeta(14)\zeta(12)}{\zeta(26)}\sigma_{-25}(n).
$$

按同余类展开与权重 12 一样。沿 FKL 的变形，把混合支写成 $x+\varepsilon n(1/2-U)$ 的 $3+\epsilon$ 次幂的期望。展开至位移四次项。四次项的系数 $\binom{3+\epsilon}{4}\sim\epsilon/4$，正好乘到 $Z_\varepsilon(1-\epsilon;\epsilon)$ 的极点。五次以上余项在 $\epsilon=0$ 逐项消失，并由可求和的 $O_n(x^{-2+\Re\epsilon})$ 控制。于是每支归一化后等于

$$
W_{n,\varepsilon}=-\frac{1140}{n^3}
\left(\sum_{j=0}^3\binom3j(\varepsilon n)^j\mu_j Z_\varepsilon(j-3)
-\frac{n^4\mu_4 A_n}{4}\right).
\tag{3}
$$

使用[校正符号的迹公式](../weight-12/s-cocycle-proof.md)，得到 $W_n=-B_n+T_n$，其中

$$
\begin{aligned}
B_n&=\frac{\Gamma(21)}{\pi}
\left(-\mathcal Z^{13,11}_{-4}\sigma_{-13}(n)+\mathcal Z^{11,13}_{-4}\sigma_{-11}(n)\right),\\
T_n&=-\frac{715}{4096}\frac{\Lambda(8)\Lambda(21)}{\langle\Delta,\Delta\rangle}\frac{a_\Delta(n)}{n^{20}}.
\end{aligned}
\tag{4}
$$

式 (35) 给出 $\mathcal Z$。尖点项的负号来自 $i^{18}=-1$。完成 $L$ 值后，系数为 $-16!/(2^{16}9!7!)$。

## 2. 三阶小高度恒等式

为便于复算，把 DKS 深度二系数写成
$(c_0,\dots,c_4)=(792792,28314,-217074,-47190,-9438/5)$。
第 $j$ 个系数分别乘
$\beta^{\rm sv}\!\left[\begin{smallmatrix}6+j&5-j\\12&14\end{smallmatrix}\right]$
和
$\beta^{\rm sv}\!\left[\begin{smallmatrix}7+j&4-j\\14&12\end{smallmatrix}\right]$。
对任一项，记内外指标为 $a,b$，权重减 2 为 $d_i,d_o$。定义
$b_{j,d,\ell}=[z^\ell]\,(1-z/2)^{d-j}(1+z/2)^j$。
嵌套指数积分 $I_{A,B}$ 沿用[权重 12 定义](../weight-12/s-cocycle-proof.md#两个正频率)。

**正正块。** 提取 $\sigma_{-13}(p)\sigma_{-11}(r)$ 后，每项的三次系数为

$$
\frac{4c_j p^{13}r^{11}}{13!11!}
\sum_{\ell=0}^3 b_{b,d_o,\ell}\,b_{a,d_i,3-\ell}\,
I_{d_o-\ell,d_i-3+\ell}(p_i,p_o).
$$

按级数权重指定 $(p_i,p_o)$。将 $I$ 写成有限阶乘和，合计十项得到 $-143n^3P_{++}/600$。该恒等式可令 $r=n-p$ 并取 $n=1$ 做有理函数化简。左侧对频率联合齐次为三次，右侧 $P_{++}$ 为零次。

**正零块。** 常数模归一化之比为 $691/65520$（$G_{12}$）及 $-1/24$（$G_{14}$）。在上式中改为一个常数模，外零项先积外层，内零项取正则化原函数。每个 $\ell$ 的阶乘和替换成 $20!/(d_o-\ell+1)$ 或 $-20!/(d_i-3+\ell+1)$，再乘对应常数模归一化之比和 $n^{d_{\rm freq}-20}\sigma_{-d_{\rm freq}-1}(n)$。得到

$$
\frac{31916599}{37800}n^{-8}\sigma_{-13}(n)
-\frac{46189}{30}n^{-10}\sigma_{-11}(n).
$$

附件的深度一项用相同 $b_{j,d,3}$ 提取，其贡献为
$ (143/5)n^3\zeta(11)\sigma_{-13}(n)+(1573/30)n^3\zeta(13)\sigma_{-11}(n)$。
两部分合计 $-143n^3B_n/600$。反全纯积分常数与纯零模只给非正频率。

**混合块。** 由 DKS I 式 (7.5)，每个深度二项给权重 $a_\varepsilon(h)e^{-Yx}$，除去此权重的多项式为

$$
\frac{4c_j}{13!11!}
\sum_{u=0}^b\sum_{v=0}^{d_i-a}
\binom bu\binom{d_i-a}{v}(d_o-u)!(d_i-v)!
\bigl(Y(x+n/2)\bigr)^u\bigl(Y(x-n/2)\bigr)^v
$$

其中外层 $G_{14}$ 属于 $\varepsilon=+1$ 支。合计后，带指数因子的三次项对固定 $h$ 等于
$\frac{143\cdot1140}{600}\,\mathbb E(x+\varepsilon n(1/2-U))^3$。
有限求和时还要叠加同余类热核的极点。对每个 $cY^j x^\ell$，只有 $j-\ell-1=3$ 的项参与这一阶，其权重是 $c\,\ell! A_n$。相加得到 $637n^4A_n/640$ 每支。加上 (3) 的正则化修正后，每支的三次系数可以写成

$$
-\frac{143n^3}{600}W_{n,\varepsilon}+\frac{4199}{4032}n^4 A_n .
$$

此处按类展开到 $Y^3$ 的余项可统一控制为 $O_n((YL)^4)$，乘除数权重后可求和，因 $\sum_{a,b}a^{-13}b^{-11}L^4<\infty$。

合并三块，取 $\alpha=\zeta(14)\zeta(12)/(42\zeta(26))$。$\alpha E_{13}$ 在 $Y^{-8}$ 的系数为 $4199n^4A_n/2016$，抵消两支混合块的余项。再用 (4) 得到

$$
[Y^{-8}][\check F-\alpha E_{13}]_n=-\frac{143}{600}n^3 T_n.
\tag{5}
$$

## 3. 补全与完整余循环

在模不变解的大高度 Laurent 多项式中减去同一个 $\alpha E_{13}$，再减适当倍数 $E_2$ 消去 $Y^2$ 项，剩下的增长有界。由模不变性，每个固定 Fourier 模在小高度为 $O(Y^{-\eta})$，任取 $0<\eta<1$。$E_2$ 没有 $Y^{-8}$ 项。零模的齐次项 $y^{-8}$ 由 DGW 式 (4.4) 给出 $-79860781\zeta(21)/1149947084800$，与附件特解一致，且两个解的 Laurent 多项式都无 $y^9$ 项。[DGW Laurent 公式及边界估计](https://arxiv.org/html/2501.03996v2#S4)

因此 $F-\check F$ 的所有非零模由 (5) 的相反数确定，零模消失。除以 $H^+_\Delta$ 的 $Y^{-8}$ 模系数 $(16!/(2\cdot8!))\,a_\Delta(n)/n^{17}$，得到
$-143\Lambda(8)\Lambda(21)/(891813888000\langle\Delta,\Delta\rangle)$。
权重 18 的函数方程给 $\Lambda(8)=-\Lambda(10)$，故为式 (1)。作为核对，DGW 式 (3.50) 给 $\widetilde C_{9;7,6}=302400/143$，代入式 (4.22) 得到同一个系数。

写 $P_{18}^+=\sum_{0\le a,b\le8,\,a+b\ {\rm 偶}}\binom8a\binom8b(-1)^{(a+b)/2}\tau^{8-a}\bar\tau^{8-b}\Lambda(a+b+1)/\Lambda(17)$。
对 DKS 的 $H_\Delta^+$ 作 $z\mapsto-1/z$，其周期积分给 $\delta_S H_\Delta^+=\pi^{17}\Lambda(17)P_{18}^+/(20160y^8)$，等于 DKS II 附录 (B.4a)。所以在本例中 $\delta_S\check F$ 等于这一完整多项式表达乘以式 (1) 补全系数的负值。[DKS II 的权重 18 周期](https://arxiv.org/html/2109.05018v3#App2.SS2)
