# 权重 12 · 非对角算例 $(6;5,3$\)

令 $Y=4y=4\pi\Im\tau\$，记 $\Delta=\Delta_{12}\$、$\Lambda(j$=\Lambda(\Delta,j)\)。DGW 的 $(s;s_1,s_2$=(6;5,3)\) 对应 DKS 的 $(s,m,k$=(6,3,5)\)。沿用[对角算例](s-cocycle-proof.md)建立的 Fourier 与齐次解方法，直接计算得

$$
F-\check F=-\frac{5}{11056}\frac{\Lambda(13)}{\Lambda(11)}H^+_\Delta,\qquad
\delta_S\check F=\frac{\pi^{11}\Lambda(13)}{132672\,y^5}p_{12}^+(\tau,\bar\tau).
\tag{1}
$$

以下计算采用 DKS I 式 (4.22) 的特解和附录 D 的 Laurent 多项式。式 (1) 的系数由 FKL 卷积与小高度展开确定，最后与 DKS II 式 (3.9)、(3.11) 对照。[DKS I](https://arxiv.org/html/2109.05017v3) · [DKS II](https://arxiv.org/html/2109.05018v3#S3.SS1)

## 1. FKL 卷积与两个混合频率方向

FKL 参数是 $(r_1,r_2,d,k$=(9,5,-2,12)\)。将其式 (33) 的正则化值乘 $\Gamma(13$/\pi\)，记作 $W_n\$。式 (32) 给正正频率 $p,r>0,\ p+r=n\$ 的归一化核

$$
P_{++}(p,r)=-2(r/n)^5P_7^{(5,-9)}(2p/n-1)-12p/n+8,
\tag{2}
$$

其权重是 $\sigma_{-9}(p$\sigma_{-5}(r)\)。与对角情形不同，右边最后两个单项在正正频率求和中保留。

对两支混合频率，用 $h\ge1\$、$x_h=h+n/2\$。以 $\varepsilon=-1\$ 表示 $\sigma_{-9}\$ 带正频 $h+n\$，以 $\varepsilon=+1\$ 表示 $\sigma_{-5}\$ 带正频 $h+n\$。相应除数权重分别为
$a_-(h$=\sigma_{-9}(h+n)\sigma_{-5}(h)\) 与
$a_+(h$=\sigma_{-5}(h+n)\sigma_{-9}(h)\)。
令 $Z_\varepsilon(s$\) 为 $\sum_{h\ge1}a_\varepsilon(h$x_h^{-s}\) 的延拓。展开除数条件后的 Hurwitz zeta 级数与对角算例同理。两支 $s=1\$ 留数相同，均为

$$
A_n=\frac{\zeta(10)\zeta(6)}{\zeta(16)}\sigma_{-15}(n).
$$

沿 FKL 的共同指数变形 $r_i\mapsto r_i+\epsilon,\ d\mapsto d-\epsilon\$，Euler 积分用 $U\sim\mathrm{Beta}(8,4$\)。两支的位移为 $\varepsilon n(U-1/2$\)，一阶矩是 $\varepsilon n/6\$，二阶矩是 $7n^2/156\$。在归一化核趋于 $-12(x_h+\varepsilon n/6$/n\) 时，二阶矩项的因子 $\epsilon(1+\epsilon$/2\) 与 $Z_\varepsilon(1-\epsilon;\epsilon$\) 相乘，留下有限项 $-7n^2 A_n/312\$。此处带 $\epsilon\$ 的 $Z_\varepsilon\$ 同时变形两个除数指数。三阶及更高余项可用 $O_n(x_h^{-2+\Re\epsilon}$\) 支配，且逐项在 $\epsilon=0\$ 消失。因而每支的正则化值是

$$
W_{n,\varepsilon}=-\frac{12}{n}
\left(Z_\varepsilon(-1)+\frac{\varepsilon n}{6}Z_\varepsilon(0)
-\frac{7n^2 A_n}{312}\right).
\tag{3}
$$

合并 (2)、(3)，以[已校正方向的 Petersson 迹公式](s-cocycle-proof.md#2-从-petersson-迹公式核定-fkl-尖点项的符号)评价 $W_n\$，得到

$$
\begin{aligned}
W_n&=-B_n+T_n,\\
B_n&=\left(\frac{22}{7n^5}-4\zeta(5)\right)\sigma_{-9}(n)
    +\left(\frac{5}{3n^9}-8\zeta(9)\right)\sigma_{-5}(n),\\
T_n&=\frac{15}{128}\frac{\Lambda(4)\Lambda(13)}{\langle\Delta,\Delta\rangle}\frac{\tau(n)}{n^{12}}.
\end{aligned}
\tag{4}
$$

边界项由 FKL 式 (35) 与 zeta 函数方程给出。尖点项常数在完成两个 $L\$ 值后为 $10!/(2^{10}7!3!$=15/128\)。留数和边界的符号可另用 (2)、(3) 的同余类展开核对。[FKL 原文](https://arxiv.org/html/2512.21413v1#S1.SS2)

## 2. DKS 第 $n>0\$ 模的小高度展开

深度二的四项用 $(w,j_i,j_o,d_i,d_o$\) 记，其中 $i,o\$ 表示内、外，$d_i,d_o\$ 等于相应 Eisenstein 权重减 2：

$$
(1260,3,3,4,8),\ (315,4,2,4,8),\
(1260,5,1,8,4),\ (315,6,0,8,4).
$$

下面计算 $Y^6[\check F]_n\$ 的 $Y^1\$ 系数，保留 $\exp(-2\pi n\Im\tau$\) 等指数因子。

**正正块。** 令 $I_{A,B}(p_i,p_o$\) 为[对角证明](s-cocycle-proof.md#两个正频率)的嵌套指数积分。固定对应 $G_{10}\$ 的频率 $p\$ 和对应 $G_6\$ 的频率 $r\$。提取 $\sigma_{-9}(p$\sigma_{-5}(r)\) 后，每个深度二项贡献

$$
\frac{4w p^9 r^5}{9!5!}\left[
(j_o-d_o/2)I_{d_o-1,d_i}(p_i,p_o)
+(j_i-d_i/2)I_{d_o,d_i-1}(p_i,p_o)
\right].
$$

两种次序用各自的 $(p_i,p_o$=(p,r)\) 或 $(r,p$\)。代入有限阶乘和后，四项总和逐项等于 $(21n/4$P_{++}(p,r)\)。这里用的是 (2) 的整个核。

**正零块。** $G_6,G_{10}\$ 的正频率系数均含 $i^k=-1\$，所以常数模与正频率归一化之比分别为 $-1/504\$ 和 $-1/264\$。对内层常数模取正则化原函数，对外层常数模先积外层变量。若 $\nu_i=j_i-d_i/2,\ \nu_o=j_o-d_o/2\$，给对应内层非零模的有理因子是 $\nu_o/d_o+\nu_i/(d_o+1$\)，给对应外层非零模的因子是 $-\nu_o/(d_i+1$-\nu_i/d_i\)。再乘 $4w\,12!/(9!5!$\) 和另一个级数的常数模比例。四项总和给
$ (33/2$n^{-4}\sigma_{-9}(n)+(35/4)n^{-8}\sigma_{-5}(n)\)。
DKS I 式 (4.22) 的深度一项给
$-21n\zeta(5$\sigma_{-9}(n)-42n\zeta(9)\sigma_{-5}(n)\)。
合计为 $(21n/4$B_n\)。

**混合块。** DKS I 式 (7.5) 的一支为 $\sum_{h\ge1} a_-(h$e^{-Yx_h} K_-(Y,x_h,n)\)，另一支将 $n\$ 在多项式中换成 $-n\$，同时用 $a_+\$ 取代 $a_-\$。其显式核为

$$
\begin{aligned}
K_-={}&140+(77x+21n/2)Y+(17x^2+13nx/2-n^2)Y^2\\
&+(11x^3/6+5nx^2/4-n^2x/8-7n^3/48)Y^3\\
&+(x^4/12+nx^3/12-n^3x/48-n^4/192)Y^4.
\end{aligned}
$$

固定 $h\$ 时，$e^{-Yx_h}K_\varepsilon\$ 的一阶项为 $-63(x_h+\varepsilon n/6$\)。按同余类先求和后展开，还要计入极点 $A_n/Y\$。对核中的 $cY^j x^\ell\$ 使用 $cY^j(-1$^\ell\partial_Y^\ell\)。满足 $j-\ell-1=1\$ 的各项合计贡献 $-9n^2 A_n/8\$，两支相同。因此第 $\varepsilon\$ 支在 $Y^1\$ 的系数为

$$
-63\left(Z_\varepsilon(-1)+\frac{\varepsilon n}{6}Z_\varepsilon(0)\right)-\frac98 n^2 A_n
=\frac{21n}{4}W_{n,\varepsilon}-\frac{33}{13}n^2 A_n.
$$

其中 $9/8+63\cdot7/312=33/13\$。有限余项的可求和估计沿用同余类热核的 $O_n((LY$^2)\) 估计，权重改为 $a^{-9}b^{-5}\$ 或其交换。[DKS I 的 Fourier 展开](https://arxiv.org/html/2109.05017v3#S7.SS1)

## 3. 从全部模到 $S\$-余循环和多重模值

合并三块，再减去 $\alpha E_8\$，其中 $\alpha=\zeta(10$\zeta(6)/(13\zeta(16))\)。$E_8\$ 的第 $n\$ 模乘 $\alpha\$ 后，$Y^{-5}\$ 项恰为 $-66n^2 A_n/13\$，抵消两支混合块留下的量。用 (4) 消去边界项，得到

$$
[Y^{-5}][\check F-\alpha E_8]_n=\frac{21n}{4}\,T_n.
\tag{5}
$$

此例的模不变解还含 $Y^3\$ 增长项。再减适当倍数的 $E_3\$，模不变解在尖点的增长至多为 $O(Y$\)，故每个 Fourier 模在小高度满足 $O(Y^{-\eta}$\)，任取 $0<\eta<1\$。$E_3\$ 本身的 $n\ne0\$ 模最低只有 $Y^{-2}\$，所以它不改动 (5)。由此 $F-\check F\$ 的 $Y^{-5}\$ 项为 $-21nT_n/4\$。齐次方程的衰减支确定全部非零模，零模仍由两解相同的 Laurent 多项式消去。[DGW 小高度估计](https://arxiv.org/html/2501.03996v2#S4.SS2)

用对角算例已算出的 $H^+_\Delta\$ 系数 $15120\,\tau(n$/n^{11}\)，补全的乘数是

$$
-\frac{21}{4}\frac{15}{128\cdot15120}
 \frac{\Lambda(4)\Lambda(13)}{\langle\Delta,\Delta\rangle}
=-\frac1{24576}\frac{\Lambda(4)\Lambda(13)}{\langle\Delta,\Delta\rangle}
=-\frac5{11056}\frac{\Lambda(13)}{\Lambda(11)}.
$$

最后一步用了 $\Lambda(4$=\Lambda(8)\) 及 $\langle\Delta,\Delta\rangle=(691/7680$\Lambda(8)\Lambda(11)\)。同一个 $H^+_\Delta\$ 的周期给出 (1) 的整个多项式。这也说明在 $W_n\$ 的归一化下，从尖点项的 Fourier 序列到齐次积分的额外乘数是 **$-1/2880\$**。

从 (1) 提取 DKS II 式 (3.9) 的未求值多重模值组合。结果为

$$
\begin{aligned}
\mathcal M\!\begin{bmatrix}0&1\\6&10\end{bmatrix}
-\mathcal M\!\begin{bmatrix}1&0\\6&10\end{bmatrix}
={}&-\frac{i\pi^{11}}{1275750}\zeta(5)-\frac{i\pi^7}{6075}\zeta(9)
 +\frac{1001i\pi^3}{41460}\zeta(13)
 +\frac{32i\pi^{14}}{1958985}\Lambda(13).
\end{aligned}
$$

按 DKS 的 shuffle 与反射公式从深度二积分独立数值复核。将两支以 $i\$ 为分点，每个 $G_6,G_{10}\$ 保留 22 个正频率，得到该组合除以 $i\$ 为 $3.00368571138371245627564168885\$，与上式在 55 位运算下的实际数值差约 $8\times10^{-52}\$。这是有限截断的核对。证明使用的是上述精确核恒等式与解析延拓。[DKS 的多重模值定义](https://arxiv.org/html/2109.05018v3#S2.SS1) · [DKS 的结果](https://arxiv.org/html/2109.05018v3#S3.SS1)
