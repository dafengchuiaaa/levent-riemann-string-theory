# 权重 12 · 从 FKL 卷积推导 DKS 的完整 $S$-余循环

令 $y=\pi\Im\tau$、$Y=4y$。用 $\check F$ 表示 DKS I 式 (4.18) 的 $\widecheck{\rm F}^{+(6)}_{4,4}$，用 $F$ 表示同一文献给出的模不变 Poincaré 解。本文从 FKL 的超几何正则化和 Petersson 迹公式出发，推出

$$
\delta_S \check F(\tau)
=\check F(\tau)-\check F(-1/\tau)
=-\frac{7\pi^{11}}{621900\,y^5}\Lambda(\Delta_{12},13)p^+_{12}(\tau,\bar\tau).
\tag{1}
$$

推导使用 DKS I 的特解和模不变解的存在性，以及 $\Delta=\Delta_{12}=\sum_{n\ge1}\tau(n)q^n$ 的临界周期比与 Petersson 范数。记 $\Lambda(j)=\Lambda(\Delta,j)$。式 (1) 和 DKS II 式 (3.11) 的补全系数均由下面的计算得出。FKL arXiv v1 的式 (37) 尖点项需要反号，见第 2 节。[DKS I](https://arxiv.org/html/2109.05017v3) · [DKS II](https://arxiv.org/html/2109.05018v3) · [FKL](https://arxiv.org/html/2512.21413v1)

## 1. FKL 的混合频率正则化

固定整数 $n>0$，写
$x_h=h+n/2$、$a(h)=\sigma_{-7}(h)\sigma_{-7}(h+n)$（$h\ge1$）。定义移位的 Dirichlet 级数和它在 $s=1$ 的留数

$$
Z_n(s;\epsilon)=\sum_{h\ge1}\sigma_{-7-\epsilon}(h)\sigma_{-7-\epsilon}(h+n)x_h^{-s},\qquad
A_n=\mathop{\rm Res}_{s=1}Z_n(s;0)=\frac{\zeta(8)^2}{\zeta(16)}\sigma_{-15}(n).
$$

将两项除数函数展开。对除数 $a,b\ge1$，仅在 $g=\gcd(a,b)\mid n$ 时有解。此时 $h\equiv0\bmod a,\ h\equiv-n\bmod b$ 的正解是 $h_0+jL$，其中 $1\le h_0\le L=\operatorname{lcm}(a,b)$。因此 $Z_n$ 延拓为

$$
\sum_{\substack{a,b\ge1\\\gcd(a,b)\mid n}}(ab)^{-7-\epsilon}L^{-s}\zeta\!\left(s,\frac{h_0+n/2}{L}\right).
\tag{2}
$$

将 Hurwitz zeta 分成首项 $u^{-s}$ 与 $\zeta(s,1+u)$，其中 $u=(h_0+n/2)/L\le1+n/2$。在 $s=-1,0,1$ 附近，除 $s=1$ 的简单极点外，(2) 局部正规收敛。留数是对同一组 $(a,b)$ 求和的 $\sum (ab)^{-7-\epsilon}/L$。置 $a=ga',b=gb'$，用互素 Euler 乘积即得 $A_n$。

取 FKL 式 (142) 的变形 $r_1=r_2=7+\epsilon,d=-2-\epsilon$。混合频率 $(n+h,-h)$ 和 $(-h,n+h)$ 的超几何自变量分别在 $(0,1)$ 和负实轴，故没有割线跳跃。Euler beta 积分把两支合为同一个期望

$$
{\mathbb E}\bigl(x_h+n(1/2-U)\bigr)^{1+\epsilon},
\qquad U\sim\mathrm{Beta}(6,6).
$$

具体地，对正支 $p=n+h$，以 $(p/n)^{1+\epsilon}$ 乘超几何函数得到 $n^{-1-\epsilon}\mathbb E(p-nU)^{1+\epsilon}$。负支得到 $n^{-1-\epsilon}\mathbb E(h+nU)^{1+\epsilon}$，由 $U\leftrightarrow1-U$ 两者相同。因为方差为 $1/52$，把期望的前两项展开得

$$
x_h^{1+\epsilon}+\frac{(1+\epsilon)\epsilon n^2}{104}\,x_h^{-1+\epsilon}
+O_n\!\left(x_h^{-3+\Re\epsilon}\right).
$$

在 $\epsilon=0$，余项逐项为零，且由可求和上界控制。这个扣项展开把绝对收敛点（例如 $\epsilon=-3$）延拓到 0。第二项因 $Z_n(1-\epsilon;\epsilon)$ 的极点留下 $-n^2A_n/104$。用式 (32) 的归一化换回式 (36) 的带符号卷积。正正项在式 (36) 中交换两下标后相消，每个混合频率对的极限权重为 $-4x_h/n$。故式 (36) 按 FKL 变形得到的值是

$$
C_n=-\frac4n\left(Z_n(-1;0)-\frac{n^2 A_n}{104}\right).
\tag{3}
$$

## 2. 从 Petersson 迹公式核定 FKL 尖点项的符号

令 $Q_\Sigma=\sum_{n_1+n_2=n,\ n_j\ne0}\sigma_{-r_1}(n_1)\sigma_{-r_2}(n_2)\mathcal Q(n_1,n_2)$，其中 $\mathcal Q$ 取 FKL 式 (8)。先在 FKL 第 3 节的绝对收敛区域计算，再在固定 $k=2d+r_1+r_2+2=12$ 上延拓。将 FKL 的 Bessel–Kloosterman 和 $\mathfrak S$ 除以

$$
B_*=\frac{(2\pi)^{r_2}n^{(r_1+r_2-1)/2}\Gamma(d+1)\Gamma(d+r_1+1)}{\zeta(1+r_2)}.
$$

**迹公式一侧。** FKL 式 (43)、(97)–(99) 的迹项带 $+i^k$，对角项带 $-i^k$。因此

$$
\frac{\mathfrak S}{B_*}
=i^k D_n
-i^k\frac{(2\pi)^{-r_2}\zeta(1+r_2)n^{-r_2}}
{\Gamma(d+1)\Gamma(d+r_1+1)}\sigma_{-r_1}(n),
$$

$$
D_n=\frac{\Gamma(k-1)n^{-d-r_1-r_2}}
{(4\pi)^{k-1}(2\pi)^{r_2}\Gamma(d+1)\Gamma(d+r_1+1)}
\frac{L(1+d+r_2,\Delta)L(1+d+r_1+r_2,\Delta)}{\langle\Delta,\Delta\rangle}.
$$

**Mellin 一侧。** FKL 式 (106) 的两条竖直积分均向上走。从左线移到右线时，左线积分等于右线积分**减**穿越的留数。这修正了式 (109) 所印的加号。代入式 (68)–(69) 后，这部分除以 $B_*$ 等于 $Z_d^{(r_2,r_1)}\sigma_{-r_2}(n)$，用的是 FKL 式 (9) 的 $Z_d^{(\alpha,\beta)}$。

对右线积分，令 $x=n_1/n>0$。在 Lemma 2.4 的 Mellin 积分中用 $w=-1\pm i0$，式 (117) 的 $\cos(\pi s+\pi r_1/2)$ 将两支合成
$ (e^{i\theta}G_-+e^{-i\theta}G_+)/2$，其中 $\theta=\pi(d+1+r_1/2)=\pi(k/2-r_2/2)$。乘上 $\Gamma(d+1)\Gamma(d+r_1+1)x^{-d-1}/(\pi\Gamma(k))$，结果正好是 $\Gamma(d+1)\Gamma(d+r_1+1)\mathcal Q/(2\pi)$。负频支用 $w=1$，得到同一正的比例系数。于是除以 $B_*$ 后，$Q_\Sigma$ 在右线积分中的系数是 $+1$。

右线积分还包含 $n_1=n$。在 $x=1$ 用 Gauss 的超几何求值，并用 $\sigma_{-r_2}(0)=\zeta(r_2)$ 的 Ramanujan 和约定。除以 $B_*$ 后，该项等于

$$
i^k\frac{(2\pi)^{r_2}\zeta(1-r_2)}
{\Gamma(d+r_1+r_2+1)\Gamma(d+r_2+1)}\sigma_{-r_1}(n).
$$

这里用到 $\zeta(1-r)=2(2\pi)^{-r}\Gamma(r)\cos(\pi r/2)\zeta(r)$。把两侧相等并解出 $Q_\Sigma$，得到

$$
Q_\Sigma=-i^k Z_d^{(r_1,r_2)}\sigma_{-r_1}(n)-Z_d^{(r_2,r_1)}\sigma_{-r_2}(n)+i^k D_n.
\tag{4}
$$

这重现了两项边界项和尖点项的大小。迹公式的符号使 (4) 的尖点项为 $+i^kD_n$，与预印本式 (11)、(34) 印出的尖点项符号相反。[FKL 原文，式 (43)、(65)、(97)–(112)、(117)、(124)](https://arxiv.org/html/2512.21413v1#S3)

按 FKL 第 3.6、4 节延拓，乘以 $\Gamma(d+1)$ 再取 $\epsilon\to0$。记

$$
\begin{aligned}
U_n&=\frac13\sum_{p=1}^{n-1}(1-p/n)^7P_5^{(7,-7)}(2p/n-1)\sigma_{-7}(p)\sigma_{-7}(n-p),\\
V_n&=\left(\frac{11}{10n^7}+2\zeta(7)\right)\sigma_{-7}(n),\\
r_n&=\frac{21}{512}\frac{\Lambda(6)\Lambda(13)}{\langle\Delta,\Delta\rangle}\frac{\tau(n)}{n^{12}}.
\end{aligned}
$$

式 (32) 给卷积一侧 $6\pi(U_n-C_n)/\Gamma(13)$。在 (4) 中 $i^k=1$，用 FKL 式 (35) 或 (147) 化边界项，得到

$$
C_n=U_n+V_n-r_n.
\tag{5}
$$

所以保留式 (36) 的带符号核时，式 (37) 中尖点项应取负号。复核时也用 (2) 的 Bernoulli 多项式值和精确的留数公式计算了 $n=1,2,4$。例如 $C_1\approx1.88774094674$，而印刷正号给出约 $4.34565616279$。这项数值核对独立于 DKS 的补全系数。

## 3. DKS 特解的 $Y^{-5}$ Fourier 系数

取关于 $x=\Re\tau$ 的第 $n>0$ 模 $[\,\cdot\,]_n(Y)$。$\check F$ 的深度二部分为 $ (98/3)(75\beta_{4,2}+24\beta_{5,1}+\beta_{6,0})$，记 $\beta_{a,b}=\beta^{\rm sv}\!\left[\begin{smallmatrix}a&b\\8&8\end{smallmatrix}\right]$。用 DKS I 式 (2.24)、(2.29)、(7.5) 直接取模，保留物理 Fourier 系数中的指数因子。下列系数都是 $Y^6[\check F]_n(Y)$ 的 $Y^1$ 项。

### 两个正频率

令 $p,r\ge1,\ p+r=n$，内层频率为 $p$。用 $2\pi$ 缩放虚轴坐标 $u,v$，外层是 $u$，积分区域 $v\ge u\ge Y/2$。在 $\beta_{6-j,j}$ 中，外层多项式为 $(u-Y/2)^{6-j}(u+Y/2)^j$，内层为 $(v-Y/2)^j(v+Y/2)^{6-j}$。下限在 $Y=0$ 的一阶变化为零。用三项权重加总后，外层和内层的一阶因子分别给 $-126$ 和 $126$。

定义
$I_{A,B}(p,r)=\int_0^\infty u^A e^{-ru}\int_u^\infty v^B e^{-pv}\,dv\,du$。则

$$
I_{A,B}(p,r)=\sum_{b=0}^B
\frac{B!(A+B-b)!}{(B-b)!\,p^{b+1}(p+r)^{A+B-b+1}}.
$$

去掉 $\sigma_{-7}(p)\sigma_{-7}(r)$，一阶核为
$336\,p^7r^7(I_{6,5}-I_{5,6})/(6!)^2$。
代入上式并用 $p+r=n$ 化简。这一核与交换 $p,r$ 后的核之和等于

$$
-\frac{28n}{3}\left[(r/n)^7P_5^{(7,-7)}(2p/n-1)
+(p/n)^7P_5^{(7,-7)}(2r/n-1)\right].
$$

因此整个正正块给出 $-28nU_n$。这里的比较对象是完整 $\check F$ 的小高度系数。

### 一个正频率和一个零频率

常数模与正频率归一化之比是
$\rho=2\zeta(8)/(2(2\pi)^8/7!)=1/480$。
对外层为零的项先对外层积分，对内层为零的项用幂函数在无穷远的正则化原函数。正频率外层乘内层零模时，式 (2.24) 中反全纯腿相减，消去同一多项式原函数的端点常数。剩下的有限端点在 $Y/2$ 与 $-Y/2$，只涉及 $Y^7$ 及更高次。两种次序各给
$-126/6+126/7=-3$，共同乘 $ (8/3)\rho\,12!/(6!)^2$ 及 $n^{-6}\sigma_{-7}(n)$。于是双层积分的单零模贡献为 $-154n^{-6}\sigma_{-7}(n)/5$。

式 (4.18) 的深度一项乘 $Y^6$ 后权重依次为 $(-28/3,-224,-700)\zeta(7)$。用式 (2.29) 的正频项取一阶，得到
$(-2\cdot5!/7!)\sum_{j=0}^2(j-3)(-28/3,-224,-700)_j=-56$，乘 $n\zeta(7)\sigma_{-7}(n)$。纯零频项和反全纯积分常数不给正 Fourier 模。合计正零块为 $-28nV_n$。

### 一正一负的频率尾

写 $p=h+n,r=h,x_h=(p+r)/2$。由 DKS I 式 (7.5)，乘 $Y^6$ 后的这一块为

$$
\sum_{h\ge1}a(h)e^{-Yx_h}K(Y,x_h,n),
$$

$$
\begin{aligned}
K(Y,x,n)={}&\frac{800}{3}+\frac{464}{3}xY+
 \left(\frac{112}{3}x^2-\frac83 n^2\right)Y^2\\
 &+\left(\frac{40}{9}x^3-\frac{10}{9}n^2x\right)Y^3+
 \left(\frac29x^4-\frac19 n^2x^2+\frac{n^4}{72}\right)Y^4 .
\end{aligned}
$$

例如上式可由 $K=(8/3)\sum_{j=0}^2 c_j\sum_{a,b=0}^j\binom ja\binom jb(6-a)!(6-b)!(Yp)^a(Yr)^b/(6!)^2$ 展开，其中 $(c_0,c_1,c_2)=(1,24,75)$。

展开时先按 (2) 的同余类求和。在每一类中，生成函数为
$e^{-(h_0+n/2)Y}/(1-e^{-LY})$。它的 Laurent 展开给出
$\sum a(h)e^{-Yx_h}=A_n/Y+Z_n(0;0)-Z_n(-1;0)Y+O_n(Y^2)$。
对于 $Y^m\sum x_h^m a(h)e^{-Yx_h}$（$0\le m\le4$），逐类使用 $(-1)^mY^m\partial_Y^m$。在 $LY\le1$ 时截断误差由 $O_n((LY)^2)$ 控制，在 $LY\ge1$ 时同一上界由 $z^m e^{-z}$ 的有界性得到。因为 $\sum (ab)^{-7}L^2$ 收敛，可将展开逐类相加。

该块的 $Y^{-1}$ 项为 $528 A_n$。它的 $Y^1$ 项为

$$
-112Z_n(-1;0)-4n^2 A_n
=28n C_n-\frac{66}{13}n^2 A_n.
\tag{6}
$$

前一个等式右侧的 $-4n^2 A_n$ 来自 $K$ 中 $n^2Y^2,-n^2 xY^3,-n^2x^2Y^4$ 三项对极点的贡献。后一个等式用了 (3) 的有限留数修正。逐项只取 $h$ 固定时的 $Y^1$ 项会漏掉极点贡献。[DKS I 式 (2.24)、(4.18)、(7.5)](https://arxiv.org/html/2109.05017v3#S7.SS1)

## 4. 从边界系数到模不变补全

置 $\alpha=\zeta(8)^2/(13\zeta(16))$。$F-\alpha E_8$ 的 $Y^8$ 项抵消，在尖点为 $O(Y)$。由模不变性和与 $E_{1+\eta}$ 的比较，它的**每个固定 Fourier 模**在 $Y\to0^+$ 时为 $O(Y^{-\eta})$，任取 $0<\eta<1$。[DGW 第 4.2 节的小高度估计](https://arxiv.org/html/2501.03996v2#S4.SS2)

对第 $n>0$ 模，$E_8$ 的 $Y^{-5}$ 系数是 $-66 n^2\sigma_{-15}(n)$，可从 DKS I 式 (2.20) 的半整数阶 Bessel 有限展开取出。乘 $\alpha$ 后恰好抵消 (6) 的第二项。于是

$$
[Y^{-5}]\,[\check F-\alpha E_8]_n
=28n(C_n-U_n-V_n)=-28n r_n .
$$

$F-\check F$ 解齐次方程 $(\Delta_\tau-30)f=0$。每个非零模在无穷高处取衰减的 $K_{11/2}$ 支，这个支的 $Y^{-5}$ 系数唯一确定整模。由上式和小高度估计，$F-\check F$ 的该系数为 $28n r_n$。负模由偶性确定。零模只可能含 $Y^6,Y^{-5}$，两者因 DKS I 附录 D.1 给出的相同 Laurent 多项式而消失。

DKS II 式 (2.33) 的 $H_\Delta^+$ 的第 $n$ 模在 $Y^{-5}$ 处给
$\Gamma(11)\tau(n)/(2\cdot5!\,n^{11})=15120\tau(n)/n^{11}$。
因此恢复整个齐次项，得到

$$
F-\check F=\frac{r_1}{540}H_\Delta^+
=\frac{7}{10365}\frac{\Lambda(13)}{\Lambda(11)}H_\Delta^+ .
\tag{7}
$$

第二个等式用 $\Lambda(6)/\Lambda(8)=4/5$ 与 DGW 的 $\langle\Delta,\Delta\rangle=(691/7680)\Lambda(8)\Lambda(11)$。这里 $r_1>0$ 是 (5) 中被减去的尖点项系数大小。把 FKL 有符号的剩余项提升到用 DKS 核写出的补全项，实际乘数是 $-1/540$，其大小来自 $28/15120$。[DKS II 的齐次积分](https://arxiv.org/html/2109.05018v3#S2.SS2) · [DGW 的 Petersson 范数](https://arxiv.org/html/2501.03996v2#S4.SS2)

## 5. 计算整个余循环，并解出多重模值

对式 (2.31) 的 $H_\Delta^+$ 作 $z\mapsto-1/z$。六次特征值所对应的积分核为 $(\tau-z)^5(\bar\tau-z)^5$。$\Delta(-1/z)=z^{12}\Delta(z)$ 与链式求导给出所需的抵消，变换后积分路径为 $\tau\to0$。所以相减得到周期积分

$$
\delta_S H_\Delta^+
=\frac{\pi^{11}}{60y^5}\operatorname{Re}
\left[i\int_0^{i\infty}(\tau-z)^5(\bar\tau-z)^5\Delta(z)\,dz\right].
$$

展开两个五次幂，$\int_0^{i\infty}z^j\Delta(z)dz=i^{j+1}\Lambda(j+1)$。取实部后只保留偶数 $j$。用 DKS II 式 (2.45) 的奇数临界周期比及反射公式，得到 $\delta_S H_\Delta^+=\pi^{11}\Lambda(11)p_{12}^+/(60y^5)$。这里完整多项式也可由下列有限和生成，$\eta_j=\Lambda(j+1)/\Lambda(11)$：

$$
p_{12}^+(\tau,\bar\tau)=
\sum_{\substack{0\le a,b\le5\\a+b\ \mathrm{偶}}}
\binom5a\binom5b(-1)^{1+(a+b)/2}\eta_{a+b}\tau^{5-a}\bar\tau^{5-b},
\quad
(\eta_0,\eta_2,\eta_4,\eta_6,\eta_8,\eta_{10})
=\left(1,\frac{691}{1620},\frac{691}{2520},\frac{691}{2520},\frac{691}{1620},1\right).
$$

将 (7) 代入 $\delta_S\check F=-\delta_S(F-\check F)$，得到式 (1) 的符号、标量与完整多项式。

要提取 $M_{01}$，使用 DKS II 式 (3.9) 中 $(m,k)=(4,4)$ 那一组尚未求值的多重模值表达，将其常数单项式系数与式 (1) 对照。以 $(\tau-\bar\tau)^{-5}$ 为公因子，该常数给

$$
M_{01}
=-\frac{4i\pi^9}{416745}\zeta(7)
+\frac{2860i\pi^3}{304731}\zeta(13)
-\frac{256i\pi^{14}}{22854825}\Lambda(\Delta_{12},13).
$$

此处推导顺序是先从正则化卷积得到完整余循环，再解出多重模值。上篇的[双重积分数值复核](dks-multiple-modular-value.md)从定义出发，可独立检查这三个常数合成后的数值。一般权重的统一转换公式仍待建立。
