# 整数指标谱的卷积、齐次补全与 S-余循环

**结论。** 对 Dorigoni–Green–Wen（DGW）式 (4.13) 的整数指标谱，以 Dorigoni–Kleinschmidt–Schlotterer（DKS）的 $\widecheck F$ 为特解，FKL 正则化卷积产生的尖点项给出的补全系数，等于 DGW 猜想 1 的系数。下面证明混合、正正、正零三个 Fourier 块的统一恒等式。证明稿尚待独立审阅。

沿用 Neon 提出的 Green–Parseval 路线。三 Bessel 积分的 Gauss 方程在下文直接推出。混合尾先减去连续密度，再用局部 Green 恒等式处理连续部分。两个有限块通过反演同一个双变量多项式计算。

[DKS I](https://arxiv.org/html/2109.05017v3) · [DKS II](https://arxiv.org/html/2109.05018v3) · [FKL](https://arxiv.org/html/2512.21413v1) · [DGW](https://arxiv.org/html/2501.03996v2)

## 1. 归一化与目标

设 $m=s_2\le s_1=b$，$w=m+b$，$\delta=b-m$，$s\in\{\delta+2,\delta+4,\ldots,w-2\}$。用

$$
Y=4\pi\Im\tau,\quad p=(w-s)/2,\quad \ell=2p-1,\quad
k=2s,\quad d=-2p,\quad q=w-2,\quad
A=s+\delta,\quad B=s-\delta .
$$

下文 $q$ 只表示整数 $w-2$，第 $n$ 模指 $e^{2\pi i n\Re\tau}$ 前的系数。FKL 的两个除数指数为 $-r_1,-r_2$，其中 $r_1=2b-1=A+\ell,\ r_2=2m-1=B+\ell$。记 $(a)^{\overline j}=\Gamma(a+j)/\Gamma(a)$，$(a)_{\underline j}=\Gamma(a+1)/\Gamma(a-j+1)$。取

$$
\Xi=(-1)^{s+m}\frac{(k)^{\overline\ell}}{\ell!},\qquad
c=(-1)^{m+s+p}\frac{(p-1)!(2s-2)!}
{(s-1)!(s+p-1)!(b-1)_{\underline p}(m-1)_{\underline p}},
\qquad C=c\Xi .
$$

$W_n$ 为 FKL 式 (33) 的值乘 $\Gamma(k+\ell)/\pi$，$T_n$ 为其中的尖点项。采用[校正后的 FKL 迹公式符号](../weight-12/s-cocycle-proof.md#2-从-petersson-迹公式核定-fkl-尖点项的符号)。设 $D_n$ 为 $W_n$ 两项除数和的相反数，所以 $W_n=-D_n+T_n$。定义 $\alpha=2\zeta(2b)\zeta(2m)/([w(w-1)-s(s-1)]\zeta(2w))$。要证

$$
[Y^{1-s}][\widecheck F-\alpha E_w]_n
=c\,n^\ell T_n,\qquad n>0 .
\tag{1}
$$

第 2–6 节的核计算对 $k\ge4$ 都成立。第 7 节只在 $k\ge6$ 时调用 FKL 迹公式。$k=4$ 的零补全由 DKS I §4.4 及 $S_4=0$ 给出。

## 2. 标记混合频率与 DKS 的系数空间

将两个 Eisenstein 因子标记为权 $2b$ 和权 $2m$，即使两权相等也保留标记。DKS I 式 (2.35)、(3.24)–(3.25) 的深度二组合使用两种次序，指标分别为

$$
(a_j,b_j)=(2m-2-j,b-m+j),\qquad
(2b-2-b_j,\,2m-2-a_j),\qquad 0\le j\le m-2 .
$$

两次序的系数相同，记为 $v_j$。为避免以 DKS 的谱数值观察作为前提，以下直接核定该谱。设 $N=2m-2$。DKS 式 (3.11) 的完整矩阵作用于 $j=0,\ldots,N$ 上的函数，可写成

$$
q(q+1)f(j)+a(j)[f(j-1)-f(j)]+b(j)[f(j+1)-f(j)],\quad
a(j)=j(b-m+j),\quad b(j)=(N-j)(q-j).
$$

端点的越界系数为零。它在次数不超过 $h$ 的多项式上三角化，在顶次 $j^h$ 上的系数是 $q(q+1)-2qh+h(h-1)=(q-h)(q-h+1)$。反射 $j\mapsto N-j$ 交换两差分，故反射奇子空间恰对应左上 $m-1$ 阶矩阵（中点值为零）。在该子空间取 $h=1,3,\ldots,2m-3$，给出全部 $s(s-1)$。相应左特征向量的末分量非零，因为三对角相邻系数均非零。故可按 DKS 式 (3.10) 将源系数归一为 1。

对权 $2b$ 的正频 $P>0$ 和权 $2m$ 的负频 $-R<0$，记除去两项除数函数的特解为 $g(Y;P,R)$。以 $u,v$ 为两条竖直积分的无量纲高度，分别对应权 $2b,2m$。DKS I 式 (2.24)、(7.5) 给

$$
g(Y;P,R)=Y^{-q} P^{r_1}R^{r_2}\!\!
\int_{Y/2}^\infty\!\int_{Y/2}^\infty e^{-Pu-Rv}\mathscr B(Y,u,v)\,dv\,du ,
\tag{2}
$$

$$
\mathscr B=\frac4{r_1!r_2!}\sum_{j=0}^{m-2}v_j
 (u-Y/2)^{r_1-1-b_j}(u+Y/2)^{b_j}
 (v+Y/2)^{r_2-1-a_j}(v-Y/2)^{a_j}.
$$

交换正负频率时，第二种指标次序给出同一个 $\mathscr B$。同权时两条标记支仍各对应一次积分，最后相加。写

$$
F_t(f,Y)=\frac{2}{\Gamma(t)} f^{t-1/2}\sqrt{Y/\pi}K_{t-1/2}(fY/2),\quad
\mathcal L_{n,s}=Y^2\partial_Y^2-\frac{n^2Y^2}{4}-s(s-1).
$$

单个标记对已经满足 $\mathcal L_{|P-R|,s}g=F_b(P,Y)F_m(R,Y)$。这个断言可在积分核上核验。对混合腿，设 $D=\tau-\bar\tau$，将两个多项式核分别记为 $U_a,V_{b'}$。有 $D\partial_\tau U_a=(2m-2-a)(U_a-U_{a+1})$，$D\partial_{\bar\tau}U_a=a(U_{a-1}-U_a)$，大权腿同理。作用于 $D^{-q}U_aV_{b'}$ 时 $a+b'=q$。先作 $\partial_\tau$，两项只升一个下标，再作 $-D^2\partial_{\bar\tau}$，给出 DKS 式 (3.11) 的三个系数。原范围内 $a\ge m$、$2b-2-b'\ge b$，导数没有积分端点项。到 $j=m-1$ 的出界项是两个对称幂积分的乘积。由 $K$ 的实积分表示，

$$
F_t(f,Y)=\frac{2f^{2t-1}Y^{1-t}}{\Gamma(t)^2}
\int_{Y/2}^\infty e^{-fu}(u^2-Y^2/4)^{t-1}\,du .
$$

故出界项的源系数为 $b m v_{m-2}\Gamma(b)^2\Gamma(m)^2/(r_1!r_2!)$，与 DKS 的源归一化一致。这个推导对任意正实 $P,R$ 有效。

把 (2) 的积分按有限阶乘和计算，可写成 $Y^{-q}e^{-(P+R)Y/2}$ 乘 $YP,YR$ 的多项式，其次数至多 $b-m+2(m-2)=q-2$。这个上界也控制下文对混合频率作连续积分时的无穷远端点。

## 3. 单对 Green 积分的 Gauss 方程

令 $n=|P-R|>0$，$\psi=\sqrt{nY}I_{s-1/2}(nY/2)$，$w_{n,s}=(2s-1)n^s/[2^{2s-1}\Gamma(s+1/2)]$。写 $\mathscr W(\psi,g)=\psi g'-\psi'g$。对一对严格混合频率，$\mathscr W(\psi,g)(\infty)=0$。由 Green 恒等式，

$$
\int_\eta^\infty\psi F_b(P,Y)F_m(R,Y)Y^{-2}\,dY
=-\mathscr W(\psi,g)(\eta).
$$

常数项一般还含 $g$ 的 $Y^{1-s-2j}$ 项。对固定的一对频率，这些低阶项全为零：源的两条奇异支的乘积只含 $Y^{2-w+2h}$，其奇偶性与 $1-s$ 相反；其余支从不低于 $Y^{1-\delta}$，且 $1-\delta\ge3-s$。由 $\mathcal L_{n,s}$ 的 Laurent 递推，低于 $1-s$ 的同奇偶系数逐项为零。因此有限部除以 $w_{n,s}$ 给 $g$ 的 $Y^{1-s}$ 系数。

为计算积分，变 $m_\epsilon=m+\epsilon/2,\ b_\epsilon=b+\epsilon/2$，取实收敛域 $m_\epsilon,b_\epsilon>1/2,\ w+\epsilon<s$，随后作亚纯延拓。令 $\nu=s-\tfrac12,\ \mu=m_\epsilon-\tfrac12,\ \rho=b_\epsilon-\tfrac12$。考虑

$$
J(a,u,v)=a^{-\nu}u^\mu v^\rho\int_0^\infty
t^{-1/2}I_\nu(at)K_\mu(ut)K_\rho(vt)\,dt .
$$

它对三参数联合齐次 $h=\mu+\rho-\nu-\tfrac12=w+\epsilon-s-1$。各 Bessel 方程给

$$
(\partial_a^2+(2\nu+1)a^{-1}\partial_a)J
=(\partial_u^2+(1-2\mu)u^{-1}\partial_u)J
=(\partial_v^2+(1-2\rho)v^{-1}\partial_v)J.
$$

在共线 $v=a+u$ 处取 $x=a/v,\ z=u/v=1-x$。写 $J=v^h f(x,z)$，$D=x\partial_x+z\partial_z$。第三个算子化为 $(h-D)(h-D-2\rho)$，前两个分别化为 $f_{xx}+(2\nu+1)f_x/x$、$f_{zz}+(1-2\mu)f_z/z$，三者都有共同因子 $v^{h-2}$。在 $x+z=1$，用 $D^2 f=x f_{xx}+z f_{zz}-xz(d/dx)^2 f(x,1-x)+Df$ 消去法向导数。函数 $j(x)=f(x,1-x)$ 满足

$$
x(1-x)j''+[2s-(2\rho-2h+1)x]j'-h(h-2\rho)j=0.
$$

因此其在 $x=0$ 的解析支是 ${}_2F_1(-h,2\rho-h;2s;x)$。归一化后的 $I_\nu$ 在此解析，小邻域的三 Bessel 积分由绝对收敛控制。原点常数由二 $K$ 矩给出。具体地，用两次 Schwinger 积分，先积 $t$，再把两个 Schwinger 变量换为和与比，得到
\(\int_0^\infty t^{s-1}K_\mu(t)K_\rho(t)dt=
2^{s-3}\prod_{\pm,\pm}\Gamma((s\pm\mu\pm\rho)/2)/\Gamma(s)\)。
恢复 $Y/2$ 的变量及 $F_t,w_{n,s}$ 的因子，得到

$$
\frac1{w_{n,s}}\int_0^\infty \psi F_{m_\epsilon}(R,Y)F_{b_\epsilon}(P,Y)Y^{-2}\,dY
=C_\epsilon P^h {}_2F_1(-h,A;k;n/P)
\quad(P>R),
\tag{3}
$$

其中

$$
C_\epsilon=\frac{2^{2s-1}}{\pi(2s-1)\Gamma(s)\Gamma(m_\epsilon)\Gamma(b_\epsilon)}
\prod_{\pm,\pm}\Gamma((s\pm\mu\pm\rho)/2).
$$

当 $R>P$ 时交换小、大参数的角色及相应阶，或用 Pfaff 变换。对两支而言，Euler beta 积分即 $\mathbb E[(1-U)P+UR]^h$，其中 $U\sim\mathrm{Beta}(A,B)$。这个权在变形下保持相同。

在整数点，小高度源的在动幂次中，奇异支乘积的 Green 被积函数幂为 $-2p-\epsilon+2j$，永不在 $\epsilon=0$ 成为 $-1$。其他三支已经可积。对固定 $P,R>0$，按两支 Bessel 级数逐项扣去有限个发散幂，(3) 在 $\epsilon=0$ 的值就是 Green 截断积分的有限部。Gamma 倍元和反射公式给 $C_0=C=c\Xi$。因此

$$
[Y^{1-s}]g(Y;P,R)=C\,\mathbb E[(1-U)P+UR]^\ell .
\tag{4}
$$

## 4. 混合尾的连续密度与有限部

固定 $n>0$，写 $x=h+n/2,\ x_0=n/2$。混合两支用 $P=h+n,R=h$ 或 $P=h,R=h+n$，相应权重是 $a_\pm(h)=\sigma_{-r_1}(P)\sigma_{-r_2}(R)$。记 $Z_\pm(v;\epsilon)=\sum_{h\ge1}a_{\pm,\epsilon}(h)x^{-v}$。两支在 $v=1,\epsilon=0$ 的留数都为

$$
A_n=\frac{\zeta(2b)\zeta(2m)}{\zeta(2w)}\sigma_{1-2w}(n).
$$

将除数条件展开为同余类 $h=h_0+jL$，$L=\operatorname{lcm}(a,b')$，$\gcd(a,b')\mid n$。用相应 Hurwitz zeta 级数延拓 $Z_\pm$，再取留数，便得到上式。记
$Q_\pm(x)=\mathbb E[x\pm n(1/2-U)]^\ell=\sum_{j=0}^\ell q_{\pm,j}x^j$，
$\mu_{\ell+1}=\mathbb E(1/2-U)^{\ell+1}$。
FKL 变形的两支混合核按 Euler 积分写为 $\mathbb E[x\pm n(1/2-U)]^{\ell+\epsilon}$，极限公共因子为 $\Xi/n^\ell$。展开到位移的 $\ell+1$ 阶。最后一项的系数
$\binom{\ell+\epsilon}{\ell+1}\sim\epsilon/(\ell+1)$
乘 $Z_\pm(1-\epsilon;\epsilon)$，留下一项 $-n^{\ell+1}\mu_{\ell+1}A_n/(\ell+1)$。剩余项在零点附近以 $O_n(x^{-2+\Re\epsilon})$ 乘除数权重控制，逐点趋零。于是

$$
W_{n,\mathrm{mix}}=\frac{\Xi}{n^\ell}\sum_\pm
\left(\sum_{j=0}^\ell q_{\pm,j} Z_\pm(-j;0)
-\frac{n^{\ell+1}\mu_{\ell+1} A_n}{\ell+1}\right).
\tag{5}
$$

将 DKS 离散和 $G_{\rm mix}$ 分解为连续密度
$\mathcal H(Y)=A_n\int_0^\infty(g_+(h,Y)+g_-(h,Y))dh$
与差 $G_{\rm mix}-\mathcal H$。

**离散减连续。** 每个同余类的指数和减去密度是
$e^{-(h_0+x_0)Y}/(1-e^{-LY})-e^{-x_0Y}/(LY)$。
逐类 Taylor 展开到 $\ell$ 阶，给第 $j$ 阶的矩
$Z_\pm(-j;0)+A_n x_0^{j+1}/(j+1)$（Taylor 项另乘 $(-Y)^j/j!$）。
对 (2) 中的 $x^uY^v$，使用 $(-1)^uY^v\partial_Y^u$。为控制换序，对每类用尺度 $LY$。在 $LY\le1$ 处，移除简单极点后 Taylor 余项（含 $Y^u\partial_Y^u$）统一为 $O_n((LY)^{\ell+\theta})$，$0<\theta<1$。在 $LY\ge1$ 处同一上界由有界的指数矩与减去的有限 Taylor 项得到。这里 $(h_0+x_0)/L$ 位于有界区间，且 $t^u e^{-t}$ 有界。$\sum a^{-r_1}(b')^{-r_2}L^{\ell+\theta}$ 收敛，因为 $r_2\ge\ell+2$。故可把取系数与求和交换。把上述矩代入 (4) 的多项式，得到

$$
[Y^{1-s}](G_{\rm mix}-\mathcal H)
=C\sum_\pm\sum_{j=0}^\ell q_{\pm,j}
\left(Z_\pm(-j;0)+\frac{A_n x_0^{j+1}}{j+1}\right).
\tag{6}
$$

**连续部分。** 两支连续源之和由 Fourier 卷积分成 $2A_nF_{w+\epsilon}(n,Y)$ 减去带密度 $A_n$ 的有限频带 $[0,n]$。这是因为 $F_t(|f|,Y)$ 是 $Y^t/(u^2+Y^2/4)^t$ 的 Fourier 变换乘 $1/\pi$。用
$\mathcal E_\epsilon=2 A_nF_{w+\epsilon}(n,Y)/[(w+\epsilon)(w+\epsilon-1)-s(s-1)]$，
则 $\mathcal E_0=[\alpha E_w]_n$。

连续源的 Green 积分分为 $[0,\eta]$ 和 $[\eta,\infty)$。在收敛域，$[0,\eta]$ 的第一项直接用局部 Green 恒等式得 $\mathscr W(\psi,\mathcal E_\epsilon)(\eta)$。有限频带的内积分在 $\epsilon\to0$ 后无常数部：其奇异支的指数为 $-2p-\epsilon+2j$，其余支在零点可积，频率端点的系数也可积（使用小参数 Bessel 的两支表达，系数是非负整数幂或次幂 $2m_\epsilon-1+2j,2b_\epsilon-1+2j$）。外积分可在 $\epsilon=0$ 直接积分。式 (2) 的多项式次数上界给 $\mathcal H=O(Y^{-3}e^{-nY/2})$ 于无穷远，故外积分是 $-\mathscr W(\psi,\mathcal H)(\eta)$。变形外积分的局部解析性也可由 $hY$ 缩放得到：小 $h$ 区间在乘 $\psi Y^{-2}$ 并对 $h$ 积分后为 $O(Y^{-2-\min(\Re m_\epsilon,\Re b_\epsilon)})$，其余有指数衰减。

所以连续 Green 值的常数部是 $-\mathscr W(\psi,\mathcal H-\mathcal E_0)$ 的归一化常数部。此时交叉项确实消失。$\mathcal H-\mathcal E_0$ 的源是负的有限频带，其到 $Y^{1-s}$ 为止只含与 $1-s$ 相反奇偶的幂。两函数有下界有限的 Laurent 展开。由 $\mathcal L$ 的系数递推，$Y^{1-s-2j}$（$j\ge1$）全部消失。故该 Green 值等于 $[Y^{1-s}](\mathcal H-\mathcal E_0)$。

另一方面，在收敛域对连续频率使用 (3) 后再延拓，同一 Green 值由 $A_nC_\epsilon\int_{x_0}^\infty\mathbb E[x\pm n(1/2-U)]^{\ell+\epsilon}\,dx$ 的两支给出。展开至 $\ell+1$ 阶；低阶幂的积分在 $\epsilon=0$ 给 $-x_0^{j+1}/(j+1)$，最后一项给 $-n^{\ell+1}\mu_{\ell+1}/(\ell+1)$。Fubini 在上述实收敛域成立：小高度连续源为 $O(Y^{1-w-\epsilon})$，大高度的估计如上。积分余项在 $\epsilon$ 的邻域可求和，零点逐项消失。与 (6) 相加，$x_0$ 项抵消，留下 (5)。得到

$$
[Y^{1-s}](G_{\rm mix}-[\alpha E_w]_n)=c\,n^\ell W_{n,\mathrm{mix}} .
\tag{7}
$$

## 5. 从混合矩反演正正块

式 (2) 的 $\ell$ 次系数可以把积分下限改为零。理由是每条边界变量在 $Y=0$ 的次数至少 $r_2-1\ge\ell+1$，对截断小区间的误差阶严格高于 $\ell$。因此 (4) 是 $\mathscr B_\ell=[Y^\ell]\mathscr B$ 的双 Laplace 变换再乘 $P^{r_1}R^{r_2}$。按 $P,R$ 逐项比较：

$$
\mathscr B_\ell(u,v)=D_0\,u^{A-1}v^{B-1}(u+v)^\ell,\qquad
D_0=\frac{C}{(k)^{\overline\ell}\Gamma(A)\Gamma(B)} .
\tag{8}
$$

具体地，$P^iR^{\ell-i}$ 的 beta 矩因子为
$\binom{\ell}{i}(B)^{\overline i}(A)^{\overline{\ell-i}}/(k)^{\overline\ell}$。
除以两次 Laplace 阶乘 $(r_1-1-i)!(r_2-1-\ell+i)!$，只剩下 $D_0\binom{\ell}{i}$。

对于 $P,R>0,\ P+R=n$ 的正正频率，改回两条全纯腿。区域 $v\ge u$ 翻转小权腿的 $Y$ 符号，$(u+v)^\ell$ 变成 $(v-u)^\ell$。区域 $u\ge v$ 翻转大权腿的符号，得到 $(u-v)^\ell$。两块合成 $|u-v|^\ell$。因此正正系数为

$$
D_0 P^{r_1}R^{r_2}
\int_0^\infty\!\int_0^\infty u^{A-1}v^{B-1}|u-v|^\ell e^{-Pu-Rv}\,dv\,du .
$$

因 $\ell$ 为奇数，把绝对值项写成 $(v-u)^\ell+2(u-v)^\ell\mathbf1_{u>v}$。第一项给 $C\mathbb E[(1-U)P-UR]^\ell$，正对 FKL 式 (32) 的全符号核。第二项用 $u=v+t$ 展开 $u^{A-1}$。去掉因子 2 后，积分（含 $P^{r_1}R^{r_2}$）等于

$$
R^{B+\ell}\sum_{j=0}^{A-1}
\binom{A-1}{j}\frac{\Gamma(B+j)\Gamma(A+\ell-j)P^j}{n^{B+j}}
= n^\ell\Gamma(A)\Gamma(B)\ell!
(R/n)^{r_2}P_{A-1}^{(r_2,-r_1)}(2P/n-1).
$$

右式即 Jacobi 的反射超几何展开，展开系数由 $(1-A)^{\overline j}/(1-A-\ell)^{\overline j}$ 给出。乘 $2D_0$ 得 $2c(-1)^{s+m}n^\ell(R/n)^{r_2}P_{A-1}^{(r_2,-r_1)}$，与式 (32) 的有限支一致。两项合计给 $c n^\ell W_{n,++}$。这一步使用整个嵌套积分的三角域，因而保留了同指数齐次支的归一化。

## 6. 正零与深度一边界块

先令权 $2m$ 那条腿取全纯常数模。式 (8) 的同一三角划分现在可取 $R=0$，并用 DKS 的无穷端有限部分积分。混合项中的反全纯常数腿仅出现一个有限端点原函数，其齐次次数为 $r_2>\ell$，所以对目标系数无贡献。于是需要的零腿积分是

$$
\operatorname{FP}\int_0^\infty v^{B-1}|u-v|^\ell dv
=2u^{B+\ell}\frac{\Gamma(B)\ell!}{\Gamma(B+\ell+1)} .
$$

证明是先减去整半轴上的多项式 $v^{B-1}(v-u)^\ell$，其切向有限部分为零，剩下 $0<v<u$ 的普通 beta 积分。权 $2m$ 的常数模与正频率前系数的比例为
$(-1)^m r_2!\zeta(2m)/(2\pi)^{2m}$。
乘 $D_0 n^{r_1}$ 后积 $u$，得到

$$
c\,n^\ell\frac{\Gamma(k+\ell)}{\pi}\,
i^k\,\frac{(2\pi)^{-r_2}\zeta(1+r_2)n^{-r_2}}{\Gamma(A)}\,\sigma_{-r_1}(n).
$$

交换两腿即另一项，$i^k$ 消失，因为 $s+m+b$ 偶。两项是 FKL 式 (35) 两个 $Z$ 的第二部分，乘 $c n^\ell\Gamma(k+\ell)/\pi$ 和各自的迹公式相位。

还要加 DKS I 式 (4.6) 的深度一修正。删除内层小权腿后，剩余系数只取 (8) 中原来该腿的 $Y^0v^{r_2-1}$ 项，即 $D_0u^{A-1}v^{r_2-1}$。理由是该修正系数只依赖被删腿的权，$Y$ 的分母与保留的深度一项合计仍为 $Y^q$。单腿正频率的 $-2/r_1!$ 与修正 $-2\zeta(r_2)/r_2$ 相乘，除以双腿归一化 $4/(r_1!r_2!)$ 给 $\Gamma(r_2)\zeta(r_2)$。积 $u$ 后，深度一贡献为

$$
c\,n^\ell(-1)^{s+m}
\frac{\Gamma(r_2)\zeta(r_2)}{\Gamma(B)\ell!}\,\sigma_{-r_1}(n).
$$

用 $\zeta'(1-r_2)=(-1)^{m-1}\Gamma(r_2)\zeta(r_2)/(2(2\pi)^{r_2-1})$，这正是 $i^kZ_d^{r_1,r_2}$ 的第一部分乘 $c n^\ell\Gamma(k+\ell)/\pi$。交换两腿得到另一项。纯零频和纯反全纯积分常数不进入 $n>0$ 模。故整个正零和深度一块给 $c n^\ell D_n$。

## 7. 一般补全及余循环

三个块相加得到 $c n^\ell(W_n+D_n)$。$k\ge6$ 时，用校正后的迹公式化成 $c n^\ell T_n$，即 (1)。在模不变解 $F$ 中减 $\alpha E_w$；若 $\delta>0$，再减适当倍数的 $E_{\delta+1}$，消去余下在尖点高于线性的增长。由模不变性与 $E_{1+\eta}$ 的比较，每个固定 Fourier 模在小高度为 $O(Y^{-\eta})$，任取 $0<\eta<1$。$E_{\delta+1}$ 的非零模最小幂次是 $-\delta>1-s$，不改变 (1)。DKS I §4.3 使 $F,\widecheck F$ 的 Laurent 多项式相同，齐次差的零模为零。非零模只剩 $K_{s-1/2}$ 衰减支，由 $Y^{1-s}$ 系数确定。偶性给出负模。$k=4$ 时，由 DKS 的尖点空间论证有 $F=\widecheck F$，上述增长估计也给 $W_n+D_n=0$。

在任一归一化 Hecke 尖点形式 $f$ 上，校正符号后的 FKL 迹项为

$$
T_{n,f}=(-1)^s2^{2-2s}\frac{\Gamma(2s-1)}{\Gamma(A)\Gamma(B)}
\frac{\Lambda_f(B)\Lambda_f(w+s-1)}{\langle f,f\rangle}
\frac{a_f(n)}{n^{w+s-2}} .
$$

这是[已有迹公式计算](../weight-12/s-cocycle-proof.md)逐项用于权 $2s$ 的本征形式，随后完成两个 $L$ 值。$H_f^+$ 的 $Y^{1-s}$ 系数为 $\Gamma(2s-1)a_f(n)/(2\Gamma(s)n^{2s-1})$。用 $\Lambda_f(B)=(-1)^s\Lambda_f(A)$，(1) 的相反数给

$$
F-\widecheck F=\sum_f \lambda_f H_f^+,\qquad
\lambda_f=-c\,2^{3-2s}\frac{\Gamma(s)}{\Gamma(A)\Gamma(B)}
\frac{\Lambda_f(A)\Lambda_f(w+s-1)}{\langle f,f\rangle}.
$$

Gamma 倍元把 $c$ 化为
$(-1)^{A/2}\Gamma(A)\Gamma(B)\Gamma(p)/(\Gamma(s+p)\widetilde C_{s;s_1,s_2})$，
其中 $\widetilde C$ 为 DGW 式 (3.50)。代入即 DGW 式 (4.21)–(4.22) 的猜想系数，同时给出一般 $S$-余循环

$$
\delta_S\widecheck F(\tau)=
-\sum_f\lambda_f\frac{2(-1)^s\pi^{2s-1}}{\Gamma(s)y^{s-1}}
\operatorname{Re}\!\left[i\int_0^{i\infty}
(\tau-z)^{s-1}(\bar\tau-z)^{s-1}f(z)\,dz\right].
$$

该周期积分由 $H_f^+$ 中代换 $z\mapsto-1/z$ 得到，可用 $i^{j+1}\Lambda_f(j+1)$ 逐项展开。证明涵盖一般整数指标的补全与整个余循环。由单一余循环反解特定约化多重模值时，还要看 DKS 的变换式辨认了哪种多重模值线性组合。

## 8. 精确复核范围

使用第 2 节的末分量 $v_{m-2}=r_1!r_2!/[bm\Gamma(b)^2\Gamma(m)^2]$，由左特征方程反向逐项求出 $v_j$。在 $2\le m\le11$、$m\le b\le16$ 的全部 495 组允许的 $(s;b,m)$ 中，检查了第一行边界方程及 (8) 的全部单项式系数。权重 12 的两个既有算例、权重 18 算例的 Jacobi 核、留数、边界项和余循环也与专稿的精确运算一致。
