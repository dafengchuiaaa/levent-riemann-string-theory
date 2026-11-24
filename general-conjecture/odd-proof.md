# 奇谱的 FKL 卷积、DKS 补全与 S 余循环

证明稿，待独立审阅。单模 Gauss 方程沿用 Neon 给出的全导数恒等式。混合尾的次极点用两个移位除数因子的联合 Dirichlet 级数计算。符号采用本档案的 [Petersson 迹公式订正](../weight-12/s-cocycle-proof.md)。

设 $b=s_1>m=s_2\ge2$，$w=b+m$，$\delta=b-m$，$s=\delta+1,\delta+3,\ldots,w-1$。记

$$
Y=4\pi\Im\tau,\quad K=2s,\quad r_1=2b-1,\quad r_2=2m-1,\quad
d_0=s-w,\quad \ell=-d_0-1=2p-2,\quad q=w-2,
\quad A=s+m-b,\quad B=s+b-m.
$$

这里 $A+B=K$，$r_2=A+\ell$，$r_1=B+\ell$。对 $n>0$，令 $T_n$ 为符号订正后 FKL (34) 的尖点项乘 $\Gamma(K+\ell)/\pi$。主恒等式为

$$
[Y^{1-s}][\widecheck F^{-(s)}_{m,b}]_n=c\,n^\ell T_n,\qquad
c=(-1)^{m+s+p+1}
\frac{2(p-1)!(2s-2)!(b-p)!(m-p)!}{(s-1)!(s+p-2)!(b-1)!(m-1)!}.
\tag{1}
$$

$K=4$ 时右边取零，其卷积恒等式由下文直接推出。底谱的混合块另有一个有限项，正零块包含抵消它的齐次 Eisenstein 项。

## 1. 单对频率与 DKS 核

用 $D=Y^2\partial_Y^2-n^2Y^2/4-s(s-1)$ 表示径向算子。两个频率 $f_m+f_b=n$ 的奇源是
$$
-\tfrac{Y^2}{2}(f_m F_m F'_b-f_b F_b F'_m),\qquad
F_j(f,Y)=\frac{2}{\sqrt{\pi}\Gamma(j)}\sqrt Y\, |f|^{j-1/2}K_{j-1/2}(|f|Y/2).
$$
此式省略了除数函数。DKS I (5.6) 的 $m$ 阶矩阵是 (3.7) 的全三对角算子在反射偶子空间上的表示。中点行的非对角元加倍。全算子保持多项式次数，次数 $h$ 的最高项特征值为 $(q-h)(q-h+1)$，所以 $h=0,2,\ldots,2m-2$ 给出上述全部 $s$。用 (5.5) 把左特征向量的末项定为 $-r_1!r_2!/(2[(m-1)!(b-1)!]^2)$。写该向量为 $v_j$，$0\le j\le m-1$。

对小权负频 $-R$ 和大权正频 $P$，从 DKS I (2.23)、(2.34)、(7.4) 得到混合核
$$
Y^{-q} R^{r_2}P^{r_1}\int_{Y/2}^{\infty}\int_{Y/2}^{\infty}
e^{-Ru-Pv}\mathcal B(Y,u,v)\,dv\,du,
$$
$$
\mathcal B=\frac4{r_1!r_2!}\sum_{j=0}^{m-1}v_j
(u+Y/2)^j(u-Y/2)^{2m-2-j}
(v-Y/2)^{q-j}(v+Y/2)^{\delta+j}.
$$
反向频率取此核的负号。DKS I (3.7)、(5.3)、(5.5) 的微分与 shuffle 恒等式逐对作用于两条非零指数腿。沿标记保留两个权的顺序时，出界差在 (5.3) 恰成上式的反对称源。深度一项和反全纯积分常数没有 $q^P\bar q^R$ 模。故混合表达式满足逐对径向方程。这些恒等式只微分指数和多项式权，对正实 $P,R$ 同样适用。

置 $\psi=\sqrt{nY}I_{s-1/2}(nY/2)$，$w_{n,s}=(2s-1)[Y^s]\psi$。对单个严格混合频率，$\mathscr W(\psi,g)=\psi g'-\psi' g$ 在无穷端归零。奇源的双奇异支幂次为 $3-w+2j$，与 $1-s$ 反奇偶，其余支至少为 $Y^{2-\delta}$。径向递推从最低幂开始消去所有 $Y^{1-s-2j}$，$j\ge1$。故 $\operatorname{fp}\int_0^\infty\psi\,\phi\,dY/Y^2$ 除以 $w_{n,s}$ 给出 $[Y^{1-s}]g$。

单模积分的计算可从确实收敛的参数域开始。用临时参数 $(s,m',b')$，$z=n/(h+n)$，将高度按 $(h+n)/2$ 缩放。记 $v=\sqrt Y I_{s-1/2}(zY)$，$u=\sqrt Y K_{m'-1/2}((1-z)Y)$，$w_*=\sqrt Y K_{b'-1/2}(Y)$，$J=v[(1-z)u w_*'+w_*u']$。对 $\mathcal P=z^{1/2-s}(1-z)^{m'-1/2}$，参数 $(s-m'-b'+1,s+b'-m',2s)$ 的 Gauss 算子作用于 $\mathcal P J$ 等于 $\mathcal P\,\partial_Y G$，其中
$$
G=\big((1-z)b'(b'-1)-m'(m'-1)+z(1-z)Y^2\big)uv w_*
+Y^2\big(zu'v w_*'-(1-z)uv'w_*'-u'v'w_*\big).
$$
代入 $u''=((1-z)^2+m'(m'-1)/Y^2)u$ 及另外两条 Bessel 方程即可核对。$\Re m',\Re b'>1/2$、$\Re(s-m'-b')>0$ 时，两端边界归零。$z\to0$ 的解析支由 $v=O(z^{s-1/2})$ 确定。对 $z=0$ 的 $\frac{h}{2}(F_{m'}F_{b'})'$ 分部积分，再用两个 $K$ 的 Gamma 矩，得到归一化常数
$$
C(s,m',b')=-\frac{2^{2s}}{\pi\Gamma(m')\Gamma(b')\Gamma(s)(2s-1)}
\Gamma\!\left(\tfrac{s+m'+b'}2\right)
\Gamma\!\left(\tfrac{s+m'-b'+1}2\right)
\Gamma\!\left(\tfrac{s-m'+b'+1}2\right)
\Gamma\!\left(\tfrac{s-m'-b'+2}2\right).
$$
两权可独立延拓。整数点的 Green 积分有限部由小高度幂次扣除给出，双奇异支在该点没有 $Y^{-1}$ Green 被积项。记 $C=C(s,m,b)=c(-1)^b(K)_\ell/\ell!$，其中 $(K)_\ell$ 为上升阶乘。令 $U\sim{\rm Beta}(B,A)$、$X=\tfrac12-U$，$x=h+n/2$。对 $(P,R)=(h+n,h)$ 和 $(h,h+n)$，两支单模系数分别是
$$
C\,\mathbb E(x+nX)^\ell,\qquad -C\,\mathbb E(x-nX)^\ell .
\tag{2}
$$

## 2. 两个移位除数因子的 Dirichlet 级数

令 $a_\tau(h)=\sigma_{-R_1}(h+n)\sigma_{-R_2}(h)$ 当 $\tau=+$，反向放置两权当 $\tau=-$。这里 $R_1,R_2$ 是复参数，整数点为 $r_1,r_2$。定义
$$
Z_\tau(v;R_1,R_2)=\sum_{h\ge1}a_\tau(h)(h+n/2)^{-v}.
$$
两支在 $v=1$ 的留数为
$\rho_n=\zeta(R_1+1)\zeta(R_2+1)\sigma_{-R_1-R_2-1}(n)/\zeta(R_1+R_2+2)$。
将两个除数写成 $g a,g b'$，其中 $g\mid n$、$(a,b')=1$，即可由同余类密度 $1/\operatorname{lcm}(ga,gb')$ 得到此式。

在底谱附近还有一个共同的极面 $R_2+v=1$，其留数沿变量 $R_2+v$ 为
$$
\mathcal R(v;R_1,n)=
\zeta(v)\frac{\zeta(R_1+1)}{\zeta(R_1+v+1)}\,
\sigma_{-R_1-v}(n).
\tag{3}
$$
证明如下。设 $q_g=n/g$。小权腿的频率是 $g b' j$，大权腿的频率是 $g b'j+\tau n$，其中 $b'j\equiv-\tau q_g\bmod a$。对 $g b'>n$，所有 $j\ge1$ 都可用。取 $j_0\in[1,a]$ 为该同余类首项。对固定 $a,b'$，$j$ 求和等于
$$
(g a b')^{-v}\zeta\!\left(v,\frac{j_0}{a}+\frac{\tau q_g}{2ab'}\right).
$$
在 $v=1-r_2$ 附近移除第二自变量中的位移，误差有额外的 $O_n((ab')^{-1})$，因为 $\partial_z\zeta(v,z)=-v\zeta(v+1,z)$ 在此实部范围对 $0<z\le 3/2$ 一致有界于紧参数集。有限个 $g b'\le n$ 的项也在该处正则。余下将 $b'$ 按模 $a$ 的单位类分组。每一类在 $R_2+v=1$ 的留数为 $1/a$，因此联合留数是
$$
\sum_{g\mid n}g^{-R_1-1}\sum_{a\ge1}a^{-R_1-1}
\sum_{u\bmod a}^{*}\sum_{\substack{j\ge1\\j\equiv-\tau q_g u^{-1}\bmod a}}j^{-v}.
$$
最后的和以 Hurwitz zeta 延拓。换 $u$ 为 $-u$ 显示两支的留数相同。设 $D_g=(q_g,a)$，内层计数等于
$\frac{\varphi(a)}{\varphi(a/D_g)}D_g^{-v}\zeta(v)\prod_{\mathfrak p\mid a/D_g}(1-\mathfrak p^{-v})$。
对 $\mathfrak p^e\Vert n$ 令 $L=\mathfrak p^{-R_1-1}$、$M=\mathfrak p^{-v}$。去掉 $\zeta(v)$ 后，按 $e$ 取生成函数，$g=\mathfrak p^t$ 给 $L^t$，$a=\mathfrak p^j$ 的两段在 $j\le e-t$ 和 $j>e-t$ 分开。合计为
$$
\frac1{1-L z}\left[
\frac{1-LMz}{(1-z)(1-\mathfrak p LMz)}
+\frac{L(1-M)}{(1-L)(1-\mathfrak p LMz)}
\right]
=\frac{1-LM}{(1-L)(1-z)(1-\mathfrak p LMz)}.
$$
这就是 (3) 的 Euler 因子。

上述分解给出实际的亚纯延拓。对 $b'$ 同余级数，扣去 $1/[a(R_2+v-1)]$ 后，每一类的 Hurwitz 余项之和至多为 $O(a^\eta(1+\log a))$，在 $|\Re(R_2+v)-1|<\eta$ 可取任意足够小的 $\eta>0$。这是将 Hurwitz 级数的首项 $(u/a)^{-R_2-v}$ 单列，再对余项使用 Euler–Maclaurin 所得。外层因 $r_1+1-r_2=K-1\ge3$ 而可求和。相同估计在固定竖条中仅多出 $|\Im v|$ 的多项式因子。离开次极面右侧，原同余类表达直接给出延拓和同类竖条界，只需单列 $v=1$ 的留数。正实部时可用首项 $h+n/2\gg_n\max(ga,gb')$，因此从 $v>1$ 可接到整数点附近；沿固定权重的实路径可选 $R_i>0$、$R_i+v>1$、$R_1+R_2+v>2$，底谱仅在终点触及 $R_2+v=1$。

## 3. 混合 Green 值与底谱有限项

先在固定整数参数处取小高度。$Z_\tau(v;r_1,r_2)$ 在 $v=-\ell$ 正则，即使在底谱也成立，因为 (3) 含 $\zeta(-\ell)=0$ 且底谱 $\ell\ge2$。对热核 $\sum a_\tau(h)e^{-xY}$ 作 Mellin 反演，将 $\Gamma(v)Z_\tau(v)$ 的积分线移过 $-\ell$。第 2 节的竖条界与 Gamma 的指数衰减允许移线，也允许有限阶 $Y\partial_Y$ 微分。主极点给 $\rho_n/Y$，其余到该阶的热矩为 $Z_\tau(-j)$，乘 $(-Y)^j/j!$。减去密度 $\rho_n\int_{x_0}^\infty e^{-xY}dx$，其中 $x_0=n/2$，热矩随之改成 $Z_\tau(-j)+\rho_n x_0^{j+1}/(j+1)$。

DKS 单对表达式乘 $Y^q$ 后是 $e^{-xY}$ 乘 $YP,YR$ 的有限多项式，次数不超过 $q$。将其单项 $Y^t x^u$ 写成 $Y^t(-\partial_Y)^u$ 作用在热核，且 $u\le t$。所以第 $\ell$ 阶系数恰是对 (2) 的两个多项式使用上述扣除密度后的热矩。Mellin 余项在操作后仍是 $O(Y^{\ell+\eta})$。由此给出离散减密度的全部端点项，不需要逐频假定无穷和的奇偶性。

连续部分从另一端计算。$F_t(f,Y)$ 是 $[Y/(u^2+Y^2/4)]^t$ 的 Fourier 变换除以 $\pi$。两个幂的梯度成比例，全频轴的反对称卷积为零，故混合连续源等于负的正正有限频带。该频带的奇异支幂次和交叉支下界与第 1 节相同，端点频率可积。连续特解由 DKS 混合核积分，其无穷端为 $O(Y^{-1}e^{-nY/2})$，因为 $hY$ 缩放后，多项式的 $Y$ 次数不超过 $q$。因此 Green 恒等式的无穷端为零，小高度没有低阶 Wronskian 交叉项。连续特解的目标系数等于连续 Green 积分有限部。

保留两个权的独立变形，从收敛域应用单模公式再延拓，连续 Green 值由下式的两端原函数给出（写 $\nu=-d-1$，beta 参数继续为 $v+R_1,v+R_2$，$v=d+1$）：
$$
\operatorname{AC}\int_{x_0}^\infty\mathbb E\big[(x+nX)^\nu-(x-nX)^\nu\big]dx
=-\frac{\mathbb E[(x_0+nX)^{\nu+1}-(x_0-nX)^{\nu+1}]}{\nu+1}.
\tag{4}
$$
有限频带的小高度扣除在整数点无共振，故延拓值确为上面的 Green 有限部。将 (4) 乘 $C\rho_n$ 并取 $\nu=\ell$。低阶矩逐项抵消离散减密度所加的端点项。令 $\mu_j=\mathbb E X^j$、$b_j=\binom{\ell}{j}n^j\mu_j$，得到实际 DKS 混合系数
$$
M_n=C\left\{
\sum_{j=0}^{\ell} b_j[Z_+(j-\ell)-(-1)^j Z_-(j-\ell)]
-\frac{2\rho_n n^{\ell+1}\mu_{\ell+1}}{\ell+1}
\right\}.
\tag{5}
$$
式 (5) 的最后一项是主极点与中心奇矩的有限乘积。

比较 FKL 时先在 $2d+R_1+R_2+2=K$ 上延拓完整带相位和，再限制到 $d=d_0-\epsilon,\ R_j=r_j+\epsilon$。例如可从 $(d,R_1,R_2)=(1/4,1/2,K-3)$ 起步；该点满足固定权重与绝对收敛，对 $K\ge6$ 也满足 FKL (66)。有限频带由 FKL (157)–(162) 延拓，Gamma 参数 $A,B,r_2+1$ 均为正，边界和尖点项用 (9)、(11) 直接取极限。这覆盖 $d_0=-1$ 及 $d_0=-r_2$。$d_0=-1$ 所需的全纯尖形式 $L$ 值没有极点。

混合支的 FKL 核，在乘 $\Gamma(d+1)$ 前，是
$$
\frac{2 n^v}{\Gamma(K)}\,\mathfrak p_\tau\,
\mathbb E (x+\tau nX)^{-v},
\quad \mathfrak p_+=i^K\cos(\pi R_2/2),\quad \mathfrak p_-=\cos(\pi R_1/2).
$$
在 $\ell+1$ 阶展开。截断余项在整数点附近以 $O_n(x^{-2+\Re\epsilon})$ 控制，且逐项趋零。除底谱次极点外，各项在等量变形上给 (5)。其中 $j=\ell+1$ 由 $\binom{\ell+\epsilon}{\ell+1}\sim\epsilon/(\ell+1)$ 乘 $Z_\pm(1-\epsilon)\sim-\rho_n/\epsilon$ 得到，符号与 (4) 一致。

在底谱只需对 $j=0$ 另算一个有限项。设 $\alpha=d-d_0$、$\beta=R_2-r_2$，$\kappa=d+R_2=\alpha+\beta$。固定权重使 $R_1-r_1=-2\alpha-\beta$，因而
$$
\mathfrak p_++\mathfrak p_-
=-2(-1)^b\sin(\pi\kappa/2)\cos(\pi\alpha/2).
$$
先在两变量中以此式消去 (3) 的 $1/\kappa$，再限制到 $\beta=-\alpha$。相位商给 $-\pi(-1)^b$，而 $\Gamma(v)\zeta(v)\to\zeta'(-\ell)/\ell!$。在 $Z_\tau$ 的正则部分，等量变形使相位相反，因此第 (5) 式的差仍按原式。两支在固定整数权时的可去极点值还各有同一个 $\partial_v\mathcal R$，在差中抵消。可见完整 FKL 混合值乘 $c n^\ell$ 等于 $M_n$ 加
$$
-2C\,\zeta'(-\ell)\frac{\zeta(2b)}{\zeta(K)}\sigma_{1-K}(n).
\tag{6}
$$
这是两移位除数因子产生的有限项。它没有小高度对数，因为 (3) 在固定整数参数的分子为零。用 zeta 函数方程及 Bernoulli 偶值公式，(6) 等于
$$
\lambda_{m,b}\,\zeta(r_2)[Y^{1-s}][E_s^{(-)}]_n,\qquad
\lambda_{m,b}=(-1)^{m-1}\frac{(s-1)!s!}{(b-1)!b!}\frac{B_{2b}}{B_{2s}} .
\tag{7}
$$
这里 $[Y^{1-s}][E_s^{(-)}]_n=2(2s-2)!\sigma_{1-K}(n)/(s-1)!^2$。例如 $\lambda_{2,7}=455/4146$，所以 $n=1$ 的 FKL 混合值相对 DKS 多 $38220\,\zeta(3)/691$。

## 4. 正正、正零和 DKS 底谱项

以 $u$ 对应小权，$v$ 对应大权。两次 Laplace 反演 (2) 给
$$
[Y^\ell]\mathcal B(Y,u,v)=D_0 u^{A-1}v^{B-1}(u+v)^\ell,\qquad
D_0=\frac C{(K)_\ell\Gamma(A)\Gamma(B)}.
$$
反演可把下限移至零。宽度为 $Y/2$ 的截断条带中，该腿的多项式是 $O(Y^{r_i-1})$，积分后为 $O(Y^{r_i})$，而 $r_i>\ell$。由 DKS (2.23) 的同向嵌套顺序，$u>v$ 属于小权在内的项，翻小权腿即交换该腿的两个线性因子。$v>u$ 属于大权在内的项，翻大权腿并取反号。核因此变成 $\operatorname{sgn}(u-v)|u-v|^\ell=(u-v)^\ell-2(v-u)^\ell\mathbf1_{v>u}$。设 $P+R=n$，$P,R$ 分别为大、小权正频。全象限第一项给 $C\mathbb E[(1-U)P-UR]^\ell$。在后三角令 $v=u+t$，其积分连同 $P^{r_1}R^{r_2}$ 为
$$
R^{r_2}\sum_{j=0}^{B-1}\binom{B-1}{j}
\frac{\Gamma(A+j)\Gamma(B+\ell-j)P^j}{n^{A+j}}
= n^\ell\Gamma(A)\Gamma(B)\ell!\,(R/n)^{r_2}
P_{B-1}^{(r_2,-r_1)}(2P/n-1).
$$
最后一步是 Jacobi 多项式的有限二项式展开。乘 $-2D_0$ 后相位为 $-2c(-1)^b=2ci^{K+r_2+1}$，两项合计等于 FKL (32) 的有限正正块乘 $c n^\ell/\pi$。

取一条腿的全纯常数模时，无穷端使用 DKS 的切向有限部。整条半轴上的多项式有限部为零，剩下普通的 beta 积分：
$$
\operatorname{fp}\int_0^\infty u^{A-1}\operatorname{sgn}(u-v)|u-v|^\ell du
=-2v^{r_2}\frac{\Gamma(A)\ell!}{\Gamma(r_2+1)},\qquad
\operatorname{fp}\int_0^\infty v^{B-1}\operatorname{sgn}(u-v)|u-v|^\ell dv
=2u^{r_1}\frac{\Gamma(B)\ell!}{\Gamma(r_1+1)}.
$$
反全纯常数腿只给齐次次数 $r_i>\ell$ 的下端点多项式，在目标阶无贡献。小权常数与其归一化正频前因子的比为 $(-1)^m r_2!\zeta(1+r_2)/(2\pi)^{1+r_2}$，大权交换指标。代入上式并积另一腿，恰给 FKL 两个 $Z$ 中的 $\zeta(1+r_i)$ 项，乘 $c n^\ell\Gamma(K+\ell)/\pi$ 及各自相位。

DKS I (4.6) 的深度一项删除内腿。小权被删时，剩余深度一的 $Y$ 分母与 (4.6) 的 $Y^{-j}$ 合为 $Y^{-q}$，所以被删腿只取 $Y^0u^{r_2-1}$，原核留下 $D_0v^{B-1}$。以单腿归一化 $-2/r_1!$ 乘 $-2\zeta(r_2)/r_2$，再除以双腿的 $4/(r_1!r_2!)$，得到删腿因子 $\Gamma(r_2)\zeta(r_2)$。积 $v$ 后结果为
$$
c n^\ell (-1)^b\frac{\Gamma(r_2)\zeta(r_2)}{\Gamma(A)\ell!}\sigma_{-r_1}(n).
$$
大权被删时交换 $A,B$ 与 $r_1,r_2$，并加负号。用
$\zeta'(1-r_2)=(-1)^{m-1}\Gamma(r_2)\zeta(r_2)/(2(2\pi)^{r_2-1})$
及大权对应式，得到两个 $Z$ 的 $\zeta'$ 部分。四个普通边界项合为 FKL (34) 中两个边界项的相反数，采用和混合块相同的归一化。

底谱 DKS I §5.5 还有一项。可用零模直接算其常数。取 $N=2m-2$，此时 DKS I (5.12) 给 $v_j=V(-1)^j\binom Nj/(1+\mathbf 1_{j=m-1})$，$V=(-1)^m r_1!r_2!/[(b-1)!^2N!]$。在未加入底谱项时，普通删小权腿的 (4.8) 产生 $Y^s$。用 (4.10) 展开，系数为
$$
-\frac{2\zeta(r_2)}{r_2}\frac{V(-1)^s B_{2b}}{(2b)!r_1!}
\sum_{j=0}^{m-1}\frac{\binom Nj(\delta+j)!(q-j)!}{1+\mathbf 1_{j=m-1}}.
$$
反射延拓后有限和是全和的一半。把阶乘化为 beta 积分，和等于 $r_1!(\delta!)^2/[2(2\delta+1)!]$。DKS 通过 $E_s$ 对应的深度一组合消去该正幂，其 $Y^s$ 系数为 $(-1)^{s+1}B_{2s}/(2s)!$。除以此系数并取相反数，得到 (7) 的 $\lambda_{m,b}\zeta(r_2)$。按 DKS (5.38) 完成后，它在非零模是 $+\lambda_{m,b}\zeta(r_2)E_s^{(-)}$。于是它与第 3 节实有混合块相对 FKL 的负偏移精确相消。

## 5. 权重 4 与全谱补全

权重 4 即 $s=2,\ b=m+1$。不调用 FKL 定理 A。直接检查 DKS 特解的二阶 Cauchy–Riemann 导数。$m=2$ 的公式为 DKS I (5.40)。当 $m\ge3$，令 $a=m-1$、$N=2a$，用 $-4\pi\nabla$ 两次作用于 (4.6) 修正的 (2.34) 组合，没有外腿端点项。按 (4.7)，在第一种次序中，新指标写作 $(N-l,3+l)$。系数由 $j=l,l+1,l+2$ 三项给出：
$$
v_l(N+1-l)(N-l)+2v_{l+1}(l+1)(N-l)+(l+2)(l+1)v_{l+2}.
$$
另一种次序交换标签后得到同一系数。代入底谱向量，只有 $l=a,a-1,a-2$ 留下，系数依次为 $v_a(a(a+1),-2a,-a(a-1))$。所以双积分合为三个 shuffle。按 (2.29) 和 (4.6)，每个 shuffle 是 Eisenstein 导数乘积加同一个纯幂项，纯幂因三个系数和为零而消去。底谱的深度一项贡献 $6\lambda_{m,b}\zeta(r_2)(\Im\tau)^4G_4$。它的反全纯修正只有 $Y^0,Y^{-1}$ 乘反全纯函数，被二阶导数消去。余下零模常数也消去。因此 $(\pi\nabla)^2\widecheck F$ 与 $(\pi\nabla)^2 F$ 都有相同的模权。

令 $H=F-\widecheck F$。齐次方程和 DKS I (5.46) 说明 $(\pi\nabla)^2H/(\Im\tau)^4$ 是权 4 全纯模形式。两解在尖点指数衰减，所以它是尖形式，且 $S_4=0$。正频衰减 Bessel 支的二阶导数非零，奇性又确定负频，零模为零。因此 $H=0$。这也推出按第 3 节的算术定义取到的权 4 卷积恒等式。

当 $K\ge6$，三块相加并用 FKL (11) 的延拓，边界项和底谱项抵消，留下 (1)。模不变的 $F$ 在尖点衰减，从而在上半平面有界，其每个固定 Fourier 模的 $Y^{1-s}$ 系数为零。$F-\widecheck F$ 的非零模式只含衰减齐次支，由该系数唯一确定。订正后的 FKL 尖点项用完成 $L$ 值写成
$$
T_{n,f}=(-1)^s\,2^{2-2s}\frac{\Gamma(2s-1)}{\Gamma(A)\Gamma(B)}
\frac{\Lambda(f,A)\Lambda(f,w+s-1)}{\langle f,f\rangle}
\frac{a_f(n)}{n^{w+s-2}}.
$$
按 DKS II (2.31)，对归一化本征形式 $f$，该支在 $H^-_f$ 中的系数为 $(2s-2)!\,a_f(n)/(2(s-1)!n^{2s-1})$。因此
$$
F-\widecheck F=\sum_f a^-_f H^-_f,\qquad
a^-_f=(-1)^{m+p}2^{4-2s}
\frac{(p-1)!(2s-2)!(m-p)!(b-p)!}
{(s+p-2)!(m-1)!(b-1)!\Gamma(A)\Gamma(B)}
\frac{\Lambda(f,A)\Lambda(f,w+s-1)}{\langle f,f\rangle}.
$$
写 $y=\pi\Im\tau$，采用 DKS II 的 $\delta_S G=G(\tau)-G(-1/\tau)$，由 (2.35) 得完整余循环
$$
\delta_S\widecheck F
=-\sum_f a^-_f\,\frac{2(-1)^s i\pi^{2s-1}y^{1-s}}{\Gamma(s)}
\sum_{\substack{0\le j,k\le s-1\\j+k\ \mathrm{odd}}}
(-1)^{(j+k-1)/2}
\binom{s-1}{j}\binom{s-1}{k}
\tau^{s-1-j}\bar\tau^{s-1-k}\Lambda(f,j+k+1).
$$
权重 4 的同一公式为空和。

[FKL 预印本](https://arxiv.org/abs/2512.21413) · [DKS I](https://arxiv.org/abs/2109.05017) · [DKS II](https://arxiv.org/abs/2109.05018)
