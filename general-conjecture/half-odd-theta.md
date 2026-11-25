# 行列式插入的奇 theta 提升及其周期消去

DGW 的四维核有一个在 $B=0$ 消失的奇分量。保留它的一阶法向导数，可以对一整族局部 Maass 函数定义奇提升。下面从该积分自身导出它选择的广义 Eisenstein 级数组合，再计算每个全纯尖形式的总系数。采用[一般半整数奇源命题](half-odd-general.md)的参考特解。原始格点核、三个模变量的 Laplace 关系及 Schwinger 积分归一化取自 [DGW §§2–3、附录 A](https://arxiv.org/html/2501.03996v2#S2)。

**命题 2。** 写 $l=s-1$，$s\ge3$。取实多项式 $P(z)=\sum_{q=0}^{2l}p_qz^q$，允许次数严格小于 $2l$。要求
$$
 z^{2l}P((z-1)/z)=P(z),\qquad P(1-z)=-P(z). \tag{2}
$$
令 $\rho=u+ih$，$\mathcal T$ 为顶点 $0,1,\infty$ 的开理想三角形。将 $P$ 极化并除以 $h^l$，即置
$$
 A_P(\rho)=h^{-l}\sum_q\frac{p_q}{\binom{2l}q}
 [T^q](1+\rho T)^l(1+\bar\rho T)^l .
$$
对 $\tau=X+iY_0$ 用 $Q(p)=|p_1+p_2\tau|^2/Y_0$，$\langle p,r\rangle_Q$ 是其极化内积。$p,r$ 遍历 $\mathbb Z^2$，$[p,r]=p_1r_2-p_2r_1$。定义
$$
 \Theta^-(\tau;u,h,V)=i\pi\sum_{p,r}[p,r]\,
  \exp\{-\pi(V/h)Q(p+ur)-\pi Vh Q(r)\},
\quad
 L_P^w(\tau)=\int_0^\infty V^{w+1}\frac{dV}{V}
 \int_{\mathcal T} A_P(\rho)\Theta^-(\tau;u,h,V)\frac{du\,dh}{h^2}.
 \tag{3}
$$
核相当于 $\partial_B\Gamma_{2,2}|_{B=0}/(2V)$。$V^{w+1}$ 使边界的两条 Eisenstein 腿权重之和为 $w$。式 (3) 给出关于 $\tau$ 模不变且在 $\tau\mapsto-\bar\tau$ 下变号的函数。令 $e_q=q-l,\ m_q=(w-e_q)/2,\ b_q=(w+e_q)/2$，$c_q=\binom l{q/2}/\binom{2l}q$（$q$ 偶）。有
$$
 L_P^w=\frac32\sum_{q\ {\rm even}}p_q c_q\Gamma(m_q)\Gamma(b_q)
 F(s;m_q,b_q). \tag{4}
$$
在 $w$ 为整数且 $w\equiv s\bmod2$ 时，(4) 中每个 $H_f^-$ 的总系数为零。该结论是对 (3) 所定义的全部提升的恒等式。允许指数超出正半整数区间时，(4) 以 Eisenstein 的亚纯延拓读之。若只需所有腿均 $\ge3/2$，只保留满足该支撑条件的 $P,w$。本构造未假定物理振幅选出这些提升，也未给出所有可能奇提升的分类。

## 格点积分与三个边

带行列式的求和在 $p=0,r=0,p+r=0$ 上自动为零。整数行变换证明两个变量 $\rho,\tau$ 各自的 $\mathrm{SL}_2(\mathbb Z)$ 不变性，$\tau$ 的反射将行列式变号。由 DGW 的 $\Delta_\tau\Gamma_{2,2}=\Delta_\rho\Gamma_{2,2}$ 对 $B$ 微分即得两 Laplacian 在 $\Theta^-$ 上相等。令 $U\rho=(\rho-1)/\rho$。极化与实线性换元相容，因此 $A_P$ 在 $U$ 下不变。它满足 $\Delta_\rho A_P=l(l+1)A_P$：对实 $T$，生成函数 $h^{-l}|1+\rho T|^{2l}$ 是 $h^{-l}$ 的实行列式一分式变换（可取下行为 $(T,1)$），逐项即得。边 $u=0$ 上 $A_P(ih)=\sum_{q\ {\rm even}}p_qc_q h^{q-l}$。

以下估计也处理了不紧积分域的端点。对于固定 $\tau$，在 $h\ge1$、有界 $u$ 的尖点区，对 $p$ 作二维 Poisson 求和。被移位的线性因子 $[p+ur,r]$ 使零对偶频率消失，$r=0$ 本来即为零。余项有 Gaussian 因子
$$
 \exp[-\pi hQ^*(k)/V-\pi VhQ(r)],\qquad k,r\ne0 ,
$$
乘 $k,r$ 的多项式、$h,V,V^{-1}$ 的幂和模一相位。有限次求导只改变多项式界（零对偶项在求导前已恒等消失）。在 $\tau$ 的紧集上，求和后由这些变量的 Laurent 多项式界乘 $\exp[-c h(V+V^{-1})]$ 控制：从两个非零格矢各留一半指数，剩余 Gaussian 矩用对应的二维积分界。$U$ 将其他两尖点送到此处。在三角形的紧核部分同样可用 $h$ 的正下界。于是 $L_P^w$ 对 $w$ 整，分部积分、求导及尖点截断极限均可交换。对 $|X|\le1/2,\ Y_0\ge1$，$Q,Q^*$ 最小特征值均 $\gg1/Y_0$，系数的变化至多是 $Y_0$ 的幂次。上述估计改用 $\exp[-c h(V+V^{-1})/Y_0]$，以 $e^{-t}\ll_a t^{-a}$ 控制积分两端，给出模尖点的多项式增长。

格林公式的边项是 $A_P\partial_{\boldsymbol n}\Theta^--\Theta^-\partial_{\boldsymbol n}A_P$ 配欧氏弧长。$U$ 循环三条边且保持双线性边积分，故只需算 $u=0$ 再乘 3。该边上 $p\mapsto-p$ 给 $\Theta^-=0$，外法向为 $-\partial_u$。写 $x=V/h,\ y=Vh$，则
$$
 -\partial_u\Theta^-|_{u=0}
 =2i\pi^2x\sum_{p,r}[p,r]\langle p,r\rangle_Q e^{-\pi xQ(p)-\pi yQ(r)}
 =y^{-1}\mathcal J(\theta_x,\theta_y),\quad
 \theta_x=\sum_p e^{-\pi xQ(p)} .
$$
最后一步可直接对 $Q$ 求导，$Y_0^2(Q(p)_X Q(r)_{Y_0}-Q(p)_{Y_0} Q(r)_X)=2[p,r]\langle p,r\rangle_Q$。由 $dV\,dh/(Vh)=dx\,dy/(2xy)$ 及 $V^{w+1}h/y=V^w$，边积分中每一项给 $3p_qc_q/2$ 乘 $\mathcal J(\Gamma(m_q)E_{m_q},\Gamma(b_q)E_{b_q})$。可先取 $\Re w$ 足够大，使所有 Epstein 积分绝对收敛；常数格点的导数为零。再沿 $w$ 延拓。由此得到 (4) 的源等式及其中的绝对常数。

两边都是奇的模函数。模不变、光滑、多项式增长的奇齐次差具有零常数模，其余模在尖点仅有衰减 $K$ 支。由固定高度处光滑 Fourier 数据的衰减，向更高处的 $K$ 模级数一致指数衰减。故差在 $L^2$ 中，$\Delta\le0$、$s(s-1)>0$ 给出差为零。这证明了 (4)。

## 每个尖形式的整体消去

有
$$
 c_q=\frac{\Gamma((q+1)/2)\Gamma((2l-q+1)/2)}
 {\sqrt\pi\,\Gamma(l+1/2)} .
$$
式 (1) 中 $A=q+1,\ B=2s-q-1$，所以两个依赖 $q$ 的 Gamma 因子恰被 $c_q$ 消去。再用 $\Lambda_f(B)=(-1)^s\Lambda_f(q+1)$ 及 $\sin(\pi m_q)=(-1)^{q/2}\sin(\pi(w+l)/2)$，(4) 中 $H_f^-$ 的系数精确等于
$$
 \frac{3(-1)^s\,2^{5-4s}\pi^{3/2}\Gamma(2s-1)\Lambda_f(N)}
 {\Gamma(l+1/2)\sin(\pi(w+l)/2)\Gamma(N/2)\Gamma(v/2)\langle f,f\rangle}
 \sum_{q\ {\rm even}}p_q(-1)^{q/2}\Lambda_f(q+1). \tag{5}
$$
在目标整数 $w$ 上，此正弦为 $\pm1$。对归一化一级 Hecke 本征形式，$f(it)$ 为实数。故 (5) 的有限和等于 $\Im\int_0^{i\infty} f(z)P(z)\,dz$。由 (2)，微分 $f(z)P(z)\,dz$ 在 $U$ 下不变。$U$ 把有向边 $0\to i\infty$ 依次送到 $i\infty\to1$ 和 $1\to0$。三条边的积分相等，且因全纯其有向和为零。尖形式在三个尖点的衰减允许取理想三角形的极限。因此该积分为零。这里的消去对每个 $f$ 成立，不用假设不同本征形式的周期线性无关。

## 两个非零单项的例子

取 $s=8,w=6$，令
$$
 P=z(z-1)(z+1)(2z-1)(z-2)\;z^2(z-1)^2(z^2-z+1).
$$
直接代换给出 (2)，偶次项的系数依 $q=4,6,8,10$ 为 $11,21,-21,-11$。由 (4) 与源的反对称性，提升为
$$
 L_P^6=-\frac{945\pi}{4576}\,[11F(8;3/2,9/2)+5F(8;5/2,7/2)] .
$$
在标准化到方括号后的组合中，$\Delta_{16}=E_4^{\rm hol}\Delta_{12}$ 的两项分别有 $K_8(3/2,9/2)=-208/(945\pi)$、$K_8(5/2,7/2)=208/(225\pi)$。式 (5) 给 $11\Lambda_f(5)=21\Lambda_f(7)$，故两项乘共同的 $\Lambda_f(13)/\langle f,f\rangle$ 后精确相抵。这里 $f(it)>0$，故三个实 Mellin 值都非零。提升本身也非零：源的首模在无穷端的 $Y^{11/2}e^{-Y/2}$ 项来自 $9/2$ 零腿的增长支和 $3/2$ 的首模，系数非零，另一对的最高幂较低。

例如将这个例子中的第一因子 $z(z-1)(z+1)(2z-1)(z-2)$ 乘上 $[z^2(z-1)^2]^i(z^2-z+1)^j$，取 $l=3+3i+j$，便得到一族 (2)。这些多项式的齐次次数按 $2l$ 计，例如上面的 $P$ 次数为 $11$，齐次化后是 $14$ 次。

## 网页 Neon 独立复审

网页 Neon 复核了提交 c7805b6 的命题 2。报告核对了极化、格点核边项的 $3/2$、Poisson 界、齐次差的 $L^2$ 唯一性，以及 Manin 三角形的逐尖形式消去，未指出实质缺口。$\Delta_{16}$ 的两个单项系数和精确关系 $11\Lambda(5)=21\Lambda(7)$ 也获独立复算。正文已补明三条有向边的轮换和齐次次数。

[后续命题](half-odd-span.md)给出临界带内消去的逆向判据、基和维数。
