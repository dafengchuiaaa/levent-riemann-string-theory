# 临界带外第一档的端点周期与格点积分

取基线 5831b20，固定 $s\ge3,\ w=s+2,\ l=s-1$。$F(s;m,b)$、参考特解 $F_{\rm p}$、$H_f^-$ 与完成 $L$ 值均取[命题 1](half-odd-general.md)的归一化。与[命题 3](half-odd-span.md)一样，系数属于两腿均至少为 $3/2$ 的形式向量，按 $m>b$ 定序；等腿的奇源为零。

**命题 4。** 同时消去每个一级 Hecke 尖形式的全核维数为
$$
 1+\left\lfloor\frac{s-1}{3}\right\rfloor .
$$
[三角形 theta 提升](half-odd-theta.md)恰张成全核中 $q=0$ 系数为零的超平面。余下的一维商有显式有理代表，由下文 (4) 的 Bernoulli 数给出。该代表的函数非零，且有下文 (8) 的四维格点积分表示。后者沿一条测地线取边界值，以横向反射配对处理端点。

## 先算 $s=6,w=8$

写 $F_q=F(s;m_q,b_q)$，取
$$
 m_q=\frac{w+l-q}{2},\qquad b_q=\frac{w-l+q}{2},\qquad
 c_q=\frac{\binom l{q/2}}{\binom{2l}q},\qquad
 C_q=3c_q\Gamma(m_q)\Gamma(b_q)\,p_q . \tag{1}
$$
$C_q$ 是 $\sum C_q F_q$ 中实际的系数，$p_q$ 是周期矩阵的坐标。本例有三列。取 $f=\Delta_{12}$，下表的 $K_6$ 是命题 1 中乘在 $\Lambda_f(11-q)\Lambda_f(13)/\langle f,f\rangle$ 前的常数。

| $q$ | $(m_q,b_q)$ | $3c_q\Gamma(m_q)\Gamma(b_q)/\pi$ | $\pi K_6(m_q,b_q)$ |
|---:|---:|---:|---:|
| 0 | $(13/2,3/2)$ | $31185/128$ | $-128/7203735$ |
| 2 | $(11/2,5/2)$ | $945/128$ | $128/218295$ |
| 4 | $(9/2,7/2)$ | $225/128$ | $-128/51975$ |

写 $r_q=\Lambda_\Delta(q+1)$。在一个共同的非零实周期单位下，$(r_0,r_2,r_4)=(192/691,16/135,8/105)$。这是 Kohnen–Zagier 对 $\Delta$ 的精确周期，亦可从其偶周期多项式按二项式系数读出。函数方程给 $\Lambda_\Delta(11-q)=r_q$。因此矩阵行与 $(r_0,-r_2,r_4)$ 成比例。由 $691r_0=1620r_2$ 和 $9r_2=14r_4$，在 $p$ 坐标有基 $(691,1620,0)$、$(0,9,14)$。转换为实际系数，各自去掉公共因子后得到
$$
 \boxed{7601 F_0+540 F_2},\qquad 27F_2+10F_4 . \tag{2}
$$
第二个方向是三角形多项式 $JQ^2$ 的提升，其最低两个偶系数是 $9,14$；第一个方向代表商空间。三列、一条非零行，故全核恰为二维。这里 $13$ 已超出权 $12$ 的临界带，所有消去都对同一个非零的 $\Lambda_\Delta(13)$ 进行。[Kohnen–Zagier，§1.1、§4.2，尤其印刷页 242 的权 12 算式](https://people.mpim-bonn.mpg.de/zagier/files/scanned/ModFormsRatPeriods/fulltext.pdf)

## 一般权重的周期核

列标正是 $q=0,2,\ldots<l$，共 $\lfloor s/2\rfloor$ 列。特别地，$q=0$ 对应 $(s+\tfrac12,\tfrac32)$。命题 1 中 $A=q+1,\ B=2s-q-1,\ N=2s+1,\ v=-1$。将 (1) 代入其 Gamma 公式，利用 $c_q=\Gamma((q+1)/2)\Gamma((2l-q+1)/2)/(\sqrt\pi\Gamma(l+1/2))$ 及函数方程，消去条件逐行成为
$$
 \sum_{\substack{0\le q<l\\q\ {\rm even}}}
 p_q(-1)^{q/2}\Lambda_f(q+1)=0 . \tag{3}
$$
删去的公共行因子非零：$\Gamma(v/2)=\Gamma(-1/2)=-2\sqrt\pi$，$\sin(\pi(w+l)/2)=(-1)^s$，$\Lambda_f(2s+1)$ 的 Euler 乘积在绝对收敛区非零。这里对 $q=0$ 使用的是命题 1 在 $w=s+2$ 的总模延拓参考解。

删去第零列即得命题 3 在 $w=s$ 的周期矩阵；改变 $w$ 引起的原系数缩放已由 (1) 吸收。该子矩阵满行秩。具体地，两端以外的偶次幂周期若全为零，受限偶周期多项式的单射性给出对应尖形式为零；函数方程补齐另一半，$s$ 为奇数时中点周期本来为零。增列后秩仍为 $d_s=\dim S_{2s}=\lfloor s/6\rfloor-\mathbf1_{s\equiv1\bmod6}$。所以全核维数是 $\lfloor s/2\rfloor-d_s=1+\lfloor l/3\rfloor$。三角形向量均满足 $p_0=0$，其多项式 $JD^jQ^{s-4-3j}$ 按最低偶次数 $2j+2$ 排成三角形，给出 $\lfloor l/3\rfloor$ 个独立向量。它们已填满 $p_0=0$ 的核。这也涵盖 $d_s=0$ 的情形。

**余维一的有理代表。** 令 $B_j$ 为 Bernoulli 数，$B_0=1$，超界的二项式系数取零。对偶数 $q$ 定义
$$
 b_{s,q}= B_{2s}\!\left[1+\binom{2s-1}{q}-\binom{2s-1}{q+1}\right]
  +2\sum_{r=1}^s\binom{2r-1}{q}\binom{2s}{2r}B_{2r}B_{2s-2r}.
$$
Kohnen–Zagier 定理 9 在 $r_q(f)=\int_0^\infty f(it)t^qdt=\Lambda_f(q+1)$ 的约定下给出
$b_{s,2l-q}=-b_{s,q}$ 及 $\sum_{q=0,\,q\ {\rm even}}^{2l}(-1)^{q/2}b_{s,q}r_q(f)=0$。由函数方程，两互补项相等，中项若存在其 $b$ 为零。端点系数 $b_{s,0}=-3(2s-1)B_{2s}\ne0$。故在 (1) 取 $p_q=b_{s,q}/b_{s,0}$ 即得 $p_0=1$ 的代表；在原系数中去掉公共 $\pi$ 后可写成
$$
 \mathcal C_s=\sum_{\substack{0\le q<l\\q\ {\rm even}}}
 a_q F_q,\qquad
 a_q=\frac{3c_q\Gamma(m_q)\Gamma(b_q)}{\pi}\,
       \frac{b_{s,q}}{b_{s,0}}\in\mathbb Q . \tag{4}
$$
它与三角形的有理基一起给出有理全核，也张成实数或复数上的全核。$s=6$ 时 (4) 与 (2) 的第一个向量可相差一个旧方向和非零比例。此处的端点关系是尖形式特有的额外周期关系，单靠 Manin 三角形给不出 $p_0\ne0$。[Kohnen–Zagier，§4.2，定理 9](https://people.mpim-bonn.mpg.de/zagier/files/scanned/ModFormsRatPeriods/fulltext.pdf)

此外，$a_0\ne0$ 也能在函数上检测。写 $Y=4\pi\Im\tau$，作用 $\Delta_\tau-s(s-1)$ 后，源的首正模在尖点含非零 $Y^{m_0+1}e^{-Y/2}$ 首项。它来自 $m_0$ 零腿的 $Y^{m_0}$ 项求导及 $b_0$ 腿的首模，外面还有 $Y^2$。其余各对的最大腿严格小于 $m_0$；首模的严格混合对还多出指数衰减。因此 (4) 的函数也不能等于旧提升的线性组合。

## 格点实现：在虚轴上指定单条边的迹

下面给出独立于周期矩阵的积分恒等式，再将 (4) 的系数放入积分。固定本档的 $w$，写 $\rho=u+ih$，$Q_\tau(p)=|p_1+p_2\tau|^2/\Im\tau$，$[p,r]=p_1r_2-p_2r_1$。取命题 2 的行列式核关于 $V$ 的矩
$$
 T_w(\tau;\rho)=\frac{i\Gamma(w+1)}{\pi^w}
 \sum_{\substack{p,r\in\mathbb Z^2\\[p,r]\ne0}}
 \frac{[p,r]}{\big[Q_\tau(p+ur)/h+hQ_\tau(r)\big]^{w+1}} . \tag{5}
$$
这正是 $\int_0^\infty V^{w+1}\Theta^-(\tau;u,h,V)\,dV/V$，由 Gamma 积分逐项得到。固定内点处四维格点和绝对收敛，因为分子为二次，分母总次数为 $2w+2$ 且 $w>2$。它在 $u$ 上是周期一的奇函数，模 $\tau$ 不变、对 $\tau$ 的反射为奇，且 $\Delta_\tau T_w=\Delta_\rho T_w$。这些性质均从同一个核继承。

对每列取 **正** 指数 $e=l-q=m_q-b_q$，于是 $1\le e\le s-1$。在 $u>0,\ h>0$ 写 $\rho=R e^{i\theta}$、$0<\theta<\pi/2$。定义
$$
 A_e(\rho)=R^e k_e(\theta),\qquad
 k_e(\theta)=
 \frac{(\sin\theta)^s\,{}_2F_1((s-e)/2,(s+e)/2;s+1/2;\sin^2\theta)}
 {{}_2F_1((s-e)/2,(s+e)/2;s+1/2;1)} . \tag{6}
$$
分母有限且为正，因两上参数为正且 Gauss 余量为 $1/2$。Gauss 方程给出 $k_e''+[e^2-s(s-1)\csc^2\theta]k_e=0$。在极坐标 $\Delta_\rho=\sin^2\theta[(R\partial_R)^2+\partial_\theta^2]$，故 $\Delta_\rho A_e=s(s-1)A_e$，$A_e(ih)=h^e$。在 $\theta=0$ 有 $k_e=O(\theta^s),\,k_e'=O(\theta^{s-1})$。在 $\pi/2$，常微分方程是常点，解及其单侧导数有限。

定义积分沿完整周期取极限
$$
 I_e(\tau)=\lim_{L\to\infty,\ L\in\mathbb N}
 \int_0^\infty\int_0^L A_e(u+ih)T_w(\tau;u+ih)\frac{du\,dh}{h^2}.
 \tag{7}
$$
则
$$
 I_{l-q}=-\frac{\Gamma(m_q)\Gamma(b_q)}2 F_q,\qquad
 \boxed{\displaystyle
 \mathcal C_s=-2\sum_q\frac{a_q}{\Gamma(m_q)\Gamma(b_q)}I_{l-q}.} \tag{8}
$$
因此一个格点积分便可实现新方向，在 (7) 的被积式中将 $A_e$ 换成 (8) 的有限线性组合即可。下面给出端点收敛和边界证明。

**端点与横向配对。** 对紧集内的 $\tau$，$T_w$ 及其固定阶 $\tau$ 导数在整个 $\rho$ 平面有界；同样有 $h|\nabla_{\!(u,h)}T_w|\ll1$。在 $h\ge1$ 这些界均可再乘 $O_N(h^{-N})$。理由是在 $\rho$ 的标准基本域，用行列式核对 $p$ 的 Poisson 展开。非零的两矢量 $k,r$ 给出 $\exp[-\pi h Q_\tau^*(k)/V-\pi Vh Q_\tau(r)]$，多项式因子与 $V$ 的矩都可被吸收；零对偶模因行列式消失。借助 $\rho$ 模不变性和双曲梯度范数的不变性即可从基本域得到全局界。对 $|\Re\tau|\le1/2,\ \Im\tau\ge1$，同一估计的常数至多按 $\Im\tau$ 的幂增长，因 $Q_\tau,Q_\tau^*$ 的最小特征值均 $\gg1/\Im\tau$。详见[行列式核的 Poisson 估计](half-odd-theta.md)。

当 $u\ge h$ 时，(6) 给 $|A_e|\ll h^su^{e-s}$，$|\partial_u A_e|\ll h^su^{e-s-1}$，$|\partial_h A_e|\ll h^{s-1}u^{e-s}$。当 $u\le h$ 时函数和一阶导数分别是 $O(h^e),O(h^{e-1})$。原点附近 $\int |A_e|du\,dh/h^2$ 收敛，因为极坐标中径向为 $R^{e-1}$，角向为 $k_e/\sin^2\theta$。在有限横区间的其余下边界用 $O(h^s)$ 即可。对 $e<s-1$，横向尾绝对可积；端点 $e=s-1$ 的单项界只有 $1/u$。在每个 $j\le u\le j+1$ 用 $T_w(j+1-t+ih)=-T_w(j+t+ih)$，将整格积分写成
$$
 \int_0^{1/2}[A_e(j+t+ih)-A_e(j+1-t+ih)]\,
 T_w(\tau;t+ih)\,dt .
$$
当 $j\ge\max(1,h)$，方括号为 $O(h^s j^{e-s-1})$。即使在端点也有 $O(h^s j^{-2})$。乘 $dh/h^2$ 后小 $h$ 可积，大 $h$ 由快速衰减控制。$1\le j<h$ 的至多 $O(h)$ 个周期也由同一大高度界吸收。因此反射配对的格积分和绝对收敛，(7) 存在；同样的估计适用于 $\tau$ 微分并给出目标尖点的多项式增长。

**边界方程与函数的唯一性。** 在 $0<u<L$ 作格林分部积分，先固定正整数 $L$。上边因快速衰减消失。在下边 $h=\varepsilon$，用 $|\partial_h T_w|\ll1/h$ 和上述两段 $A_e$ 界，有边通量绝对值 $O_L(\varepsilon^e[1+|\log\varepsilon|])$，趋于零。两条竖边上 $T_w=0$。右边留下 $A_e(L+ih)\partial_uT_w(L+ih)$，其 $h$ 积分为 $O(L^{e-s})+o(1)$，因为在 $h\le L$ 可用 $h^sL^{e-s}$ 及 $|\partial_u T_w|\ll_N h^{-1}(1+h)^{-N}$，$h>L$ 则快速衰减。故极限中只剩 $u=0$ 的外法向导数。

在这条边上，令 $x=V/h,\ y=Vh$。取 $\theta_x(\tau)=\sum_p e^{-\pi xQ_\tau(p)}$ 和命题 1 的反对称梯度积 $\mathcal J$。对原 Gaussian 求导得到
$-\partial_u\Theta^-|_{u=0}=y^{-1}\mathcal J(\theta_x,\theta_y)$。用 $A_e(ih)=h^e$、$dV\,dh/(Vh)=dx\,dy/(2xy)$，于是
$$
 (\Delta_\tau-s(s-1))I_e=
 \tfrac12\Gamma((w-e)/2)\Gamma((w+e)/2)
 \mathcal J(E_{(w-e)/2},E_{(w+e)/2}).
$$
两指数均 $\ge3/2$，Eisenstein 的 Mellin 积分及其梯度在此可绝对积分。这一边的次序是 $(b_q,m_q)$，正是 (8) 负号的来源。$I_e$ 模不变、为奇并具有多项式增长；与 (8) 右边的齐次差零模为零，各非零模只含衰减的 $K$ 支。因此该差在 $L^2$ 中，$\Delta\le0<s(s-1)$ 给出差为零。至此 (8) 在函数层面成立。

**结论的范围。** (3)–(4) 的维数、余维及逐本征形式的消去仅依赖周期矩阵、参考系数和 Kohnen–Zagier 定理。(5)–(8) 另给出了带明确收敛次序的格点实现，其输入是该周期关系确定的边界权重。三角形多项式的反射条件始终排除 $q=0$。若进一步要求一个几何对称性自行选出 (4) 的 Bernoulli 权重，或要求该代表来自某个指定物理振幅，本文的积分尚未给出这种选择机制。
