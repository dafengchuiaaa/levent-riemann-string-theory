# 半整数奇源的全模公式

取 DGW (5.1) 的源次序 $(m,b)$，$w=m+b$，$Y=4\pi\Im\tau$。记 $\mathcal J(f,g)=(\nabla f\,\bar\nabla g-\nabla g\,\bar\nabla f)/(2(\Im\tau)^2)$，$F=(\Delta-s(s-1))^{-1}\mathcal J(E_m,E_b)$。用 DKS II (2.31) 的 $H_f^-$，正模为正，$a_f(1)=1$。以下给出一种确定参考特解的频率延拓。参考解的选择影响单个齐次系数。[DGW §5.1](https://arxiv.org/html/2501.03996v2#S5.SS1) · [DKS II](https://arxiv.org/abs/2109.05018)

**命题 1。** 设 $s\ge3$ 为整数，$m,b\in\frac12+\mathbb Z_{\ge0}$，$s\equiv w\bmod2$。置
$$
 A=s+b-m,\qquad B=s+m-b,\qquad N=s+w-1,\qquad v=s-w+1 .
$$
以下定义的 $F_{\rm p}$ 在这些参数处正则，且
$$
 F-F_{\rm p}=\sum_{f\in\mathcal B_{2s}}
 K_s(m,b)\frac{\Lambda_f(B)\Lambda_f(N)}{\langle f,f\rangle}\,H_f^- ,
\qquad
 K_s(m,b)=
 \frac{2^{6-4s}\pi^2\Gamma(2s-1)}
 {\sin(\pi m)\Gamma(m)\Gamma(b)\Gamma(N/2)\Gamma(v/2)
 \Gamma(A/2)\Gamma(B/2)} . \tag{1}
$$
$\mathcal B_{2s}$ 是归一化 Hecke 本征基。在 $w\le s$ 时，$1\le A,B,N\le2s-1$，两个 $L$ 值都是奇临界值。$w>s$ 时 (1) 仍给延拓参考解的系数，$N$ 已超出临界带。交换两腿时用 $\Lambda_f(2s-B)=(-1)^s\Lambda_f(B)$ 即得反号。等腿的源为零，许可的 $s$ 为奇数，中心周期为零。DGW §5.1 只提出了半整数腿和奇偶条件。这里的临界带由上述不等式确定。

## 参考解与带符号的频率延拓

先在实开集 $m,b>1,\ m+b<s$ 计算，避开整数腿及下文有限部的共振点，再作亚纯延拓。写 $M=m-1/2,\ r_2=2m-1,\ r_1=2b-1,\ e=(-1)^s\sin\pi m$。每对省略 $\sigma_{-r_i}$。单腿是 $F_j(p,Y)=2p^{j-1/2}\sqrt Y K_{j-1/2}(pY/2)/(\sqrt\pi\Gamma(j))$。频率 $f_m+f_b=n>0$ 的源为 $-Y^2(f_mF_m F'_b-f_bF_b F'_m)/2$。衰减解的迹指 $\lim_{Y\to0}Y^{s-1}g(Y)$。两个解相差衰减齐次支时，迹唯一确定其差。

* 小权标签 $m$ 带负频 $-R$、$b$ 带正频 $P$ 时，取在无穷端无较慢齐次尾的 Green 解。这里「小权」只是标签，无须 $m<b$。反向混合对先取同样的 Green 解。
* 正正对从前一支出发。固定 $P>0$，从 $0<R<P$ 延拓到 $R=-R_0$ 的两岸，$n=P-R$。延拓对象是局部解模掉零点正则齐次支 $Y^s(1+O(Y^2))$ 的等价类。若 $G_+,G_-$ 依次对应 $n/P$ 的上、下岸，取 $(e^{i\pi M}G_++e^{-i\pi M}G_-)/(2\cos\pi M)$，再选取无增长支的代表。每个等价类由源及其迹指定。从 $W'=\psi_n\phi/Y^2$ 在零点积分，$W(0)$ 由迹确定。起始开集上 $\phi=O(Y^{3-w})$，故 $\psi_n\phi/Y^2=O(Y^v)$ 且 $\Re v>1$，零端积分沿避开 $R=0,n=0$ 的路径在紧集上局部一致收敛。再从固定正高度 $Y_*>0$ 积分 $(g/\psi_n)'=W/\psi_n^2$，其中 $\psi_n=\sqrt{nY}I_{s-1/2}(nY/2)$。第二次积分的常数只改变被取商的 $\psi_n$ 支。由于 $n/P=1-R/P$，上、下岸 $n/P\pm i0$ 分别对应 $R=R_0e^{\mp i\pi}$。切到负实轴时，$K_M(xe^{\pm i\pi})=e^{\mp i\pi M}K_M(x)\mp i\pi I_M(x)$。连同 $R^M$，上述权重恰消掉 $I_M$，留下 $F_m(R_0,Y)$，带符号的频率 $-R$ 同时变为 $R_0$。因此所得源就是实际的正正源。它指数衰减，变参数公式允许唯一去掉增长的齐次支。
* 零腿先用无穷端 Wronskian 的零有限部，按原高度 $Y$ 扣除发散幂。[首模基准](half-odd-benchmark.md)给出这一定义和 Mellin 积分。为了同时定义二复变量的延拓，对大权标签 $b$ 为负频、或为零频的块，将上述初始迹乘 $\eta=-\sin(\pi b)/e$，做法是加上相应倍数的单位迹衰减齐次支，源保持原值。其余块不动。这个离开物理参数面的约定使两种余弦相位各自延拓；在固定合宜整数差 $b-m$ 的整条共同平移线上，$\eta=1$。在目标半整数点，正正的两岸各占一半。

先合并同一 $n$ 的全部块，再延拓总模。此取序也是对边界 $s=w$、$m=1/2$ 及非绝对收敛区的定义。它不要求在极点逐腿取常数项。在绝对收敛物理区该定义直接给出未变形的参考解；在原基准 $(6;3/2,5/2)$，由[有限乘积的迹](half-odd-positive.md)及衰减唯一性，正正块也逐对等于先前的有限乘积。一般参数的构造不依赖有限 Laurent–Bessel 表示的存在。

## 混合迹、零腿与迹公式

令 $U\sim{\rm Beta}(A,B)$、$D=C\Gamma(2s)/(2e)$，其中
$$
 C=-\frac{2^{2s}\Gamma((s+w)/2)\Gamma((B+1)/2)
 \Gamma((A+1)/2)\Gamma((v+1)/2)}
 {\pi\Gamma(m)\Gamma(b)\Gamma(s)(2s-1)} .
$$
[单对 Green 恒等式](odd-proof.md)第 1 节在上述收敛开集直接给出正向混合迹 $C\,\mathbb E[(1-U)P+UR]^{-v}$，反向在调整前取负号。其证明将归一化的三 Bessel 积分化为 Gauss 方程，在 $P=R$ 用两 $K$ 的 Mellin 矩固定 $C$。所以正正迹直接取 $CP^{-v}{}_2F_1(v,A;2s;n/P)$ 的上述加权割线值。令 $G_\pm={}_2F_1(v,A;2s;n/P\pm i0)$。有
$$
 e\frac{e^{i\pi M}G_++e^{-i\pi M}G_-}{\cos\pi M}
 =(-1)^s[\cos(\pi M)(G_++G_-)+i\sin(\pi M)(G_+-G_-)] .
$$
反向混合迹乘上 $\eta$ 后的相位为 $\sin(\pi b)$。与 FKL (8) 比较，三种非零腿的共同乘子正是 $Dn^{-v}$。固定 $n$ 时混合尾由 $O_n(h^{-\Re v})$ 支配。在实起始区，$F'_j<0$ 对 $j>1$ 成立，每支的绝对 Green 矩即上述矩的绝对值，利用 Wronskian 的二次积分可将 $Y^{s-1}g$ 一致控制于该矩的常数倍。额外的单位迹齐次支也满足此界。两除数权在该区有界，故可取总迹。

实起始区的定模 Green 界也排除了求和后产生增长支的可能。具体地，写 $J_h=\int_0^\infty|\psi_n\phi_h|\,dY/Y^2$，则未调整的每个混合块满足 $|g_h(Y)|\le J_h\,\psi_n(Y)\int_Y^\infty\psi_n(t)^{-2}dt$。这来自先从无穷端积分 Wronskian，再积分 $g_h/\psi_n$。右侧公共高度因子在无穷端为 $O(e^{-nY/2})$，在零端为 $O(Y^{1-s})$。加上的单位迹衰减支也有此界。固定 $Y>0$ 时原源的频率尾局部一致指数收敛，故定模和满足相应微分方程（也可先按分布取和再用常微分方程的正则性）。

对 $n=1$ 写 $\omega=(2s-1)/(2^{2s-1}\Gamma(s+1/2))$。零腿 $Y^t$、另一腿 $j$ 的初始迹（不计次序符号）等于
$$
 B_s(t,j)=\frac{t\,2^{2t}}{\omega\sqrt\pi\Gamma(j)}
 \frac{\Gamma((s+j+t)/2)\Gamma((s-j+t+1)/2)\Gamma(-t)}
 {\Gamma((s+1-t-j)/2)\Gamma((s+j-t)/2)}.
$$
这里按 $I K$ 的 Mellin 矩从 $t<0$ 延拓。记 $a_m,b_m$ 为 $E_m$ 在 $Y^m,Y^{1-m}$ 的系数。分别有
$$
\begin{aligned}
 a_m B_s(m,b)&=D(-1)^s
 \frac{\zeta(2m)}{(2\pi)^{r_2}\Gamma(v)\Gamma(A)},\\
 b_m B_s(1-m,b)&=D(-1)^s
 \frac{(2\pi)^{r_2}\zeta(1-r_2)}{\Gamma(N)\Gamma(B)}.
\end{aligned}
$$
第一行代 $a_m=2^{1-2m}\pi^{-2m}\zeta(2m)$，用 $m\Gamma(-m)\Gamma(m)=-\pi/\sin(\pi m)$，对 $v,A$ 用二倍公式。第二行代 $b_m=2^{2m-1}\Gamma(m-1/2)\zeta(r_2)/(\sqrt\pi\Gamma(m))$，对 $N,B$ 用二倍公式及 $\zeta(1-r)=2(2\pi)^{-r}\Gamma(r)\cos(\pi r/2)\zeta(r)$。大权零腿交换 $m,b$ 再乘 $-\eta$。对任意 $n$，每个 $Y^t$ 项另乘 $n^{j-t-s}\sigma_{1-2j}(n)$。故两零腿恰给 $Dn^{-v}[(-1)^s Z_d^{(r_1,r_2)}\sigma_{-r_1}+ Z_d^{(r_2,r_1)}\sigma_{-r_2}]$，其中 $d=v-1$，$Z$ 含 FKL (9) 的 $n$ 因子。

应用 [FKL 定理 A](https://arxiv.org/html/2512.21413v1#S1) 及[本档案从 Petersson 公式推导的迹号](../weight-12/s-cocycle-proof.md)。两个边界抵消，每个 $f$ 的剩余迹为
$$
 D\frac{(-1)^s2^{2-2s}\pi\Gamma(2s-1)}
 {\Gamma(v)\Gamma(A)\Gamma(B)\Gamma(N)}
 \frac{\Lambda_f(B)\Lambda_f(N)}{\langle f,f\rangle}
 \frac{a_f(n)}{n^{2s-1}} .
$$
源是光滑且在尖点指数衰减的奇模函数，模曲面上的热预解式给唯一有界解，其定模小高度迹为零。负的剩余迹除以 $H_f^-$ 的共同齐次归一化 $\Gamma(2s-1)/(2\Gamma(s))$（保留 $a_f(n)/n^{2s-1}$）即得 (1)，四次 Gamma 二倍公式给出其简式。负模由奇性确定，零模恒为零。

最后说明延拓的存在性。开集上逐模恒等式已成立，无须把割线积分或逐对无限和延拓到发散点。等式右侧只有有限多个全纯尖形式的完整 $L$ 函数、Gamma 因子及固定 Eichler 级数。$F$ 关于两腿也亚纯：Eisenstein 级数在紧集亚纯，尖点处梯度的反对称乘积消掉两常数项之间的乘积，其余各项在紧参数集上连同有限阶微分指数衰减，乘固定多项式界即可。热预解式遂保持此参数依赖。故 $F$ 减去该有限和给出开集上已独立构造的总特解的唯一亚纯延拓。所有目标点 $\sin(\pi m)=\pm1$，其余分母均为整函数 $1/\Gamma$ 的取值，$E_{1/2}$ 取两常数支合并的正则值。由此得到总模的正则性与函数恒等式。此证明使用 FKL 的 $2s\ge6$。原文没有规定奇 SMF 的物理选择原则；其他规定齐次分拆的方式需要另行比较。

## 网页 Neon 独立复审

网页 Neon 复核了提交 c7805b6 的命题 1。报告核对了 Gamma 归一化、上下岸的 $I_M$ 消去、零腿配对及目标点正则性，未指出实质缺口。权 16 两个源 $(8;3/2,9/2)$、$(8;5/2,7/2)$ 的首模总迹在 $h\le400$ 截断下之比与 $-5/11$ 一致，报告的相对差约为 $6\times10^{-13}$。正文已补明绕向、第二次积分的基点和零端收敛指数。
