# 临界带内奇 theta 提升的消去空间

**命题 3。** 取整数 $s\ge3$、$3\le w\le s$、$w\equiv s\pmod 2$。用[命题 1](half-odd-general.md)的参考特解，考虑所有 $m+b=w$、$m,b\in\frac12+\mathbb Z_{\ge1}$ 的奇源 $F(s;m,b)$ 的形式系数向量，每对不同的腿按 $m>b$ 定序。逐个 Hecke 尖形式的补全系数同时为零，当且仅当该向量属于下列 [命题 2](half-odd-theta.md) 的提升系数向量所张成的空间：
$$
 P_j(z)=J(z)D(z)^jQ(z)^{s-4-3j},\qquad
 \frac{s-w}{2}\le j\le\left\lfloor\frac{s-4}{3}\right\rfloor ,
$$
其中
$$
 J(z)=z(z-1)(z+1)(2z-1)(z-2),\quad
 D(z)=z^2(z-1)^2,\quad Q(z)=z^2-z+1 .
$$
这些向量构成一组基，故空间维数为
$$
 \max\left(0,\left\lfloor\frac{s-1}{3}\right\rfloor-\frac{s-w}{2}\right).
$$
将提升除以 $\pi$ 后可取有理系数基。这里在系数空间内分类，所需的函数等式由命题 2 给出。

## 三角形对称多项式

记 $l=s-1$，齐次化为 $\widetilde P(X,Y)=Y^{2l}P(X/Y)$。命题 2 的条件 (2) 等价于 $\widetilde P\circ U=\widetilde P$ 和 $\widetilde P\circ R=-\widetilde P$，其中 $U(X,Y)=(X-Y,X)$，$R(X,Y)=(Y-X,Y)$。因 $U^3=-I$ 且 $RUR=U^{-1}$，$\widetilde P$ 在六个反射 $U^iR$ 的固定直线上均消失。直线的定义因子为 $X,Y,X-Y,X+Y,2X-Y,X-2Y$。乘积 $\widetilde J$ 在 $U$ 下不变，在 $R$ 下变号。于是 $\widetilde P/\widetilde J$ 是次数 $2l-6$ 的齐次不变量。当 $l<3$ 时只有零多项式。

用 $u=X+\omega Y$、$v=X+\bar\omega Y$，其中 $\omega=\exp(2\pi i/3)$。$U$ 在这两个坐标上乘互逆的六次单位根，$R$ 将其送为 $(-v,-u)$。偶次齐次不变量由 $uv$ 和 $u^6+v^6$ 生成。具体地，允许的单项式满足指数差为 $6$ 的倍数，$R$ 将它们对称配对，所得 $u^{6i}+v^{6i}$ 可用 $u^6+v^6$ 和 $(uv)^6$ 递推。又有 $uv=\widetilde Q$ 和 $u^6+v^6=2\widetilde Q^3-27\widetilde D$，其中 $\widetilde D=X^2Y^2(X-Y)^2$。因此全部满足 (2) 的多项式由 $P_j$ 张成，此时 $0\le j\le\lfloor(l-3)/3\rfloor$。

由 $RU(X,Y)=(Y,X)$ 得 $p_{2l-q}=-p_q$。对 $P_j$，写 $a=l-3-3j$，最低的非零偶次项为
$$
 [z^{2j+2}]P_j=5+2a+4j.
 \tag{6}
$$
这是由 $J=-2z+5z^2+O(z^3)$、$D^j=z^{2j}(1-2jz+O(z^2))$、$Q^a=1-az+O(z^2)$ 直接相乘得到。于是偶次投影在线性组合上是单射，基向量数为 $\lfloor l/3\rfloor$。

## 由周期矩阵得到必要性

先取 $w=s$。有序腿与偶数 $q$ 一一对应：
$$
 2\le q<l,\qquad m_q=(w+l-q)/2,\qquad b_q=(w-l+q)/2 .
$$
给一个形式系数向量，对其 $q$ 坐标除以 $3c_q\Gamma(m_q)\Gamma(b_q)$，记结果为 $p_q$。若该向量来自 (4)，按 $p_{2l-q}=-p_q$ 和源的反对称性配对，恰好恢复多项式的偶次系数。对任意形式系数向量套用命题 1，式 (5) 中相同的 Gamma 化简表明，$H_f^-$ 的系数是一个非零因子乘
$$
 \sum_{\substack{2\le q<l\\q\ {\rm even}}}p_q(-1)^{q/2}\Lambda_f(q+1).
 \tag{7}
$$
这里的非零因子是式 (5) 的公共前因子乘 $2$。配对的两项同号，因为 $\Lambda_f(2s-q-1)=(-1)^s\Lambda_f(q+1)$，且 $l=s-1$。公共因子的 $\Lambda_f(N)$ 在 $N=s+w-1\ge s+2$ 处由绝对收敛 Euler 乘积保证非零，其余因子亦非零。

取所有归一化 Hecke 本征形式作为 (7) 的行。若有行关系，则相应的尖形式线性组合的奇周期 $\Lambda(q+1)$ 在所有偶数 $2\le q<l$ 处同时为零。由函数方程，互补的偶 $q>l$ 处亦然。若 $q=l$ 为偶数，$s$ 为奇数，中心周期自动为零。这些周期恰是去掉常数项和顶次项的偶周期多项式的系数，差别只有非零二项式因子和符号。受限偶周期多项式映射由 Eichler–Shimura 对应为单射，故该行关系为零。[Baumard–Schneps §1.1](https://arxiv.org/pdf/1109.3786#page=1)

周期矩阵因此具有满行秩 $\dim S_{2s}=\lfloor s/6\rfloor-\mathbf 1_{s\equiv1\bmod6}$。在 $w=s$，列数为 $\lfloor(s-2)/2\rfloor$，核维数为 $\lfloor l/3\rfloor$。命题 2 已将每个 $P_j$ 的偶次投影置于此核内，而 (6) 给出恰好这么多个线性无关向量。所以它们张成全核。

对 $3\le w\le s$，条件 $b_q\ge3/2$ 恰好截为 $q\ge s+2-w$。把缺少的低 $q$ 坐标补零后，(7) 仍是相同的周期矩阵。由 (6) 的三角结构，基展开中 $j<(s-w)/2$ 的系数必须依次为零，剩余基向量全满足支撑条件。转回原系数时 $\Gamma(m_q)\Gamma(b_q)/\pi$ 为非零有理数，故 $L_{P_j}^w/\pi$ 给出有理系数基。这也证明有理消去向量可展开为该基的有理线性组合。

## 例子与待复审范围

$s=8,w=6$ 时只有 $j=1$，故命题 2 中的 $11F(8;3/2,9/2)+5F(8;5/2,7/2)$ 给出唯一的消去方向。$s=8,w=8$ 时有 $j=0,1$ 两个方向。

本命题为新增推导，待独立复审。临界带外可先看 $w=s+2$，此时 $q=0$ 首次对应两条均不低于 $3/2$ 的腿，而三角形反射条件迫使 $p_0=0$。物理方面仍需找到具体的半整数奇源候选并比较其系数。DGW §5.1 提出的是奇 SMF 的研究问题。[DGW §5.1](https://arxiv.org/html/2501.03996v2#S5.SS1)
