# 权重 12 的系数比较

取 $s=6$、$s_1=s_2=4$。尖点形式的模权为 $k=2s=12$，源项指标之和为 $w=s_1+s_2=8$。记 $\Delta=\Delta_{12}=\sum_{n\ge1}\tau(n)q^n$，首项系数为 1。下文使用完成的 $L$ 值及周期因子

$$
\Lambda(t)=(2\pi)^{-t}\Gamma(t)L(t,\Delta),\qquad
A=\frac{\Lambda(6)\Lambda(13)}{\langle\Delta,\Delta\rangle}.
$$

FKL 指 Fedosova 与 Klinger-Logan，DGW 指 Dorigoni、Green 与 Wen，DKS 指 Dorigoni、Kleinschmidt 与 Schlotterer。下列 DGW 编号采用期刊版。[版本与公式索引](../references/README.md)

## 1. 特殊化 FKL 定理 B · 已核验

由 $r_i=2s_i-1$ 及 $k=2d+r_1+r_2+2$ 得 $(r_1,r_2,d)=(7,7,-2)$。此时 $-6\le d\le-2$，也避开式 (31) 的例外值。式 (34) 的两个自变量为 $1+d+r_2=6$ 和 $1+d+r_1+r_2=13$，分母中 $n$ 的指数为 $d+r_1+r_2=12$。权重 12 尖点空间由 $\Delta$ 张成。

对 $n>0$，令 $n_1+n_2=n$、$n_1n_2\ne0$、$x=n_1/n$。式 (32) 涉及以下两个多项式。第一式可由二项式展开逐项复核，第二式取 $t=n/n_1$。

$$
\begin{aligned}
P_5^{(7,-7)}(2x-1)
&=\sum_{j=0}^5 {12\choose j}{-2\choose 5-j}(x-1)^{5-j}x^j\\
&=6+30x+84x^2+168x^3+252x^4+252x^5,\\
P_1^{(11,-7)}(1-2t)&=12-6t.
\end{aligned}
$$

记第一式为 $P(x)$。代入式 (32)，含 $\mathbf1_{0<x<1}$ 的项只出现在有限和中；符号项延伸到区间外。

$$
\frac{\Gamma(13)}{\pi}\mathcal P(n_1,n_2)
=2(1-x)^7P(x)\mathbf1_{0<x<1}
-6\,\mathrm{sgn}(n_1)(1-2x).
$$

在以 $\sigma_{-7}(n_1)\sigma_{-7}(n_2)$ 加权的求和中移出有限项，剩余的带符号卷积即 FKL 式 (36)。这里采用定理 B 的正则化。下面只提取式 (37) 的尖点项，完整恒等式还含有限和及除数和项。[FKL 定理 B 与算例](https://arxiv.org/html/2512.21413v1#S1.SS2)

## 2. 同一周期因子前的系数 · 已核验

FKL 式 (37) 的尖点项为

$$
\frac{\Gamma(11)\Gamma(13)}{3\cdot2^{30}\pi^{19}\Gamma(6)}
\frac{L(6,\Delta)L(13,\Delta)}{\langle\Delta,\Delta\rangle}
\frac{\tau(n)}{n^{12}}
=\frac{21}{512}\,A\,\frac{\tau(n)}{n^{12}}.
$$

完成两个 $L$ 值后，$\pi^{19}$ 抵消，有理系数化成 $10!/(3\cdot2^{11}(5!)^2)=21/512$。

DGW 期刊版式 (4.22) 给出 $\mathcal E(6;4,4)$ 中乘在 $H_\Delta$ 前的系数。由 DKS 式 (2.42)，$\Lambda(6)/\Lambda(10)=5/12$。因此

$$
\lambda_\Delta
=\frac{35}{6\cdot184320}\,
\frac{\Lambda(10)\Lambda(13)}{\langle\Delta,\Delta\rangle}
=\frac{35}{6\cdot184320}\frac{12}{5}\,A
=\frac{7}{92160}\,A .
$$

取 **FKL 除以 DGW**，两处有理系数之比为

$$
\frac{21/512}{7/92160}=540 .
$$

FKL 一侧的系数乘 $\tau(n)/n^{12}$，DGW 一侧的系数乘函数 $H_\Delta$。这里核验的是提取共同周期后两个数的比值。[DGW 预印本 v2 的对应公式 (4.23)](https://arxiv.org/html/2501.03996v2#S4.SS2)

## 3. 从 DGW 换回 DKS 的周期比 · 已核验

DGW 期刊版式 (4.19) 给出 $\langle\Delta,\Delta\rangle=(691/7680)\Lambda(8)\Lambda(11)$。DKS 式 (2.42) 给出 $\Lambda(8)=(25/48)\Lambda(10)$。把两式代入上一节首个 $\lambda_\Delta$ 表达式可得

$$
\lambda_\Delta
=\left(\frac{35}{6\cdot184320}\frac{7680}{691}\frac{48}{25}\right)
\frac{\Lambda(13)}{\Lambda(11)}
=\frac{7}{10365}\frac{\Lambda(13)}{\Lambda(11)}.
$$

这正是 DKS 式 (3.11) 中 $(4,4)$ 一行的正号系数。另从 DKS 式 (3.9) 和 (2.43a) 比较同一多项式的两项，抵消所需的有理因子为 $60\cdot7/621900=7/10365$。本节复核已发表公式间的换算，采用论文给出的周期恒等式。[DKS 期刊版，式 (2.42)、(2.43a) 与 (3.9)–(3.11)](https://doi.org/10.1007/JHEP01(2022)134)

## 4. 从卷积到模变换 · 待证明

DKS 已给出这个算例的模不变补全。若要从 FKL 独立导出该结果，需要把正则化卷积识别为 DKS 同一特解的 Fourier 数据，控制重组级数及正则化产生的附加项，并在 $S:\tau\mapsto-1/\tau$ 下得到式 (3.9) 的完整变换差。上述两次分数换算只确定系数，尚未建立这一步。[一般指标的覆盖问题](../general-conjecture/README.md)

本算例的后续计算从 DKS 的双重 Eisenstein 积分入手，给出完整卷积、两侧核的比较及 540 对尖点项的作用。[DKS 多重模值续算](dks-multiple-modular-value.md)
