# 权重 12 · 从双重 Eisenstein 积分到模变换

取 $(s,m,k)=(6,4,4)$，记 $M_{01}=\mathcal M\!\left[\begin{smallmatrix}0&1\\8&8\end{smallmatrix}\right]$。下文用 DKS 预印本 v3 的编号。期刊版式 (3.9)–(3.11) 的编号相同。[DKS 原文](https://arxiv.org/html/2109.05018v3)

后续已从卷积推出完整 $S$-余循环，并订正 FKL 式 (37) 所印的尖点项符号。本页保留迭代积分与数值复核。[完整证明](s-cocycle-proof.md)

## 1. 从定义展开 $M_{01}$ · 已核验

记 $m_{ab}$ 为 DKS 式 (2.24b) 的原始二重模值，内层积分指标为 $a$，外层为 $b$，两个 Eisenstein 级数均为 $G_8$。记 $m_a$ 为相应的深度一周期。式 (2.50b) 和 shuffle 恒等式给出 $M_{01}=2(m_0m_1-m_{01})=2m_{10}$。因此实际需要展开的积分内层带 $z_1$，外层不带 $z_2$。

在 $S$ 的不动点 $i$ 使用式 (2.25b)。将 $\mathcal G_{ab}(i)$ 记为从 $i$ 出发的正则化积分，将 $\mathcal G_a(i)$ 记为深度一积分，则有

$$
M_{01}=2\bigl(\mathcal G_{10}(i)+\mathcal G_{56}(i)-\mathcal G_6(i)m_1\bigr),
\qquad m_1=\mathcal G_1(i)+\mathcal G_5(i).
$$

令 $\lambda=2\pi$、$c=2\zeta(8)$、$B=2\lambda^8/7!$，设 $s_r=\sigma_7(r)$。沿虚轴有 $G_8(it)=c+B\sum_{r\ge1}s_r e^{-\lambda rt}$。记

$$
J_j(v)=\int_1^\infty t^j e^{-\lambda v t}\,dt
=e^{-\lambda v}\sum_{h=0}^{j}\frac{j!}{(j-h)!(\lambda v)^{h+1}} .
$$

代入深度一和深度二积分。切向基点正则化给常数项的内层积分赋值 $-c\,t^{a+1}/(a+1)$。从 1 积到无穷时，纯多项式部分取原函数在上端的有限部分。

$$
\begin{aligned}
i^{-a-1}\mathcal G_a(i)&=-\frac{c}{a+1}+B\sum_{r\ge1}s_r J_a(r),\\
i^{-a-b-2}\mathcal G_{ab}(i)
&=\frac{c^2}{(a+1)(a+b+2)}
 -\frac{cB}{a+1}\sum_{r\ge1}s_r J_{a+b+1}(r)\\
&\quad+B\sum_{p\ge1}s_p
  \sum_{h=0}^{a}\frac{a!}{(a-h)!(\lambda p)^{h+1}}
   \left[cJ_{a-h+b}(p)+B\sum_{r\ge1}s_r J_{a-h+b}(p+r)\right].
\end{aligned}
$$

这一步同时保留了双零模、单零模和双正频项。下限为 1 的指数衰减允许逐项积分。直接把下限改为 0 会遇到尖点发散。

## 2. 对照 DKS 式 (3.9) · 已核验

取式 (3.9) 中 $(\tau-\bar\tau)^{-5}$ 后常数单项式的系数。利用 $y=\pi\Im\tau$、$p_{12}^+(0,0)=1$，从文献给出的两种写法解出

$$
M_{01}
=-\frac{4i\pi^9}{416745}\zeta(7)
+\frac{2860i\pi^3}{304731}\zeta(13)
-\frac{256i\pi^{14}}{22854825}\Lambda(\Delta_{12},13).
$$

这里取的是多项式系数，未在边界点 $\tau=0$ 求值。数值复核分别用上一节的积分级数和 Ramanujan $\Delta$ 的乘积展开。完成的 $L$ 值在 $t=1$ 拆分 Mellin 积分，用 $\Delta(i/t)=t^{12}\Delta(it)$ 得
$\Lambda(13)=\int_1^\infty(t^{12}+t^{-2})\Delta(it)\,dt$。

| 每层正频截断 | $M_{01}/i$ | 与式 (3.9) 反解的绝对差 |
|---:|---:|---:|
| 12 | $-2.0493122666895251592150634632$ | $1.54\times10^{-27}$ 以内 |
| 22 | $-2.0493122666895251592150634647$ | $2.38\times10^{-53}$ 以内 |

差值来自 70 位计算，属于数值复核。公式的精确来源是 DKS 式 (3.9)，计算未给出它的独立解析证明。[DKS 式 (3.9)–(3.10)](https://arxiv.org/html/2109.05018v3#S3.SS1)

## 3. 对照 FKL 的核 · 已核验

记 $C_n$ 为 FKL 式 (36) 按超几何变形正则化后的值，$P(x)=P_5^{(7,-7)}(2x-1)$。采用重新计算 Petersson 迹公式所得的尖点项符号，有限和取 $1\le p<n$。

$$
C_n=\frac13\sum_{p=1}^{n-1}(1-p/n)^7 P(p/n)
         \sigma_{-7}(p)\sigma_{-7}(n-p)
+\left(\frac{11}{10 n^7}-\frac{16}{45}\pi^6\zeta'(-6)\right)\sigma_{-7}(n)
-\frac{21}{512}\,A\,\frac{\tau(n)}{n^{12}},
$$

其中 $A=\Lambda(6)\Lambda(13)/\langle\Delta,\Delta\rangle$。双正频部分的比较可取固定的下限 $T>0$。在 $\mathcal G_{10}(iT)$ 的双正频项中，令 $p$ 为内层频率、$r$ 为外层频率、$n=p+r$。去掉公共因子 $-iB^2 e^{-\lambda nT}$ 和 $\sigma_7(p)\sigma_7(r)$ 后，积分核为

$$
\frac{T}{\lambda^2 p n}+\frac{1}{\lambda^3}
   \left(\frac{1}{p n^2}+\frac{1}{p^2n}\right).
$$

比较该核的 $T$ 无关部分。FKL 有限和用 $\sigma_{-7}(p)=p^{-7}\sigma_7(p)$ 改写后，核为 $P(p/n)/(3n^7p^7)$。把两个方向 $(p,r)$ 与 $(r,p)$ 相加，去掉 DKS 核中的公共因子 $\lambda^{-3}$，精确计算得

| 固定 $n=4$ | FKL 核之和 / DKS 核之和 |
|---|---:|
| $(p,r)=(1,3)$ | $2121571/1656225792$ |
| $(p,r)=(2,2)$ | $77/524288$ |

两个比例不同，排除了对这两个孤立核逐项做常数缩放。DKS 完整的 $M_{01}$ 还包含反射得到的 $\mathcal G_{56}$ 及深度一乘积。FKL 式 (36) 的正正贡献交换 $p,r$ 后相消，其原始带符号和来自一正一负的频率对，正正有限项出现在正则化结果中。比较完整 DKS 特解的小高度系数时，正正块恰好给出 Jacobi 核。[完整比较](s-cocycle-proof.md) · [FKL 定理 B，式 (32)–(37)](https://arxiv.org/html/2512.21413v1#S1.SS2)

## 4. 540 在尖点项上的归一化 · 已核验

从 $C_n$ 中减去上节有限和与除数和项，称有符号的剩余系数为 $R_n$。订正后的卷积公式给出

$$
R_n=R_1\frac{\tau(n)}{n^{12}},\qquad
R_1=-\frac{252}{691}\frac{\Lambda(13)}{\Lambda(11)}.
$$

第二式用 DKS 的周期比和 DGW 的 Petersson 内积式换算。令 $D=q\,d/dq$，则 $D^{12}\sum_{n\ge1}R_n q^n=R_1\Delta_{12}$。把这个尖点形式放入 DKS 式 (2.31) 的偶齐次积分，得到 $R_1 H_\Delta^+$。换到 DKS 式 (3.11) 的补全项需再乘 $-1/540$，因为

$$
\lambda_\Delta^{\rm DKS}
=\frac{7}{10365}\frac{\Lambda(13)}{\Lambda(11)}
=-\frac{R_1}{540}.
$$

这个标量来自卷积与 DKS 齐次核在小高度的系数之比。已发表的 DKS 补全系数可用于对照。[DKS 齐次积分](https://arxiv.org/html/2109.05018v3#S2.SS2) · [标量推导](s-cocycle-proof.md)

## 5. 与完整余循环的关系

完整推导先从 FKL 的混合符号正则化得到各模小高度系数，再计算尖点形式的周期，最终得到整个余循环。将所得余循环代回 DKS 中含 $M_{01}$ 的表达式，可精确求出这里通过双重积分复核过的模值。一般权重的统一转换仍待建立。[完整证明](s-cocycle-proof.md)
