# 一个 Rallis 核内积引理能否计算 DGW 的 $(6;4,4)$ 系数？

2026-11-25 归档。下文的源码核对及迁移判断出自 Astra 的测试报告。

**判断。** Alpöge–Claude 源码中的 $\,\texttt{GlobalTheory.T.rallis_inner_product_kernel_identity}\,$ 不能直接计算 DKS 特解的 $\Delta_{12}$ 补全系数。两项必要假设在自然代入中失败：正定稳定范围、两条输入的尖性。即使先取尖形式 $\Delta_{12}$ 作测试，结论也是两个提升的内积，还需一条把该内积识别为 DKS 有符号边界迹的独立恒等式。本测试没有替代原证明中的步骤。

## 原说明与实际引理

Alpöge 的说明第 4 页第 8 项给出全定酉群之间的迭代 theta 提升，第 12 项说可用 seesaw 将迭代提升的周期写为与 Hecke 特征无关的 Bessel 周期，并提到 Rallis。该页没有给出 DGW 核或局部常数。这里测试的是发布源码中有明确定义的一条 Rallis **核内积**引理。第 12 项的整条 Bessel 周期结论与这条引理的结论应分别读。[作者原说明，第 4 页](https://alpo.ge/remark.pdf#page=4)

源码定位到发布档案的镜像提交 **7cbf7a0c**。记 $G=\mathrm{GA}(D,W),\,H=\mathrm{GA}(D,V)$，$\mathcal G,\mathcal H$ 为相应有理点作用的基本域，$\omega$ 为源码的 Weil 作用。用 $f_i$ 表示 theta 测试函数，$\phi_i$ 表示自守输入。令
$$
 U_i(h)=\int_{\mathcal G}\overline{\theta_\omega(h,g;f_i)}\,\phi_i(g)\,d\nu_W(g),
 \qquad f_i^{g_i}=\omega(1,g_i)f_i .
$$
引理在给出第一提升被积函数对每个 $h$ 的可积性、以及 $\overline{U_1}U_2$ 的可积性后，断言
$$
 \int_{\mathcal H}\overline{U_1(h)}U_2(h)\,d\nu_V(h)
 =\nu_V(\mathcal H)_{\mathbb R}
 \int_{\mathcal G^2}\overline{\phi_1(g_1)}\phi_2(g_2)
 \sum_{Y\in\mathrm{Lag}(W)}
 \mathrm{lagTerm}_Y(f_1^{g_1},f_2^{g_2})\,
 d\nu_W(g_1)d\nu_W(g_2).
$$
$\mathrm{Lag}(W)$ 是 $W\oplus W$ 对差厄米型的最大各向同性子空间集。若选定框架 $\mathrm{fr}(Y)=(C_1,C_2)$，末项是 $\int_{\mathcal F}f_1^{g_1}(C_1u)\overline{f_2^{g_2}(C_2u)}\,d\mu(u)$。左右两处复共轭的方向及基本域体积均保留。声明在下列第一个文件第 14 行，实际证明在第二个文件第 17–108 行。证明先用 $\texttt{KerAbs}$ 控制绝对积分，第 92 行调用输入的 $\texttt{SWPair}$，第 105 行交换积分。此次核对读取了声明、定义和证明源码，未重新编译 Lean 项目。[Lean 声明](https://github.com/CocoSgt/riemann-hypothesis-lean/blob/7cbf7a0c32fca1badf05931879a37c3f8a2b03b4/src/Development/AbelianVarietyTheory/Proofs/InertiaEigenvalues/RallisInnerProductKernelIdentityStatement.lean#L14-L15) · [Lean 证明](https://github.com/CocoSgt/riemann-hypothesis-lean/blob/7cbf7a0c32fca1badf05931879a37c3f8a2b03b4/src/Development/AbelianVarietyTheory/Proofs/InertiaEigenvalues/RallisInnerProductKernelIdentity.lean#L17-L108) · [各向同性空间与积分项](https://github.com/CocoSgt/riemann-hypothesis-lean/blob/7cbf7a0c32fca1badf05931879a37c3f8a2b03b4/src/Development/AbelianVarietyTheory/Definitions/NonGaloisDescent/IsWeilDatumOOneOfRankZero.lean#L25-L39)

## 假设逐项对照

| Lean 的输入 | 在 DGW 格点核及 DKS 特解上的测试 |
|---|---|
| $D:\texttt{D7}$，加性特征 $\psi$，两个厄米空间 $V,W$，酉特征 $\chi_V,\chi_W$ 及 $\texttt{IsWeilDatumO}$ 的作用 $\omega$ | DGW 采用分裂签名 $(2,2)$ 的 Narain 格，在 $B=0$ 取 Gaussian 核。其正定 Gaussian 随模参数变化；底层格型仍是 $(2,2)$。DGW 未给出与该酉 Weil 数据的等距识别。 |
| $\texttt{V.SigOK}$、$\texttt{V.IsPosDef}$、$2\,\mathrm{rank}(W)<\mathrm{rank}(V)$ | 源码从前两项推出 $V$ 在每个复嵌入的 Gram 矩阵正定。直接把 DGW 的实四维格当成 $V$ 的底层实空间，所需复厄米秩为 2，实签名 $(2,2)$ 要求厄米签名 $(1,1)$。正定性已经失败；若 $W$ 非零，严格秩不等式也失败。这只排除这个直接识别。 |
| 两群的 Haar 测度和有理点基本域；$\texttt{ratL W V}$ 的余体积为一；每个 $Y$ 的框架 $\mathrm{fr}(Y)$；两个 $\texttt{IsThetaFunO}$ 测试函数 | DGW §2 的 $\Gamma'_{2,2}$ 对 $p_1+p_2+p_3=0$ 求和并删去任一 $p_i=0$ 的项。提升对局部 Maass 输入和体积参数做积分。将此删项核写成上述两个测试函数的核配对，需要给出有限处格、测度、删项和无穷处向量的字典。四个格点整数本身不能确定这些归一化。 |
| 无穷型族 $\theta$ 满足 $\texttt{IsInfFam}$；$\phi_1,\phi_2$ 各自满足 $\texttt{IsAutFunSm}$ 和 $\texttt{IsCuspFun}$ | DGW 的源是 $E(4;\tau)^2$。若将两腿放入自守输入槽，$E(4)$ 的常数模含 $2\zeta(8)\tau_2^4/\pi^4$，直接违反尖性。若改放 DKS 的 $\widecheck F$，它有非零 $S$ 余循环，不能作为所需自守输入。$\Delta_{12}$ 本身是尖形式，但即使另建酉群上的实现，把它放进一个槽也没有产生另一个槽所需的 DKS 线性泛函。 |
| $\texttt{KerAbs D}$ 与 $\texttt{SWPair D}$ | 前者给有限体积、带尖形式的 theta 核绝对积分及关于 $h$ 的一致界。后者给 theta 核乘积的 Siegel–Weil 求和等式。两者都是此引理的前提。DGW/DKS 本例所需的混合尾使用移位除数和的有限部，还要保留其留数项；上述绝对积分换序没有给出该有限部。 |

[Lean 的 SWSetup、SWPair、KerAbs](https://github.com/CocoSgt/riemann-hypothesis-lean/blob/7cbf7a0c32fca1badf05931879a37c3f8a2b03b4/src/Development/AbelianVarietyTheory/Definitions/NonGaloisDescent/IsWeilDatumOOneOfRankZero.lean#L41-L73) · [复嵌入处正定性的推论](https://github.com/CocoSgt/riemann-hypothesis-lean/blob/7cbf7a0c32fca1badf05931879a37c3f8a2b03b4/src/Development/AbelianVarietyTheory/Definitions/NonGaloisDescent/IsWeilDatumOOneOfRankZero.lean#L104-L108) · [DGW §2 的格点核](https://arxiv.org/html/2501.03996v2#S2) · [DGW 的 Eisenstein 归一化](https://arxiv.org/html/2501.03996v2#A1)

即便只从 DGW 已有格点提升取数据，式 (4.47) 中本测试对象也是以 $112\,\mathcal E(6;3,5)+75\,\mathcal E(6;4,4)$ 的组合出现，另有 $E(8)$ 项。对组合的尖形式消去只给一条齐次关系，不能固定 $(6;4,4)$ 的单项大小。这里没有使用另一组参数的补全系数。[DGW §4.3](https://arxiv.org/html/2501.03996v2#S4.SS3)

## 直接代入要算出的数

取 $\Delta=\Delta_{12}$、$a_\Delta(1)=1$，$Y=4\pi\Im\tau$，$\Lambda(j)=(2\pi)^{-j}\Gamma(j)L(\Delta,j)$。只取第一正模即可固定标量。以 $n=1$ 的 DKS 边界迹作为待算量，原证明给
$$
 \alpha=\frac{\zeta(8)^2}{13\zeta(16)},\qquad
 r_1=\frac{21}{512}\frac{\Lambda(6)\Lambda(13)}{\langle\Delta,\Delta\rangle},
 \qquad [Y^{-5}][\widecheck F-\alpha E(8)]_1=-28r_1 .
$$
最后一式来自 FKL 的**带负号**的剩余项 $C_1-U_1-V_1=-r_1$ 与 DKS 三块的匹配。模不变解 $F-\alpha E(8)$ 的对应迹为零，$[Y^{-5}][H_\Delta^+]_1=15120$，因此在 $F=\widecheck F+\lambda H_\Delta^+$ 的约定中
$$
 \lambda=+\frac{28r_1}{15120}
 =+\frac{7}{92160}\frac{\Lambda(6)\Lambda(13)}{\langle\Delta,\Delta\rangle}
 =+\frac{7}{10365}\frac{\Lambda(13)}{\Lambda(11)}.
$$
从 FKL 的有符号剩余项换到 DKS 补全系数的乘数为 $-1/540$。这些数来自基线证明。DKS 的特解归一化在这里必不可少：若将特解平移 $tH_\Delta^+$，源保持相同，而待求的 $\lambda$ 变为 $\lambda-t$。[整数偶谱证明 §1、§7](../general-conjecture/proof.md) · [权 12 基线 §2–4](../weight-12/s-cocycle-proof.md) · [DKS 的齐次项](https://arxiv.org/html/2109.05018v3#S2.SS2)

**迁移决策。** 若想用该引理替代 FKL 的尖点迹求值，必须从其右侧实际导出 $-r_1$；若想连同 DKS 的块匹配一起替代，则必须直接识别出 $-28r_1$。该引理给的是双提升内积；取对角数据只得到范数平方，对提升改号不敏感。混合内积原则上保留相对符号，但这里没有将第二向量校准为 DKS 的 $Y^{-5}$ 迹，也没有产生 $\Lambda(6)\Lambda(13)$ 的局部积分计算。**不建议继续直接套用这条引理。** 只有已有适用于分裂核和非尖输入的正则化恒等式，且能明确识别 DKS 边界泛函，才有理由重新测试这条路线。
