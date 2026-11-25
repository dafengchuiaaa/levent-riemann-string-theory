# Levent–Riemann 弦理论研究档案

本档案的主线是 Dorigoni–Green–Wen（DGW）的广义 Eisenstein 级数猜想。工作从权重 12 的 $\Delta_{12}$ 算例出发，比较 FKL 的正则化除数卷积与 DKS 的模不变补全，再推广到 DGW 式 (4.13) 的整数偶谱。[DGW 猜想 1、2](https://arxiv.org/html/2501.03996v2#S1.SS2)

## 原猜想

| 问题 | 档案状态 |
|---|---|
| **猜想 1**。求 DKS 特解在整数偶谱上的尖形式齐次补全 | 一般证明稿已给出全谱的系数及 $S$-余循环并经修订。已有的网页复审覆盖较早的主稿。待对照 FKL 与 DKS 原文复核的核心是全谱上混合、正正和正零三个 Fourier 块如何合成 DKS 特解的有符号边界迹。[一般证明](general-conjecture/proof.md) · [复审范围](research-log/2026-11-25-even-audit.md) |
| **猜想 2**。DGW 的格点提升逐尖形式消去 | 将猜想 1 的系数代入 DGW 的周期关系，得到整数二圈模图函数（MGF）部分的推导。DGW 原文已处理其半整数 SMF 范围。[整数 MGF 消去](general-conjecture/integer-manin.md) · [DGW §4.3](https://arxiv.org/html/2501.03996v2#S4.SS3) |

[权重 12 的算例与符号订正](weight-12/s-cocycle-proof.md)记录了完整 $S$-余循环及其与多重模值的联系。

## 莱文特–Claude 方法迁移

已在权重 12 测试一条发布源码中的 Rallis 核内积引理。分裂格的正定与稳定秩条件、Eisenstein 输入的尖性均无法直接满足，测试未导出 DKS 的有符号边界迹。[迁移测试](method-transfer/rallis-dgw-6-4-4.md)

## 奇谱扩展

奇谱源来自 DGW §5.1。档案包含整数奇谱的补全及权重 24 底谱算例，也包含半整数奇谱的全模系数和四维 theta 提升。临界带消去空间与带外第一档 $w=s+2$ 的端点周期已有证明稿，待独立复审。[整数奇谱](general-conjecture/odd-proof.md) · [权重 24 算例](general-conjecture/weight-24-odd-bottom.md) · [半整数全模](general-conjecture/half-odd-general.md) · [奇 theta 提升](general-conjecture/half-odd-theta.md) · [临界带分类](general-conjecture/half-odd-span.md) · [带外第一档](general-conjecture/half-odd-first-exterior.md) · [DGW §5.1](https://arxiv.org/html/2501.03996v2#S5.SS1)

[研究状态](STATUS.md) · [研究交接规则](AGENTS.md) · [一般整数指标索引](general-conjecture/README.md) · [权重 12 研究日志](research-log/2026-11-24-weight-12.md) · [文献与公式](references/README.md) · [撤回与订正](corrections/README.md)
