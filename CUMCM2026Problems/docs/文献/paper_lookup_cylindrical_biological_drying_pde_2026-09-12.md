# 圆柱形生物多孔材料热风干燥：内部温度场与含水率场定向检索（第一轮）

访问日期：2026-09-12（Asia/Hong_Kong）  
范围：对流/热风干燥；圆柱、根/根茎及可迁移的多孔食品模型；重点核查内部温湿场 PDE、界面传质闭合、有效扩散率与蒸发潜热。  
本文件是检索与证据清单，不是文献综述。

## 1. 先给结论

1. **不能把固体干基含水率与空气湿度直接相减后套入通常气相定义的 `h_m`。** `X_s` 的单位基准是 kg-water/kg-dry-solid；空气湿度可能是 kg-water/kg-dry-air、相对湿度或 kg-water/m³-air。即使数值上都写成“kg/kg”，两者也不是同一热力学变量。`h_m(C_s-C_∞)` 只有在两端都是同一种气相水蒸气浓度（质量浓度或摩尔浓度）时才是标准形式。
2. **严谨边界通常需要吸附/脱附等温线或水活度。** 推荐链条为 `X_s,T_s → a_w(X_s,T_s) → p_{v,s}=a_w p_sat(T_s) → c_{v,s}=M_w p_{v,s}/(RT_s)`，再写 `j_m=h_m(c_{v,s}-c_{v,∞})`。较简化的固相写法 `-ρ_d D_eff ∂X/∂n = k_X(X_s-X_e)` 也可用，但 `X_e=X_e(T_∞,RH_∞)` 必须由等温线给出，且 `k_X` 不是未经换算的气相膜系数。
3. **把 `D_eff` 用在干基含水率方程中非常常见，但它是宏观有效参数。** 固定域、干固体体积密度近似常数时，可简化为 `∂X/∂t=∇·(D_eff∇X)`；更一般的守恒式是 `∂(ρ_d X)/∂t=∇·(ρ_d D_eff∇X)`。收缩或 `ρ_d` 变化时，应保留密度/Jacobian/固相速度项，不能无条件约去 `ρ_d`。
4. **忽略潜热没有一个仅由“温度低”决定的通用条件。** 应验证表面或体内蒸发吸热相对于对流、辐射、导热和显热储存确实很小。可作误差判据：

   `Λ_L = |j_m L_v| / (|h_T(T_∞-T_s)| + |q_rad|) ≪ 1`。

   对体内蒸发模型还应比较 `|S_evap L_v|` 与 `|∇·(k∇T)|`、`|ρc_p ∂T/∂t|`。只有在晚期低失水通量、外加温度被规定而不求解能量方程，或灵敏度计算证明删项误差可接受时，才有充分理由忽略。小热 Biot 数只说明内部温度可能近似均匀，**不等于**蒸发潜热可忽略。

### 1.1 两种 Robin 传质边界不要混写

气相膜控制形式：

`-ρ_d D_eff ∂X/∂n = h_m [c_v*(X_s,T_s)-c_{v,∞}]`,

其中 `c_v*=M_w a_w(X_s,T_s)p_sat(T_s)/(RT_s)`。

若在平衡点附近线性化，

`c_v*(X_s,T_s)-c_{v,∞} ≈ (∂c_v*/∂X)_e (X_s-X_e)`，

于是固相形式中的等效系数是

`k_X = (h_m/ρ_d)(∂c_v*/∂X)_e`，

相应质量 Biot 数为

`Bi_m = k_X L/D_eff = h_m L(∂c_v*/∂X)_e/(ρ_d D_eff)`。

因此，只有在变量、基准和系数定义相容时，才可使用简写 `Bi_m=h_mL/D_eff`。这是由界面守恒与线性化得到的推论；Defraeye & Verboven (2017) 对扩散率、渗透率和由等温线得到的湿容量之间的关系提供了直接模型证据。

## 2. 英文术语扩展与检索块

| 主题 | 主词 | 扩展词 |
|---|---|---|
| 内部场 | internal temperature field; internal moisture field | temperature distribution; moisture profile; spatially resolved; transient hygrothermal field |
| 热湿耦合 | coupled heat and mass transfer | simultaneous heat and moisture transfer; hygrothermal model; thermo-hydro model; Luikov model |
| 材料与几何 | cylindrical biological material | finite cylinder; cylindrical body; potato cylinder; carrot core/cortex; root; rhizome; porous food; biomaterial |
| 内部传质 | Fickian diffusion | effective moisture diffusivity; moisture permeability; liquid-water diffusion; vapor diffusion; capillary-porous transport |
| 表面传质 | Robin boundary condition | convective mass-transfer boundary; external mass-transfer resistance; air–solid interface; conjugate mass transfer |
| 界面热力学 | sorption isotherm | desorption isotherm; equilibrium moisture content; water activity; GAB; BET; vapor pressure; surface relative humidity |
| 含水率基准 | dry-basis moisture content | kg water/kg dry solid; moisture ratio; dry-solid density; moisture capacity |
| 无量纲数 | mass Biot number | Bi–Di correlation; Sherwood number; external/internal resistance ratio |
| 热效应 | latent heat of evaporation | evaporative cooling; enthalpy of vaporization; phase-change source; non-isothermal drying |
| 变形 | shrinkage | moving boundary; solid-matrix velocity; ALE; Lagrangian coordinate; thermo-hydro-mechanical coupling |

本轮主概念式：

`("coupled heat and mass transfer" OR "simultaneous heat and moisture transfer" OR hygrothermal) AND ("convective drying" OR "hot air drying") AND (cylinder* OR potato OR carrot OR root OR rhizome OR "porous food" OR biomaterial) AND ("moisture profile" OR "temperature distribution" OR Fick* OR "effective moisture diffusivity")`

界面闭合式：

`drying AND ("sorption isotherm" OR "equilibrium moisture content" OR "water activity") AND ("mass transfer coefficient" OR "mass Biot number" OR "vapor concentration") AND (model OR diffusion)`

潜热式：

`("latent heat" OR "evaporative cooling" OR "enthalpy of vaporization") AND ("convective drying" OR "hot air drying") AND (food OR biological) AND (model OR simulation)`

## 3. 实际数据库、查询和返回量

| 数据库/端点 | 实际查询或用途 | 总命中/取回 | 状态与筛选说明 |
|---|---|---:|---|
| Crossref | `cylindrical biological material convective drying heat mass transfer Fickian diffusion Robin boundary`；`type=journal-article`；20 条 | 259,882/20 | 成功；噪声较高，按题名、对象、内部场人工筛选 |
| Crossref | `cylindrical root rhizome potato carrot cassava convective drying temperature moisture distribution finite element` | 162,115/20 | 成功但“root”引入统计学 unit-root 噪声，仅作负面检索记录 |
| Crossref 单 DOI | `GET /works/{doi}` | 15/15 | 成功；用于题名、作者、期刊、年份核验 |
| OpenAlex 主题搜索 | `search=cylindrical biological material convective drying heat mass transfer Fickian diffusion Robin boundary` | 未取得 | 初次及一次重试均 HTTP 429；未将失败解释为零命中 |
| OpenAlex 单 DOI | `GET /works/https://doi.org/{doi}` | 15/15 | 成功；与 Crossref 交叉核对，并记录 OA 状态、引用数 |
| Semantic Scholar | `cylindrical biological material convective drying coupled heat mass transfer moisture profile` | 未取得 | 初次及一次重试均 HTTP 429；未继续撞限流 |
| Europe PMC | 内部场式：`("internal moisture distribution" OR "moisture profile" OR "temperature distribution") AND ("convective drying" OR "hot air drying") AND (...)` | 132/20 | 成功；发现 Teleken、Hu、Chen 等 |
| Europe PMC | 潜热式 | 87/15 | 成功 |
| Europe PMC | 等温线/传质/Biot 式 | 46/15 | 成功 |
| Europe PMC | 根/根茎材料式：`(potato OR carrot OR cassava OR ginger OR rhizome OR root) ...` | 266/20 | 成功；大量品质、混合干燥、薄层模型被排除 |
| Europe PMC 引用端点 | 3 个开放全文种子的 references | 未取得 | 服务维护，HTTP 503；改从 JATS `<ref>` 安全提取引用网络 |
| Europe PMC/PMC 全文 | PMC11817341、PMC7692062、PMC12861252、PMC10178041 | 4/4 | JATS 含正文；用于逐式核查边界、潜热和数值方法 |
| doi.org | 15 个 DOI 重定向 | 15/15 | 每个 DOI 均发生至少一次有效重定向到出版社；部分最终页返回 403/202 是出版社机器人策略，不是 DOI 解析失败 |
| Unpaywall | 未执行 | — | API 强制发送真实邮箱；本轮未获将本机邮箱外发的授权。OA 链接仅采用 PMC、OpenAlex 与机构库可核验链接 |

未声称检索 Scopus、Web of Science、CNKI 或万方。

筛选：排除只有平均含水率的经验薄层拟合、设备侧流场-only、微波/冻干/喷雾-only、机器学习-only，以及只用总干燥时间接近作验证的论文。按 DOI 去重。卷期年与 online-first 年冲突时采用 Crossref 的正式卷期年，并在必要处注明 online 年。

## 4. 15 篇核心论文证据卡

证据等级：A＝DOI+Crossref/OpenAlex+开放全文方法核查；B＝DOI+双库元数据+出版社摘要/正文摘录；C＝仅元数据或题名层，本清单没有用 C 级信息支撑关键公式。相关性满分 5。

### 1. Teleken et al. (2025) — 虾体内部温湿场与气相浓度边界

- 题名/作者/期刊：*Heat and Mass Transfer in Shrimp Hot-Air Drying: Experimental Evaluation and Numerical Simulation*；Jhony T. Teleken, Suélen M. Amorim, Sarah S. S. Rodrigues, Thailla W. P. de Souza, João P. Ferreira, Bruno A. M. Carciofi；*Foods*。
- 标识符/OA：[DOI 10.3390/foods14030428](https://doi.org/10.3390/foods14030428)；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11817341/)。
- 对象/几何/PDE：虾；真实不规则 3D，并与同体积有限圆柱比较；内部 Fourier 导热 + Fick 水分扩散。
- 边界/物性/收缩：质量通量用表面 `a_w p_sat(T_s)/(RT_s)` 与环境 `RH·p_sat(T_∞)/(RT_∞)` 的**气相蒸汽浓度差**；热边界含对流与 `j_mL_v`；密度、比热、导热率随局部组成/温度；未考虑变形。
- 数值/验证：COMSOL FEM、BDF、网格独立性；60/70 °C 实验温度与含水率，R²≥0.95。
- 本题对应/论断：最直接支持“`X_s` 先经水活度映射到蒸汽浓度”、潜热不可默认删除、有限圆柱仅是局部场近似。
- 局限/评分：虾而非药材，且不收缩；5.0/5，证据 A。

### 2. Yang, Sakai & Watanabe (2001) — 圆柱马铃薯、热湿–非各向同性收缩

- 题名/期刊：[ *Drying Model with Non-Isotropic Shrinkage Deformation Undergoing Simultaneous Heat and Mass Transfer*](https://doi.org/10.1081/DRT-100105299)；Hui Yang, Noboru Sakai, Manabu Watanabe；*Drying Technology*。
- 对象/几何/PDE：圆柱马铃薯，二维轴–径向；热传导、湿扩散与虚功固体力学耦合。
- 边界/物性/收缩：摘要确认初边值问题、非恒定物理/热物性和轴/径不同收缩；边界的具体热力学写法未从可得摘要核验。
- 数值/验证：Galerkin FEM；中心温度、平均含水率、轴向和径向收缩实验。
- 本题对应/论断：精确支持圆柱二维基准与非各向同性收缩；可用于检验“一维径向+端面等效源”的降阶误差。
- 局限/评分：单工况、边界公式与网格验证未核验；4.9/5，证据 B。

### 3. Srikiatden & Roberts (2008) — 马铃薯/胡萝卜实测内部含水率剖面

- 题名/期刊：[ *Predicting Moisture Profiles in Potato and Carrot during Convective Hot Air Drying Using Isothermally Measured Effective Diffusivity*](https://doi.org/10.1016/j.jfoodeng.2007.06.009)；Jaruk Srikiatden, John S. Roberts；*Journal of Food Engineering*。
- 对象/几何/PDE：直径 1.4/2.8 cm 马铃薯、胡萝卜芯和皮层圆柱样品；Fick 第二定律与同步热–湿模型；输出位置相关的温度和含水率。
- 边界/物性/收缩：`D_eff(T)` 由真正等温实验独立测得，再用于非等温干燥；具体 Robin 公式和收缩处理未从摘要核验。
- 数值/验证：Crank–Nicolson 有限差分；70 °C、1.5 m/s 下直接比较实测温度与水分剖面。
- 本题对应/论断：强支持在干基含水率/水分场中使用宏观 `D_eff`，并强调不可用非等温整条干燥曲线不加辨识地反演扩散率。
- 局限/评分：作者称样品为低孔隙“hygroscopic non-porous”；不是药材；4.9/5，证据 B。

### 4. Curcio et al. (2008) — 圆柱蔬菜共轭热–质传递 FEM

- 题名/期刊：[ *Simulation of Food Drying: FEM Analysis and Experimental Validation*](https://doi.org/10.1016/j.jfoodeng.2008.01.016)；Stefano Curcio, Maria Aversa, Vincenza Calabrò, Gabriele Iorio；*Journal of Food Engineering*。
- 对象/几何/PDE：圆柱蔬菜；空气动量、热量和水分与固体内部热湿场的非线性非稳态 PDE 共轭求解。
- 边界/物性/收缩：界面采用热/质通量连续，从而不必预设经验 `h_T,h_m`；具体物性函数和收缩处理在本轮可得摘要中未核验。
- 数值/验证：FEM；实验与预测吻合良好。
- 本题对应/论断：提供“固定 Robin 系数”之外的高保真对照；也说明外部阻力控制时气流特征才强烈影响干燥。
- 局限/评分：材料名称与内部剖面测量细节未从摘要核验；4.8/5，证据 B。

### 5. Niamnuy et al. (2008) — 三维热湿–力学收缩耦合

- 题名/期刊：[ *Modeling Coupled Transport Phenomena and Mechanical Deformation of Shrimp during Drying in a Jet Spouted Bed Dryer*](https://doi.org/10.1016/j.ces.2008.07.031)；Chalida Niamnuy, Sakamon Devahastin, Somchart Soponronnarit, G. S. Vijaya Raghavan；*Chemical Engineering Science*。
- 对象/几何/PDE：高度收缩、不规则虾体；3D 热传导、质量扩散、弹性固体力学。
- 边界/物性/收缩：明确求解初边值问题并预测应力/变形；边界具体浓度变量与物性函数未从出版社摘要核验。
- 数值/验证：COMSOL FEM；比较含水率、中层温度、收缩实验，并对比有/无变形模型。
- 本题对应/论断：支持移动固体/收缩会反馈热湿场，不能只事后修改半径。
- 局限/评分：喷动床且温度较高，非圆柱；4.6/5，证据 B。

### 6. Lamnatou et al. (2010) — Luikov 多孔体与有限体积共轭/非共轭对照

- 题名/期刊：[ *Finite-Volume Modelling of Heat and Mass Transfer during Convective Drying of Porous Bodies—Non-Conjugate and Conjugate Formulations Involving the Aerodynamic Effects*](https://doi.org/10.1016/j.renene.2009.11.008)；Chr. Lamnatou, E. Papanicolaou, V. Belessiotis, N. Kyriakis；*Renewable Energy*。
- 对象/几何/PDE：建筑多孔体和农产品；Luikov 热–湿方程；应用示例为矩形板及干燥室。
- 边界/物性/收缩：可选固体-only 的传递系数边界，或空气+固体共轭界面；两者均含蒸发相变；收缩未见摘要说明。
- 数值/验证：有限体积法；用文献实验/数值结果验证。
- 本题对应/论断：支持节点/控制体守恒框架，并可把 Robin 边界与共轭气相模型交叉验证；相变项是标准能量耦合之一。
- 局限/评分：不是圆柱，且不能直接证明本题半控制体离散公式；4.5/5，证据 B。

### 7. Zhu et al. (2021) — 多相多孔介质 THM 双向耦合

- 题名/期刊：[ *Multiphase Porous Media Model with Thermo-Hydro and Mechanical Bidirectional Coupling for Food Convective Drying*](https://doi.org/10.1016/j.ijheatmasstransfer.2021.121356)；Yueqiang Zhu, Peng Wang, Dongliang Sun, Zhiguo Qu, Bo Yu；*International Journal of Heat and Mass Transfer*。
- 对象/几何/PDE：香菇；多相多孔介质热–水–力学模型，内部蒸发率、温度和含水率均为空间场。
- 边界/物性/收缩：质量与能量方程显式考虑固体基质移动速度；双向收缩耦合。精确界面公式和几何未从可得摘要核验。
- 数值/验证：数值法名称未核验；实验验证 THM，并与不含力学反馈的 TH 模型比较；峰值蒸发率差异约 29.5%，内部温差最高约 3 °C。
- 本题对应/论断：强烈反驳“收缩只改几何、无需速度/Jacobian 耦合”的一般化说法。
- 局限/评分：香菇而非圆柱药材，闭源全文；4.8/5，证据 B。

### 8. Adrover, Venditti & Brasiello (2020) — 非等温移动边界与等温线闭合

- 题名/期刊：[ *A Non-Isothermal Moving-Boundary Model for Continuous and Intermittent Drying of Pears*](https://doi.org/10.3390/foods9111577)；Alessandra Adrover, Claudia Venditti, Antonio Brasiello；*Foods*；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7692062/)。
- 对象/几何/PDE：球形梨；含收缩速度的水分与能量守恒型对流–扩散方程。
- 边界/物性/收缩：表面水蒸气浓度由脱附等温线/RH 和饱和蒸气压得到；热边界含对流和潜热；`ρ,c_p,k` 随局部水分，`D_eff(T)`；点态收缩速度预测移动边界。
- 数值/验证：FEM + ALE moving mesh + BDF；连续/间歇干燥，R²>0.99。
- 本题对应/论断：支持 `X_s→RH_eq/a_w→c_v` 的边界链条，也说明坐标映射后必须在速度项或 Jacobian 中保持守恒。
- 局限/评分：球体且理想各向同性收缩；5.0/5，证据 A。

### 9. Nguyen, Ngo & Le (2019) — 干基含水率、Bi–Di 关联与内部场

- 题名/期刊：[ *Experimental and Numerical Investigation of Transport Phenomena and Kinetics for Convective Shrimp Drying*](https://doi.org/10.1016/j.csite.2019.100465)；Minh Phu Nguyen, Tu Thien Ngo, Thanh Danh Le；*Case Studies in Thermal Engineering*；出版社开放获取。
- 对象/几何/PDE：真实虾形；ANSYS 求内部温度和含水率分布。
- 边界/物性/收缩：由多温度/风速实验和 Bi–Di 关联获得干燥常数、`h_m`、`D_eff`；具体 Robin 变量定义、潜热和收缩处理未从摘要核验。
- 数值/验证：输运参数来自实验；报告中心与尾部水分差异、升温过程和达到 0.25 d.b. 的时间。
- 本题对应/论断：说明 `D_eff`、`h_m` 与质量 Biot 数常用于连接内部扩散和外部阻力，但不能据此把固体 d.b. 与空气湿度直接相减。
- 局限/评分：边界公式和局部场实验验证不足；4.2/5，证据 B。

### 10. Hu et al. (正式卷期 2026；online 2025) — 多相、多尺度、水活度与持续蒸发冷却

- 题名/期刊：[ *A Multiphase and Multiscale Mechanistic Model for Hot Air Drying of Shiitake Mushroom*](https://doi.org/10.1016/j.crfs.2025.101296)；Lina Hu, Xin Jin, Yitong Xie, Jinfeng Bi, Ruud G. M. van der Sman；*Current Research in Food Science*；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC12861252/)。
- 对象/几何/PDE：香菇帽半球+菌褶双孔隙层；菌丝液态水 Fick 扩散与孔隙蒸汽扩散，总水分守恒；能量式含液/汽焓输运和汽化焓。
- 边界/物性/收缩：`j_m=h_ext(c_{v,s}-c_air)`，其中 `c_{v,s}` 与表面组织 `a_w` 平衡；热边界含对流/辐射与蒸发冷却。Flory–Huggins/Flory–Rehner、双 Maxwell 黏弹模型把 `a_w`、收缩和细胞结构耦合起来。
- 数值/验证：控制体离散框架；具体时间积分器未核验。用多个温度下平均水分与中心温度拟合，35 °C 独立 RH 数据验证，并与在线 NMR/动态吸附数据互证。
- 本题对应/论断：是四个特别核查点最完整的近期证据；直接显示温和低温干燥仍可能存在持续内部蒸发冷却。
- 局限/评分：复杂度远高于一维 Fick 模型，几何不是圆柱；5.0/5，证据 A。

### 11. Chen et al. (2023) — 花生壳–仁二维轴对称内部场

- 题名/期刊：[ *A Heat and Mass Transfer Model of Peanut Convective Drying Based on a Two-Component Structure*](https://doi.org/10.3390/foods12091823)；Pengxiao Chen, Nan Chen, Wenxue Zhu, Dianxuan Wang, Mengmeng Jiang, Chenling Qu, Yu Li, Zhuoyun Zou；*Foods*；[PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10178041/)。
- 对象/几何/PDE：单粒花生，壳+仁二维轴对称嵌套椭圆；瞬态导热与含水率扩散。
- 边界/物性/收缩：水分写 `-D_i∂M/∂n=h_m(M-M_e)`，`M,M_e` 都是固体 d.b. 含水率；`h_m` 由 Sherwood 关联估计。能量方程含 `ρh_g∂M/∂t` 体积潜热耦合，表面热边界为对流；`k,c_p` 随含水率实测；忽略收缩。
- 数值/验证：COMSOL，网格核验，1 min 步长；平均干燥曲线、壳/仁对比和低场 NMR 水分状态验证。
- 本题对应/论断：证明 `D_eff` 用于 d.b. 方程很常见；也展示 `M-M_e` 型边界的合法前提是 `M_e` 为同基准平衡含水率，而不是空气湿度。
- 局限/评分：正文写法省略了显式 `ρ_d`，因此其 `h_m` 是固相等效系数还是气相膜系数存在定义含混；不宜原样照搬。4.7/5，证据 A。

### 12. Defraeye & Radu (2017) — 空气–果实界面共轭模型

- 题名/期刊：[ *Convective Drying of Fruit: A Deeper Look at the Air-Material Interface by Conjugate Modeling*](https://doi.org/10.1016/j.ijheatmasstransfer.2017.01.002)；Thijs Defraeye, Andrea I. Radu；*International Journal of Heat and Mass Transfer*；[机构库开放稿](https://www.dora.lib4ri.ch/empa/dload/empa:13545/PDF2/view)。
- 对象/几何/PDE：半圆苹果片；空气湍流热湿输运与果实组织热湿场同时求解。
- 边界/物性/收缩：界面通量由共轭解连续得到，不预设常数 `h_m,h_T`；具体组织物性函数和收缩处理未从摘要核验。
- 数值/验证：连续体共轭模型，经实验验证；可解析局部、时变界面传递系数和内部含水率分布。
- 本题对应/论断：发现局部高湿微环境可造成负“传质系数”表象和局部再吸湿；说明单一常数 Robin 边界可能掩盖空间/时间变化。论文还明确认为即便近等温低温工况也不宜轻率采用等温模型。
- 局限/评分：计算代价高、几何非圆柱；4.8/5，证据 B。

### 13. Defraeye & Verboven (2017；online 2016) — 扩散率、渗透率与吸附湿容量

- 题名/期刊：[ *Convective Drying of Fruit: Role and Impact of Moisture Transport Properties in Modelling*](https://doi.org/10.1016/j.jfoodeng.2016.08.013)；Thijs Defraeye, Pieter Verboven；*Journal of Food Engineering*；[开放作者稿入口](https://www.sciencedirect.com/science/article/am/pii/S0260877416303004)。
- 对象/几何/PDE：苹果组织验证过的湿热连续体模型，输出内部温度/含水率分布；几何细节本轮未从摘要确认。
- 边界/物性/收缩：比较 moisture permeability `K_m` 与 diffusivity `D_m`；二者通过由吸附等温线导出的 moisture capacity `C_m` 相关：`K_m=D_m C_m`。边界具体式和收缩未核验。
- 数值/验证：参数敏感性/连续体模拟；基于已验证模型比较内部场与干燥曲线。
- 本题对应/论断：支持 `D_eff` 是把微观机制汇总后的宏观参数；也解释了为什么不同状态变量下的 `h_m`、Biot 数不能机械通用。
- 局限/评分：主要是参数形式敏感性而非新实验圆柱模型；4.8/5，证据 B。

### 14. Rahman & Kumar (2011) — 有限圆柱、外部阻力与收缩的解析扩散模型

- 题名/期刊：[ *Evaluation of Moisture Diffusion Coefficient of Cylindrical Bodies Considering Shrinkage During Natural Convection Drying*](https://doi.org/10.2202/1556-3758.1352)；Najmur Rahman, Subodh Kumar；*International Journal of Food Engineering*。
- 对象/几何/PDE：长 0.05 m、直径 0.01 m 马铃薯有限圆柱；Fick 扩散解析模型，比较有限/无限圆柱及有/无收缩。
- 边界/物性/收缩：显式考虑固–气界面外部传质阻力；`D_eff(T,X)` 反演；收缩进入有限几何。边界驱动力的精确状态变量未从摘要核验。
- 数值/验证：Levenberg–Marquardt 反演，自然对流干燥曲线；`D_eff` 约 3.93–8.63×10⁻¹⁰ m²/s（40–60 °C）。
- 本题对应/论断：无限圆柱假设会高估 `D_eff`，考虑收缩会降低反演值；端面效应不能在未验证时直接忽略或均匀体源化。
- 局限/评分：无内部温度场，验证主要是整体含水率；4.7/5，证据 B。

### 15. Prachayawarakorn et al. (2002) — GAB 脱附等温线的界面闭合证据

- 题名/期刊：[ *Desorption Isotherms and Drying Characteristics of Shrimp in Superheated Steam and Hot Air*](https://doi.org/10.1081/DRT-120002823)；Somkiat Prachayawarakorn, Somchart Soponronnarit, Somboon Wetchacama, Donrudee Jaisut；*Drying Technology*。
- 对象/几何/模型：虾；50–80 °C 脱附等温线，70–140 °C 热风干燥；GAB 比 BET 更好描述平衡含水率，参数用 Arrhenius 温度关联。
- 物性/收缩/验证：比较干燥速率、`D_eff` 与收缩；实验等温线和热风数据。不是内部空间 PDE 论文。
- 本题对应/论断：为 `X_e(T,RH)`/`a_w(X,T)` 提供原始实验依据，说明表面固体含水率必须经过等温线才能与空气蒸汽状态连接。
- 局限/评分：温区较高且对象非药材；只作为边界本构证据，不用于验证内部场。4.3/5，证据 B。

## 5. 对四个特别问题的证据矩阵

| 问题 | 直接支持 | 修正或反驳 | 本题建议 |
|---|---|---|---|
| `h_m(C_surface-C_air)` 能否直接用 | Teleken 2025、Adrover 2020、Hu 2026 都在气相蒸汽浓度/分压层面写驱动力 | Chen 2023 的 `M-M_e` 说明固相形式也存在，但两项必须同为 d.b.，系数需重新定义 | 不要用 `X_s - RH_air` 或 `X_s - Y_air`；选择气相形式，或明确 `X_e(T,RH)` 和固相等效 `k_X` |
| 是否需要平衡含水率/蒸汽浓度/等温线 | Teleken、Adrover、Hu、Prachayawarakorn；Defraeye & Verboven 给出 sorption capacity 关系 | 高湿边界会局部再吸湿，Defraeye & Radu 表明固定常数边界可能失真 | 优先用 `a_w(X_s,T_s)`；数据不足时用 GAB/Henderson 等拟合的 `X_e`，并报告外推范围 |
| `D_eff` 用于干基含水率是否常见 | Srikiatden & Roberts、Chen、Rahman、Nguyen | Defraeye & Verboven 指出 `D_eff` 与 permeability 不是同一参数，数值随含水率、温度与辨识方法变化 | 可用，但把它称为“宏观有效扩散率”；收缩时写守恒式并做独立辨识/敏感性 |
| 何时可忽略潜热 | 没有找到支持“热风温度不高即可忽略”的核心论文 | Teleken、Adrover、Lamnatou、Hu、Chen 都保留相变吸热；Hu 显示温和干燥仍有持续蒸发冷却；Defraeye & Radu 反对近等温即等温化 | 先计算 `Λ_L` 或做有/无潜热灵敏度；全程最大温差、平均含水率和终止时间误差均在容许范围后再删项 |

## 6. 可直接用于当前模型的最小修改建议

1. 若状态变量是干基含水率 `X`，内部式优先写成守恒形式：`∂(ρ_dX)/∂t=(1/r)∂[rρ_dD_eff∂X/∂r]/∂r + S_X`。固定、恒 `ρ_d` 时再化简。
2. 表面质量边界首选气相形式：`-ρ_dD_eff∂X/∂r=h_m[c_v*(X_s,T_s)-c_{v,∞}]`，并用实验吸附/脱附等温线给 `a_w(X_s,T_s)`。
3. 若改用 `-D_eff∂X/∂r=k_X(X_s-X_e)`，必须说明 `k_X` 的定义/单位、`X_e(T,RH)` 模型，以及质量 Biot 数采用哪个系数。
4. 能量边界保留 `j_mL_v` 作为基线模型；只有在无潜热对照的场温度、干燥终止时间和峰值误差都足够小时才删项。
5. 用 Yang 2001、Curcio 2008 或有限圆柱模型建立二维轴对称基准，再检验一维径向模型的端面等效源。Rahman 2011 已表明有限/无限圆柱会系统改变反演的 `D_eff`。

## 7. 开放问题与未核验项

- 本轮没有找到一篇同时满足“中药根茎 + 精确有限圆柱 + 内部温度/含水率实测剖面 + 收缩 + 水活度边界”的单篇论文；证据需要由圆柱根类模型和多孔界面模型拼接。
- Yang、Niamnuy、Zhu、Rahman 的边界具体公式未获得合法开放全文，已明确标注，未凭题名猜写。
- Unpaywall 尚未查询；若用户授权提供一个可发送给 Unpaywall 的邮箱，可继续逐 DOI 核查合法 OA 位置。
- Semantic Scholar 与 OpenAlex 的主题搜索接口本轮均受匿名限流；这不影响 15 篇 DOI 的 Crossref/OpenAlex 单条双库核验。

## 8. 检索方法技能引用

本次检索流程使用 paper-lookup 技能所提供的数据库路由、JATS 正文抽取、OpenAlex 摘要重建和可复现记录方法。按技能要求引用：Timothy Kassis, Vinayak Agarwal, Yuhuan He, Darshil Patel, Aubrey M. Brueckner (2026), *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*, arXiv:2609.00065, current version v2；[arXiv](https://arxiv.org/abs/2609.00065)，[DOI 形式](https://doi.org/10.48550/arXiv.2609.00065)。
