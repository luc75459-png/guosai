# 2026 数学建模 A 题：可复现的定向论文检索（第一轮）

访问日期：2026-09-11（Asia/Hong_Kong）  
范围：候选发现、标识符核验和证据分级；**不构成文献综述**。  
去重规则：优先按小写 DOI 去重；无 DOI 时按去标点、折叠空格后的标准化英文题名去重。

## 1. 英文术语扩展表

| 主题 | 核心英文术语 | 扩展/同义词 | 容易误检、需结合上下文的词 |
|---|---|---|---|
| 热风对流干燥 | hot-air drying; convective drying | air drying; forced-convection drying; dehydration of biological materials | air-drying 也用于涂层、污泥、木材 |
| 内部耦合温湿场 | coupled heat and mass transfer; simultaneous heat and moisture transfer | conjugate heat and mass transfer; non-isothermal moisture diffusion; internal temperature/moisture distribution | coupled model 可能只是设备空气侧 |
| 水分输运 | Fickian diffusion; effective moisture diffusivity | liquid diffusion model; moisture concentration equation; water activity; sorption/desorption isotherm | drying kinetics 常仅指薄层经验式 |
| 热输运 | transient heat conduction; Fourier heat equation | energy conservation; thermal diffusivity; latent heat of evaporation | heat transfer 可能仅为空气流场 |
| 圆柱几何 | finite cylinder; cylindrical body; axisymmetric cylinder | one-dimensional radial model; radial diffusion; finite-length cylinder | infinite cylinder 会忽略端面 |
| 表面对流边界 | Robin boundary condition; third-kind boundary condition | mixed boundary; convective heat/mass boundary; evaporative boundary; external transfer resistance | constant surface concentration/temperature 是 Dirichlet，不等价 |
| 端面/降维 | end effects; axial diffusion; finite-length correction | radial–axial model; aspect-ratio reduction; volume averaging; effective volumetric source | “source term”也可能指化学反应或内热源 |
| 变物性 | temperature-dependent properties; moisture-dependent properties | composition-dependent density/specific heat/thermal conductivity; nonlinear diffusivity | apparent/effective property 可能是纯拟合参数 |
| 收缩 | shrinkage; volume change; deforming biological material | isotropic/non-isotropic shrinkage; ideal shrinkage; shrinkage velocity | shrinkage kinetics 可能没有内部场方程 |
| 移动域 | moving boundary; free boundary | ALE (Arbitrary Lagrangian–Eulerian); moving mesh; fixed-domain transformation; front fixing | moving boundary 也常指相变 Stefan 问题 |
| 材料坐标 | material coordinate; Lagrangian coordinate | referential coordinate; ξ=r/R(t); normalized radial coordinate; grid velocity | “convective term”需区分物理速度与坐标网格速度 |
| 有限体积 | cell-centered finite volume method | control-volume method; conservative flux; half control volume; axisymmetric/radial finite volume | node-centered 与 cell-centered 不可混称 |
| 面扩散系数 | harmonic mean/average; arithmetic mean/average | interface transmissibility; harmonic averaging point; integral/path average | 对连续平滑系数与跳跃系数的结论不相同 |
| 时间推进 | Heun method; SSPRK(2,2); explicit trapezoidal rule | strong-stability-preserving Runge–Kutta; BDF; backward differentiation formula | “二阶”不代表无扩散 CFL 限制 |
| 验证 | code verification; solution verification | method of manufactured solutions (MMS); grid refinement; observed order; analytical benchmark; mesh independence | validation（对实验）与 verification（对方程/代码）应分开 |
| 半径插值 | PCHIP; shape-preserving piecewise cubic Hermite interpolation | local monotone piecewise cubic interpolation; monotonicity-preserving interpolation | PCHIP 不会把本来非单调的含噪数据强制变为单调 |
| 终止准则 | maximum local dry-basis moisture content | worst-point moisture criterion; max-norm stopping condition | average moisture threshold 不能保证全场阈值 |

## 2. 分主题检索式（可直接复用）

以下为第一轮采用或建议继续使用的英文检索式。数据库字段语法不同，保留了可迁移的概念块。

### T1 内部耦合温湿场、Robin 边界与变物性

- 通用概念式：`("coupled heat and mass transfer" OR "simultaneous heat and moisture transfer") AND ("hot air drying" OR "convective drying") AND (food OR fruit OR vegetable OR biomass OR "biological material") AND (diffusion OR "internal distribution")`
- PubMed：`(("hot air drying"[Title/Abstract] OR "convective drying"[Title/Abstract]) AND (model*[Title/Abstract] OR simulation[Title/Abstract]) AND (plant*[Title/Abstract] OR food[Title/Abstract] OR fruit*[Title/Abstract]))`
- Europe PMC：`("heat and mass transfer" AND drying AND (food OR fruit OR vegetable OR biomass) AND (model OR simulation))`

### T2 圆柱有限长度、端面效应与一维降维

- `("finite cylinder" OR "cylindrical body") AND drying AND ("radial diffusion" OR "axial diffusion" OR "moisture diffusion") AND ("convective boundary" OR "external mass transfer resistance")`
- `("finite cylinder" AND drying AND ("end effects" OR "finite-length correction" OR "radial-axial"))`
- 对比检索：`("finite cylinder" OR "infinite cylinder") AND drying AND shrinkage AND diffusivity`

### T3 收缩、移动边界与材料/ALE 坐标

- `(shrinkage AND drying AND ("moving boundary" OR Lagrangian OR "material coordinate" OR ALE))`
- `("convective drying" AND shrinkage AND ("Arbitrary Lagrangian Eulerian" OR "moving mesh" OR "fixed domain transformation"))`
- `((isotropic OR anisotropic OR "non-isotropic") AND shrinkage AND "heat and mass transfer" AND drying)`

### T4 圆柱有限体积、变系数面通量与平均方式

- `("cell-centered finite volume" OR "control volume") AND (cylindrical OR radial OR axisymmetric) AND (diffusion OR conduction)`
- `("finite volume" AND "nonlinear diffusion" AND ("harmonic averaging" OR "arithmetic average" OR "interface coefficient" OR transmissibility))`
- arXiv：`all:"finite volume" AND all:"nonlinear diffusion"`

### T5 时间推进与数值验证

- `("SSPRK(2,2)" OR Heun OR "strong stability preserving Runge Kutta") AND (diffusion OR parabolic)`
- `("backward differentiation formula" OR BDF) AND (nonlinear diffusion OR drying)`
- `("method of manufactured solutions" OR MMS) AND (diffusion OR heat conduction) AND ("grid refinement" OR "observed order")`

### T6 PCHIP 半径数据

- `("monotone piecewise cubic interpolation" OR "shape-preserving piecewise cubic Hermite" OR PCHIP) AND (algorithm OR slopes OR monotonicity)`

### T7 中药材/药食两用植物的物理场模型

- `(("medicinal plant" OR herbal OR rhizome OR root) AND ("hot air drying" OR "convective drying") AND (diffusion OR "heat transfer" OR "mass transfer") AND (model OR simulation))`
- 排除块（需人工复筛，不建议直接全部写成 NOT 以免漏检）：`spray drying; freeze drying; microwave-only; machine-learning-only; thin-layer-only; dryer-airflow-only`。

## 3. 实际数据库与查询日志

### 3.1 数据库状态

| 数据库/服务 | 本轮实际状态 | 说明 |
|---|---|---|
| Crossref | 成功 | 用于候选发现与 DOI 单条元数据核验 |
| OpenAlex | 搜索失败；DOI 单条接口成功 | 搜索先 HTTP 429，按要求重试一次后 HTTP 503；18 个 DOI 单条查询均成功 |
| Semantic Scholar | 失败 | 无 API key 的搜索返回 429；一次批量 DOI 重试未返回可用响应，不再继续重试 |
| PubMed | 成功 | 两个有界主题检索；返回 PMID 后用 eSummary 取元数据 |
| Europe PMC | 成功 | 主题检索，并通过 `fullTextXML` 读取 3 篇开放全文 |
| arXiv | 部分成功 | 一个精确式 0 条；一个移动边界式无可用响应；一个较宽数值式 9 条 |
| Unpaywall | **未执行（待授权）** | API 强制发送真实邮箱；系统不允许未经明确授权将本机 Git 邮箱外发。没有伪称已查 |
| doi.org | 成功 | 18 个入选 DOI 均返回 HTTP 302 有效重定向 |
| Scopus/WoS/CNKI/万方 | 未检索 | 当前没有相应可访问接口，故不作覆盖声明 |

### 3.2 发现性查询与返回数量

“总命中/本次取回”；所有检索均为有界首轮，未声称穷尽。

| 数据库 | 实际查询 | 总命中/取回 | 筛选备注 |
|---|---|---:|---|
| OpenAlex | `search=coupled+heat+mass+transfer+hot+air+drying+cylindrical+food&sort=cited_by_count:desc` | 429；重试 503 | 搜索服务不可用；改用 DOI 单条核验 |
| Crossref | `query.bibliographic=coupled+heat+mass+transfer+convective+drying+cylinder; type=journal-article; sort=is-referenced-by-count` | 1,801,764/12 | 过宽、噪声高，仅作试探 |
| Crossref | `coupled heat mass transfer convective drying food finite element` | 3,175,708/8 | 保留 Murugesan、Koukouch 等题名强匹配项 |
| Crossref | `shrinkage moving boundary convective drying food` | 1,490,357/8 | 保留 Adrover 系列、Mayor–Sereno 综述作术语入口 |
| Crossref | `cylindrical food drying moisture diffusion temperature distribution` | 2,970,778/8 | 排除微波、设备侧等；保留 Rahman–Kumar |
| Crossref | `finite volume cylindrical diffusion variable coefficient harmonic mean` | 1,711,961/8 | 保留 Zhang–Shan；其余多为分数阶/无关问题 |
| Crossref | `strong stability preserving Runge Kutta method second order` | 2,732,178/8 | 再以精确题名锁定 Gottlieb–Shu–Tadmor |
| Crossref | `monotone piecewise cubic interpolation Fritsch Carlson` | 146,564/8 | 锁定 Fritsch–Carlson 与 Fritsch–Butland |
| Crossref | `Lagrangian coordinate shrinkage drying heat mass transfer food` | 2,653,405/10 | 保留 Yang、Ning、Aprajeeta、Das 等 |
| Crossref | `material coordinate shrinking cylinder diffusion moving boundary` | 1,046,369/10 | 食品相关性低；未进入核心证据 |
| Crossref | `Arbitrary Lagrangian Eulerian convective drying shrinkage` | 164,138/10 | 保留 2025 banana ALE 与 2017 banana ALE 候选 |
| Crossref | `fixed domain transformation moving boundary diffusion shrinking sphere` | 1,253,621/10 | 多为非干燥数学问题；作方法词网络 |
| Crossref | `finite volume cylindrical coordinates radial diffusion origin boundary control volume` | 788,770/10 | 多为一般圆柱守恒格式；未发现直接支持“半控制体+本题方程”的单篇证据 |
| Crossref | `cell centered finite volume radial heat conduction cylinder axis` | 4,171,575/10 | 噪声高；未进入核心证据 |
| Crossref | `finite volume diffusion axisymmetric cylindrical coordinate singularity` | 1,368,112/10 | 多为电磁/等离子体/流固；未进入核心证据 |
| Crossref | `finite cylinder drying axial radial moisture diffusion convective boundary` | 1,066,628/10 | 发现 2026 apricot moving-boundary 论文 |
| Crossref | `finite cylinder heat mass transfer drying analytical solution end effects` | 337,108/10 | 发现 Koukouch 与 1998 cylinder 应用 |
| Crossref | `finite cylinder diffusion convective boundary average moisture content` | 1,257,356/10 | 噪声高；未直接支持均匀端面体源 |
| Europe PMC | `(drying AND (heat transfer OR moisture diffusion) AND (cylinder OR cylindrical) AND (food OR plant OR biological))` | 3,912/12 | 过宽，含冻干/喷雾等排除项 |
| Europe PMC | `("heat and mass transfer" AND drying AND (food OR fruit OR vegetable OR biomass) AND (model OR simulation))` | 1,114/15 | 按引用数排序仍有高引无关项，说明需题名/摘要人工筛选 |
| Europe PMC | `(shrinkage AND drying AND ("moving boundary" OR Lagrangian OR "material coordinate" OR ALE))` | 152/1（计数复核） | 查询回显与发送式一致 |
| Europe PMC | `((medicinal plant OR herbal OR rhizome OR root) AND ("hot air drying" OR "convective drying") AND (diffusion OR "heat transfer" OR "mass transfer"))` | 299/1（计数复核） | 大量仅品质/动力学论文，核心物理场项少 |
| PubMed | `(("hot air drying"[Title/Abstract] OR "convective drying"[Title/Abstract]) AND (model*[Title/Abstract] OR simulation[Title/Abstract]) AND (plant*[Title/Abstract] OR food[Title/Abstract] OR fruit*[Title/Abstract]))` | 110/20 | 发现 Teleken 2025、Zambra 2022 等 |
| PubMed | `(drying[Title/Abstract] AND shrinkage[Title/Abstract] AND ("moving boundary"[Title/Abstract] OR diffus*[Title/Abstract]))` | 106/20 | 发现 Adrover 2020 等 |
| arXiv | `all:"nonlinear diffusion" AND all:"finite volume" AND (all:"harmonic averaging" OR all:"variable coefficient")` | 0/0 | 合法零命中 |
| arXiv | `all:drying AND (all:"moving boundary" OR all:Lagrangian OR all:"material coordinate")` | 无可用响应 | 不解释为零命中 |
| arXiv | `all:"finite volume" AND all:"nonlinear diffusion"` | 9/9 | 发现 Nishikawa–Diskin 2022（算术平均+MMS）等 |
| Semantic Scholar | `coupled heat mass transfer convective drying cylinder food` | 429/0 | 共享池限流 |
| Semantic Scholar | 16 个 DOI 的 `/paper/batch` | 无可用响应 | 429 后唯一一次重试；停止 |

### 3.3 单条核验端点

- DOI：`HEAD https://doi.org/{doi}`；入选 18 个 DOI 均为 302 并给出出版社跳转。
- Crossref：`GET https://api.crossref.org/works/{doi}`；每篇 1/1。
- OpenAlex：`GET https://api.openalex.org/works/doi:{doi}?select=id,title,authorships,publication_year,primary_location,is_retracted,open_access`；每篇 1/1。
- Europe PMC 全文：`GET /PMC11817341/fullTextXML`、`/PMC7692062/fullTextXML`、`/PMC9265362/fullTextXML`；三次均有 JATS `<body>`，并经技能脚本抽取模型/方法/结果段。
- arXiv 技能引用核验：`id_list=2609.00065`，1/1，当前 v2（引用时不附版本号）。

## 4. 第一轮入选论文证据表

证据等级：A＝DOI+双库元数据+开放全文方法/结果；B＝DOI+双库元数据+摘要或权威开放章节；C＝DOI+双库元数据，但模型细节仅题名层或未获摘要。相关性 1–5 分针对本题模型，而非论文质量。

### 4.1 热风干燥、收缩、有限几何

#### A1. Teleken et al. (2025), *Heat and Mass Transfer in Shrimp Hot-Air Drying: Experimental Evaluation and Numerical Simulation*

- 作者/期刊/标识符：Jhony T. Teleken et al.; *Foods* 14(3), 428; DOI [10.3390/foods14030428](https://doi.org/10.3390/foods14030428); [PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC11817341/)。Crossref/OpenAlex 题名、作者、年份一致。
- 对象/几何：热风干燥虾；真实不规则 3D 几何，并与同体积有限圆柱比较。
- 方程/边界：内部 Fick 水分扩散 + Fourier 导热；质量 Robin 通量由水活度、表面饱和蒸气压、环境 RH 决定；热边界含对流换热和蒸发潜热。
- 物性/收缩：密度、比热、导热系数按局部水/蛋白组成与温度计算；`D_eff` 与 `h_m` 按工况拟合；介质假设各向同性、**不变形**。
- 数值/验证：COMSOL 一阶四面体 FEM，PARDISO，BDF，最大步长 10 s；网格独立性；60/70 °C 的中心附近温度与平均水分实验，R²>0.95，RMSE<1.12 °C 和 0.22 kg/kg。
- 对应论断：直接支持内部温湿场、Robin 换热/传质、蒸发潜热耦合、变物性、BDF 与网格独立性；也说明有限圆柱对总体曲线可近似，但局部场不如真实几何。
- 局限：3D 虾体、参数拟合且无收缩；不能证明一维端面均匀体源或 Heun 稳定性。相关性 5.0/5；证据 A。

#### A2. Adrover, Venditti & Brasiello (2020), *A Non-Isothermal Moving-Boundary Model for Continuous and Intermittent Drying of Pears*

- 期刊/标识符：*Foods* 9(11), 1577; DOI [10.3390/foods9111577](https://doi.org/10.3390/foods9111577); [PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC7692062/)。双库元数据一致。
- 对象/几何：Rocha 梨，球形一维径向移动边界。
- 方程/边界：水分和热量均写为含收缩速度 `v_s` 的守恒型对流–扩散；表面质量 Robin 条件由脱附等温线/RH 和 `h_m` 控制；热 Robin 条件包含 `h_T` 和蒸发潜热。
- 物性/收缩：产品密度、比热、导热系数为局部含水量函数；`D_eff(T)`；理想收缩系数 α=1，边界速度与局部扩散水通量耦合。
- 数值/验证：FEM + ALE moving mesh；二次 Lagrange 元；UMFPACK；BDF；边界附近加密。连续/间歇干燥全曲线，R²>0.99；部分参数来自独立实验或关联式。
- 对应论断：强支持移动域中必须一致处理收缩速度、Jacobian/体积变化与表面 Robin 通量；支持 BDF 对切换边界工况。论文还显示等温模型会误判时变扩散率，非等温模型更合适。
- 局限：球体、理想各向同性收缩；其 `v_s` 是物理收缩速度。映射到 `ξ=r/R(t)` 后“显式拖曳项是否出现”取决于 Eulerian/ALE/严格材料坐标的变量定义，不能机械照抄。相关性 5.0/5；证据 A。

#### A3. Zambra et al. (2022), *Experimental and Numerical Study of a Turbulent Air-Drying Process for an Ellipsoidal Fruit with Volume Changes*

- 期刊/标识符：*Foods* 11(13), 1880; DOI [10.3390/foods11131880](https://doi.org/10.3390/foods11131880); [PMC 全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC9265362/)。双库元数据一致。
- 对象/几何：灯笼果（*Physalis peruviana*），3D 椭球，收缩主要沿一个方向；同时求解干燥室空气流场。
- 方程/边界：空气侧 RANS `k–ε`、热量和湿分；果体 3D Fourier/Fick 扩散；固–气共轭传递。设备壁绝热、不可渗，出口零梯度。
- 物性/收缩：空气物性随温度；果体 `D_eff` 为入口温度和平均含水率函数；基于图像/体积实验的各向异性离散收缩。
- 数值/验证：自编 Fortran 有限体积，隐式 Euler，SIMPLE，Gauss–Seidel/TDMA；3 套网格和 3 个时间步比较；验证干燥曲线、体积/形状，含收缩后相对误差显著降低。
- 对应论断：反驳“收缩必然各向同性”；提示端面/外部非均匀流会破坏纯径向对称；支持 FVM、网格/时间步检查和实验半径/体积驱动几何更新。
- 局限：`D_eff` 的预指数因子逐时用实验平均含水率校准，预测独立性有限；未采用 `ξ=r/R(t)`。相关性 4.7/5；证据 A。

#### A4. Xanthopoulos, Lentzou & Papadakis (2026), *Moving-boundary finite element modeling and inverse identification of moisture-dependent diffusivity in convective drying of apricot hemispheres: Uncertainty quantification*

- 期刊/标识符：*Drying Technology*; DOI [10.1080/07373937.2026.2686407](https://doi.org/10.1080/07373937.2026.2686407)。Crossref/OpenAlex 的题名、作者、2026 年一致；OpenAlex 记 closed OA。
- 对象/几何：杏半球，变形域。
- 方程/边界/物性：摘要确认内部水分扩散、有效界面传质阻力、ALE 收缩；`D_eff(T,MR)` 为 Arrhenius 型温度–含水率函数；未核验热方程和边界公式细节。
- 数值/验证：moving-boundary FEM + 非线性逆优化 + Monte Carlo 不确定性；拟合完整干燥曲线并复现实验收缩。
- 对应论断：是 2026 年直接支持“移动边界 + 含水率依赖扩散率 + 参数不确定性”的近期原始研究。
- 局限：摘要层证据；等温模型、半球而非圆柱；全文与 Unpaywall 尚未核验。相关性 4.8/5；证据 B。

#### B1. Yang, Sakai & Watanabe (2001), *Drying Model with Non-Isotropic Shrinkage Deformation Undergoing Simultaneous Heat and Mass Transfer*

- 期刊/标识符：*Drying Technology*; DOI [10.1081/DRT-100105299](https://doi.org/10.1081/DRT-100105299)。双库元数据一致。
- 对象/几何：圆柱马铃薯，二维轴–径向；中心温度、平均含水率、轴向/径向收缩实验。
- 方程/物性/收缩：摘要确认同时热湿传递与虚功原理耦合，非恒定物理/热物性，非各向同性收缩。
- 数值/验证：Galerkin FEM；单一热风工况实验，预测与实验吻合。边界公式未在可得摘要中核验。
- 对应论断：高度支持本题圆柱、轴向端部与径向收缩不能先验等同；为“仅一维径向+端面体源”提供必要的二维对照基准。
- 局限：老论文、摘要层；定量网格验证未核验。相关性 4.9/5；证据 B。

#### B2. Rahman & Kumar (2011), *Evaluation of Moisture Diffusion Coefficient of Cylindrical Bodies Considering Shrinkage During Natural Convection Drying*

- 期刊/标识符：*International Journal of Food Engineering*; DOI [10.2202/1556-3758.1352](https://doi.org/10.2202/1556-3758.1352)。双库元数据一致。
- 对象/几何：长 0.05 m、直径 0.01 m 的马铃薯有限圆柱；同时比较有限/无限圆柱及有/无收缩。
- 方程/边界/物性：解析 Fick 扩散模型；明确考虑固–气界面外部传质阻力；`D_eff` 拟合为温度与含水率函数；收缩纳入几何。
- 数值/验证：Levenberg–Marquardt 参数反演；自然对流干燥曲线验证。
- 对应论断：直接表明用无限圆柱代替有限圆柱会高估 `D_eff`，考虑收缩则推得更低 `D_eff`；支持端面效应不可无验证忽略。
- 局限：无内部温度场、自然对流、主要验证整体干燥曲线；**不直接证明把两个端面通量均匀折算为体源是正确的**。相关性 4.8/5；证据 B。

#### C1. Ning et al. (Crossref 卷期年 2026；OpenAlex 在线年 2025), *Modeling of heat and moisture transfer during convective drying of in-hull almonds considering non-isotropic shrinkage and variable moisture diffusivity*

- 期刊/标识符：*International Journal of Heat and Mass Transfer*; DOI [10.1016/j.ijheatmasstransfer.2025.128286](https://doi.org/10.1016/j.ijheatmasstransfer.2025.128286)。双库题名/作者一致；年份口径差异已保留。
- 对象/几何/模型：题名层确认带壳杏仁、热湿传递、非各向同性收缩、变水分扩散率。
- 边界、物性具体式、数值法和验证：**未核验**（OpenAlex/Crossref 无摘要；OA closed；Unpaywall 待授权）。
- 对应论断：近期且主题高度贴合，建议第二轮优先获取摘要/全文。
- 局限：本轮只能作为高优先级候选，不能据题名填充方程细节。相关性 4.9/5；证据 C。

#### C2. Das et al. (2025), *Finite element modeling of banana slices in convective drying, studies on isotropic shrinkage kinetics using Arbitrary Lagrangian–Eulerian (ALE) approach*

- 期刊/标识符：*Thermochimica Acta*; DOI [10.1016/j.tca.2025.180024](https://doi.org/10.1016/j.tca.2025.180024)。双库元数据一致，OA closed。
- 题名层确认：香蕉片、对流干燥、FEM、ALE、各向同性收缩。
- 方程、边界、物性、验证细节：**未核验**。
- 对应论断：用于补齐 2015–2026 的 ALE 应用；不能仅凭题名判断有无坐标拖曳项。相关性 4.4/5；证据 C。

#### C3. Aprajeeta, Gopirajah & Anandharamakrishnan (Crossref 卷期年 2015；OpenAlex 在线年 2014), *Shrinkage and porosity effects on heat and mass transfer during potato drying*

- 期刊/标识符：*Journal of Food Engineering*; DOI [10.1016/j.jfoodeng.2014.08.004](https://doi.org/10.1016/j.jfoodeng.2014.08.004)。双库题名/作者一致。
- 题名层确认：马铃薯干燥中收缩、孔隙率对热湿传递的影响。
- 几何、边界、物性具体式、数值法和验证：**未核验**。
- 对应论断：提醒仅改变半径而不更新孔隙率/有效物性可能不闭合。相关性 4.3/5；证据 C。

#### C4. Koukouch et al. (2020), *Analytical solution of coupled heat and mass transfer equations during convective drying of biomass: experimental validation*

- 期刊/标识符：*Heat and Mass Transfer*; DOI [10.1007/s00231-020-02817-w](https://doi.org/10.1007/s00231-020-02817-w)。双库元数据一致。
- 题名层确认：生物质对流干燥的耦合热湿方程、解析解、实验验证。
- 几何、Robin 形式、物性和解析假设：**未核验**。
- 对应论断：第二轮可作为本题解析基准候选；在未读方法前不能宣称适用于圆柱变系数问题。相关性 4.3/5；证据 C。

#### C5. Murugesan et al. (2007), *Convective drying analysis of three-dimensional porous solid by mass lumping finite element technique*

- 期刊/标识符：*Heat and Mass Transfer*; DOI [10.1007/s00231-007-0260-9](https://doi.org/10.1007/s00231-007-0260-9)。双库元数据一致。
- 题名层确认：3D 多孔固体对流干燥、质量集中的 FEM。
- 方程、边界、物性、收缩与验证：**未核验**。
- 对应论断：可用于检查三维/降维和质量矩阵处理，但不能据题名作为当前 FVM 细节依据。相关性 4.0/5；证据 C。

### 4.2 有限体积、时间推进、验证与插值

#### N1. Eymard, Gallouët & Herbin (2000), *Finite Volume Methods*

- 来源/标识符：*Handbook of Numerical Analysis*; DOI [10.1016/S1570-8659(00)07005-8](https://doi.org/10.1016/S1570-8659(00)07005-8); [开放稿](https://hal.science/hal-02100732v2/file/bookevol.pdf)。双库元数据一致。
- 内容/对应点：权威方法章节，支持以控制体积分建立守恒离散、面通量配对和离散一致性/稳定性分析；适合作为节点/单元中心 FVM 的总方法依据。
- 局限：本轮未逐页核验其圆柱中心半控制体公式；不能单独证明端面等效体源或某一种非线性面平均最优。相关性 4.6/5；证据 B。

#### N2. Zhang & Shan (2021), *Finite volume schemes of evolutionary diffusion equations based on harmonic averaging interpolation*

- 期刊/标识符：*Scientia Sinica Mathematica*; DOI [10.1360/SSM-2020-0262](https://doi.org/10.1360/SSM-2020-0262)。双库元数据一致。
- 方程/方法/验证：cell-centered FVM；调和平均点插值；backward Euler；对连续、间断、异质扩散张量以及非线性扩散做稳定性、误差估计和多网格数值实验；报告部分格式 L² 二阶、H¹ 一阶。
- 对应论断：支持把调和型处理纳入变系数扩散面通量对比，并要求用误差阶/鲁棒性而非“曲线看起来接近”选择。
- 局限：一般非结构网格与张量问题，不等于“均匀径向网格上的所有平滑 `D(X,T)` 都应使用调和平均”。相关性 4.5/5；证据 B。

#### N3. Gottlieb, Shu & Tadmor (2001), *Strong Stability-Preserving High-Order Time Discretization Methods*

- 期刊/标识符：*SIAM Review*; DOI [10.1137/S003614450036757X](https://doi.org/10.1137/S003614450036757X)。双库元数据一致。
- 内容/对应点：SSP 方法把 forward Euler 在给定步长下的强稳定性质传递到高阶 Runge–Kutta；为 Heun/SSPRK2 的方法定位与稳定性测试提供依据。
- 局限：重点是方法线/双曲问题的一般 SSP 理论；对扩散半离散仍需满足显式抛物型步长限制，不能用“SSP”替代 `Δt~Δr²/D_max` 检查。相关性 4.8/5；证据 B。

#### N4. Roache (2001/卷期 2002), *Code Verification by the Method of Manufactured Solutions*

- 期刊/标识符：*Journal of Fluids Engineering* 124(1); DOI [10.1115/1.1436090](https://doi.org/10.1115/1.1436090)。Crossref/OpenAlex 均记 2001（卷期通常列 2002，需按引用样式保留 online/issue 差别）。
- 内容/对应点：明确区分代码验证与计算误差估计；MMS 与系统网格加密结合，可检测离散和实现错误。
- 局限：MMS 验证“是否正确求解所写方程”，不验证方程是否适合中药材，也不替代实验 validation。相关性 5.0/5；证据 B。

#### N5. Curtiss & Hirschfelder (1952), *Integration of Stiff Equations*

- 期刊/标识符：*PNAS* 38(3); DOI [10.1073/pnas.38.3.235](https://doi.org/10.1073/pnas.38.3.235); [开放 PDF](https://pmc.ncbi.nlm.nih.gov/articles/PMC1063538/pdf/pnas01576-0089.pdf)。双库元数据一致。
- 内容/对应点：刚性常微分方程隐式多步积分的奠基文献，可作为 BDF 思路的历史方法依据。
- 局限：OpenAlex 的所谓摘要实际是期刊说明，不是论文摘要；现代 BDF 自适应实现、误差控制和 PDE 半离散比较需另找近期求解器文献。相关性 3.9/5；证据 C。

#### N6. Fritsch & Carlson (1980), *Monotone Piecewise Cubic Interpolation*

- 期刊/标识符：*SIAM Journal on Numerical Analysis*; DOI [10.1137/0717021](https://doi.org/10.1137/0717021)。双库元数据一致。
- 内容/对应点：给出三次多项式区间单调的充要条件并构造单调分段三次插值；支持避免普通 cubic spline 对半径数据产生过冲。
- 局限：若实验半径序列本身因噪声非单调，该算法不等于先验强制整体收缩单调；需先说明是否做物理约束/去噪。相关性 4.5/5；证据 B。

#### N7. Fritsch & Butland (1984), *A Method for Constructing Local Monotone Piecewise Cubic Interpolants*

- 期刊/标识符：*SIAM Journal on Scientific and Statistical Computing*; DOI [10.1137/0905021](https://doi.org/10.1137/0905021)。双库元数据一致。
- 内容/对应点：完全局部、简单实现的单调分段三次斜率构造；与常用 PCHIP 的局部斜率思想更直接。
- 局限：同样只对单调输入保持单调；插值误差应通过留点或半径测量不确定度评估。相关性 4.8/5；证据 B。

### 4.3 单库待复核候选（不进入核心依据）

- Hiroaki Nishikawa & Boris Diskin (2022), *Arithmetic Averages of Viscosity Coefficient are Sufficient for Second-Order Finite-Volume Viscous Discretization on Unstructured Grids*, arXiv:[2203.08334](https://arxiv.org/abs/2203.08334)。arXiv 摘要给出：算术平均在其 cell-centered viscous discretization 中可保持二阶，并用 1D 非线性扩散和 3D MMS 展示。Crossref 精确题名未找到同一记录，Semantic Scholar 又限流，因此作者/年份第二库核验未完成。它可用来**反驳“调和平均普遍优于算术平均”这种过强说法**，但问题设置不是圆柱干燥。

## 5. 对当前模型的第一轮证据标记（不是综述结论）

1. **控制方程与 Robin 边界：强支持。** Teleken 2025 和 Adrover 2020 直接给出内部 Fourier/Fick 方程、表面传热/传质阻力及蒸发潜热耦合。
2. **收缩坐标中的“拖曳项”：不能简单删或加。** Adrover 2020 在空间坐标/移动域中明确保留收缩速度导致的守恒型对流项。推论是：若用 `ξ=r/R(t)` 的固定计算坐标，会出现网格速度/Jacobian 项；若用严格材料坐标并把未知量定义为单位参考体积（或参考质量）量，显式拖曳项可能被 Jacobian 吸收。**省略显式项同时也不更新 Jacobian/体积测度通常不守恒。** 这是基于守恒变换的推论，不是某篇论文对本题公式的逐字结论。
3. **端面均匀等效体源：本轮未找到直接文献背书。** Rahman–Kumar 2011 证明有限/无限圆柱会导致不同的反演扩散率；Yang 2001 提供轴–径向二维基准。当前“两个端面通量轴向平均后均匀加到一维径向方程”的做法应作为待验证的降阶闭合：至少与二维轴对称解比较局部最大误差和全局守恒误差。
4. **节点中心 FVM 与半控制体：一般守恒思想有依据，具体圆柱公式仍需第二轮。** Eymard et al. 支持控制体守恒框架；本轮针对“圆柱中心/表面半控制体”的检索没有得到一篇能逐式核对的核心论文。
5. **面扩散系数平均没有普遍赢家。** Zhang–Shan 支持调和型插值在异质/间断与非线性扩散中的分析和测试；Nishikawa–Diskin 单库记录显示算术平均在另一类 cell-centered 粘性离散中也可二阶。应按本题连续性、非线性与保守通量推导，再用 MMS/解析解比较误差阶，不能仅看某次曲线接近。
6. **Heun/SSPRK2 + BDF 的验证组合方向正确，但需独立性。** SSPRK2 仍受扩散 CFL 限制；BDF 与 Heun 相互接近只说明两种实现一致，强代码验证仍应包括 MMS、解析特例和系统网格/时间步阶次。实验拟合属于模型 validation，不等同 code verification。
7. **PCHIP 合理但不自动施加物理单调性。** 对已单调的半径数据可抑制过冲；若测量噪声造成反弹，应单独说明是否先做单调回归/物理约束，并传播半径不确定度。
8. **最大局部干基含水率阈值是更保守的终止量。** 本轮未检索到支持 `0.15 kg/kg` 这一具体值的通用论文；它应被表述为题设/工艺准则，而非由上述模型论文推出。

## 6. 排除与保留记录

- 排除：只有平均含水率的薄层经验拟合；微波-only、冻干、喷雾干燥；只做干燥器空气流场；机器学习-only；仅以总时长接近作为正确性证据。
- 保留但降级：含内部扩散模型、有限几何或收缩反演，即使只验证平均干燥曲线（如 Rahman–Kumar），因为它们能检验端面/收缩导致的参数偏差；证据等级相应降低。
- 综述（如 Mayor & Sereno 2004）仅用于术语和引用网络，本轮未作为核心控制方程依据。

## 7. 本轮缺口与第二轮优先级

1. 经用户明确授权后，用真实邮箱逐 DOI 查询 Unpaywall；当前 OA 链接仅来自 Europe PMC/OpenAlex。
2. 获取 Ning 2025/2026、Das 2025、Aprajeeta 2014/2015、Koukouch 2020、Murugesan 2007 的摘要或合法全文，补齐边界式、几何、离散和验证字段。
3. 定向检索“finite-length cylinder → 1D radial source closure / aspect-ratio asymptotics / volume averaging”，并寻找可逐式验证端面体源的数学论文。
4. 定向检索圆柱坐标 cell-centered FVM 在 `r=0` 与 Robin 表面的半控制体离散及收敛证明。
5. 对 `ξ=r/R(t)` 分别按 Eulerian ALE 与材料守恒变量推导，建立同一物理问题的守恒等价形式，再用 MMS/二维轴对称解验证。

## 8. 检索方法引用

本报告使用了 Paper Lookup 技能的数据库路由、静默失败检查和可复现记录流程：

Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065. [https://doi.org/10.48550/arXiv.2609.00065](https://doi.org/10.48550/arXiv.2609.00065)
