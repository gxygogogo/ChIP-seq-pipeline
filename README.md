# ChIP-seq流程
## 染色质免疫沉淀-测序
染色质免疫沉淀-测序（英语：ChIP-sequencing，简称为ChIP-seq）被用于分析蛋白质与DNA的交互作用。该技术将染色质免疫沉淀（ChIP）与大规模并行DNA测序结合起来以鉴定与DNA相关蛋白的结合部位。其可被用于精确绘制任意目的蛋白在全基因组上的结合位点。在此之前，ChIP-on-chip是研究这些蛋白-DNA联系的最常用的技术。
![概况](https://github.com/gxygogogo/ChIP-seq-pipeline/blob/main/img/overview.png)
## ChIP-seq应用
ChIP-seq主要用于确定转录因子和其他染色质相关蛋白质如何影响表型影响机制。确定蛋白质如何与DNA相互作用来调控基因表达对于充分了解许多生物过程和疾病状态至关重要。这种表观遗传信息与基因型和表达分析是互补的。 ChIP-seq技术目前主要被看作是需要杂交阵列的ChIP-on-Chip的替代品。 这必然会引入一些偏差，因为一个阵列仅限于有固定数量的探针。相反，尽管不同测序技术的测序偏差尚未完全了解，测序被认为偏差较小。
与转录因子和其他蛋白质直接物理相互作用的特定DNA位点可通过染色质免疫沉淀被分离出来。 ChIP产生在活体内In vivo兴趣蛋白质结合的目标DNA位点文库。 大规模平行序列分析与全基因组序列数据库结合使用来分析任何蛋白质与DNA的相互作用模式或任何表观遗传染色质修饰的模式。 这可以应用于一系列ChIP蛋白和修饰，如转录因子，聚合酶和转录机制，结构蛋白质，蛋白修饰和DNA修饰。 作为对特异性抗体依赖性的替代方法，已经开发了不同的方法来发现基因组中所有核小体耗尽的(nucleosome-depleted)或核小体扰动(nucleosome-disrupted)的活性调控区的超集，如DNase-Seq和FAIRE-Seq。
## 分析流程
### Step1. Raw data质控
### Step2. 比对分析
### Step3. 比对统计
### Step4. peak calling以及观察
### Step5. peak基因组注释
### Step6. peak下游分析
