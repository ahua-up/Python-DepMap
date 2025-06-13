# Sustech-生物医学Python编程入门期末项目报告  

## project  

基于DepMap数据初步探究KRAS突变在结直肠癌治疗中的重要意义  

## 1. Background  

KRAS蛋白是一类GTP/GDP结合蛋白，KRAS基因能够产生形成两个剪接变体 (KRAS 4A and KRAS 4B)。其中，KRAS 4B由于在人类癌症中广泛存在和高表达被认为是主要亚型。并且近年来，KRAS 4A也逐渐被证明在各种癌症中普遍表达，并能增加肿瘤细胞在应激条件下的适应性。KRAS蛋白的分子量为21kDa，由六条β链和五个α螺旋组成，形成两个主要结构域：G结构域和C末端 。G结构域高度保守，包含开关I和开关II环，负责GDP-GTP交换。C末端是一个可变区域，包括CAAX(C = 半胱氨酸，A = 任何脂肪氨基酸，X = 任何氨基酸) 基序，是各种翻译后修饰的靶点，在新合成和处理的 KRAS 运输以及最终质膜锚定中起着至关重要的作用。  

KRAS蛋白作为一个 “分子开关”, 在GDP结合的非活性状态和GTP结合的活性状态之间循环。当GTP与KRAS的结合，会触发几个下游途径，例如快速激活的RAF-MEK-ERK和PI3K-AKT-mTOR途径，这些途径促进细胞生长和存活。而当GDP结合KRAS，KRAS会失去其活性并阻止其持续的信号转导激活。  

因此，一旦KRAS发生突变，GTP的水解被破坏或核苷酸交换被增强，KRAS就会以活性状态积累，活性KRAS持续激活下游信号通路，从而促进肿瘤细胞增殖与生存。临床研究表明，结直肠癌携带KRAS突变，与晚期疾病状态差、肿瘤分化不良、远端转移、患者生存率低有关。  

DepMap数据库，是一个关于癌症细胞的大型开放式数据库，它提供了大量关于癌症细胞系、基因突变、基因表达、药物敏感性等相关信息。可以帮助我们更好地理解癌症的机制并寻找潜在的治疗靶点。因此，通过DepMap的数据，我们可以初步探究KRAS突变对与结直肠癌研究的意义，为后续结直肠癌的研究提供可能的方向。  


## 2. Method  

首先，在DepMap数据库上，我们需要下载好Model（细胞模型信息）、OmicsCNGene（基因拷贝信息）、OmicsSomaticMutations（细胞系突变信息）、CRISPRGeneEffect（基因CRISPR效应评分）和OmicsExpressionProteinCodingGenesTPMLogp1（标准化的基因表达数据）。  

而后在该项目中，我们会使用相应的python代码来处理depmap下载的数据，分别得到相应的结果。  
方法简单介绍如下：  

* 结合CRISPRGeneEffect和OmicsCNGene的数据，我们可以利用py文件中的代码①处理数据，分析得到KRAS在不同器官细胞系的表达情况和结直肠癌细胞对KRAS的依赖情况。  

* OmicsSomaticMutations数据中储存有各个细胞系的突变数据，我们可以利用py文件的代码②，进行处理，可以得到结直肠癌中不同的KRAS突变类型与占比情况。  

* OmicsSomaticMutations和CRISPRGeneEffect，同样可以利用代码③进行结合分析，从而得知KRAS的敲除是否对于KRAS突变的结直肠癌细胞产生显著影响。  

* 最后利用OmicsExpressionProteinCodingGenesTPMLogp1这个文件，我们可以选择某类KRAS突变类型的结直肠癌细胞和KRAS野生型的结直肠癌细胞进行差异基因分析  


## 3. Result  

根据不同处理方法，通过python代码，我们可以得到以下类似结果：  
① 相较于其它肿瘤细胞，CRC细胞KRAS敲除的效应评分：  
[![pVAPMUe.png](https://s21.ax1x.com/2025/06/13/pVAPMUe.png)](https://imgse.com/i/pVAPMUe)  

② 结直肠癌细胞的KRAS突变类型TOP10，KRAS突变的结直肠癌细胞占比：  
[![pVAP1Cd.png](https://s21.ax1x.com/2025/06/13/pVAP1Cd.png)](https://imgse.com/i/pVAP1Cd)  

[![pVAP8gI.png](https://s21.ax1x.com/2025/06/13/pVAP8gI.png)](https://imgse.com/i/pVAP8gI)  

③ KRASG12D的结直肠癌细胞系CRISPR效应评分 VS KRAS野生型的结直肠癌细胞系CRISPR效应评分：  
[![pVAPsvq.png](https://s21.ax1x.com/2025/06/13/pVAPsvq.png)](https://imgse.com/i/pVAPsvq)  

④ 差异基因分析：KRAS G12D结直肠癌细胞 VS KRAS野生型结直肠癌细胞：  
[![pVAPgbT.png](https://s21.ax1x.com/2025/06/13/pVAPgbT.png)](https://imgse.com/i/pVAPgbT)  


## 4. Discussion  

本项目基于DepMap数据库，分析得到了KRAS基因在结直肠癌中的突变情况、结直肠癌细胞对KRAS的功能依赖性和其它可能受到KRAS影响的基因。  

我进行了初步分析，发现：尽管KRAS在CRC细胞中的拷贝数没有显著扩增，但其敲除对细胞生存影响显著，提示突变型KRAS在CRC中可能是驱动基因。而在结直肠癌细胞中，G12D是其最常见的突变类型，携带该突变的细胞对KRAS敲除更加敏感。  

并且，利用相应的数据与代码，进一步进行差异表达分析可以得到一些可能参与KRAS G12D突变功能通路的基因。这为后续研究KRAS的功能提供了相应的依据。  

欢迎大家使用或参考该项目中的分析流程，并结合自己的数据进行扩展。




