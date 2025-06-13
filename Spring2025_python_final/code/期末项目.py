# 以下代码主要利用的是从DepMap下载的数据，初步探索结直肠癌中KRAS突变的重要性

# ①：在结直肠癌中，检测KRAS的敲除效应分数情况
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

data_dir = r"C:\Users\31904\Desktop\shangke\python\qimo"

crispr = pd.read_csv(os.path.join(data_dir, "CRISPRGeneEffect.csv"))
cnv = pd.read_csv(os.path.join(data_dir, "OmicsCNGene.csv"))
models = pd.read_csv(os.path.join(data_dir, "Model.csv"))

crispr['ModelID'] = crispr['Unnamed: 0']
cnv['ModelID'] = cnv['Unnamed: 0']

kras_crispr_col = [col for col in crispr.columns if col.startswith('KRAS')][0]
kras_cnv_col = [col for col in cnv.columns if col.startswith('KRAS')][0]

df_crispr = crispr[['ModelID', kras_crispr_col]].rename(columns={kras_crispr_col: 'KRAS_GeneEffect'})
df_cnv = cnv[['ModelID', kras_cnv_col]].rename(columns={kras_cnv_col: 'KRAS_CNV'})
df_models = models[['ModelID', 'OncotreePrimaryDisease']]

merged = df_crispr.merge(df_cnv, on='ModelID').merge(df_models, on='ModelID')
merged = merged.dropna()

merged['Group'] = merged['OncotreePrimaryDisease'].apply(
    lambda x: 'Colorectal Adenocarcinoma' if x == 'Colorectal Adenocarcinoma' else 'Other')

# 图一：CRC中的KRAS基因敲除影响
plt.figure(figsize=(10, 7))
sns.scatterplot(
    data=merged, 
    x='KRAS_CNV', 
    y='KRAS_GeneEffect', 
    hue='Group', 
    palette={'Colorectal Adenocarcinoma': 'red', 'Other': 'blue'},
    alpha=0.7
)
plt.axhline(-0.5, linestyle='--', color='gray', label=' -0.5')
plt.xlabel(' Copy Number of KRAS')
plt.ylabel('KRAS CRISPR GeneEffect Score')
plt.title('KRAS GeneEffect — CRC vs Other Cancers')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


# ②：结直肠癌中KRAS突变类型分布TOP，KRASG12D突变占比
import pandas as pd
import os
import matplotlib.pyplot as plt

data_dir = r"C:\Users\31904\Desktop\shangke\python\qimo"

mutations = pd.read_csv(os.path.join(data_dir, "OmicsSomaticMutations.csv"))
models = pd.read_csv(os.path.join(data_dir, "Model.csv"))

crc_models = models[models['OncotreePrimaryDisease'] == 'Colorectal Adenocarcinoma']
crc_model_ids = crc_models['ModelID'].tolist()
print(f"结直肠癌细胞系总数: {len(crc_model_ids)}")

kras_mutations = mutations[mutations['HugoSymbol'] == 'KRAS']

crc_kras_mutations = kras_mutations[kras_mutations['ModelID'].isin(crc_model_ids)]

kras_mut_counts = crc_kras_mutations['ProteinChange'].value_counts()
print("\nKRAS突变类型分布:")
print(kras_mut_counts.head(10)) 

g12d_mutations = crc_kras_mutations[crc_kras_mutations['ProteinChange'] == 'p.G12D']
g12d_models = g12d_mutations['ModelID'].unique()
print(f"\nKRAS G12D突变的结直肠癌细胞系数量: {len(g12d_models)}")
print(f"KRAS G12D突变在结直肠癌细胞系中的占比: {len(g12d_models)/len(crc_model_ids):.2%}")

kras_mutated_models = crc_kras_mutations['ModelID'].unique()
print(f"\n携带KRAS突变的结直肠癌细胞系数量: {len(kras_mutated_models)}")
print(f"KRAS突变在结直肠癌细胞系中的总占比: {len(kras_mutated_models)/len(crc_model_ids):.2%}")
print(f"KRAS G12D在所有KRAS突变结直肠癌细胞系中的占比: {len(g12d_models)/len(kras_mutated_models):.2%}")

# 图二：结直肠癌细胞系,KRAS突变类型TOP10
plt.figure(figsize=(12, 7))
top10_mutations = kras_mut_counts.head(10)

bars = plt.bar(top10_mutations.index, top10_mutations.values, color='skyblue')
plt.title('Colorectal Cancer Cell Lines: Top 10 KRAS Mutations', fontsize=15)  
plt.xlabel('Mutation Type', fontsize=12)  
plt.ylabel('Number of Cell Lines', fontsize=12)  
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.grid(axis='y', linestyle='--', alpha=0.7)

for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 0.1,
             f'{height:.0f}',
             ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.show()

# 图三：展示G12D vs 其他KRAS突变 vs 野生型的占比
labels = ['KRAS G12D', 'other KRAS mutation', 'KRAS WT']
other_kras = len(kras_mutated_models) - len(g12d_models)
wt = len(crc_model_ids) - len(kras_mutated_models)
sizes = [len(g12d_models), other_kras, wt]

plt.figure(figsize=(8, 8))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#ff9999','#66b3ff','#99ff99'])
plt.axis('equal')
plt.title('Distribution of KRAS mutation status in colorectal cancer cell lines')
plt.tight_layout()
plt.show()


# ③：KRAS敲除分数：在结直肠癌中，KRAS G12D突变细胞对比KRAS野生型细胞
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

data_dir = r"C:\Users\31904\Desktop\shangke\python\qimo"

mutations = pd.read_csv(os.path.join(data_dir, "OmicsSomaticMutations.csv"))
models = pd.read_csv(os.path.join(data_dir, "Model.csv"))
crispr = pd.read_csv(os.path.join(data_dir, "CRISPRGeneEffect.csv"))

crc_models = models[models['OncotreePrimaryDisease'] == 'Colorectal Adenocarcinoma']
crc_model_ids = crc_models['ModelID'].tolist()
print(f"结直肠癌细胞系总数: {len(crc_model_ids)}")

kras_mutations = mutations[mutations['HugoSymbol'] == 'KRAS']
crc_kras_mutations = kras_mutations[kras_mutations['ModelID'].isin(crc_model_ids)]

g12d_mutations = crc_kras_mutations[crc_kras_mutations['ProteinChange'] == 'p.G12D']
g12d_models = g12d_mutations['ModelID'].unique()
print(f"\nKRAS G12D突变的结直肠癌细胞系数量: {len(g12d_models)}")

crc_kras_mut_models = crc_kras_mutations['ModelID'].unique()

crc_kras_wt_models = list(set(crc_model_ids) - set(crc_kras_mut_models))

kras_col = [col for col in crispr.columns if col.startswith("KRAS")][0]
crispr['ModelID'] = crispr['Unnamed: 0']  

g12d_df = crispr[crispr['ModelID'].isin(g12d_models)][['ModelID', kras_col]].copy()
g12d_df['Group'] = 'KRAS G12D'

wt_df = crispr[crispr['ModelID'].isin(crc_kras_wt_models)][['ModelID', kras_col]].copy()
wt_df['Group'] = 'KRAS WT'

compare_df = pd.concat([g12d_df, wt_df], ignore_index=True)
compare_df = compare_df.rename(columns={kras_col: 'GeneEffect'})
compare_df = compare_df.dropna()

# 图4：结直肠癌中KRAS敲除分数在KRASG12D细胞中相较于KRAS野生型细胞的情况
plt.figure(figsize=(8, 6))
sns.boxplot(x='Group', y='GeneEffect', data=compare_df, palette='Set2')
plt.axhline(-0.5, color='gray', linestyle='--', label=' -0.5')
plt.title('KRAS G12D vs WT (In CRC)')
plt.ylabel('Gene Effect Score (KRAS)')
plt.legend()
plt.grid(True)
plt.show()


# ④：对KRAS G12D CRC细胞与 KRAS WT CRC细胞进行分组，差异基因分析
import pandas as pd
import numpy as np
import os
from scipy.stats import ttest_ind
import statsmodels.stats.multitest as smm
import matplotlib.pyplot as plt

data_dir = r"C:\Users\31904\Desktop\shangke\python\qimo"

expr = pd.read_csv(os.path.join(data_dir, "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), index_col=0)

models = pd.read_csv(os.path.join(data_dir, "Model.csv"))
mutations = pd.read_csv(os.path.join(data_dir, "OmicsSomaticMutations.csv"))

crc_models = models[models['OncotreePrimaryDisease'] == 'Colorectal Adenocarcinoma']['ModelID']

kras_mut = mutations[(mutations['HugoSymbol'] == 'KRAS') & (mutations['ModelID'].isin(crc_models))]

g12d_models = kras_mut[kras_mut['ProteinChange'] == 'p.G12D']['ModelID'].unique()
wt_models = crc_models[~crc_models.isin(kras_mut['ModelID'])]  
group_g12d = [m for m in g12d_models if m in expr.index]
group_wt = [m for m in wt_models if m in expr.index]

print(f"G12D组细胞系数量: {len(group_g12d)}")
print(f"WT组细胞系数量: {len(group_wt)}")

results = []
for gene in expr.columns:
    g1 = expr.loc[group_g12d, gene].dropna()
    g2 = expr.loc[group_wt, gene].dropna()
    if len(g1) >= 3 and len(g2) >= 3:  
        stat, p = ttest_ind(g1, g2, equal_var=False)  
        log2fc = g1.mean() - g2.mean()  
        results.append([gene, log2fc, p, len(g1), len(g2)])

deg = pd.DataFrame(results, columns=["Gene", "log2FC", "pvalue", "n_g12d", "n_wt"])
deg['FDR'] = smm.multipletests(deg['pvalue'], method='fdr_bh')[1]
deg = deg.sort_values("FDR")

significant_genes = deg[deg['FDR'] < 0.05]
print(f"显著基因数量：{len(significant_genes)}")
print(significant_genes[['Gene', 'log2FC', 'pvalue', 'FDR']])

# 保存结果
#deg.to_csv(os.path.join(data_dir, "DEG_KRAS_G12D_vs_WT.csv"), index=False)
#print("差异表达分析结果已保存")

# 由于FDR过滤后没有显著的基因，所以用pvalue作为初步筛选标准看一下
top_candidates = deg.sort_values("pvalue").head(20)
print("\nTop candidate genes based on raw p-value:")
print(top_candidates[['Gene', 'log2FC', 'pvalue']])

# 图5：用火山图展示差异分析结果
def get_color_and_size(row):
    if row['pvalue'] < 0.01 and row['log2FC'] > 1:
        return 'red', 10  
    elif row['pvalue'] < 0.01 and row['log2FC'] < -1:
        return 'blue', 10  
    else:
        return 'gray', 5  

colors, sizes = zip(*deg.apply(get_color_and_size, axis=1))

plt.figure(figsize=(7,6))
plt.scatter(deg["log2FC"], -np.log10(deg["pvalue"]), s=sizes, c=colors)

plt.axvline(x=1, linestyle='--', color='black')
plt.axvline(x=-1, linestyle='--', color='black')
plt.axhline(y=-np.log10(0.05), linestyle='--', color='black')
plt.scatter([], [], color='red', s=10, label='Upregulated (p<0.01, FC>2)')
plt.scatter([], [], color='blue', s=10, label='Downregulated (p<0.01, FC<0.5)')
plt.scatter([], [], color='gray', s=5, label='Not significant')

plt.xlabel("log2 Fold Change")
plt.ylabel("-log10(p-value)")
plt.title("Volcano Plot (KRAS G12D VS KRAS WT)")
plt.legend(
    loc='center left',
    bbox_to_anchor=(-0.002, 0.925),
    frameon=True,
    fontsize=9
)
plt.tight_layout()
plt.show()


