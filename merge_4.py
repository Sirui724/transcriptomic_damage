import pandas as pd
import numpy as np
import re

# ======================================================
# 步骤1：读取并处理每个文件
# ======================================================

# 读取junc_gene_ratio_result.txt
junc = pd.read_csv('junc_gene_ratio_result.txt', sep='\t', index_col=0)
junc.index = junc.index.astype(str) + '_junc'  # 添加_junc后缀

# 读取stopsite_ratio_result.txt
stop = pd.read_csv('stopsite_ratio_result.txt', sep='\t', index_col=0)
stop.index = stop.index.astype(str) + '_stop'  # 添加_stop后缀

# 读取normalized_merge_re.txt
merge_re = pd.read_csv('normalized_merge_re.txt', sep='\t', index_col='name')

# ======================================================
# 步骤2：统一列名格式（可选）
# 如果列名存在类似"626-1a"需要转为"626-1"的情况，执行：
# ======================================================
def clean_colname(col):
    """将类似626-1a转为626-1"""
    return re.sub(r'([0-9]+-[0-9]+)[a-zA-Z]*', r'\1', col)

for df in [junc, stop, merge_re]:
    df.columns = [clean_colname(col) for col in df.columns]

# ======================================================
# 步骤3：合并所有数据
# ======================================================
merged = pd.concat([junc, stop, merge_re], axis=0)

# ======================================================
# 步骤4：保存结果
# ======================================================
merged.to_csv('merged_features.txt', sep='\t', na_rep='NA')

print("合并完成！结果文件：merged_features.txt")
print("最终数据维度：", merged.shape)
