# tcga-information-extraction

# 注释

可以提取TCGA数据库中的表达矩阵和临床信息, 其中表达矩阵可以含counts, tpm等值

# Description

The expression matrix and clinical information in TCGA database can be extracted, and the expression matrix can contain counts, tpm, etc

# 使用指南

## prepare.R

主要用来设置root_dir, output_dir, imput_dir三个工作目录,
其中,
root_dir为根目录, 一般设置为本地账户登入目录的"RStudio", 可自行更改但不推荐,
output_dir为输出路径, 可自行编辑文件夹名及目标目录, 默认创建至root_dir下,
imput_dir为数据读入目录, 建议使用全路径.


如缺少obgetDEGs包可使用`devtools::install_github('sandy9707/obgetDEGs')`命令从github安装(科学上网最佳).