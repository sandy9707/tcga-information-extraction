# !准备
# 清除当前环境中的所有对象
rm(list = ls())

# 设置文件路径
root_dir <- paste0("~", "/Rstudio")

# 加载 librarian 包
library(librarian)

# 安装并载入 obgetDEGs 包
# devtools::install_github('sandy9707/obgetDEGs')
shelf(obgetDEGs, cran_repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

# 设置输出目录
output_dir <- "TCGA-BRCA/output"
# shelfEnvironment(output_dir, path = root_dir)

# 设置输入目录
imput_dir <- "/Users/sandy/Downloads/tcga/bcca/exp"

# 载入 openxlsx、dplyr 和 jsonlite 包
shelf(openxlsx, dplyr, jsonlite, cran_repo = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
