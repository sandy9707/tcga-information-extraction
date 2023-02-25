source("/Users/sandy/Rstudio/TCGA-BRCA/code/prepare.R") # 准备配置

# !1.从文件获取总矩阵-------------------------------------------------------
# 读取 GDC 数据下载清单
shelfEnvironment(imput_dir, path = "")
manifest_matrix <- read.table("gdc_manifest_20230217_131521.txt",
  sep = "\t", header = TURE, stringsAsFactors = FALSE
)

# 在查看器中查看清单
# View(manifest_matrix)

# 获取当前目录中的文件名
dirnames <- list.files() %>% data.frame()

# 查找清单中存在但当前目录中不存在的文件名
setdiff(manifest_matrix[, 1], dirnames[, 1])

# 查找具有 .json 扩展名的文件
list.files(pattern = ".json")

# 读取 json 文件
json <- fromJSON(txt = "metadata.cart.2023-02-17.json")
# View(json)
# 创建空列表
MMRF_ID <- list() # 用于存储 MMRF_ID
matrix_MMRF <- list() # 用于存储表达矩阵
filedir_in_json <- list() # 用于存储 file_id

# 提取 MMRF_ID 和 file_id
for (i in seq_len(nrow(json))) {
  MMRF_ID[i] <- json[[3]][[i]] # MMRF_ID 存储到列表中
}
MMRF_ID <- do.call(rbind, MMRF_ID) # 将列表转换为矩阵

for (i in seq_len(nrow(json))) {
  filedir_in_json[i] <- json[i, ]["file_id"] # file_id 存储到列表中
}
filedir_in_json <- do.call(rbind, filedir_in_json) %>% as.character() # 将列表转换为字符向量

# 提取表达量至一个数据框(以tibble格式),counts值选4,fpkm选8,tpm选7
exctract_MMRF <- function(filedir_in_json, exct_num) {
  lapply(filedir_in_json, function(x) {
    setwd(paste(imput_dir, x, sep = "/")) # 切换目录到 file_id 所在的目录
    read.table(file = list.files(pattern = ".tsv"), header = F, fill = TRUE)[, exct_num] # 读取文件中的第exct_num列作为表达矩阵，存储到列表中
  })
}
matrix_MMRF <- exctract_MMRF(filedir_in_json = filedir_in_json, exct_num = 7)

# 输出待处理数据
shelfEnvironment(output_dir, path = root_dir)
save(list = c("filedir_in_json", "matrix_MMRF", "MMRF_ID"), file = "matrix_MMRF_tpm.RData") # 保存数据

# counts值, counts值选4,fpkm选8,tpm选7
matrix_MMRF <- exctract_MMRF(filedir_in_json = filedir_in_json, exct_num = 4)
shelfEnvironment(output_dir, path = root_dir)
save(list = c("filedir_in_json", "matrix_MMRF", "MMRF_ID"), file = "matrix_MMRF_count.RData") # 保存数据


# !从matrix_MMRF总提取信息-------------------------------------------------------
## !tpm
source("/Users/sandy/Rstudio/TCGA-BRCA/code/prepare.R") # 准备配置
load("matrix_MMRF_tpm.RData")
View(matrix_MMRF[[1]]) # 显示matrix_MMRF中的第一个样本的信息，确认数据是否读取正确

paste(imput_dir, filedir_in_json[1], sep = "/") %>% setwd() # 进入json文件中的第一个样本文件夹
# 读取该样本文件夹中包含基因信息的tsv文件，将其转换为数据框probematrix，选择第2列和第3列作为基因名称和基因类型, 第一列作为基因ensb_id
probematrix <- read.table(list.files(pattern = ".tsv"), header = F, fill = TRUE)[, 1:3]

matrix_MMRF <- do.call(cbind, matrix_MMRF) # 将之前读取的所有样本基因表达量数据拼接为一个矩阵
matrix_MMRF <- cbind(probematrix, matrix_MMRF) # 将当前样本的基因信息和基因表达量数据拼接为一个数据框，并加入到之前的矩阵中
tibble_MMRF <- tibble::as_tibble(matrix_MMRF) # 将矩阵转换为tibble格式的数据框
colnames(tibble_MMRF) <- c("gene_id", "gene_name", "gene_type", MMRF_ID[, 1]) # 给数据框中的列命名，第一列是基因名称，第二列是基因类型，后面的列名对应每个样本的样本ID

# 提取蛋白编码基因
pcg <- c("protein_coding", "IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_C_gene", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene")
mRNA_exprSet <- tibble_MMRF %>%
  dplyr::filter(gene_type %in% pcg) %>%
  dplyr::select(-gene_type) %>%
  tidyr::separate(gene_id, into = c("gene_id"), sep = "\\.") %>%
  tidyr::separate(gene_name, into = c("gene_name"), sep = "\\.") %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)


# 重新排序,将癌旁排在前面便于下一步筛选,0-9为癌数据,排在后面
mRNA_exprSet <- mRNA_exprSet %>%
  select(-matches(".-[0-9]{2}[^A]$")) %>%
  select(-matches(".-0[1-9][A]$"), everything())
# if have some repeat
# geneNames <- as.data.frame(mRNA_exprSet[,1])
# mRNA_exprSet[duplicated(geneNames),]

# 写出表达矩阵_tpm
shelfEnvironment(output_dir, path = root_dir)
write.table(mRNA_exprSet, "expSet_TPM_ProCoding.txt", sep = "\t", row.names = F, col.names = T, quote = F)

## !counts
source("/Users/sandy/Rstudio/TCGA-BRCA/code/prepare.R") # 准备配置
load("matrix_MMRF_count.RData")
View(matrix_MMRF[[1]]) # 显示matrix_MMRF中的第一个样本的信息，确认数据是否读取正确

paste(imput_dir, filedir_in_json[1], sep = "/") %>% setwd() # 进入json文件中的第一个样本文件夹
# 读取该样本文件夹中包含基因信息的tsv文件，将其转换为数据框probematrix，选择第2列和第3列作为基因名称和基因类型, 第一列作为基因ensb_id
probematrix <- read.table(list.files(pattern = ".tsv"), header = F, fill = TRUE)[, 1:3]

matrix_MMRF <- do.call(cbind, matrix_MMRF) # 将之前读取的所有样本基因表达量数据拼接为一个矩阵
matrix_MMRF <- cbind(probematrix, matrix_MMRF) # 将当前样本的基因信息和基因表达量数据拼接为一个数据框，并加入到之前的矩阵中
tibble_MMRF <- tibble::as_tibble(matrix_MMRF) # 将矩阵转换为tibble格式的数据框
colnames(tibble_MMRF) <- c("gene_id", "gene_name", "gene_type", MMRF_ID[, 1]) # 给数据框中的列命名，第一列是基因名称，第二列是基因类型，后面的列名对应每个样本的样本ID

# 提取蛋白编码基因
pcg <- c("protein_coding", "IG_V_gene", "IG_D_gene", "IG_J_gene", "IG_C_gene", "TR_V_gene", "TR_D_gene", "TR_J_gene", "TR_C_gene")
mRNA_exprSet <- tibble_MMRF %>%
  dplyr::filter(gene_type %in% pcg) %>%
  dplyr::select(-gene_type) %>%
  tidyr::separate(gene_id, into = c("gene_id"), sep = "\\.") %>%
  tidyr::separate(gene_name, into = c("gene_name"), sep = "\\.") %>%
  dplyr::distinct(gene_id, .keep_all = TRUE)


# 重新排序,将癌旁排在前面便于下一步筛选,0-9为癌数据,排在后面
mRNA_exprSet <- mRNA_exprSet %>%
  select(-matches(".-[0-9]{2}[^A]$")) %>%
  select(-matches(".-0[1-9][A]$"), everything())
# if have some repeat
# geneNames <- as.data.frame(mRNA_exprSet[,1])
# mRNA_exprSet[duplicated(geneNames),]

# 写出表达矩阵_tpm
shelfEnvironment(output_dir, path = root_dir)
write.table(mRNA_exprSet, "expSet_COUNTS_ProCoding.txt", sep = "\t", row.names = FALSE, col.names = T, quote = FALSE)
