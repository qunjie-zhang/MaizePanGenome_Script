# ==========================================
# Version=2.0 
# ========================================== 

# v2
# 主要解决filter_vcf 函数中 type 为 custom 时参数失效问题
# https://github.com/ZhangRenL/geneHapR/issues/2
# 2024.10.15

suppressPackageStartupMessages(library(getopt))
# 接受关键词参数
spec <- matrix(
  c(
    "geneid",  "gene", 2, "character", "显示基因ID信息",
    "chr", "ch", 2, "character",  "染色体名称",
    "start", "s", 2, "integer",  "开始位置",
    "end", "e", 2, "integer",  "结束位置",
    "happrefix", "hp", 2, "character",  "单倍型前缀",
    "vcf", "v", 2, "character",  "VCF文件路径",
    "gff", "f", 2, "character",  "gff文件路径",
    "pheno", "p", 2, "character",  "pheno表型文件路径",
    "accinfo", "ai", 2, "character",  "AccINFO关联信息文件路径",
    "filter_vcf_mode",'fvm',2,"character","filter_vcf过滤模式(POS / both)",
    "filter_vcf_type",'fvt',2,"character","filter_vcf过滤项目( CDS,exon.. 多个使用','进行分割)",
    "hetero_remove", "hr", 2, "logical",  "是否移除包含杂合位点的样本",
    "na_drop", "nd", 2, "logical",  "是否移除包含基因型缺失的样本",
    "output",'o',2,"character","输出路径"
    ),
  byrow=TRUE, ncol=5)
opt <- getopt(spec=spec)
if( 
    !is.null(opt$help) || 
    is.null(opt$geneid) || 
    is.null(opt$chr) ||
    is.null(opt$end) ||
    is.null(opt$happrefix) ||
    is.null(opt$vcf) ||
    is.null(opt$gff) ||
    is.null(opt$pheno) ||
    is.null(opt$accinfo) ||
    is.null(opt$hetero_remove) ||
    is.null(opt$na_drop)
    )
{
    # ... 这里你也可以自定义一些东放在里面
    cat(paste(getopt(spec=spec, usage = T), "\n"))
    quit()
}

# ------------------------------------------------------------------------------------------

# 统一日志输出
printLog <- function(message) {
  current_time <- sprintf("[%s]", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  res <- paste(current_time, message, sep = " - ")
  print(res)
}

#  filter_vcf Version2.0 fix issue#2
POS2GRanges <- function(Chr, POS) {
  POSRange <- GenomicRanges::GRanges(seqnames = Chr,IRanges::IRanges(start = POS,width = rep(1, length(POS))))
  return(POSRange)
}
filter_vcf <- function(vcf,
                       gff = gff,
                       mode = c("POS", "type", "both"),
                       Chr = Chr,
                       start = start,
                       end = end,
                       type = c("CDS", "exon", "gene", "genome", "custom"),
                       cusTyp = c("CDS", "five_prime_UTR", "three_prime_UTR"),
                       geneID = geneID) {
  if (mode == "POS" | mode == "both") {
    POS <- vcfR::getPOS(vcf) %>% as.numeric()
    Chrs <- vcfR::getCHROM(vcf)
    probe <- rep(TRUE, length(POS))
    if(!missing(Chr))
      probe <- probe & Chrs == Chr
    if(!missing(start))
      probe <- probe & (POS >= min(start, end))
    if(!missing(end))
      probe <- probe & (POS <= max(start, end))
    
    nr <- sum(probe)
    if(nr > 1){
      vcf@fix <- vcf@fix[probe,]
      vcf@gt <- vcf@gt[probe,]
    } else if(nr == 1){
      nf <- colnames(vcf@fix)
      ng <- colnames(vcf@gt)
      vcf@fix <- matrix(vcf@fix[probe,], byrow = TRUE, nrow = nr)
      colnames(vcf@fix) <- nf
      vcf@gt <- matrix(vcf@gt[probe,], byrow = TRUE, nrow = nr)
      colnames(vcf@gt) <- ng
    } else stop(" No variants after filter by position")
  }
  
  if (mode == "type" | mode == "both") {
    if (missing(gff))
      stop("gff is missing!")
    
    if (type == "custom")
      type <- cusTyp
    
    p <- type %in% unique(gff$type)
    m <- paste(unique(gff$type), collapse = "','")
    if (FALSE %in% p)
      stop("type should in c('",m,"')")
    
    if ("genome" %in% type) {
      gff <- gff
    } else {
      if(!missing(geneID)){
        ids <- tolower(gff$ID)
        nms <- tolower(gff$Name)
        id <- tolower(geneID)
        p <- (stringr::str_detect(ids,id)) | (stringr::str_detect(nms,id))
      }
      gff <- gff[(gff$type %in% type) & p]
    }
    
    POS <- vcfR::getPOS(vcf)
    POS <- as.numeric(POS)
    POSRange <- POS2GRanges(Chr = vcfR::getCHROM(vcf),
                            POS = POS)
    POSRange_rm <- POSRange[!(POSRange %over% gff)]
    
    POS_rm <- IRanges::start(POSRange_rm)
    probe <- !(POS %in% POS_rm)
    nr <- sum(probe)
    if(nr > 1){
      vcf@fix <- vcf@fix[probe,]
      vcf@gt <- vcf@gt[probe,]
    } else if(nr == 1){
      nf <- colnames(vcf@fix)
      ng <- colnames(vcf@gt)
      vcf@fix <- matrix(vcf@fix[probe,], byrow = TRUE, nrow = nr)
      colnames(vcf@fix) <- nf
      vcf@gt <- matrix(vcf@gt[probe,], byrow = TRUE, nrow = nr)
      colnames(vcf@gt) <- ng
    } else stop(" No variants after filter by ", type)
  }
  
  return(vcf)
}


# ------------------------------------------------------------------------------------------
printLog('The application is starting')
printLog('Start loading the program runtime environment')
# 屏蔽警告信息
options(warn = -1)
# 加载所需包
suppressPackageStartupMessages(library(geneHapR))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(muscle))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(trackViewer))
suppressPackageStartupMessages(library(vcfR))

printLog('Program runtime environment loaded successfully')


# 定义基本参数变量
geneID <- opt$geneid        # 基因ID
Chr <- opt$chr              # 基因所处的染色体名称
start <- opt$start          # 基因的起始位置（染色体坐标）
end <- opt$end              # 基因的终止位置（染色体坐标）
hapPrefix <- opt$happrefix  # 单倍型名称的前缀   

output <-  opt$output  # 输出路径

filter_vcf_mode <- opt$filter_vcf_mode
filter_vcf_type <- unlist(strsplit(opt$filter_vcf_type, ","))

# 用户输入字符串转化为布尔值
hetero_remove <- as.logical(opt$hetero_remove)
na_drop <- as.logical(opt$na_drop)

printLog('-----------------------------------------')
printLog('User-defined parameters:')
print(paste("geneID:", geneID))
print(paste("Chr:", Chr))
print(paste("start:", start))
print(paste("end:", end))
print(paste("hapPrefix:", hapPrefix))
print(paste("vcf:", opt$vcf))
print(paste("gff:", opt$gff))
print(paste("pheno:", opt$pheno))
print(paste("accinfo:", opt$accinfo))
print(paste("filter_vcf_mode:", filter_vcf_mode))
print(paste("filter_vcf_type:", filter_vcf_type))
print(paste("hetero_remove:", hetero_remove))
print(paste("na_drop:", na_drop))
print(paste("output:", output))
printLog('-----------------------------------------')

printLog('Start importing the data required for this analysis')
# 设置工作路径

setwd(output)
printLog(paste("Switched to output directory:", getwd()))
# setwd("D:/OneDrive - stu.scau.edu.cn/Project/haplotype")

# GFF3 单个基因
printLog('Reading gff file')
gff <- import_gff(opt$gff,format='gff3')
# gff <- read.delim(opt$gff, header=FALSE, comment.char="#")

# 表型文件 全部
printLog('Reading pheno file')
pheno <- read.csv(opt$pheno ,header=T, row.names = 1)

# 分类数据 全部
printLog('Reading AccINFO file')
AccINFO <- import_AccINFO(opt$accinfo, 
                          sep = ",",                      # 分隔符号，默认为制表符"\t"
                          na.strings = "NA")              # 导入其他样本信息


printLog('Loading the required VCF file data into the workspace.')
# 导入指定VCF文件

vcf <- import_vcf(opt$vcf)


# 相关源代码
# https://github.com/ZhangRenL/geneHapR/blob/1390bf269effd745509cfd3b70e73c5d97dd6d03/R/filter.R#L45
# 本次这里使用重写方法。
# mode使用POS或Both。POS为仅通过位置过滤数据，Both则加入gff文件共同过滤。用户仅输入起始结束位置时采用POS，输入基因ID时采用both
# type统一使用custom。type为custom时，cusType值可以为多个。
# 使用gff 时使用both,即指定范围内指定类型的变异信息
printLog('Start filtering VCF data')

# filterLargeVCF(VCFin = vcf,
#                VCFout = '/tmp/vcf.vcf',
#                Chr = Chr,
#                POS = c(start,end),
#                override = TRUE)


vcf <- filter_vcf(
                  vcf, gff,
                  mode = filter_vcf_mode,
                  start = start,
                  end = end, 
                  Chr = Chr,
                  type = 'custom',
                  #cusTyp = c("CDS",'exon','gene','five_prime_UTR','three_prime_UTR'),
                  cusTyp = filter_vcf_type
                  )
printLog('VCF data filtering complete')

# 切换到输出路径
if (!dir.exists(output)) {dir.create(output, recursive = TRUE)}
setwd(output)

# 从VCF开始单倍型鉴定
printLog('Haplotype identification in progress')
hapResult <- vcf2hap(vcf, hapPrefix = hapPrefix,
                     hetero_remove = hetero_remove, # 移除包含杂合位点的样本
                     na_drop = na_drop ) # 移除包含基因型缺失的样本
write.hap(hapResult, file = "hapResult.txt")


# 对单倍型结果进行汇总整理
printLog('Data summarization of haplotype identification results is in progress')
hapSummary <- hap_summary(hapResult, hapPrefix = hapPrefix)
write.hap(hapSummary, file = "hapSummary.txt")


# # 将单倍型分析结果输出保存
printLog('The results of the identification of relevant haplotype data are being saved')

# # # 导入之前的单倍型分析结果
# hapResult <- import_hap(file = "hapResult")
# hapSummary <- import_hap(file = "hapSummary")

plot_width <- 16
plot_height <- 9


# 2 单倍型结果可视化
tryCatch({
    printLog('Haplotype visualization output in progress')
    plot <- plotHapTable(hapSummary,
                     hapPrefix = hapPrefix,
                     angle = 45,
                     displayIndelSize = 0,
                     title = geneID,
                     INFO_tag = "ANN", tag_field = 11, geneName = geneID) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

  pdf(
  'plotHapTable.pdf', 
    width= plot_width,
    height = plot_height, 
  )
  print(plot)
  while (!is.null(dev.list()))  dev.off()
}, error = function(e) {
  printLog('An error occurred while generating the haplotype visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
})


# 3 变异位点信息可视化
# 这个功能调用了lolliplot，只能先打开画布后绘图，无法将图像数据保存到变量
tryCatch({
printLog('Variant On GeneModel visualization output in progress')
pdf(
  'displayVarOnGeneModel.pdf',
  width= plot_width,
  height = plot_height, 
  )
displayVarOnGeneModel(gff = gff, hapSummary = hapSummary,
                      startPOS = start-10,
                      endPOS = end+10,
                      cex = 0.8) # size of variants
while (!is.null(dev.list()))  dev.off()
} , error = function(e) {
  printLog('An error occurred while generating the variant visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
})

# 4 单倍型网络分析可视化
# 这个功能调用了lolliplot，只能先打开画布后绘图，无法将图像数据保存到变量

tryCatch({
printLog('Haplotype network visualization output in progress')
pdf(
  'plotHapNet.pdf',
  width= plot_width,
  height = plot_height, 
  )
hapSummary[hapSummary == "DEL"] = "N"
hapnet <- get_hapNet(hapSummary,                  # 单倍型结果
                     AccINFO = AccINFO,           # 包含样本分类信息的数据框(data.frame)
                     groupName = "group", # 含有样本分类信息的列名称
                     na.label = "Unknown")        # 未知分类样本的类别
plotHapNet(hapnet,                          # 单倍型网络
           scale = "log2",                  # 标准化方法"log10"或"log2"或"none"
           show.mutation = 2,               # 是否展示变异位点数量, 0,1,2,3
           col.link = 2, link.width = 2,    # 单倍型之间连线的颜色和宽度
           main = geneID,                   # 主标题
           pie.lim = c(0.5, 2),             # 圆圈的大小
           legend_version = 0,              # 图例形式（0或1）
           labels = T,                      # 是否在单倍型上添加label
           # legend = FALSE)                # 不添加图例
           # legend = TRUE)                 # 添加图例,但需要单击添加的位置
           legend = c(8,1),                # 图例的坐标
           cex.legend = 0.6)                # 图例中文字的大小
while (!is.null(dev.list()))  dev.off()
} , error = function(e) {
  printLog('An error occurred while generating the haplotype network visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
})


# 6 连锁不平衡分析
tryCatch({
printLog('Linkage disequilibrium(LD) visualization output in progress')
pdf(
  'plot_LDheatmap.pdf', 
  width= plot_width,
  height = plot_height, 
  )
plot_LDheatmap(hap = hapResult, # 单倍型结果
               add.map = TRUE,  # 是否添加基因模式图
               gff = gff,       # 注释信息
               Chr = Chr,       # 染色体名称
               start = start,   # 基因的起始位置
               end = end)       # 基因的终止位置（更多参数参见帮助文档）
while (!is.null(dev.list()))  dev.off()
} , error = function(e) {
  printLog('An error occurred while generating the LD visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
})

# 7 表型关联分析
tryCatch({
printLog('Phenotypic correlation analysis visualization output in progress')
pdf(
  'hapVsPheno.pdf',
  width= plot_width,
  height = plot_height, 
  units='px', 
  res=plot_ppi
  )
hapVsPheno(hap = hapResult,       # 单倍型分析结果
           pheno = pheno,         # 表型
           # phenoName = '',      
           hapPrefix = hapPrefix, # 单倍型名称的前缀
           title = geneID,        # 主标题
           minAcc = 4,            # 参与p值计算所需的最小样本数
           symnum.args = list(    # 定义显著性标注方式
             cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
             symbols = c("***", "**", "*", "ns")),
           mergeFigs = TRUE)     # 结果包括两个图，是否融合成一张图
while (!is.null(dev.list()))  dev.off()
}, error = function(e) {
  printLog('An error occurred while generating the phenotypic correlation analysis visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
}) 

# 8 位点效应
tryCatch({
printLog('Variant effect estimation visualization output in progress')
pdf(
  'plotEFF.pdf',
  width= plot_width,
  height = plot_height, 
  units='px', 
  res=plot_ppi
  )
EFF <- siteEFF(hapResult, pheno)
plotEFF(EFF, gff = gff,
        Chr = Chr, start = start, end = end,
        showType = c("five_prime_UTR", "CDS", "three_prime_UTR"), # see help(plotEFF)
        y = "effect",                      # the means of y axis, one of effect or pvalue
        ylab = "effect",                  # label of y axis
        cex = 0.5,                         # Cex
        legend.cex = 0.8,                  # legend size
        main = geneID,                     # main title
        CDS.height = 1,                    # controls the height of CDS, heights of others will be half of that
        markMutants = TRUE,                # mark mutants by short lines
        mutants.col = 1, mutants.type = 1, # parameters for appearance of mutants
        pch = 20)                          # points type

while (!is.null(dev.list()))  dev.off()
} , error = function(e) {
  printLog('An error occurred while generating the variant effect estimation visualization output:')
  printLog(e)
  while (!is.null(dev.list()))  dev.off()
  NULL
})

printLog('Done.')