create_design_matrix <- function(groups) {
  # 替换空格和加号，并转为因子
  group <- as.factor(gsub(" |\\+", "_", groups))
  
  # 创建不含截距的设计矩阵
  design <- model.matrix(~ 0 + group)
  
  return(design)
}