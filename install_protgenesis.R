# ==============================================================================
# install_ProtGenesis.R - 安装 ProtGenesis R 包
# ==============================================================================

cat("==============================================================================\n")
cat("Installing ProtGenesis Package\n")
cat("==============================================================================\n\n")

# 设置工作目录
pkg_dir <- "c:/Users/刘传扬/Desktop/Exploration/1.how_protein_space_lookslike/ProtGenesis"

# 先尝试更新关键依赖包
cat("Step 1: Updating key packages...\n")
update.packages(ask = FALSE, repos = "https://cloud.r-project.org/")

# 加载 devtools
cat("\nStep 2: Loading devtools...\n")
library(devtools)

# 安装 ProtGenesis 包
cat("\nStep 3: Installing ProtGenesis package...\n")
install(pkg_dir, quiet = FALSE, dependencies = TRUE)

# 测试加载
cat("\nStep 4: Testing package loading...\n")
library(ProtGenesis)

cat("\n==============================================================================\n")
cat("Package installation completed!\n")
cat("==============================================================================\n\n")

# 显示所有导出的函数
cat("Available functions in ProtGenesis:\n")
ls("package:ProtGenesis")
