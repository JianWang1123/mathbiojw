#!/bin/bash
# 在本地获取三物种模型并生成图片的完整步骤

echo "================================"
echo "获取三物种模型并生成图片"
echo "================================"
echo ""

# 步骤1：拉取最新代码
echo "步骤1：拉取最新代码..."
git fetch origin
git checkout claude/three-strategist-model-hs4rK
git pull origin claude/three-strategist-model-hs4rK

# 步骤2：进入模型目录
echo ""
echo "步骤2：进入模型目录..."
cd models/three_species_crossfeeding

# 步骤3：安装依赖
echo ""
echo "步骤3：安装Python依赖..."
pip install numpy scipy matplotlib seaborn

# 步骤4：运行脚本生成图片
echo ""
echo "步骤4：生成图片..."
python3 quick_start.py

# 完成
echo ""
echo "================================"
echo "完成！"
echo "================================"
echo ""
echo "图片位置: models/three_species_crossfeeding/figures/"
echo ""
echo "生成的图片："
ls -lh figures/*.png

echo ""
echo "打开图片文件夹："
echo "cd figures"
echo ""
