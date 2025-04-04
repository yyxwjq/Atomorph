#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
演示如何使用Atomorph的固定元素功能。
"""

import os
import sys

# 添加上级目录到路径，以便导入atomorph
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from atomorph import StructureConverter

def main():
    """演示固定元素功能"""
    print("演示Atomorph的固定元素功能")
    
    # 创建转换器
    converter = StructureConverter()
    
    # 设置输入和输出文件路径
    input_path = "../tests/data/IGRSI-ts.extxyz"  # 使用测试数据文件
    output_path = "fixed_elements_output.vasp"
    
    # 定义要固定的元素
    fixed_elements = ["elements", "Ru"]  # 固定所有Ru原子
    
    print(f"输入文件: {input_path}")
    print(f"输出文件: {output_path}")
    print(f"固定元素: {fixed_elements[1]}")
    
    # 进行转换，固定指定元素
    converter.convert(
        input_path=input_path,
        output_path=output_path,
        constraints=fixed_elements
    )
    
    print(f"转换完成！结果已保存到 {output_path}")
    print("提示: 使用VESTA或其他分子可视化软件查看结果")

if __name__ == "__main__":
    main() 