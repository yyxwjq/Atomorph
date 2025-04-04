#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
演示如何使用Atomorph处理多帧结构文件。
"""

import os
import sys

# 添加上级目录到路径，以便导入atomorph
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from atomorph import StructureConverter

def main():
    """演示多帧处理功能"""
    print("演示Atomorph的多帧处理功能")
    
    # 创建转换器
    converter = StructureConverter()
    
    # 设置输入和输出文件路径
    input_path = "../tests/data/IGRSI-ts.extxyz"  # 使用测试数据文件
    output_path = "multi_frame_output"
    
    print(f"输入文件: {input_path}")
    print(f"输出目录: {output_path}")
    
    # 进行多帧转换，使用并行处理
    converter.convert(
        input_path=input_path,
        output_path=output_path,
        multi_frame=True,    # 处理多帧
        parallel=True,       # 使用并行处理
        separate_dirs=True   # 将每帧保存在单独的目录中
    )
    
    print(f"转换完成！结果已保存到 {output_path} 目录")
    print("提示: 检查输出目录中的 frame_X 子目录，每个都包含一个帧的VASP文件")

if __name__ == "__main__":
    main() 