#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
测试Atomorph的固定层功能。
"""

import os
import sys
import argparse
from pathlib import Path

# 添加上级目录到路径，以便导入atomorph
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from atomorph import StructureConverter
from ase.io import read

def test_fixed_layers(input_path, output_path, layer_range, use_fractional=True):
    """
    测试固定层功能。
    
    参数:
        input_path: 输入文件路径
        output_path: 输出文件路径
        layer_range: 层范围，格式为"start,end"
        use_fractional: 是否使用分数坐标
    """
    print(f"测试固定层功能...")
    print(f"输入文件: {input_path}")
    print(f"输出文件: {output_path}")
    print(f"层范围: {layer_range}")
    print(f"使用分数坐标: {use_fractional}")
    
    # 创建转换器
    converter = StructureConverter()
    
    # 构建约束
    constraints = ["layers", layer_range]
    
    # 读取结构文件
    structure = read(input_path)
    print(f"结构包含 {len(structure)} 个原子")
    
    # 显示分数坐标范围（如果使用分数坐标）
    if use_fractional:
        cell_inv = structure.get_reciprocal_cell().T
        positions = structure.positions
        frac_positions = positions @ cell_inv
        frac_z_values = [pos[2] for pos in frac_positions]
        print(f"分数Z坐标范围: {min(frac_z_values):.5f} - {max(frac_z_values):.5f}")
        
        # 解析层范围
        start, end = map(float, layer_range.split(','))
        
        # 查找在该范围内的原子
        atoms_in_range = []
        for i, frac_pos in enumerate(frac_positions):
            if start <= frac_pos[2] <= end:
                atoms_in_range.append(i)
        
        print(f"层范围 {start}-{end} 内的原子数: {len(atoms_in_range)}")
        if atoms_in_range:
            print(f"前几个原子: {atoms_in_range[:5]} ...")
    
    # 执行转换
    converter.convert(
        input_path=input_path,
        output_path=output_path,
        constraints=constraints,
        use_fractional=use_fractional
    )
    print(f"转换完成，输出文件: {output_path}")
    
    # 验证固定原子
    output_path = Path(output_path)
    
    if output_path.is_dir():
        # 处理输出为目录的情况
        # 查找第一个frame文件
        frame_files = list(output_path.glob("frame_*.vasp"))
        if not frame_files:
            frame_files = list(output_path.glob("frame_*/POSCAR"))
            
        if not frame_files:
            print(f"警告: 在输出目录 {output_path} 中未找到VASP文件")
            return
            
        # 使用第一个找到的文件
        check_path = frame_files[0]
        print(f"检查文件: {check_path}")
    else:
        # 输出是单个文件的情况
        check_path = output_path
    
    # 验证固定原子
    if str(check_path).endswith(".vasp") or str(check_path).endswith("POSCAR"):
        # 读取VASP文件
        with open(check_path, 'r') as f:
            lines = f.readlines()
            
        # 检查是否包含selective dynamics
        if "Selective dynamics" in ''.join(lines[:10]):
            # 找到坐标开始的行
            for i, line in enumerate(lines):
                if "Selective" in line:
                    coord_line = i + 2  # 坐标通常在Selective dynamics下两行
                    break
            else:
                coord_line = 7  # 默认情况
            
            # 计算固定的原子
            fixed_count = 0
            total_atoms = 0
            
            for line in lines[coord_line:]:
                if not line.strip():
                    continue
                    
                parts = line.split()
                if len(parts) >= 6:  # 有足够的列来包含约束标记
                    total_atoms += 1
                    if 'F' in parts[3:6]:
                        fixed_count += 1
            
            print(f"VASP文件中的总原子数: {total_atoms}")
            print(f"固定的原子数: {fixed_count}")
            print(f"固定比例: {fixed_count / total_atoms:.2%}")
        else:
            print("警告: VASP文件中没有发现Selective dynamics标记")
    else:
        print(f"警告: 输出文件格式不是VASP，无法验证固定原子")

def main():
    parser = argparse.ArgumentParser(description="测试Atomorph的固定层功能")
    parser.add_argument("input", help="输入文件路径")
    parser.add_argument("output", help="输出文件路径")
    parser.add_argument("--layer-range", default="0.31,0.4", help="层范围，格式为'start,end'")
    parser.add_argument("--use-fractional", action="store_true", help="使用分数坐标")
    
    args = parser.parse_args()
    
    test_fixed_layers(args.input, args.output, args.layer_range, args.use_fractional)

if __name__ == "__main__":
    main() 