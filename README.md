# Atomorph

[English](#english) | [中文](#chinese)

<a name="english"></a>
## Atomorph

Atomorph is a tool for converting and manipulating atomic structure files. It supports conversions between various file formats, as well as handling atomic constraints.

**Current Version: 1.0.3**

### Features

- Support for multiple file formats (VASP, XYZ, EXTXYZ, etc.)
- Multi-frame structure file processing
- Atomic constraints (fixed layers, specific elements, specific indices)
- Atom sorting
- Parallel processing
- Batch conversion of multiple files
- Directory traversal with pattern matching
- Both command-line and Python API interfaces
- Simplified command syntax with aliases
- Detailed logging and error handling
- Comprehensive documentation with bilingual support

### Installation

```bash
pip install atomorph
```

### Usage

#### Command-line Usage

Basic conversion (with short command):
```bash
atomorph c input.xyz output.vasp
```

Specify formats:
```bash
atomorph convert input.xyz output.vasp -if xyz -of vasp
```

Fix specific elements (with short parameters):
```bash
atomorph c input.xyz output.vasp -c "elements:C,H"
```

Sort atoms ascending (with short parameters):
```bash
atomorph c input.xyz output.vasp -s asc
```

Custom element order (with short parameters):
```bash
atomorph c input.xyz output.vasp -e "Si,O,H"
```

#### Batch Conversion

Convert all POSCAR files in a directory:
```bash
atomorph b input_dir/ output_dir/ -p POSCAR -of xyz
```

Convert with recursive search:
```bash
atomorph batch input_dir/ output_dir/ -p "*.vasp" -of xyz
```

Merge all structures into one file:
```bash
atomorph b input_dir/ output_dir/ -p POSCAR -of xyz -m
```

#### Python API Usage

```python
from atomorph import StructureConverter, batch_convert

# Basic conversion
converter = StructureConverter()
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp"
)

# Batch conversion
batch_convert(
    input_dir="input_dir/",
    output_dir="output_dir/",
    pattern="POSCAR",
    output_format="xyz"
)

# Merge all structures into one file
batch_convert(
    input_dir="input_dir/",
    output_dir="output_dir/",
    pattern="*.vasp",
    output_format="xyz",
    merge_output=True
)
```

### Version History

See [CHANGELOG.md](CHANGELOG.md) for a complete version history.

### License

MIT

---

<a name="chinese"></a>
## Atomorph 中文说明

Atomorph是一个用于转换和操作原子结构文件的工具，支持多种文件格式间的转换，以及原子约束的处理。

**当前版本: 1.0.3**

### 特性

- 支持多种文件格式转换(VASP, XYZ, EXTXYZ等)
- 支持多帧结构文件处理
- 支持原子约束(固定层、指定元素、指定索引)
- 支持原子排序
- 支持并行处理
- 支持批量转换多个文件
- 支持使用模式匹配遍历目录
- 命令行和Python API两种使用方式
- 简化的命令语法和别名
- 详细的日志和错误处理
- 全面的双语文档支持

### 安装

```bash
pip install atomorph
```

### 使用方法

#### 命令行使用

基本转换（使用简短命令）：
```bash
atomorph c input.xyz output.vasp
```

指定格式：
```bash
atomorph convert input.xyz output.vasp -if xyz -of vasp
```

固定特定元素（使用简短参数）：
```bash
atomorph c input.xyz output.vasp -c "elements:C,H"
```

按升序排序原子（使用简短参数）：
```bash
atomorph c input.xyz output.vasp -s asc
```

自定义元素顺序（使用简短参数）：
```bash
atomorph c input.xyz output.vasp -e "Si,O,H"
```

#### 批量转换

转换目录中的所有POSCAR文件：
```bash
atomorph b input_dir/ output_dir/ -p POSCAR -of xyz
```

递归搜索并转换文件：
```bash
atomorph batch input_dir/ output_dir/ -p "*.vasp" -of xyz
```

将所有结构合并为一个文件：
```bash
atomorph b input_dir/ output_dir/ -p POSCAR -of xyz -m
```

#### Python API使用

```python
from atomorph import StructureConverter, batch_convert

# 基本转换
converter = StructureConverter()
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp"
)

# 批量转换
batch_convert(
    input_dir="input_dir/",
    output_dir="output_dir/",
    pattern="POSCAR",
    output_format="xyz"
)

# 合并所有结构为一个文件
batch_convert(
    input_dir="input_dir/",
    output_dir="output_dir/",
    pattern="*.vasp",
    output_format="xyz",
    merge_output=True
)
```

### 版本历史

完整的版本历史请参见[CHANGELOG.md](CHANGELOG.md)。

### 许可证

MIT 