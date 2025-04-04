# Atomorph

[English](#english) | [中文](#chinese)

<a name="english"></a>
## Atomorph

Atomorph is a tool for converting and manipulating atomic structure files. It supports conversions between various file formats, as well as handling atomic constraints.

**Current Version: 1.0.0**

### Features

- Support for multiple file formats (VASP, XYZ, EXTXYZ, etc.)
- Multi-frame structure file processing
- Atomic constraints (fixed layers, specific elements, specific indices)
- Atom sorting
- Parallel processing
- Both command-line and Python API interfaces
- Detailed logging and error handling
- Comprehensive documentation with bilingual support

### Installation

```bash
pip install atomorph
```

### Usage

#### Command-line Usage

Basic conversion:
```bash
atomorph input.xyz output.vasp
```

Specify formats:
```bash
atomorph input.xyz output.vasp --input-format xyz --output-format vasp
```

Fix specific elements:
```bash
atomorph input.xyz output.vasp --constraints "elements:C,H"
```

Fix specific layers:
```bash
atomorph input.xyz output.vasp --constraints "layers:0.3,0.5" --use-fractional
```

Fix specific atoms:
```bash
atomorph input.xyz output.vasp --constraints "indices:1,2,5-10"
```

#### Python API Usage

```python
from atomorph import StructureConverter

converter = StructureConverter()

# Basic conversion
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp"
)

# Fix specific layers
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp",
    constraints=["layers", "0.3,0.5"],
    use_fractional=True
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

**当前版本: 1.0.0**

### 特性

- 支持多种文件格式转换(VASP, XYZ, EXTXYZ等)
- 支持多帧结构文件处理
- 支持原子约束(固定层、指定元素、指定索引)
- 支持原子排序
- 支持并行处理
- 命令行和Python API两种使用方式
- 详细的日志和错误处理
- 全面的双语文档支持

### 安装

```bash
pip install atomorph
```

### 使用方法

#### 命令行使用

基本转换：
```bash
atomorph input.xyz output.vasp
```

指定格式：
```bash
atomorph input.xyz output.vasp --input-format xyz --output-format vasp
```

固定特定元素：
```bash
atomorph input.xyz output.vasp --constraints "elements:C,H"
```

固定特定层：
```bash
atomorph input.xyz output.vasp --constraints "layers:0.3,0.5" --use-fractional
```

固定特定原子：
```bash
atomorph input.xyz output.vasp --constraints "indices:1,2,5-10"
```

#### Python API使用

```python
from atomorph import StructureConverter

converter = StructureConverter()

# 基本转换
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp"
)

# 固定特定层
converter.convert(
    input_path="input.xyz",
    output_path="output.vasp",
    constraints=["layers", "0.3,0.5"],
    use_fractional=True
)
```

### 版本历史

完整的版本历史请参见[CHANGELOG.md](CHANGELOG.md)。

### 许可证

MIT 