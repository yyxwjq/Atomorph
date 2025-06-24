# Changelog / 变更日志

All notable changes to this project will be documented in this file.

此文件记录项目的所有重要更改。

## [1.0.5] - 2024-06-XX

### Fixed / 修复
- Fixed documentation content, including README.md / 修正文档内容，包括README.md

## [1.0.4] - 2024-06-XX

### Improved / 改进
- Optimized several functions related to atom sorting / 优化了原子排序相关的几个函数

## [1.0.3] - 2024-04-06

### Added / 新增
- Batch conversion functionality to process multiple files in directories / 批量转换功能，可处理目录中的多个文件
- Shortened command aliases for better usability (e.g., 'c' for 'convert', 'b' for 'batch') / 更短的命令别名以提高可用性（如'c'代表'convert'，'b'代表'batch'）
- Support for merging multiple structures into one multi-frame file / 支持将多个结构合并为一个多帧文件
- Shorthand parameters for commonly used options (e.g., '-s' for '--sort', '-e' for '--element-order') / 常用选项的简写参数（如'-s'代表'--sort'，'-e'代表'--element-order'）

### Improved / 改进
- Command-line interface with shorter aliases and better help documentation / 命令行界面使用更短的别名和更好的帮助文档
- Error handling and debugging for batch operations / 批处理操作的错误处理和调试
- Progress reporting for batch conversion operations / 批量转换操作的进度报告
- File pattern matching for finding structure files in directories / 用于在目录中查找结构文件的文件模式匹配

### Fixed / 修复
- Sort order parameter handling for custom element ordering / 自定义元素排序的排序参数处理
- Directory structure preservation when batch converting files / 批量转换文件时保留目录结构
- Command-line parameter parsing for abbreviated options / 缩写选项的命令行参数解析

## [1.0.2] - 2024-04-05

### Changed / 变更
- Version number update for PyPI compatibility / 更新版本号以兼容PyPI
- Maintained all functionality from 1.0.0 / 保持与1.0.0版本相同的功能

## [1.0.0] - 2024-04-05

### Added / 新增
- First official release with complete functionality / 首个官方发布版本，具有完整功能
- Comprehensive Python API for structural conversions / 全面的Python API用于结构转换
- Enhanced documentation with bilingual support / 增强的双语文档支持
- Complete test suite with unit and integration tests / 完整的测试套件，包含单元和集成测试
- Utility functions for file handling and logging / 用于文件处理和日志记录的实用功能

### Improved / 改进
- Optimized code structure with better organization / 优化的代码结构，组织更合理
- Enhanced error handling and debugging capabilities / 增强的错误处理和调试功能
- Better support for atomic constraints (fixed elements, layers, indices) / 更好地支持原子约束（固定元素、层和索引）
- Improved command-line interface with more options / 改进的命令行界面，提供更多选项

### Changed / 变更
- Standardized project structure following best practices / 标准化的项目结构，遵循最佳实践
- English debug output for better international support / 英文调试输出，以提供更好的国际支持