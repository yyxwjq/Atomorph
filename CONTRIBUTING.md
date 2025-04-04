# Contributing to Atomorph / Atomorph贡献指南

Thank you for your interest in contributing to Atomorph! This document provides guidelines and instructions for contributing.

感谢您有兴趣为Atomorph做出贡献！本文档提供了贡献的指南和说明。

## Code of Conduct / 行为准则

Please be respectful and considerate of others when contributing to this project.

请在为本项目做贡献时尊重并考虑他人。

## Getting Started / 入门指南

1. **Fork the repository / 复刻仓库**

2. **Clone your fork / 克隆您的复刻**
   ```bash
   git clone https://github.com/YOUR_USERNAME/atomorph.git
   cd atomorph
   ```

3. **Set up development environment / 设置开发环境**
   ```bash
   pip install -e ".[dev]"
   ```

4. **Create a branch / 创建分支**
   ```bash
   git checkout -b feature/your-feature-name
   ```

## Development Guidelines / 开发指南

### Code Style / 代码风格

- Follow PEP 8 style guide / 遵循PEP 8风格指南
- Use docstrings for all functions, classes, and modules / 为所有函数、类和模块使用文档字符串
- Write clear, concise comments in English / 用英文编写清晰、简洁的注释

### Testing / 测试

- Add unit tests for new features / 为新功能添加单元测试
- Run tests before submitting / 提交前运行测试
  ```bash
  pytest
  ```

### Commit Messages / 提交消息

- Use clear and descriptive commit messages / 使用清晰、描述性的提交消息
- Format: `[type]: Brief description` / 格式：`[类型]: 简短描述`
- Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore` / 类型：`feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

## Pull Request Process / 拉取请求流程

1. Update the README.md or documentation with details of changes / 更新README.md或文档，详细说明更改
2. Update the CHANGELOG.md file / 更新CHANGELOG.md文件
3. The PR should work for Python 3.7 and above / PR应适用于Python 3.7及以上版本
4. Wait for review and address any comments / 等待审查并处理任何评论

## Release Process / 发布流程

1. Update version in `__init__.py` (following semantic versioning, starting from 1.0.0) / 在`__init__.py`中更新版本（遵循语义化版本，从1.0.0开始）
2. Update CHANGELOG.md with the new version details / 用新版本详情更新CHANGELOG.md
3. Create a new GitHub release with appropriate tag / 用适当的标签创建新的GitHub发布
4. Publish to PyPI / 发布到PyPI

Thank you for contributing! / 感谢您的贡献！ 