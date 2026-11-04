# Validation record

- 16 个整理后的 R 分析/示例脚本和 1 个 R 检查工具通过 Tree-sitter R 语法解析，未发现 ERROR 或 missing 节点。
- Python velocity 脚本通过 Python 标准库 `ast.parse` 语法检查。
- 五份原始正文按来源行范围完整分配到 17 个脚本；来源文本与拆分映射一并保留。
- README 的本地链接和 ZIP 内文件完整性已检查。

本机没有可用的 Rscript，因此未执行 R 自带的 `parse()`。Tree-sitter 检查不能替代 R 解释器和包 API 验证。在配置好的 R 环境中，可从仓库根目录运行：

```bash
Rscript tools/check_syntax.R
```

没有提供表达数据，也没有运行 Seurat、muscat、WGCNA、GENIE3 或 scVelo 分析。语法通过不表示依赖、数据格式、对象状态、统计设计或生物学结论已验证。
