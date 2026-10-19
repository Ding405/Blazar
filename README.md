# Blazar

用于耀变体光学图像中中心活动核与寄主星系二维分解的研究项目。
当前提供的是**模拟数据验证脚本**，尚不能直接分析真实 HST 图像，亦未得到 S5 0716+714 的实测结果。

## 安装

已在 Python 3.12 中验证。建议使用独立虚拟环境：

```bash
python -m venv .venv
```

激活环境后安装依赖（Linux/macOS：`source .venv/bin/activate`；Windows：`.venv\Scripts\activate`）：

```bash
python -m pip install -r requirements.txt
python validate.py
```

## Jupyter Notebook

在 Notebook 单元格中安装依赖：

```python
%pip install -r requirements.txt
```

然后运行：

```python
%run validate.py
```

代码在 Notebook 内请求显示保存的对比图；命令行运行时使用无界面的 Agg 后端。
也可以将脚本内容粘贴到单元格运行，没有 `__file__` 时会使用当前工作目录。

## 输出

每次运行创建独立的 `outputs/validation-*/` 目录，不覆盖旧结果，也不在退出时自动删除：

- `synthetic.fits`：模拟观测图像，带 FITS 校验和。
- `comparison.png`：模拟观测、最佳拟合模型、残差三联图。
- `residual.png`：残差图。
- `summary.json`：拟合参数、模拟真值、残差统计及使用范围。

脚本打印输出目录。上述文件生成在运行代码的机器上；云端路径不是公开下载链接。

## 验证与科学限制

脚本检查四个参数（包括背景）的恢复误差、残差均值、约化 χ²、FITS 无损读写和 PNG 完整性。
验证阈值仅适用于当前固定随机种子的模拟案例。

模拟与拟合使用相同模型：固定高斯 PSF、中心、Sérsic 指数 n=2、椭率及位置角，噪声为独立高斯噪声。
核与寄主通量均按有限图像范围定义，寄主通量不是延伸到无限远的 Sérsic 总通量。
尚未评估 PSF 失配、相关噪声、参数退化和科学不确定性。
真实 HST 拟合还需要科学图像、适配 PSF、误差模型、掩膜、光度标定及注入恢复实验。
