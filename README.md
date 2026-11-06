# saliency-cpp

[English](#english) | 简体中文

Itti–Koch 视觉显著性模型的纯 C++17 重实现，移植自 Dirk B. Walther（Caltech，2006–2013）的 MATLAB **SaliencyToolbox**。逐行复刻原始工具箱的计算路径，无 MATLAB 依赖，单文件可执行程序即可完成显著图计算、注意焦点预测与返回抑制（IOR）。

## 特性

- 完整的自下而上显著性流水线：高斯金字塔 → 强度/颜色（RG、BY、色调）/方向（Gabor）特征 → 中央-周边差 → 跨尺度合成显著图。
- LIF（leaky integrate-and-fire）神经元 WTA 网络逐个产生注意焦点，支持三种返回抑制模式（`disk` / `shape` / `none`）。
- 形状估计（proto-object segmentation），可导出每个注视点的形状 mask。
- **Top-down 调制**：可混入外部注意图（`--topdown`），或基于 NCC 的模板匹配热力图（`--template`）。
- 批处理整个目录，输出显著图 PNG、叠加可视化 PNG、注视点坐标 CSV。
- 可复现：`--norandom` + `--seed` 固定全部随机性；`--selftest` 内置数值自检。
- 零第三方依赖：图像 IO 使用单头文件 [stb_image](https://github.com/nothings/stb)，Windows 下 MSVC 一条命令即可构建。

## 构建

需要 Windows + MSVC（Visual Studio Build Tools）：

```bat
cd src
build.bat
```

`build.bat` 内部调用 `vcvarsall.bat x64` 与 `cl /O2 /std:c++17 /EHsc /W3 /utf-8`，产出 `saliency.exe`。

## 使用方法

```
saliency.exe <图片或目录> [选项]
  --fixations N    注视点个数(默认 5,0 = 只算显著图)
  --out DIR        输出目录(默认 ./out)
  --ior MODE       disk|shape|none(默认 shape)
  --norandom       关闭坐标 jitter 与噪声(完全可复现)
  --seed N         随机种子(默认 42)
  --topdown IMG    外部 top-down 注意图(加权混入显著图)
  --template PATH  模板文件或目录,NCC 匹配热力图作为 top-down(与 --topdown 互斥)
  --tdweight W     top-down 权重(默认 1)
  --masks          额外保存每个注视点的形状 mask PNG
  --selftest       运行内置自检
```

### 示例

```bat
:: 单张图片,预测 5 个注视点
saliency.exe examples\cl1.jpg --fixations 5 --out out

:: 批量处理 examples 目录,输出形状 mask
saliency.exe examples --fixations 3 --masks --out out_batch

:: 模板匹配 top-down:在零件图中寻找模板目标
saliency.exe examples\lg2.jpg --template examples\lg3.jpg --fixations 3

:: 完全可复现运行
saliency.exe examples --norandom --seed 42 --out out_repro
```

### 输出文件（每张输入图）

| 文件 | 内容 |
|---|---|
| `*_salmap.png` | 灰度显著图 |
| `*_overlay.png` | 原图叠加注视点顺序的可视化 |
| `*_fixations.csv` | 注视点坐标（原图坐标系）与对应显著值 |
| `*_mask_*.png` | （`--masks` 时）每个注视点的形状 mask |

## 项目结构

```
src/
  main.cpp       CLI、图像 IO、叠加绘制、NCC 模板匹配、selftest
  salmap.cpp     高斯金字塔、特征金字塔、中央-周边差、归一化合成
  wta.cpp        LIF 神经元、WTA 演化、IOR
  segment.cpp    fastSegmentMap、estimateShape、形态学操作
  saliency.hpp   公共数据结构与声明
  stb_image*.h   图像读写(第三方单头文件)
examples/        5 张测试图片(工业零件图像)
```

## 与 MATLAB 版的一致性说明

本实现以逐行复刻为目标，包括原始工具箱的若干固有行为，**不是 bug**：

- `conv2PreserveEnergy` 只补偿上下边界（MATLAB 原版即如此）。
- `normalizeImage` 会放大 1e-16 级浮点尾差，在完美均匀图上显著图出现条纹属预期现象。
- 降采样核 `[1,5,10,10,5,1,1]/33` 本身不对称，与原工具箱一致。

若需与原 MATLAB 结果对比，请使用 `--norandom --seed 42` 固定随机性。

## 许可与致谢

本项目是 [SaliencyToolbox](https://www.salilab.org)（Dirk B. Walther / California Institute of Technology，Copyright © 2006–2013）的衍生重实现。原始工具箱采用 **仅限非盈利研究目的** 的许可（见 Caltech 的 LICENSE 条款），本项目同样**仅供学术研究使用，禁止商业用途与再分发授权**。使用本项目即表示同意遵守原始工具箱的许可条件。

`stb_image.h` / `stb_image_write.h` 由 Sean Barrett 发布，公有领域 / MIT 双许可。

---

<a name="english"></a>

# saliency-cpp (English)

A dependency-free C++17 reimplementation of the Itti–Koch visual saliency model, ported line-by-line from Dirk B. Walther's MATLAB **SaliencyToolbox** (Caltech, 2006–2013). A single self-contained executable computes saliency maps, predicts scanpaths via a winner-take-all network, and applies inhibition of return — no MATLAB required.

## Features

- Full bottom-up pipeline: Gaussian pyramids → intensity / color (RG, BY, hue) / orientation (Gabor) features → center-surround differences → cross-scale combination into a single saliency map.
- LIF (leaky integrate-and-fire) winner-take-all network producing successive fixations, with three IOR modes (`disk` / `shape` / `none`).
- Shape estimation (proto-object segmentation) with optional per-fixation mask export.
- **Top-down modulation**: blend an external attention map (`--topdown`), or an NCC-based template-matching heat map (`--template`).
- Batch processing of whole directories; outputs saliency PNG, overlay PNG, and fixation CSV per image.
- Reproducible: `--norandom` + `--seed` pin all randomness; `--selftest` runs built-in numeric checks.
- Zero third-party dependencies beyond the single-header [stb_image](https://github.com/nothings/stb) libraries; builds on Windows with one MSVC command.

## Build

Requires Windows + MSVC (Visual Studio Build Tools):

```bat
cd src
build.bat
```

This invokes `vcvarsall.bat x64` and `cl /O2 /std:c++17 /EHsc /W3 /utf-8`, producing `saliency.exe`.

## Usage

```
saliency.exe <image-or-directory> [options]
  --fixations N    number of fixations (default 5; 0 = saliency map only)
  --out DIR        output directory (default ./out)
  --ior MODE       disk|shape|none (default shape)
  --norandom       disable coordinate jitter and noise (fully reproducible)
  --seed N         random seed (default 42)
  --topdown IMG    external top-down attention map (blended into saliency)
  --template PATH  template image or directory; NCC match heat map as top-down
                   (mutually exclusive with --topdown)
  --tdweight W     top-down weight (default 1)
  --masks          additionally save a shape-mask PNG per fixation
  --selftest       run the built-in self-test
```

### Examples

```bat
:: Single image, 5 fixations
saliency.exe examples\cl1.jpg --fixations 5 --out out

:: Batch over the examples directory, with shape masks
saliency.exe examples --fixations 3 --masks --out out_batch

:: Template-matching top-down: find the template part in a scene
saliency.exe examples\lg2.jpg --template examples\lg3.jpg --fixations 3

:: Fully reproducible run
saliency.exe examples --norandom --seed 42 --out out_repro
```

### Output files (per input image)

| File | Content |
|---|---|
| `*_salmap.png` | grayscale saliency map |
| `*_overlay.png` | input image overlaid with the fixation scanpath |
| `*_fixations.csv` | fixation coordinates (original-image space) and saliency values |
| `*_mask_*.png` | (with `--masks`) shape mask per fixation |

## Project layout

```
src/
  main.cpp       CLI, image IO, overlay rendering, NCC template matching, selftest
  salmap.cpp     Gaussian pyramids, feature pyramids, center-surround, normalization
  wta.cpp        LIF neurons, WTA evolution, IOR
  segment.cpp    fastSegmentMap, estimateShape, morphology
  saliency.hpp   shared data structures and declarations
  stb_image*.h   image loading/saving (third-party single-header)
examples/        5 test images (industrial part photographs)
```

## Fidelity notes vs. the MATLAB version

The port intentionally reproduces several quirks of the original toolbox — these are **not bugs**:

- `conv2PreserveEnergy` compensates only the top/bottom borders (as in the MATLAB original).
- `normalizeImage` amplifies 1e-16-level floating-point residue; banding in the saliency map of a perfectly uniform image is expected.
- The decimation kernel `[1,5,10,10,5,1,1]/33` is asymmetric, matching the toolbox.

For comparison against MATLAB results, pin randomness with `--norandom --seed 42`.

## License & acknowledgements

This project is a derivative reimplementation of the [SaliencyToolbox](https://www.salilab.org) by Dirk B. Walther / California Institute of Technology (Copyright © 2006–2013). The original toolbox is licensed for **non-profit research purposes only** (see Caltech's LICENSE terms); accordingly, this project is likewise **for academic research use only — no commercial use, no sublicensing**. Using this project implies acceptance of the original toolbox license conditions.

`stb_image.h` / `stb_image_write.h` by Sean Barrett, public domain / MIT dual-licensed.
