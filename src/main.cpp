// main.cpp - SaliencyToolbox C++ 复刻版命令行入口.
// 用法: saliency.exe <图片或目录> [--fixations N] [--out DIR] [--ior disk|shape|none]
//                    [--norandom] [--seed N] [--topdown IMG] [--tdweight W]
//                    [--template IMG_OR_DIR] [--masks] [--selftest]
#include "saliency.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <stdexcept>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#define STBI_FAILURE_USERMSG
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace fs = std::filesystem;

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
// 让中文输出在 Windows 控制台正常显示.
static struct ConsoleUtf8 {
    ConsoleUtf8() { SetConsoleOutputCP(CP_UTF8); }
} g_consoleUtf8;
#endif

// ---- 图像 IO ----

static Image loadImage(const std::string& path) {
    int w = 0, h = 0, n = 0;
    unsigned char* data = stbi_load(path.c_str(), &w, &h, &n, 0);
    if (!data) throw std::runtime_error("cannot read image: " + path + " (" + stbi_failure_reason() + ")");
    Image img;
    img.w = w;
    img.h = h;
    img.isColor = (n >= 3);
    img.filename = path;
    img.r = Map(h, w);
    img.g = Map(h, w);
    img.b = Map(h, w);
    for (int i = 0; i < w * h; ++i) {
        if (n >= 3) {
            img.r.d[i] = data[3 * i] / 255.0;
            img.g.d[i] = data[3 * i + 1] / 255.0;
            img.b.d[i] = data[3 * i + 2] / 255.0;
        } else if (n == 2) {  // 灰度+alpha.
            img.r.d[i] = img.g.d[i] = img.b.d[i] = data[2 * i] / 255.0;
        } else {
            img.r.d[i] = img.g.d[i] = img.b.d[i] = data[i] / 255.0;
        }
    }
    stbi_image_free(data);
    return img;
}

// 灰度图加载(top-down 注意图/模板用),任意尺寸,值域 [0,1].
static Map loadGrayMap(const std::string& path) {
    int w = 0, h = 0, n = 0;
    unsigned char* data = stbi_load(path.c_str(), &w, &h, &n, 0);
    if (!data) throw std::runtime_error("cannot read grayscale image: " + path);
    Map m(h, w);
    for (int i = 0; i < w * h; ++i) {
        if (n >= 3)
            m.d[i] = (0.2989 * data[3 * i] + 0.5870 * data[3 * i + 1] + 0.1140 * data[3 * i + 2]) / 255.0;
        else if (n == 2)
            m.d[i] = data[2 * i] / 255.0;
        else
            m.d[i] = data[i] / 255.0;
    }
    stbi_image_free(data);
    return m;
}

// Map 归一化到 [0,255] 存灰度 PNG.
static void saveGrayPng(const std::string& path, const Map& m) {
    Map norm = normalizeImage(m, 0, 255);
    std::vector<unsigned char> px(norm.d.size());
    for (size_t i = 0; i < px.size(); ++i)
        px[i] = (unsigned char)std::clamp(std::lround(norm.d[i]), 0l, 255l);
    stbi_write_png(path.c_str(), norm.cols, norm.rows, 1, px.data(), norm.cols);
}

// ---- 模板匹配(NCC)----

// 彩色图转灰度,值域 [0,1].
static Map toGray(const Image& img) {
    if (!img.isColor) return img.r;
    Map m(img.h, img.w);
    for (int i = 0; i < img.w * img.h; ++i)
        m.d[i] = 0.2989 * img.r.d[i] + 0.5870 * img.g.d[i] + 0.1140 * img.b.d[i];
    return m;
}

// 零均值归一化互相关(NCC),有效窗口的相关系数写在窗口中心,其余为 0,负相关裁为 0.
// 积分图算窗口方差,互相关项直接滑窗;86 张零件图 x 8 模板约几分钟,可接受.
static Map nccMatch(const Map& img, const Map& tpl) {
    const int h = img.rows, w = img.cols, th = tpl.rows, tw = tpl.cols;
    Map out(h, w, 0.0);
    if (th > h || tw > w) return out;
    const int sw = w + 1;
    auto at = [&](int r, int c) { return (size_t)r * sw + c; };
    std::vector<double> S((size_t)(h + 1) * sw, 0.0), Q((size_t)(h + 1) * sw, 0.0);
    for (int r = 0; r < h; ++r)
        for (int c = 0; c < w; ++c) {
            double v = img(r, c);
            S[at(r + 1, c + 1)] = S[at(r, c + 1)] + S[at(r + 1, c)] - S[at(r, c)] + v;
            Q[at(r + 1, c + 1)] = Q[at(r, c + 1)] + Q[at(r + 1, c)] - Q[at(r, c)] + v * v;
        }
    const double n = (double)th * tw;
    double tm = 0;
    for (double v : tpl.d) tm += v;
    tm /= n;
    std::vector<double> tz(tpl.d.size());
    double tnorm = 0;
    for (size_t i = 0; i < tz.size(); ++i) {
        tz[i] = tpl.d[i] - tm;
        tnorm += tz[i] * tz[i];
    }
    tnorm = std::sqrt(tnorm);
    if (tnorm == 0) return out;  // 常数模板无区分度.
    for (int y = 0; y <= h - th; ++y)
        for (int x = 0; x <= w - tw; ++x) {
            double ws = S[at(y + th, x + tw)] - S[at(y, x + tw)] - S[at(y + th, x)] + S[at(y, x)];
            double wq = Q[at(y + th, x + tw)] - Q[at(y, x + tw)] - Q[at(y + th, x)] + Q[at(y, x)];
            double wvar = wq - ws * ws / n;
            if (wvar <= 1e-12) continue;  // 常数窗口无区分度.
            double num = 0;
            for (int i = 0; i < th; ++i) {
                const double* irow = &img.d[(size_t)(y + i) * w + x];
                const double* trow = &tz[(size_t)i * tw];
                for (int j = 0; j < tw; ++j) num += irow[j] * trow[j];
            }
            double r = num / (std::sqrt(wvar) * tnorm);
            // 负相关不是注意目标,裁掉.
            out(y + th / 2, x + tw / 2) = r > 0 ? r : 0;
        }
    return out;
}

// 多模板逐像素取 max 融合成一张热力图.
static Map matchTemplatesMax(const Map& gray, const std::vector<std::string>& tplNames,
                             const std::vector<Map>& tpls) {
    Map best(gray.rows, gray.cols, 0.0);
    for (size_t i = 0; i < tpls.size(); ++i) {
        if (tpls[i].rows > gray.rows || tpls[i].cols > gray.cols) {
            std::fprintf(stderr, "  模板 %s (%dx%d) 大于图像,跳过\n", tplNames[i].c_str(),
                         tpls[i].cols, tpls[i].rows);
            continue;
        }
        Map m = nccMatch(gray, tpls[i]);
        for (size_t k = 0; k < m.d.size(); ++k)
            if (m.d[k] > best.d[k]) best.d[k] = m.d[k];
    }
    return best;
}

// ---- overlay 绘制 ----

struct Canvas {
    int w, h;
    std::vector<unsigned char> rgb;  // 3 通道.
    void set(int x, int y, unsigned char r, unsigned char g, unsigned char b) {
        if (x < 0 || x >= w || y < 0 || y >= h) return;
        size_t i = 3 * ((size_t)y * w + x);
        rgb[i] = r; rgb[i + 1] = g; rgb[i + 2] = b;
    }
};

static void drawRing(Canvas& cv, int cx, int cy, int radius, unsigned char r, unsigned char g, unsigned char b, int thick = 2) {
    for (int t = 0; t < thick; ++t) {
        int rad = radius - t;
        if (rad <= 0) break;
        int steps = std::max(64, rad * 8);
        for (int i = 0; i < steps; ++i) {
            double a = 2 * 3.14159265358979 * i / steps;
            cv.set(cx + (int)std::lround(rad * std::cos(a)), cy + (int)std::lround(rad * std::sin(a)), r, g, b);
        }
    }
}

static void drawLine(Canvas& cv, int x0, int y0, int x1, int y1, unsigned char r, unsigned char g, unsigned char b) {
    int dx = std::abs(x1 - x0), dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1, sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    while (true) {
        cv.set(x0, y0, r, g, b);
        if (x0 == x1 && y0 == y1) break;
        int e2 = 2 * err;
        if (e2 >= dy) { err += dy; x0 += sx; }
        if (e2 <= dx) { err += dx; y0 += sy; }
    }
}

// 3x5 数字点阵.
static const uint8_t DIGITS[10][5] = {
    {0b111, 0b101, 0b101, 0b101, 0b111}, {0b010, 0b110, 0b010, 0b010, 0b111},
    {0b111, 0b001, 0b111, 0b100, 0b111}, {0b111, 0b001, 0b111, 0b001, 0b111},
    {0b101, 0b101, 0b111, 0b001, 0b001}, {0b111, 0b100, 0b111, 0b001, 0b111},
    {0b111, 0b100, 0b111, 0b101, 0b111}, {0b111, 0b001, 0b001, 0b010, 0b010},
    {0b111, 0b101, 0b111, 0b101, 0b111}, {0b111, 0b101, 0b111, 0b001, 0b111}};

static void drawDigit(Canvas& cv, int x, int y, int num, unsigned char r, unsigned char g, unsigned char b, int scale = 2) {
    char buf[8];
    std::snprintf(buf, sizeof(buf), "%d", num);
    for (int k = 0; buf[k]; ++k) {
        const uint8_t* glyph = DIGITS[buf[k] - '0'];
        for (int row = 0; row < 5; ++row)
            for (int col = 0; col < 3; ++col)
                if (glyph[row] & (1 << (2 - col)))
                    for (int sy = 0; sy < scale; ++sy)
                        for (int sx = 0; sx < scale; ++sx)
                            cv.set(x + k * 4 * scale + col * scale + sx, y + row * scale + sy, r, g, b);
    }
}

// ---- 单图处理 ----

struct Options {
    int fixations = 5;
    std::string outDir = "out";
    IORType iorType = IORType::Shape;
    bool useRandom = true;
    uint64_t seed = 42;
    std::string topdownPath;
    std::string templatePath;  // 模板文件或目录,与 --topdown 互斥.
    double tdweight = 1.0;
    bool saveMasks = false;
};

static void processImage(const std::string& path, const Options& opt, Rng& rng,
                         const std::vector<Map>& topdownPyr,
                         const std::vector<std::string>& tplNames,
                         const std::vector<Map>& tpls) {
    Image img = loadImage(path);
    std::fprintf(stderr, "处理 %s (%dx%d%s) ...\n", path.c_str(), img.w, img.h,
                 img.isColor ? "" : ", 灰度");

    std::string stem = fs::path(path).stem().string();
    std::string prefix = (fs::path(opt.outDir) / stem).string();

    // 模板匹配热力图作为 top-down 输入.
    std::vector<Map> tdPyr = topdownPyr;
    if (!tpls.empty()) {
        Map tmap = matchTemplatesMax(toGray(img), tplNames, tpls);
        saveGrayPng(prefix + "_tmap.png", tmap);
        tdPyr = makeDyadicPyramid(tmap);
    }

    SalParams params;
    params.foaSize = (int)std::lround(std::max(img.w, img.h) / 6.0);
    params.iorType = opt.iorType;
    params.useRandom = opt.useRandom;
    params.topdownWeight = opt.tdweight;
    if (!opt.useRandom) { params.noiseAmpl = 0; params.noiseConst = 0; }

    SaliencyResult sr = makeSaliencyMap(img, params, tdPyr);

    // 显著图(放大回原图尺寸便于查看).
    saveGrayPng(prefix + "_salmap.png", resizeNearest(sr.salmap, img.h, img.w));

    // WTA 注视点序列.
    WTA wta;
    if (opt.fixations > 0) wta = initializeWTA(sr.salmap, params, rng);

    std::vector<std::pair<int, int>> fixPts;   // 原图坐标 (row,col),1-based.
    std::vector<double> fixTimeMs;
    std::vector<ShapeData> shapes;
    const long maxSteps = 1000000;             // MATLAB 原版无上限,这里加保护.
    for (int f = 0; f < opt.fixations; ++f) {
        std::pair<int, int> winner{-1, -1};
        long steps = 0;
        while (winner.first == -1) {
            winner = evolveWTA(wta);
            if (++steps > maxSteps)
                throw std::runtime_error("WTA 超过 1e6 步仍未产生胜者,疑似数值问题");
        }
        ShapeData sd = estimateShape(sr, winner.first, winner.second, params);
        applyIOR(wta, winner.first, winner.second, params, sd);
        auto imgCo = winnerToImgCoords(winner.first, winner.second, params, rng);
        fixPts.push_back(imgCo);
        fixTimeMs.push_back(wta.exc.time * 1000.0);
        shapes.push_back(sd);
        std::fprintf(stderr, "  注视点 %d: (%d,%d) @ %.1f ms\n", f + 1, imgCo.second, imgCo.first,
                     fixTimeMs.back());
    }

    // CSV.
    {
        std::string csvPath = prefix + "_fixations.csv";
        FILE* fp = std::fopen(csvPath.c_str(), "w");
        if (!fp) throw std::runtime_error("cannot write " + csvPath);
        std::fprintf(fp, "fix,x,y,time_ms\n");
        for (size_t i = 0; i < fixPts.size(); ++i)
            std::fprintf(fp, "%zu,%d,%d,%.1f\n", i + 1, fixPts[i].second, fixPts[i].first, fixTimeMs[i]);
        std::fclose(fp);
    }

    // overlay:原图 + 形状轮廓(绿)+ 轨迹(黄)+ 末点(红).
    Canvas cv{img.w, img.h, std::vector<unsigned char>((size_t)img.w * img.h * 3)};
    for (int i = 0; i < img.w * img.h; ++i) {
        cv.rgb[3 * i] = (unsigned char)std::lround(img.r.d[i] * 255);
        cv.rgb[3 * i + 1] = (unsigned char)std::lround(img.g.d[i] * 255);
        cv.rgb[3 * i + 2] = (unsigned char)std::lround(img.b.d[i] * 255);
    }
    // 形状轮廓.
    for (const ShapeData& sd : shapes) {
        if (!sd.valid) continue;
        Map mask = resizeNearest(sd.binaryMap, img.h, img.w);
        for (int r = 0; r < img.h; ++r)
            for (int c = 0; c < img.w; ++c) {
                if (mask(r, c) == 0) continue;
                bool boundary = (r == 0 || r == img.h - 1 || c == 0 || c == img.w - 1) ||
                                mask(r - 1, c) == 0 || mask(r + 1, c) == 0 ||
                                mask(r, c - 1) == 0 || mask(r, c + 1) == 0;
                if (boundary) cv.set(c, r, 0, 255, 0);
            }
    }
    // 轨迹与注视点.
    int foaRadius = params.foaSize / 2;
    for (size_t i = 0; i < fixPts.size(); ++i) {
        int x = fixPts[i].second - 1, y = fixPts[i].first - 1;  // 转 0-based.
        if (i > 0)
            drawLine(cv, fixPts[i - 1].second - 1, fixPts[i - 1].first - 1, x, y, 255, 255, 0);
        bool last = (i + 1 == fixPts.size());
        drawRing(cv, x, y, foaRadius, last ? 255 : 255, last ? 0 : 255, 0);
        drawDigit(cv, x + 4, y + 4, (int)i + 1, 255, 255, 255);
    }
    stbi_write_png((prefix + "_overlay.png").c_str(), cv.w, cv.h, 3, cv.rgb.data(), cv.w * 3);

    // 可选:每个注视点的形状 mask.
    if (opt.saveMasks) {
        for (size_t i = 0; i < shapes.size(); ++i) {
            if (!shapes[i].valid) continue;
            char suffix[64];
            std::snprintf(suffix, sizeof(suffix), "_mask_%zu.png", i + 1);
            saveGrayPng(prefix + suffix, shapes[i].binaryMap);
        }
    }
}

// ---- selftest ----

static Image syntheticImage(int w, int h, bool color) {
    Image img;
    img.w = w; img.h = h; img.isColor = color;
    img.r = Map(h, w, 0.0); img.g = Map(h, w, 0.0); img.b = Map(h, w, 0.0);
    return img;
}

static int selftest() {
    int fails = 0;
    // 用 256x256(mapLevel 16x16)避免小图上迭代归一化的固有爆炸,贴近真实使用尺度.

    // 1. 黑底中央白方块:显著图峰值应落在方块区域,首个注视点也在方块内.
    {
        Image img = syntheticImage(256, 256, false);
        for (int r = 96; r < 160; ++r)
            for (int c = 96; c < 160; ++c) {
                img.r(r, c) = img.g(r, c) = img.b(r, c) = 1.0;
            }
        SalParams params;
        params.foaSize = 43;
        params.useRandom = false;
        params.noiseAmpl = 0; params.noiseConst = 0;
        SaliencyResult sr = makeSaliencyMap(img, params);
        // 显著图 16x16,方块对应 map 行列约 6..10(0-based).
        double mx = -1;
        size_t mi = 0;
        for (size_t i = 0; i < sr.salmap.d.size(); ++i)
            if (sr.salmap.d[i] > mx) { mx = sr.salmap.d[i]; mi = i; }
        int mr = (int)(mi / sr.salmap.cols), mc = (int)(mi % sr.salmap.cols);
        bool ok = (mr >= 5 && mr <= 10 && mc >= 5 && mc <= 10 && mx > 0);
        std::fprintf(stderr, "[selftest] 白方块显著性峰值在 map(%d,%d)=%.3f ... %s\n",
                     mr, mc, mx, ok ? "PASS" : "FAIL");
        if (!ok) ++fails;

        Rng rng(1);
        WTA wta = initializeWTA(sr.salmap, params, rng);
        std::pair<int, int> winner{-1, -1};
        long steps = 0;
        while (winner.first == -1 && ++steps < 1000000) winner = evolveWTA(wta);
        auto imgCo = winnerToImgCoords(winner.first, winner.second, params, rng);
        bool ok2 = imgCo.first >= 80 && imgCo.first <= 176 && imgCo.second >= 80 && imgCo.second <= 176;
        std::fprintf(stderr, "[selftest] 首个注视点 img(%d,%d),应落在方块附近 ... %s\n",
                     imgCo.first, imgCo.second, ok2 ? "PASS" : "FAIL");
        if (!ok2) ++fails;
    }

    // 2a. 纯灰图仅强度通道:无对比,center-surround 输出为 1e-16 量级浮点尾差.
    // normalizeImage 会把该尾差放大到 [0,10] 并被迭代归一化继续放大,MATLAB 原版同样如此.
    // (完美均匀图是病态输入,真实图不会出现).此处只验证输出有限、无 NaN/Inf.
    {
        Image img = syntheticImage(256, 256, false);
        for (double& v : img.r.d) v = 0.5;
        img.g = img.r; img.b = img.r;
        SalParams params;
        params.foaSize = 43;
        params.features = {"Intensities"};
        params.weights = {1};
        params.useRandom = false;
        params.noiseAmpl = 0; params.noiseConst = 0;
        SaliencyResult sr = makeSaliencyMap(img, params);
        bool ok = true;
        for (double v : sr.salmap.d) if (!std::isfinite(v)) { ok = false; break; }
        std::fprintf(stderr, "[selftest] 纯灰图(仅强度)显著图全部有限 ... %s\n",
                     ok ? "PASS" : "FAIL");
        if (!ok) ++fails;
    }

    // 2b. 纯灰图全特征:Gabor 边界效应与归一化噪声放大都会产生响应(与 MATLAB 一致),只验证有限.
    {
        Image img = syntheticImage(256, 256, false);
        for (double& v : img.r.d) v = 0.5;
        img.g = img.r; img.b = img.r;
        SalParams params;
        params.foaSize = 43;
        params.useRandom = false;
        params.noiseAmpl = 0; params.noiseConst = 0;
        SaliencyResult sr = makeSaliencyMap(img, params);
        bool ok = true;
        for (double v : sr.salmap.d) if (!std::isfinite(v)) { ok = false; break; }
        std::fprintf(stderr, "[selftest] 纯灰图(全特征)显著图全部有限 ... %s\n",
                     ok ? "PASS" : "FAIL");
        if (!ok) ++fails;
    }

    // 3. 绿底红方块:颜色对立通道应产生显著性,峰值在方块区域.
    {
        Image img = syntheticImage(256, 256, true);
        for (int i = 0; i < 256 * 256; ++i) { img.r.d[i] = 0; img.g.d[i] = 0.8; img.b.d[i] = 0; }
        for (int r = 96; r < 160; ++r)
            for (int c = 96; c < 160; ++c) { img.r(r, c) = 0.8; img.g(r, c) = 0; img.b(r, c) = 0; }
        SalParams params;
        params.foaSize = 43;
        params.useRandom = false;
        params.noiseAmpl = 0; params.noiseConst = 0;
        SaliencyResult sr = makeSaliencyMap(img, params);
        double mx = -1;
        size_t mi = 0;
        for (size_t i = 0; i < sr.salmap.d.size(); ++i)
            if (sr.salmap.d[i] > mx) { mx = sr.salmap.d[i]; mi = i; }
        int mr = (int)(mi / sr.salmap.cols), mc = (int)(mi % sr.salmap.cols);
        bool ok = (mr >= 5 && mr <= 10 && mc >= 5 && mc <= 10 && mx > 0);
        std::fprintf(stderr, "[selftest] 红方块(绿底)显著性峰值在 map(%d,%d)=%.3f ... %s\n",
                     mr, mc, mx, ok ? "PASS" : "FAIL");
        if (!ok) ++fails;
    }

    // 4. NCC 模板匹配:图案贴到已知位置,峰值应落在贴入窗口中心且接近 1.
    {
        Map img(180, 200, 0.0);
        Map pat(24, 30);
        for (int r = 0; r < 24; ++r)
            for (int c = 0; c < 30; ++c)
                pat(r, c) = 0.5 + 0.5 * std::sin(r * 1.7) * std::cos(c * 2.3);
        for (int r = 0; r < 24; ++r)
            for (int c = 0; c < 30; ++c) img(50 + r, 60 + c) = pat(r, c);
        Map m = nccMatch(img, pat);
        double mx = -1;
        size_t mi = 0;
        for (size_t i = 0; i < m.d.size(); ++i)
            if (m.d[i] > mx) { mx = m.d[i]; mi = i; }
        int mr = (int)(mi / m.cols), mc = (int)(mi % m.cols);
        bool ok = (mr == 50 + 12 && mc == 60 + 15 && mx > 0.999);
        std::fprintf(stderr, "[selftest] NCC 峰值在 (%d,%d)=%.4f,期望 (%d,%d) ... %s\n",
                     mr, mc, mx, 50 + 12, 60 + 15, ok ? "PASS" : "FAIL");
        if (!ok) ++fails;
    }

    std::fprintf(stderr, "[selftest] %s\n", fails == 0 ? "全部通过" : "存在失败项");
    return fails == 0 ? 0 : 1;
}

// ---- CLI ----

// 列出路径下的图片;单文件则直接返回该文件.
static std::vector<std::string> listImages(const std::string& path) {
    std::vector<std::string> files;
    if (!fs::is_directory(path)) {
        files.push_back(path);
        return files;
    }
    for (const auto& e : fs::directory_iterator(path)) {
        if (!e.is_regular_file()) continue;
        std::string ext = e.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".jpg" || ext == ".jpeg" || ext == ".png" || ext == ".bmp")
            files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

static void usage() {
    std::fprintf(stderr,
        "用法: saliency.exe <图片或目录> [选项]\n"
        "  --fixations N    注视点个数(默认 5,0 = 只算显著图)\n"
        "  --out DIR        输出目录(默认 ./out)\n"
        "  --ior MODE       disk|shape|none(默认 shape)\n"
        "  --norandom       关闭坐标 jitter 与噪声(完全可复现)\n"
        "  --seed N         随机种子(默认 42)\n"
        "  --topdown IMG    外部 top-down 注意图(加权混入显著图)\n"
        "  --template PATH  模板文件或目录,NCC 匹配热力图作为 top-down(与 --topdown 互斥)\n"
        "  --tdweight W     top-down 权重(默认 1)\n"
        "  --masks          额外保存每个注视点的形状 mask PNG\n"
        "  --selftest       运行内置自检\n");
}

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);
    for (const std::string& a : args) {
        if (a == "--selftest") return selftest();
    }
    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        usage();
        return args.empty() ? 1 : 0;
    }

    std::string input = args[0];
    Options opt;
    for (size_t i = 1; i < args.size(); ++i) {
        const std::string& a = args[i];
        auto needVal = [&](const char* name) -> std::string {
            if (++i >= args.size()) throw std::runtime_error(std::string("缺参数值: ") + name);
            return args[i];
        };
        if (a == "--fixations") opt.fixations = std::stoi(needVal("--fixations"));
        else if (a == "--out") opt.outDir = needVal("--out");
        else if (a == "--ior") {
            std::string v = needVal("--ior");
            if (v == "disk") opt.iorType = IORType::Disk;
            else if (v == "shape") opt.iorType = IORType::Shape;
            else if (v == "none") opt.iorType = IORType::None;
            else throw std::runtime_error("未知 --ior: " + v);
        }
        else if (a == "--norandom") opt.useRandom = false;
        else if (a == "--seed") opt.seed = std::stoull(needVal("--seed"));
        else if (a == "--topdown") opt.topdownPath = needVal("--topdown");
        else if (a == "--template") opt.templatePath = needVal("--template");
        else if (a == "--tdweight") opt.tdweight = std::stod(needVal("--tdweight"));
        else if (a == "--masks") opt.saveMasks = true;
        else throw std::runtime_error("未知选项: " + a);
    }

    try {
        fs::create_directories(opt.outDir);
        if (!opt.topdownPath.empty() && !opt.templatePath.empty())
            throw std::runtime_error("--topdown 与 --template 互斥,只能选一个");

        std::vector<Map> topdownPyr;
        if (!opt.topdownPath.empty())
            topdownPyr = makeDyadicPyramid(loadGrayMap(opt.topdownPath));

        std::vector<std::string> tplNames;
        std::vector<Map> tpls;
        if (!opt.templatePath.empty()) {
            for (const std::string& f : listImages(opt.templatePath)) {
                tplNames.push_back(fs::path(f).filename().string());
                tpls.push_back(loadGrayMap(f));
            }
            if (tpls.empty())
                throw std::runtime_error("模板路径下没有可用图片: " + opt.templatePath);
            std::fprintf(stderr, "加载 %zu 个模板\n", tpls.size());
        }

        Rng rng(opt.seed);

        if (fs::is_directory(input)) {
            std::vector<std::string> files = listImages(input);
            std::fprintf(stderr, "发现 %zu 张图片\n", files.size());
            for (const std::string& f : files) processImage(f, opt, rng, topdownPyr, tplNames, tpls);
        } else {
            processImage(input, opt, rng, topdownPyr, tplNames, tpls);
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "错误: %s\n", e.what());
        return 1;
    }
    return 0;
}
