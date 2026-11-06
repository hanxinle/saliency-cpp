// segment.cpp - 形状估计:fastSegmentMap + estimateShape + 二值形态学.
// 复刻 fastSegmentMap.m 与 estimateShape.m(默认 shapeMode='shapeFM', segmentComputeType='Fast').
#include "saliency.hpp"

#include <algorithm>
#include <cmath>
#include <queue>

// 4 连通域标记,等价于 bwlabel(bw,4)(标号顺序无关紧要,只用相等性).
static std::vector<int> bwlabel4(const Map& bw, int& numLabels) {
    const int rows = bw.rows, cols = bw.cols;
    std::vector<int> lab((size_t)rows * cols, 0);
    numLabels = 0;
    std::queue<std::pair<int, int>> q;
    for (int r0 = 0; r0 < rows; ++r0) {
        for (int c0 = 0; c0 < cols; ++c0) {
            if (bw(r0, c0) == 0 || lab[(size_t)r0 * cols + c0] != 0) continue;
            ++numLabels;
            q.push({r0, c0});
            lab[(size_t)r0 * cols + c0] = numLabels;
            while (!q.empty()) {
                auto [r, c] = q.front();
                q.pop();
                const int dr[4] = {-1, 1, 0, 0}, dc[4] = {0, 0, -1, 1};
                for (int k = 0; k < 4; ++k) {
                    int nr = r + dr[k], nc = c + dc[k];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                    if (bw(nr, nc) == 0 || lab[(size_t)nr * cols + nc] != 0) continue;
                    lab[(size_t)nr * cols + nc] = numLabels;
                    q.push({nr, nc});
                }
            }
        }
    }
    return lab;
}

// 复刻 fastSegmentMap.m:阈值 = 0.05*seedVal,保留种子所在 4 连通域.
// seedR/seedC 为 1-based.
Map fastSegmentMap(const Map& map, int seedR, int seedC, double thresh) {
    Map zeros(map.rows, map.cols, 0.0);
    double seedVal = map(seedR - 1, seedC - 1);
    if (seedVal < 0.001) return zeros;

    Map bw(map.rows, map.cols);
    for (size_t i = 0; i < bw.d.size(); ++i)
        bw.d[i] = (map.d[i] / seedVal >= thresh) ? 1.0 : 0.0;

    int nl = 0;
    std::vector<int> lab = bwlabel4(bw, nl);
    int sVal = lab[(size_t)(seedR - 1) * map.cols + (seedC - 1)];
    if (sVal == 0) return zeros;
    Map out(map.rows, map.cols);
    for (size_t i = 0; i < out.d.size(); ++i)
        out.d[i] = (lab[i] == sVal) ? 1.0 : 0.0;
    return out;
}

// 5x5 菱形结构元(与 strel('disk',2) 的 13 像素邻域相同),供 imdilate/imopen/imclose.
static const int SE5[5][5] = {
    {0, 0, 1, 0, 0},
    {0, 1, 1, 1, 0},
    {1, 1, 1, 1, 1},
    {0, 1, 1, 1, 0},
    {0, 0, 1, 0, 0}};
// 3x3 十字.
static const int SE3[3][3] = {
    {0, 1, 0},
    {1, 1, 1},
    {0, 1, 0}};

static std::vector<std::pair<int, int>> seOffsets(bool smallCross) {
    std::vector<std::pair<int, int>> off;
    if (smallCross) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                if (SE3[r][c]) off.push_back({r - 1, c - 1});
    } else {
        for (int r = 0; r < 5; ++r)
            for (int c = 0; c < 5; ++c)
                if (SE5[r][c]) off.push_back({r - 2, c - 2});
    }
    return off;
}

// 二值膨胀:越界按 0(不贡献).
static Map imdilate(const Map& bin, const std::vector<std::pair<int, int>>& off) {
    Map out(bin.rows, bin.cols, 0.0);
    for (int r = 0; r < bin.rows; ++r) {
        for (int c = 0; c < bin.cols; ++c) {
            for (auto [dr, dc] : off) {
                int nr = r + dr, nc = c + dc;
                if (nr < 0 || nr >= bin.rows || nc < 0 || nc >= bin.cols) continue;
                if (bin(nr, nc) != 0) { out(r, c) = 1.0; break; }
            }
        }
    }
    return out;
}

// 二值腐蚀:越界按 1(MATLAB imerode 的边界填充行为).
static Map imerode(const Map& bin, const std::vector<std::pair<int, int>>& off) {
    Map out(bin.rows, bin.cols, 0.0);
    for (int r = 0; r < bin.rows; ++r) {
        for (int c = 0; c < bin.cols; ++c) {
            bool all = true;
            for (auto [dr, dc] : off) {
                int nr = r + dr, nc = c + dc;
                if (nr < 0 || nr >= bin.rows || nc < 0 || nc >= bin.cols) continue;  // 越界=1,不影响.
                if (bin(nr, nc) == 0) { all = false; break; }
            }
            if (all) out(r, c) = 1.0;
        }
    }
    return out;
}

// imdilate(bin, strel('disk',2)).
Map imdilateDisk2(const Map& bin) {
    return imdilate(bin, seOffsets(false));
}

// imclose(imopen(bin,se),se);smallCross=true 时用 3x3 十字,否则 5x5 菱形.
Map imopenClose(const Map& bin, bool smallCross) {
    auto off = seOffsets(smallCross);
    Map opened = imdilate(imerode(bin, off), off);
    return imerode(imdilate(opened, off), off);
}

// 复刻 estimateShape.m(默认 shapeFM;从最深的图向浅层尝试,接受第一个 0<areaRatio<0.1 的分割).
ShapeData estimateShape(const SaliencyResult& sr, int winR, int winC, const SalParams& params) {
    ShapeData sd;
    const size_t widx = (size_t)(winR - 1) * sr.salmap.cols + (winC - 1);

    // winMap 层级:salmap -> 响应最大的 CM -> 该 CM 响应最大的 FM(MATLAB max 取首个最大值,确定性).
    std::vector<Map> winMaps;
    winMaps.push_back(sr.salmap);

    double cmMax = -1;
    int cmIdx = -1;
    for (size_t f = 0; f < sr.feats.size(); ++f) {
        double v = sr.feats[f].CM.d[widx];
        if (v > cmMax) { cmMax = v; cmIdx = (int)f; }
    }
    if (cmMax > 0 && cmIdx >= 0) {
        winMaps.push_back(sr.feats[cmIdx].CM);
        // 该特征所有 FM(列主序摊平:FM(pyr,fm) -> idx = fm*numPyr + pyr).
        const auto& FM = sr.feats[cmIdx].FM;
        int numPyr = (int)FM.size();
        double fmMax = -1;
        int fmPyr = -1, fmMap = -1;
        int nFM = numPyr ? (int)FM[0].size() : 0;
        for (int fm = 0; fm < nFM; ++fm) {
            for (int p = 0; p < numPyr; ++p) {
                double v = FM[p][fm].d[widx];
                if (v > fmMax) { fmMax = v; fmPyr = p; fmMap = fm; }
            }
        }
        if (fmMax > 0 && fmPyr >= 0)
            winMaps.push_back(FM[fmPyr][fmMap]);
    }

    // 从最深一层开始尝试分割.
    for (int idx = (int)winMaps.size() - 1; idx >= 0; --idx) {
        Map binMap = fastSegmentMap(winMaps[idx], winR, winC, 0.05);
        double area = 0;
        for (double v : binMap.d) area += v;
        double areaRatio = area / binMap.d.size();
        if (areaRatio <= 0 || areaRatio >= 0.1) continue;

        sd.valid = true;
        sd.winningMap = winMaps[idx];
        sd.winR = winR;
        sd.winC = winC;
        sd.iorMask = imdilateDisk2(binMap);

        // binaryMap 后处理:5x5 菱形开闭运算,失败改 3x3 十字;最后保留种子连通域.
        Map newMap;
        Map tmp = imopenClose(binMap, false);
        if (tmp(winR - 1, winC - 1) > 0) {
            double s = 0;
            for (double v : tmp.d) s += v;
            if (s > 0) newMap = tmp;
        } else {
            tmp = imopenClose(binMap, true);
            double s = 0;
            for (double v : tmp.d) s += v;
            if (tmp(winR - 1, winC - 1) > 0 && s > 0) newMap = tmp;
        }
        if (!newMap.empty()) {
            int nl = 0;
            std::vector<int> lab = bwlabel4(newMap, nl);
            int sVal = lab[(size_t)(winR - 1) * newMap.cols + (winC - 1)];
            Map cleaned(newMap.rows, newMap.cols, 0.0);
            if (sVal > 0)
                for (size_t i = 0; i < cleaned.d.size(); ++i)
                    cleaned.d[i] = (lab[i] == sVal) ? 1.0 : 0.0;
            sd.binaryMap = cleaned;
        } else {
            sd.binaryMap = binMap;
        }
        return sd;
    }
    return sd;  // 全部失败:valid=false.
}
