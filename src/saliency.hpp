// saliency.hpp - SaliencyToolbox (Itti-Koch) C++ 复刻版公共声明.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

// 2d double 矩阵,行主序,等价于 MATLAB 的 map.data.
struct Map {
    int rows = 0;
    int cols = 0;
    std::vector<double> d;

    Map() = default;
    Map(int r, int c, double v = 0.0) : rows(r), cols(c), d((size_t)r * c, v) {}

    double& operator()(int r, int c) { return d[(size_t)r * cols + c]; }
    double operator()(int r, int c) const { return d[(size_t)r * cols + c]; }
    bool empty() const { return d.empty(); }
};

// 输入图像,三通道平面存储,值域 [0,1].灰度图时 g/b 与 r 相同且 isColor=false.
struct Image {
    int w = 0, h = 0;
    bool isColor = false;
    Map r, g, b;
    std::string filename;
};

struct GaborParams {
    double filterPeriod = 7;
    double elongation = 1;
    int filterSize = 9;
    double stddev = 2.3333;
    std::vector<double> phases = {0, 90};
};

// 金字塔层级参数,保持 MATLAB 的 1-based 层号.
struct LevelParams {
    int minLevel = 3;
    int maxLevel = 5;
    int minDelta = 3;
    int maxDelta = 4;
    int mapLevel = 5;
};

enum class IORType { None, Disk, Shape };

struct SalParams {
    int foaSize = -1;                                    // 传图像尺寸后 = round(max(w,h)/6).
    std::vector<std::string> features = {"Color", "Intensities", "Orientations"};
    std::vector<double> weights = {1, 1, 1};
    IORType iorType = IORType::Shape;
    LevelParams lp;
    int numIter = 3;
    bool useRandom = true;
    double IORdecay = 0.9999;
    double smOutputRange = 1e-9;
    double noiseAmpl = 1e-17;
    double noiseConst = 1e-14;
    GaborParams gabor;
    std::vector<double> oriAngles = {0, 45, 90, 135};
    double topdownWeight = 1.0;  // 外部 top-down 注意图的权重(--tdweight).
};

// 单个特征的处理中间数据(对应 MATLAB 的 saliencyData(f)).
struct FeatureData {
    std::string label;
    struct Pyramid {
        std::string label;
        std::vector<Map> levels;  // levels[i] 对应 1-based 层号 i+1.
    };
    std::vector<Pyramid> pyr;
    std::vector<std::vector<Map>> FM;  // FM[p][fm]:第 p 个金字塔的第 fm 张特征图(已归一化).
    Map CM;                            // 该特征的 conspicuity map(加权后).
};

struct SaliencyResult {
    Map salmap;
    std::vector<FeatureData> feats;
};

// 形状估计结果(对应 MATLAB 的 shapeData),valid=false 表示分割失败.
struct ShapeData {
    bool valid = false;
    Map winningMap;
    Map iorMask;    // 膨胀后的二值 mask,供 shapeIOR.
    Map binaryMap;  // 形态学清理后的二值形状.
    int winR = 0, winC = 0;  // winner 在 winningMap 中的 1-based 坐标.
};

// ---- salmap.cpp ----
Map resizeNearest(const Map& m, int newRows, int newCols);
Map sepConv2PreserveEnergy(const std::vector<double>& fy, const std::vector<double>& fx, const Map& data);
Map conv2PreserveEnergy(const Map& src, const Map& f);
void attenuateBordersInPlace(Map& m, int borderSize);
std::vector<double> gaussian1d(double peak, double sigma, int maxhw);
Map normalizeImage(Map m, double lo, double hi);
Map maxNormalizeIterative(const Map& data, int numIter, double lo, double hi);
Map gaussianSubsample(const Map& img);
std::vector<Map> makeDyadicPyramid(const Map& m);
SaliencyResult makeSaliencyMap(const Image& img, SalParams& params,
                               const std::vector<Map>& topdownPyramid = {});

// ---- wta.cpp ----
struct LIF {
    double timeStep = 1e-4;
    double Eleak = 0, Eexc = 100e-3, Einh = -20e-3;
    double Gleak = 1e-8, Gexc = 0, Ginh = 0, GinhDecay = 1;
    double Ginput = 5e-8, Vthresh = 1e-3, C = 1e-9, time = 0;
    bool DoesFire = true;
    bool ginhIsMap = false;
    Map V, I, GinhMap;
};

struct WTA {
    LIF sm, exc, inhib;
};

// rngState:确定性随机源(复现用).
struct Rng {
    uint64_t s;
    explicit Rng(uint64_t seed = 42) : s(seed ? seed : 0x9E3779B97F4A7C15ULL) {}
    uint64_t next() {  // splitmix64.
        uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
    double uniform() { return (next() >> 11) * (1.0 / 9007199254740992.0); }  // [0,1).
};

WTA initializeWTA(const Map& salmap, const SalParams& params, Rng& rng);
// 推进一步;有胜者时返回其 1-based (row,col),否则返回 (-1,-1).
std::pair<int, int> evolveWTA(WTA& wta);
void diskIOR(WTA& wta, int winR, int winC, const SalParams& params);
void shapeIOR(WTA& wta, int winR, int winC, const SalParams& params, const ShapeData& sd);
void applyIOR(WTA& wta, int winR, int winC, const SalParams& params, const ShapeData& sd);
// winner(1-based map 坐标) -> 原图 1-based (row,col).
std::pair<int, int> winnerToImgCoords(int winR, int winC, const SalParams& params, Rng& rng);

// ---- segment.cpp ----
Map fastSegmentMap(const Map& map, int seedR, int seedC, double thresh = 0.05);
ShapeData estimateShape(const SaliencyResult& sr, int winR, int winC, const SalParams& params);
Map imdilateDisk2(const Map& bin);
Map imopenClose(const Map& bin, bool smallCross);
