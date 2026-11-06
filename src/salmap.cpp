// salmap.cpp - 显著图计算:金字塔、特征、center-surround、归一化、合成.
// 逐行复刻 SaliencyToolbox 的 MATLAB 实现(默认 dyadic + Iterative 归一化).
#include "saliency.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

static const double PI = 3.14159265358979323846;

// MATLAB imresize(...,'nearest'):0-based src = floor((dst+0.5)*in/out),截到 in-1.
Map resizeNearest(const Map& m, int newRows, int newCols) {
    Map out(newRows, newCols);
    for (int r = 0; r < newRows; ++r) {
        int sr = (int)std::floor((r + 0.5) * m.rows / newRows);
        sr = std::min(sr, m.rows - 1);
        for (int c = 0; c < newCols; ++c) {
            int sc = (int)std::floor((c + 0.5) * m.cols / newCols);
            sc = std::min(sc, m.cols - 1);
            out(r, c) = m(sr, sc);
        }
    }
    return out;
}

// 可分离卷积,零填充 same,边界按"实际参与权重和"做能量补偿.
// 复刻 sepConv2PreserveEnergy.m(含小数据的整图补偿分支).
Map sepConv2PreserveEnergy(const std::vector<double>& fy, const std::vector<double>& fx, const Map& data) {
    const int sd1 = data.rows, sd2 = data.cols;
    // conv2(fy,fx,data,'same'):先零填充卷积(对称核,相关=卷积;为非对称核仍按翻转实现).
    Map result(sd1, sd2, 0.0);
    const int fl1 = (int)fy.size(), fl2 = (int)fx.size();
    const int h1 = (fl1 - 1) / 2, h2 = (fl2 - 1) / 2;
    for (int r = 0; r < sd1; ++r) {
        for (int c = 0; c < sd2; ++c) {
            double acc = 0.0;
            for (int i = 0; i < fl1; ++i) {
                int sr = r + h1 - i;  // 卷积翻转.
                if (sr < 0 || sr >= sd1) continue;
                for (int j = 0; j < fl2; ++j) {
                    int sc = c + h2 - j;
                    if (sc < 0 || sc >= sd2) continue;
                    acc += fy[i] * fx[j] * data(sr, sc);
                }
            }
            result(r, c) = acc;
        }
    }

    // 垂直方向(沿行)边界补偿.
    auto cumsum = [](const std::vector<double>& v, bool rev) {
        std::vector<double> cs(v.size());
        double s = 0;
        if (!rev) {
            for (size_t i = 0; i < v.size(); ++i) { s += v[i]; cs[i] = s; }
        } else {
            // 等价于 MATLAB 的 cumsum(v(end:-1:1)):cs[j] = 末尾 j+1 个元素之和.
            for (size_t i = v.size(); i-- > 0;) { s += v[i]; cs[v.size() - 1 - i] = s; }
        }
        return cs;
    };
    const std::vector<double> fc1a = cumsum(fy, false), fc1b = cumsum(fy, true);
    const std::vector<double> fc2a = cumsum(fx, false), fc2b = cumsum(fx, true);
    const double fs1 = fc1a.back(), fs2 = fc2a.back();
    const int fl1b = (fl1 - 1) / 2, fl2b = (fl2 - 1) / 2;

    if (fl1 > sd1) {
        // tmp = [zeros(1,fl1) fc1a fs1*ones(1,fl1)],1-based;ff1(row)=tmp(fl1+fl1b+row)-tmp(fl1+fl1b+row-sd1).
        auto tmpAt = [&](int t) {  // t 为 1-based.
            if (t <= fl1) return 0.0;
            if (t <= 2 * fl1) return fc1a[t - fl1 - 1];
            return fs1;
        };
        for (int r = 0; r < sd1; ++r) {
            int range = fl1 + fl1b + (r + 1);
            double ff1 = tmpAt(range) - tmpAt(range - sd1);
            double scale = fs1 / ff1;
            for (int c = 0; c < sd2; ++c) result(r, c) *= scale;
        }
    } else {
        for (int k = 1; k <= fl1b; ++k) {  // k 为 1-based 行序.
            double st = fs1 / fc1a[fl1b + k - 1];
            double sb = fs1 / fc1b[fl1b + k - 1];
            for (int c = 0; c < sd2; ++c) {
                result(k - 1, c) *= st;
                result(sd1 - k, c) *= sb;
            }
        }
    }

    if (fl2 > sd2) {
        auto tmpAt = [&](int t) {
            if (t <= fl2) return 0.0;
            if (t <= 2 * fl2) return fc2a[t - fl2 - 1];
            return fs2;
        };
        for (int c = 0; c < sd2; ++c) {
            int range = fl2 + fl2b + (c + 1);
            double ff2 = tmpAt(range) - tmpAt(range - sd2);
            double scale = fs2 / ff2;
            for (int r = 0; r < sd1; ++r) result(r, c) *= scale;
        }
    } else {
        for (int k = 1; k <= fl2b; ++k) {
            double sl = fs2 / fc2a[fl2b + k - 1];
            double srr = fs2 / fc2b[fl2b + k - 1];
            for (int r = 0; r < sd1; ++r) {
                result(r, k - 1) *= sl;
                result(r, sd2 - k) *= srr;
            }
        }
    }
    return result;
}

// 2d 卷积 + 边界能量补偿.复刻 conv2PreserveEnergy.m,包括其 bug:
// 左/右/角点的补偿分支条件恒为 false,实际只补偿上/下边界(且不含左右 fw2 列).
Map conv2PreserveEnergy(const Map& src, const Map& f) {
    const int sh = src.rows, sw = src.cols;
    const int fh = f.rows, fw = f.cols;
    const int fh2 = (fh - 1) / 2, fw2 = (fw - 1) / 2;
    // conv2(src,f,'full') 再裁 same 区域:result(r,c) = sum f(i,j)*src(r+fh2-i, c+fw2-j).
    Map result(sh, sw, 0.0);
    for (int r = 0; r < sh; ++r) {
        for (int c = 0; c < sw; ++c) {
            double acc = 0.0;
            for (int i = 0; i < fh; ++i) {
                int sr = r + fh2 - i;
                if (sr < 0 || sr >= sh) continue;
                for (int j = 0; j < fw; ++j) {
                    int sc = c + fw2 - j;
                    if (sc < 0 || sc >= sw) continue;
                    acc += f(i, j) * src(sr, sc);
                }
            }
            result(r, c) = acc;
        }
    }

    double fsum = 0.0;
    for (double v : f.d) fsum += v;

    if (fh2 > 0 && fh2 <= sh) {
        // 上边界第 k 行(1-based k=1..fh2):用到滤波器前 fh2+k 行(1-based)的权和.
        // 下边界倒数第 k 行:用到滤波器后 fh2+k 行的权和.
        std::vector<double> rowSum(fh, 0.0);
        for (int i = 0; i < fh; ++i)
            for (int j = 0; j < fw; ++j) rowSum[i] += f(i, j);
        std::vector<double> pre(fh + 1, 0.0);  // 前缀和.
        for (int i = 0; i < fh; ++i) pre[i + 1] = pre[i] + rowSum[i];
        for (int k = 1; k <= fh2; ++k) {
            double sT = fsum / pre[fh2 + k];  // 前 fh2+k 行(1-based).
            // 下边界倒数第 k 行用到的是最后 fh2+k 行.
            double sBottom = fsum / (pre[fh] - pre[fh - (fh2 + k)]);
            for (int c = fw2; c <= sw - 1 - fw2; ++c) {
                result(k - 1, c) *= sT;
                result(sh - k, c) *= sBottom;
            }
        }
    }
    // 左/右/角点:MATLAB 原代码条件恒 false,不补偿(零填充),复刻此行为.
    return result;
}

// 复刻 attenuateBorders.m:borderSize 先按行数截断再按列数截断,系数 k/(bs+1).
void attenuateBordersInPlace(Map& m, int borderSize) {
    if (borderSize * 2 > m.rows) borderSize = m.rows / 2;
    if (borderSize * 2 > m.cols) borderSize = m.cols / 2;
    if (borderSize < 1) return;
    for (int k = 1; k <= borderSize; ++k) {
        double coef = (double)k / (borderSize + 1);
        for (int c = 0; c < m.cols; ++c) {
            m(k - 1, c) *= coef;
            m(m.rows - k, c) *= coef;
        }
    }
    for (int k = 1; k <= borderSize; ++k) {
        double coef = (double)k / (borderSize + 1);
        for (int r = 0; r < m.rows; ++r) {
            m(r, k - 1) *= coef;
            m(r, m.cols - k) *= coef;
        }
    }
}

// 复刻 gaussian.m:hw = floor(sigma*sqrt(-2*ln(0.01))),被 maxhw(>0 时)截断;peak=0 时取 1/(sigma*sqrt(2pi)).
std::vector<double> gaussian1d(double peak, double sigma, int maxhw) {
    int hw = (int)std::floor(sigma * std::sqrt(-2.0 * std::log(0.01)));
    if (maxhw > 0 && hw > maxhw) hw = maxhw;
    if (peak == 0) peak = 1.0 / (sigma * std::sqrt(2 * PI));
    std::vector<double> k(2 * hw + 1);
    k[hw] = peak;
    for (int i = 1; i <= hw; ++i) {
        double v = peak * std::exp(-(double)(i * i) / (2 * sigma * sigma));
        k[hw - i] = v;
        k[hw + i] = v;
    }
    return k;
}

// 复刻 normalizeImage.m:[0,0] 不处理;mx==mn 时平移到 0.5*(lo+hi).
Map normalizeImage(Map m, double lo, double hi) {
    if (lo == 0 && hi == 0) return m;
    double mn = m.d[0], mx = m.d[0];
    for (double v : m.d) { mn = std::min(mn, v); mx = std::max(mx, v); }
    if (mx == mn) {
        double base = -mx + 0.5 * (lo + hi);
        for (double& v : m.d) v += base;
    } else {
        double scale = std::fabs(hi - lo) / (mx - mn);
        double off = std::min(lo, hi);
        for (double& v : m.d) v = (v - mn) * scale + off;
    }
    return m;
}

static void clampLowInPlace(Map& m, double lo) {
    for (double& v : m.d) if (v < lo) v = lo;
}

// 复刻 maxNormalizeIterative.m(默认参数已内联:iterInhi=2, coEx=0.5, coIn=1.5, exSig=2%, inSig=25%).
Map maxNormalizeIterative(const Map& data, int numIter, double lo, double hi) {
    Map result = data;
    clampLowInPlace(result, 0);
    result = normalizeImage(result, lo, hi);

    int sz = std::max(result.rows, result.cols);
    int maxhw = std::max(0, std::min(result.rows, result.cols) / 2 - 1);
    double esig = sz * 2.0 * 0.01;
    double isig = sz * 25.0 * 0.01;
    std::vector<double> gExc = gaussian1d(0.5 / (esig * std::sqrt(2 * PI)), esig, maxhw);
    std::vector<double> gInh = gaussian1d(1.5 / (isig * std::sqrt(2 * PI)), isig, maxhw);

    for (int it = 0; it < numIter; ++it) {
        Map excit = sepConv2PreserveEnergy(gExc, gExc, result);
        Map inhib = sepConv2PreserveEnergy(gInh, gInh, result);
        double mx = 0.0;
        for (double v : result.d) mx = std::max(mx, v);
        double globinhi = 0.01 * 2.0 * mx;
        for (size_t i = 0; i < result.d.size(); ++i) {
            double v = result.d[i] + excit.d[i] - inhib.d[i] - globinhi;
            result.d[i] = v < 0 ? 0 : v;
        }
    }
    return result;
}

// 复刻 gaussianSubsample.m:核 [1,5,10,10,5,1,1]/33,0-based 从 0 起隔点抽取;小尺寸走特判分支.
Map gaussianSubsample(const Map& img) {
    static const std::vector<double> filt = [] {
        std::vector<double> f = {1, 5, 10, 10, 5, 1, 1};
        double s = 0;
        for (double v : f) s += v;
        for (double& v : f) v /= s;
        return f;
    }();
    const int h = img.rows, w = img.cols;

    if (std::min(h, w) > 3) {
        Map conv = sepConv2PreserveEnergy(filt, filt, img);
        Map out((h + 1) / 2, (w + 1) / 2);
        for (int r = 0; r < out.rows; ++r)
            for (int c = 0; c < out.cols; ++c)
                out(r, c) = conv(2 * r, 2 * c);
        return out;
    }

    // 小图:水平、垂直分开处理.
    Map hor;
    if (w <= 1) {
        hor = img;
    } else if (w == 2) {  // 核 [1,1]:两列平均.
        hor = Map(h, 1);
        for (int r = 0; r < h; ++r) hor(r, 0) = 0.5 * (img(r, 0) + img(r, 1));
    } else if (w == 3) {  // 核 [1,2,1]:mean(img(:,[1,2,2,3])).
        hor = Map(h, 1);
        for (int r = 0; r < h; ++r) hor(r, 0) = 0.25 * (img(r, 0) + 2 * img(r, 1) + img(r, 2));
    } else {
        static const std::vector<double> one = {1.0};
        Map conv = sepConv2PreserveEnergy(one, filt, img);
        int outW = (w - 1 + 1) / 2;  // 1:2:end-1 (1-based) -> 0-based 0,2,...,w-2.
        hor = Map(h, outW);
        for (int r = 0; r < h; ++r)
            for (int c = 0; c < outW; ++c)
                hor(r, c) = conv(r, 2 * c);
    }

    if (h <= 1) return hor;
    if (h == 2) {
        Map out(1, hor.cols);
        for (int c = 0; c < hor.cols; ++c) out(0, c) = 0.5 * (hor(0, c) + hor(1, c));
        return out;
    }
    if (h == 3) {
        Map out(1, hor.cols);
        for (int c = 0; c < hor.cols; ++c) out(0, c) = 0.25 * (hor(0, c) + 2 * hor(1, c) + hor(2, c));
        return out;
    }
    static const std::vector<double> one = {1.0};
    Map conv = sepConv2PreserveEnergy(filt, one, hor);
    int outH = (h - 1 + 1) / 2;
    Map out(outH, hor.cols);
    for (int r = 0; r < outH; ++r)
        for (int c = 0; c < out.cols; ++c)
            out(r, c) = conv(2 * r, c);
    return out;
}

// 复刻 makeDyadicPyramid.m:建到 min(h,w)<=1 为止.levels[i] 对应 1-based 层号 i+1.
std::vector<Map> makeDyadicPyramid(const Map& m) {
    std::vector<Map> pyr;
    pyr.push_back(m);
    while (std::min(pyr.back().rows, pyr.back().cols) > 1)
        pyr.push_back(gaussianSubsample(pyr.back()));
    return pyr;
}

// ---- 特征金字塔 ----

static Map toGray(const Image& img) {
    Map g(img.h, img.w);
    for (int r = 0; r < img.h; ++r)
        for (int c = 0; c < img.w; ++c)
            g(r, c) = 0.2989 * img.r(r, c) + 0.5870 * img.g(r, c) + 0.1140 * img.b(r, c);
    return g;
}

// 复刻 getRGB.m:lumThresh=0.1,in = max(r,g,b),低于阈值处三通道清零.
static void getRGB(const Image& img, Map& r, Map& g, Map& b, Map& in) {
    r = img.r; g = img.g; b = img.b;
    in = Map(img.h, img.w);
    for (int i = 0; i < img.h * img.w; ++i) {
        double mx = std::max({r.d[i], g.d[i], b.d[i]});
        in.d[i] = mx;
        if (mx < 0.1) { r.d[i] = 0; g.d[i] = 0; b.d[i] = 0; }
    }
}

static Map safeDivideMap(const Map& a, const Map& b) {
    Map out(a.rows, a.cols);
    for (size_t i = 0; i < a.d.size(); ++i)
        out.d[i] = (b.d[i] == 0) ? 0.0 : a.d[i] / b.d[i];
    return out;
}

// 复刻 makeGaborFilter.m(单个相位):9x9,去均值 + L2 归一化;si = -sin(angle).
static Map makeGaborFilter(const GaborParams& gp, double angleDeg, double phaseDeg) {
    double majorStd = gp.stddev, minorStd = majorStd * gp.elongation;
    int sz = gp.filterSize / 2;  // filterSize=9 -> sz=4.
    double rt = PI / 180 * angleDeg;
    double co = std::cos(rt), si = -std::sin(rt);
    double omega = 2 * PI / gp.filterPeriod;
    double majorSigq = 2 * majorStd * majorStd, minorSigq = 2 * minorStd * minorStd;
    double psi = PI / 180 * phaseDeg;

    int len = 2 * sz + 1;
    Map f(len, len);
    double mean = 0;
    for (int y = 0; y < len; ++y) {
        for (int x = 0; x < len; ++x) {
            double vx = x - sz, vy = y - sz;
            double major = vx * co + vy * si;   // 行主序下 x 为列方向.
            double minor = vy * si - vx * co;
            double env = std::exp(-major * major / majorSigq - minor * minor / minorSigq);
            double v = std::cos(omega * major + psi) * env;
            f(y, x) = v;
            mean += v;
        }
    }
    mean /= f.d.size();
    double nrm = 0;
    for (double& v : f.d) { v -= mean; nrm += v * v; }
    nrm = std::sqrt(nrm);
    for (double& v : f.d) v /= nrm;
    return f;
}

// 复刻 gaborFilterMap.m:两相位响应绝对值之和.
static Map gaborFilterMap(const Map& m, const GaborParams& gp, double angleDeg) {
    Map out(m.rows, m.cols, 0.0);
    for (double ph : gp.phases) {
        Map f = makeGaborFilter(gp, angleDeg, ph);
        Map c = conv2PreserveEnergy(m, f);
        for (size_t i = 0; i < out.d.size(); ++i) out.d[i] += std::fabs(c.d[i]);
    }
    return out;
}

// ---- center-surround ----

struct CSResult {
    std::vector<Map> featureMaps;
};

// 复刻 centerSurround.m:全部层最近邻缩放到 mapLevel 尺寸,|c-s| 差分 + attenuateBorders.
static CSResult centerSurround(const std::vector<Map>& pyr, const SalParams& params) {
    const LevelParams& lp = params.lp;
    int numLevels = (int)pyr.size();
    if (lp.mapLevel > numLevels)
        throw std::runtime_error("image too small: pyramid has no mapLevel");
    int mr = pyr[lp.mapLevel - 1].rows, mc = pyr[lp.mapLevel - 1].cols;

    // 预缩放所有要用到的层.
    std::vector<Map> maps(numLevels + 1);  // 1-based 访问.
    for (int l = lp.minLevel; l <= lp.maxLevel + lp.maxDelta && l <= numLevels; ++l)
        maps[l] = resizeNearest(pyr[l - 1], mr, mc);

    CSResult res;
    int borderSize = (int)std::lround(std::max(mr, mc) / 20.0);
    for (int l = lp.minLevel; l <= lp.maxLevel; ++l) {
        for (int d = lp.minDelta; d <= lp.maxDelta; ++d) {
            int l2 = l + d;
            if (l2 > numLevels) continue;
            Map fm(mr, mc);
            for (size_t i = 0; i < fm.d.size(); ++i)
                fm.d[i] = std::fabs(maps[l].d[i] - maps[l2].d[i]);
            attenuateBordersInPlace(fm, borderSize);
            res.featureMaps.push_back(std::move(fm));
        }
    }
    return res;
}

// 复刻 centerSurroundTopDown.m:各层只做 resize + attenuateBorders,无差分.
static CSResult centerSurroundTopDown(const std::vector<Map>& pyr, int targetRows, int targetCols) {
    CSResult res;
    int borderSize = (int)std::lround(std::max(targetRows, targetCols) / 20.0);
    for (const Map& lv : pyr) {
        Map fm = resizeNearest(lv, targetRows, targetCols);
        attenuateBordersInPlace(fm, borderSize);
        res.featureMaps.push_back(std::move(fm));
    }
    return res;
}

// ---- 显著图主流程,复刻 makeSaliencyMap.m ----

SaliencyResult makeSaliencyMap(const Image& img, SalParams& params,
                               const std::vector<Map>& topdownPyramid) {
    // 灰度图去除颜色特征(复刻 removeColorFeatures.m).
    if (!img.isColor) {
        std::vector<std::string> feats;
        std::vector<double> wts;
        for (size_t i = 0; i < params.features.size(); ++i) {
            if (params.features[i] == "Color") continue;
            feats.push_back(params.features[i]);
            wts.push_back(params.weights[i]);
        }
        params.features = feats;
        params.weights = wts;
    }
    // top-down 特征追加.
    if (!topdownPyramid.empty()) {
        params.features.push_back("TopDown");
        if (params.weights.size() < params.features.size())
            params.weights.push_back(params.topdownWeight);
    }

    const int numFeat = (int)params.features.size();
    SaliencyResult sr;
    sr.feats.resize(numFeat);

    Map gray;  // 延迟计算,供 Intensity 与 Orientation 复用.
    std::vector<Map> intPyr;

    for (int f = 0; f < numFeat; ++f) {
        const std::string& feat = params.features[f];
        FeatureData& fd = sr.feats[f];
        fd.label = feat;

        // 建金字塔.
        if (feat == "Intensities" || feat == "Intensity") {
            if (gray.empty()) gray = toGray(img);
            intPyr = makeDyadicPyramid(gray);
            fd.pyr.push_back({"Intensity", intPyr});
        } else if (feat == "Color") {
            Map r, g, b, in;
            getRGB(img, r, g, b, in);
            Map rg = safeDivideMap([&] { Map t(r.rows, r.cols); for (size_t i = 0; i < t.d.size(); ++i) t.d[i] = r.d[i] - g.d[i]; return t; }(), in);
            Map by = safeDivideMap([&] { Map t(r.rows, r.cols); for (size_t i = 0; i < t.d.size(); ++i) t.d[i] = b.d[i] - std::min(r.d[i], g.d[i]); return t; }(), in);
            fd.pyr.push_back({"Red/Green", makeDyadicPyramid(rg)});
            fd.pyr.push_back({"Blue/Yellow", makeDyadicPyramid(by)});
        } else if (feat == "Orientations" || feat == "Orientation") {
            if (intPyr.empty()) {
                if (gray.empty()) gray = toGray(img);
                intPyr = makeDyadicPyramid(gray);
            }
            // efficient 模式:只对 minLevel..min(maxLevel+maxDelta, 层数) 滤波.
            int maxL = std::min(params.lp.maxLevel + params.lp.maxDelta, (int)intPyr.size());
            for (double ang : params.oriAngles) {
                FeatureData::Pyramid p;
                p.label = "Gabor";
                p.levels.resize(intPyr.size());
                for (int l = params.lp.minLevel; l <= maxL; ++l)
                    p.levels[l - 1] = gaborFilterMap(intPyr[l - 1], params.gabor, ang);
                fd.pyr.push_back(std::move(p));
            }
        } else if (feat == "TopDown") {
            fd.pyr.push_back({"TopDown", topdownPyramid});
        } else {
            throw std::runtime_error("unknown feature: " + feat);
        }

        // 每个金字塔:center-surround -> 逐图归一化([0,10]) -> 跨尺度求和.
        int numPyr = (int)fd.pyr.size();
        std::vector<Map> combFM(numPyr);
        fd.FM.resize(numPyr);
        for (int p = 0; p < numPyr; ++p) {
            CSResult cs;
            if (feat == "TopDown") {
                // 与其他特征统一尺寸:用强度金字塔 mapLevel 的尺寸(对奇数尺寸比 MATLAB 更稳健).
                int mr = intPyr.empty() ? fd.pyr[p].levels[params.lp.mapLevel - 1].rows
                                        : intPyr[params.lp.mapLevel - 1].rows;
                int mc = intPyr.empty() ? fd.pyr[p].levels[params.lp.mapLevel - 1].cols
                                        : intPyr[params.lp.mapLevel - 1].cols;
                cs = centerSurroundTopDown(fd.pyr[p].levels, mr, mc);
            } else {
                cs = centerSurround(fd.pyr[p].levels, params);
            }
            for (Map& fm : cs.featureMaps) {
                fd.FM[p].push_back(maxNormalizeIterative(fm, params.numIter, 0, 10));
            }
            // combineMaps:逐像素求和.
            Map sum = fd.FM[p][0];
            for (size_t k = 1; k < fd.FM[p].size(); ++k)
                for (size_t i = 0; i < sum.d.size(); ++i) sum.d[i] += fd.FM[p][k].d[i];
            combFM[p] = std::move(sum);
        }

        // 子特征合图各自再归一化([0,0]:跳过初始归一化).
        for (Map& m : combFM) m = maxNormalizeIterative(m, params.numIter, 0, 0);

        // CM.
        if (numPyr == 1) {
            fd.CM = combFM[0];
        } else {
            Map sum = combFM[0];
            for (int p = 1; p < numPyr; ++p)
                for (size_t i = 0; i < sum.d.size(); ++i) sum.d[i] += combFM[p].d[i];
            fd.CM = maxNormalizeIterative(sum, params.numIter, 0, 0);
        }

        // 加权:weights(f)/numPyr/numFeat.
        double w = params.weights[f] / numPyr / numFeat;
        for (double& v : fd.CM.d) v *= w;
    }

    // 总显著图:求和 + 归一化到 [0,2].
    Map salmap = sr.feats[0].CM;
    for (int f = 1; f < numFeat; ++f)
        for (size_t i = 0; i < salmap.d.size(); ++i) salmap.d[i] += sr.feats[f].CM.d[i];
    sr.salmap = maxNormalizeIterative(salmap, params.numIter, 0, 2);
    return sr;
}
