// wta.cpp - winner-take-all 网络演化与返回抑制.
// 复刻 initializeWTA / evolveWTA / evolveLeakyIntFire / diskIOR / shapeIOR / winnerToImgCoords.
#include "saliency.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

WTA initializeWTA(const Map& salmap, const SalParams& params, Rng& rng) {
    WTA wta;
    const int n = salmap.rows * salmap.cols;

    // sm 层:被动积分层,永不发放,Einh=Eexc=0(抑制是把 V 拉向 0,不是超极化).
    wta.sm.C = 5e-8;
    wta.sm.Gleak = 1e-7;
    wta.sm.Einh = 0;
    wta.sm.Eexc = 0;
    wta.sm.DoesFire = false;
    wta.sm.GinhDecay = params.IORdecay;
    wta.sm.ginhIsMap = true;
    wta.sm.GinhMap = Map(salmap.rows, salmap.cols, 0.0);
    wta.sm.V = Map(salmap.rows, salmap.cols, 0.0);
    wta.sm.I = Map(salmap.rows, salmap.cols);
    for (int i = 0; i < n; ++i)
        wta.sm.I.d[i] = salmap.d[i] * params.smOutputRange
                      + params.noiseAmpl * rng.uniform() + params.noiseConst;

    // exc 层:默认值,初始 Ginh=1e-2 防止 t=0 虚假发放.
    wta.exc.I = Map(salmap.rows, salmap.cols, 0.0);
    wta.exc.V = Map(salmap.rows, salmap.cols, 0.0);
    wta.exc.Ginh = 1e-2;

    // inhib:单个抑制性中间神经元(1x1).
    wta.inhib.V = Map(1, 1, 0.0);
    wta.inhib.I = Map(1, 1, 0.0);
    return wta;
}

// 复刻 evolveLeakyIntFire.m:积分 -> 钳位 -> Ginh 衰减 -> 判发放 -> 复位.
// 返回是否有神经元发放;spikes(可选)按行主序标记.
static bool evolveLIF(LIF& L, double t, std::vector<char>* spikes = nullptr) {
    double dt = t - L.time;
    bool any = false;
    const int n = (int)L.V.d.size();
    if (spikes) spikes->assign(n, 0);
    for (int i = 0; i < n; ++i) {
        double ginh = L.ginhIsMap ? L.GinhMap.d[i] : L.Ginh;
        double v = L.V.d[i] + dt / L.C * (L.I.d[i] - L.Gleak * (L.V.d[i] - L.Eleak)
                                        - L.Gexc * (L.V.d[i] - L.Eexc)
                                        - ginh * (L.V.d[i] - L.Einh));
        if (v < L.Einh) v = L.Einh;
        L.V.d[i] = v;
    }
    if (L.ginhIsMap) {
        for (double& g : L.GinhMap.d) g *= L.GinhDecay;
    } else {
        L.Ginh *= L.GinhDecay;
    }
    for (int i = 0; i < n; ++i) {
        bool s = (L.V.d[i] > L.Vthresh) && L.DoesFire;
        if (s) {
            L.V.d[i] = 0;
            any = true;
            if (spikes) (*spikes)[i] = 1;
        }
    }
    L.time = t;
    return any;
}

// 复刻 evolveWTA.m.胜者返回 1-based (row,col);无胜者返回 (-1,-1).
// 平局按 MATLAB find 的列主序取第一个.
std::pair<int, int> evolveWTA(WTA& wta) {
    double time = wta.exc.time + wta.exc.timeStep;
    std::pair<int, int> winner{-1, -1};

    evolveLIF(wta.sm, time);

    for (size_t i = 0; i < wta.exc.I.d.size(); ++i)
        wta.exc.I.d[i] = wta.sm.V.d[i] * wta.exc.Ginput;

    std::vector<char> spikes;
    evolveLIF(wta.exc, time, &spikes);

    // 全局抑制只持续一步.
    wta.exc.Ginh = 0;

    bool anySpike = false;
    for (char s : spikes) anySpike |= (s != 0);
    if (anySpike) {
        // 列主序第一个.
        const int rows = wta.exc.V.rows, cols = wta.exc.V.cols;
        for (int c = 0; c < cols && winner.first == -1; ++c)
            for (int r = 0; r < rows; ++r)
                if (spikes[(size_t)r * cols + c]) { winner = {r + 1, c + 1}; break; }
        wta.inhib.Gexc = wta.inhib.Gleak * 10;
    }

    if (evolveLIF(wta.inhib, time)) {
        wta.exc.Ginh = 1e-2;
        wta.inhib.Gexc = 0;
    }
    return winner;
}

// 复刻 diskIOR.m.注意 psdev 用的是 2^mapLevel(不是 2^(mapLevel-1)).
void diskIOR(WTA& wta, int winR, int winC, const SalParams& params) {
    const int rows = wta.sm.V.rows, cols = wta.sm.V.cols;
    double pampl = 0.1 * wta.sm.V(winR - 1, winC - 1);
    double mampl = 1e-4 * pampl;
    double psdev = 0.3 * params.foaSize / std::pow(2.0, params.lp.mapLevel);
    double msdev = 4.0 * psdev;
    double p2 = psdev * psdev, m2 = msdev * msdev;
    for (int r = 0; r < rows; ++r) {
        double dy = (r + 1) - winR;
        for (int c = 0; c < cols; ++c) {
            double dx = (c + 1) - winC;
            double d = dx * dx + dy * dy;
            wta.sm.GinhMap(r, c) += pampl * std::exp(-0.5 * d / p2)
                                  - mampl * std::exp(-0.5 * d / m2);
        }
    }
}

// 复刻 shapeIOR.m:mask 全零回退 diskIOR;mask 最近邻缩放到 sm 尺寸后平坦注入.
void shapeIOR(WTA& wta, int winR, int winC, const SalParams& params, const ShapeData& sd) {
    double mx = 0;
    for (double v : sd.binaryMap.d) mx = std::max(mx, v);
    if (mx == 0) {
        diskIOR(wta, winR, winC, params);
        return;
    }
    double ampl = 0.1 * wta.sm.V(winR - 1, winC - 1);
    Map bin = sd.iorMask;
    if (bin.rows != wta.sm.V.rows || bin.cols != wta.sm.V.cols)
        bin = resizeNearest(bin, wta.sm.V.rows, wta.sm.V.cols);
    for (size_t i = 0; i < wta.sm.GinhMap.d.size(); ++i)
        wta.sm.GinhMap.d[i] += ampl * bin.d[i];
}

// 复刻 applyIOR.m.
void applyIOR(WTA& wta, int winR, int winC, const SalParams& params, const ShapeData& sd) {
    switch (params.iorType) {
        case IORType::None:
            break;
        case IORType::Disk:
            diskIOR(wta, winR, winC, params);
            break;
        case IORType::Shape:
            if (!sd.valid) diskIOR(wta, winR, winC, params);
            else shapeIOR(wta, winR, winC, params, sd);
            break;
    }
}

// 复刻 winnerToImgCoords.m(仅 dyadic):useRandom 时先加 U(0,1) jitter,再 round((w-1)*2^(mapLevel-1)+1).
std::pair<int, int> winnerToImgCoords(int winR, int winC, const SalParams& params, Rng& rng) {
    double r = winR, c = winC;
    if (params.useRandom) {
        r += rng.uniform();
        c += rng.uniform();
    }
    double scale = std::pow(2.0, params.lp.mapLevel - 1);
    return {(int)std::lround((r - 1) * scale + 1), (int)std::lround((c - 1) * scale + 1)};
}
