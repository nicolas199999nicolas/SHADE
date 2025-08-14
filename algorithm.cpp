#include "algorithm.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <functional>
#include <cmath>
#include "functions.h"

using namespace std;

vector<double> differential_evolution(
    int D, int NP, int G, double pb, double c, double minVal, double maxVal,
    function<double(const vector<double>&)> f
) {
    // Initialization phase
    int H = NP; // 固定 H 為 NP
    vector<vector<double>> P(NP); // Population
    vector<double> MCR(H, 0.5), MF(H, 0.5); // 記錄交叉率與縮放因子
    vector<vector<double>> Archive; // 歷史存檔
    int k = 0; // Index counter
    std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    std::uniform_real_distribution<> rand01(0.0, 1.0);
    std::uniform_int_distribution<> randH(0, H-1);
    std::uniform_int_distribution<> randNP(0, NP-1);
    
    // 初始化族群
    vector<double> fitness(NP);
    for (int i = 0; i < NP; ++i) {
        P[i] = generateRandomIndividual(D, minVal, maxVal, gen);
        fitness[i] = f(P[i]);
    }
    
    // Main loop
    for (int genIdx = 0; genIdx < G; ++genIdx) {
        vector<double> CR(NP), F(NP);
        //成功的CR和F
        vector<double> S_CR, S_F, S_delta; // 新增 S_delta 用於 weighted mean
        vector<vector<double>> newP(NP);
        vector<double> newFitness(NP);
        cout << genIdx << ": "<< endl;
        // mutate & crossover
        for (int i = 0; i < NP; ++i) {
            int r = randH(gen); // 隨機選取歷史索引

            std::normal_distribution<> randCR(MCR[r], 0.1);
            std::cauchy_distribution<> randF(MF[r], 0.1);

            CR[i] = std::min(1.0, std::max(0.0, randCR(gen)));
            F[i] = std::min(1.0, std::max(0.0, randF(gen)));

            // p_i uniform [pb, 0.2]
            double pmin = 2.0 / NP;
            double pmax = 0.2;
            double p_i = pmin + (pmax - pmin) * rand01(gen);

            // p-best selection
            vector<int> sortedIdx(NP);
            iota(sortedIdx.begin(), sortedIdx.end(), 0);
            sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b){ return fitness[a] < fitness[b]; });//排序
            int num_p = max(2, static_cast<int>(NP * p_i)); // 至少 2 個
            int pBestIdx = sortedIdx[static_cast<int>(rand01(gen) * num_p)];//隨機選一個

            // 隨機選取 r1 from P excluding i
            int r1;
            do { r1 = randNP(gen); } while (r1 == i);

            // 選取 xr2 from P U A excluding i and r1
            vector<int> candidates;
            //把族群 P 中除了自己（i）和 r1 以外的所有個體索引都加入 candidates
            for (int c = 0; c < NP; ++c) {
                if (c != i && c != r1) candidates.push_back(c);
            }
            //把 Archive（歷史存檔）裡的所有個體也加入 candidates
            for (int a = 0; a < Archive.size(); ++a) {
                candidates.push_back(NP + a); // 偏移索引
            }
            // 隨機選取 xr2 from candidates
            int xr2_idx = candidates[static_cast<int>(rand01(gen) * candidates.size())];
            vector<double> xr2(D);
            if (xr2_idx < NP) {
                xr2 = P[xr2_idx];
            } else {
                xr2 = Archive[xr2_idx - NP];
            }

            // 變異
            vector<double> vi(D);
            for (int d = 0; d < D; ++d) {
                vi[d] = P[i][d] + F[i] * (P[pBestIdx][d] - P[i][d]) + F[i] * (P[r1][d] - xr2[d]);
                // 邊界處理：反射
                if (vi[d] < minVal) vi[d] = 2 * minVal - vi[d];
                if (vi[d] > maxVal) vi[d] = 2 * maxVal - vi[d];
            }
                                
            // 交叉
            vector<double> ui(D);
            int jrand = static_cast<int>(rand01(gen) * D);
            for (int d = 0; d < D; ++d) {
                if (rand01(gen) < CR[i] || d == jrand) ui[d] = vi[d];
                else ui[d] = P[i][d];
            }
            newP[i] = ui;
            newFitness[i] = f(ui);
        }
        // selection
        for (int i = 0; i < NP; ++i) {
            if (newFitness[i] <= fitness[i]) {
                if (newFitness[i] < fitness[i]) {
                    Archive.push_back(P[i]); // 加入 archive
                    S_delta.push_back(abs(fitness[i] - newFitness[i])); // 記錄 delta f
                    S_CR.push_back(CR[i]); // 正確記錄 CR
                    S_F.push_back(F[i]);   // 正確記錄 F
                }
                P[i] = newP[i];
                fitness[i] = newFitness[i];
            }
        }
        // 控制 archive 大小
        while (Archive.size() > NP) {
            int eraseIdx = static_cast<int>(rand01(gen) * Archive.size());
            Archive.erase(Archive.begin() + eraseIdx);
        }
        // 更新 MCR, MF (使用 weighted mean)
        if (!S_CR.empty() && !S_F.empty()) {
            double sum_delta = 0.0;
            double sum_delta_cr = 0.0;
            double sum_delta_f = 0.0;
            double sum_delta_f2 = 0.0;
            for (size_t s = 0; s < S_CR.size(); ++s) {
                double delta = S_delta[s];
                sum_delta += delta;
                sum_delta_cr += delta * S_CR[s];
                sum_delta_f += delta * S_F[s];
                sum_delta_f2 += delta * S_F[s] * S_F[s];
            }
            if (sum_delta > 0.0) {
                MCR[k] = sum_delta_cr / sum_delta; // weighted arithmetic mean for CR
                MF[k] = sum_delta_f2 / sum_delta_f; // weighted lehmer mean for F
            }
            k = (k + 1) % H;
        }
    }
    // 回傳最佳解
    int bestIdx = min_element(fitness.begin(), fitness.end()) - fitness.begin();
    return P[bestIdx];
}