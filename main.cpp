#include <iostream>
#include <vector>
#include "algorithm.h"
#include "functions.h"
using namespace std;

int main(int argc, char *argv[]) {
    // 選擇目標函數
    //auto targetFunction = F2; // 只需更改這一行即可切換目標函數

    // 1. 預計的測試函數列表
    std::vector<std::function<double(const std::vector<double>&)>> funcs = {
        F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13
    };
    // 2. 對應的 G (gen) 數值
    std::vector<int> gens = {
        1500, 2000, 5000, 5000, 3000, 100, 3000, 1000, 1000, 500, 500,500,500
    };
    // 參數設置
    const int D = 30; // 維度
    const int NP = 100; // 種群大小
    const int G = 500; // 迭代次數
    const double pb = 0.05; //取前幾%的個體
    const double c = 0.1; //自適應參數
    
    //函數邊界
    const double minVal = -1.28;
    const double maxVal = 1.28;

    cout << "Initializing parameters:\n";
    cout << "D: " << D << ", NP: " << NP << ", G: " << G << ", p: " << pb << ", c: " << c << "\n";

    // 呼叫演算法
    //vector<double> best = differential_evolution(D, NP, G, pb, c, minVal, maxVal, targetFunction);
    
    //指定開頭和結尾函式
    int st = 1  ,ed = 13;
    min(int(funcs.size()),ed);
    // 3. 每個函數的邊界最大值 F1 ~ F13
    std::vector<double> bounds = {
        100, 10, 100, 100, 30, 100, 1.28, 500, 5.12, 32, 600, 50, 50
    };

    for (int i = st; i <= ed; ++i) {
        double bound = bounds[i-1];
        double minVal = -bound;
        double maxVal = bound;
        std::cout << "==== F" << (i) << " (Gen=" << gens[i-1] << endl;
        std::vector<double> best = differential_evolution(D, NP, gens[i-1], pb, c, minVal, maxVal, funcs[i-1]);
        std::cout << "Best solution: ";
        for (double val : best) std::cout << val << " ";
        std::cout << "\nFitness: " << funcs[i-1](best) << "\n";
        cout<<"F"<<i<<" result:"<<"\n";
        cout << "========================\n";
    }
    /*// 輸出最佳解
    cout << "Best solution: ";
    for (double val : best) cout << val << " ";
    cout << "\nFitness: " << targetFunction(best) << "\n";
    */
    system("pause");
    return 0;
}