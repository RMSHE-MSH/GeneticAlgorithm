#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <utility>
using namespace std;

constexpr auto PI = 3.14159265358979323846;

// 使用随机数引擎生成随机数
random_device rd;
mt19937 gen(rd());

typedef struct Genetic_Struct {
    vector<int> rna;            //个体遗传物质
    int protein;                //个体生物性状(蛋白质)
    double fitness;             //个体适应度
}VECGS;

class GeneticAlgorithm {
public:
    int rna_size = 16;                   // rna长度(碱基数)
    int population_size = 16;           // 种群规模
    int generation_total = 1000;        // 种群最大迭代次数
    double mutate_rate = 0.02;          // 碱基突变率
    double selection_rate = 0.5;          //自然选择中的幸存率

    vector<VECGS> population;           //定义一个种群

    vector<vector<VECGS>> generation;   //记录每一代的种群信息

private:
    //argsort函数，接受一个fitness值数组，并返回一个按fitness值排序后的索引数组
    vector<int> argsort(const vector<double> fitness) {
        // 将fitness值和其索引放入一个pair中
        vector<pair<int, double>> items;
        for (int i = 0; i < fitness.size(); ++i)
            items.emplace_back(i, fitness[i]);

        // 使用sort函数按fitness值从大到小排序
        sort(items.begin(), items.end(), [this](const pair<int, double> &a, const pair<int, double> &b)->bool {return a.second > b.second; });

        // 返回排序后的索引数组
        vector<int> result;
        for (const auto &item : items)
            result.push_back(item.first);

        return result;
    }

    //翻译个体遗传物质(将二进制遗传信息转为十进制)
    int translation_rna(vector<int> rna) {
        // 定义一个整型变量来存储十进制数
        int _protein = 0;
        // 循环遍历二进制数vector
        for (int i = 0; i < rna.size(); ++i) {
            // 根据位数计算十进制数
            _protein += rna[i] * pow(2, rna.size() - i - 1);
        }

        return _protein;
    }

    //计算个体适应度
    double calculate_fitness(int protein) { return target_function(protein); }

    //创建个体随机RNA片段
    vector<int> init_rna() {
        vector<int>rna;

        // 使用std::uniform_int_distribution生成介于0和1之间的随机整数
        uniform_int_distribution<> dis(0, 1);

        // 生成随机rna序列
        for (int i = 0; i < rna_size; ++i) rna.push_back(dis(gen));

        return rna;
    }

    //创建初始种群(初始化种群)
    vector<VECGS> init_population() {
        //将"population"的大小改为"population_size"指定的大小;
        population.resize(population_size);

        // 对于"population"数组中的每个元素，它会生成一个随机的RNA片段，翻译为蛋白质，并计算适应度
        for (auto &i : population) {
            i.rna = init_rna();//创建个体随机RNA片段
            i.protein = translation_rna(i.rna);//翻译个体的RNA为蛋白质
            i.fitness = calculate_fitness(i.protein);//计算个体的适应度
        }

        return population;
    }

    //自然选择
    vector<VECGS> select(vector<VECGS> population) {
        //按比例选择适应度较高的个体(适应度从大到小排序)
        vector<VECGS> population_ng;

        //创建一个向量来存放指向"population"中的所有"fitness"
        vector<double> fitness;
        for (auto &i : population) fitness.push_back(i.fitness);

        //返回一个按fitness值排序后的索引数组(大到小)
        vector<int> population_index = argsort(fitness);

        //使最适应的个体排到列表前面,然后取前"%selection_rate"最适应环境的个体
        int num_indices = population_index.size() * selection_rate;
        for (auto i = population_index.begin(); i < population_index.begin() + num_indices; ++i)
            population_ng.push_back(population[*i]);

        return population_ng;
    }

    //基因重组
    vector<VECGS> crossover(vector<VECGS> population) {
        vector<VECGS> population_ng;

        //预统计自然选择后种群中剩余的个体数下标[0,population.size() - 1]
        int selected_population_size = population.size() - 1;

        for (int i = 0; i < population_size; ++i) {
            //生成随机数列表;
            uniform_int_distribution<int> dis(0, selected_population_size);

            //如果种群中的个体数大于等于2, 那么就会进行自由交配; 如果种群中仅剩下一个个体, 那么个体将会进行无性繁殖;
            if (population.size() >= 2) {
                // 随机选择两个不同的个体
                int j, k; do {
                    j = dis(gen);
                    k = dis(gen);
                } while (j == k);

                // 随机选择两个个体的rna交换位置
                uniform_int_distribution<int> disk_rs(0, rna_size - 1);
                int rs = disk_rs(gen);

                // 交换两个个体的遗传信息
                vector<int> offspring;
                offspring.insert(offspring.end(), population[j].rna.begin(), population[j].rna.begin() + rs);
                offspring.insert(offspring.end(), population[k].rna.begin() + rs, population[k].rna.end());

                population_ng.push_back({ offspring ,NULL, NULL });
            } else {
                population.push_back(population[i]);
            }
        }

        return population_ng;
    }

    //基因突变
    vector<VECGS> mutate(vector<VECGS> population) {
        //生成随机数列表;
        uniform_int_distribution<int> dis(0, 100);

        // 遍历所有碱基,每个碱基都有概率会发生突变
        for (auto &vecgs : population) {
            for (auto &rna : vecgs.rna) {
                if (dis(gen) < (mutate_rate * 100)) rna = rna == 0 ? 1 : 0;
            }
        }

        return population;
    }

public:
    //目标函数
    double target_function(double input) {
        double x = input * (PI / 180);
        double output = -pow(x, 2) + 20 * sin(x) + 100;
        return output;
    }

    //进行演化
    vector<VECGS> iteration(bool Record_generation = false) {
        //创建初始种群
        init_population();

        //开始迭代
        for (int i = 0; i < generation_total; ++i) {
            //自然选择
            population = select(population);

            //基因重组
            population = crossover(population);

            //基因突变
            population = mutate(population);

            //翻译RNA和计算种群适应度
            for (auto &j : population) {
                j.protein = translation_rna(j.rna);//翻译个体的RNA为蛋白质
                j.fitness = calculate_fitness(j.protein);//计算个体的适应度
            }

            if (Record_generation == true) generation.push_back(population);
        }

        //返回最优解种群
        return population;
    }
}GA;

int main() {
    GeneticAlgorithm GA;

    vector<VECGS> best_population = GA.iteration();

    for (auto i = best_population.begin(); i < best_population.end(); ++i) {
        cout << (*i).fitness << "    \t\t" << (*i).protein * (PI / 180) << "    \t\t";
        for (auto j = (*i).rna.begin(); j < (*i).rna.end(); ++j)
            cout << (*j);

        cout << endl;
    }

    system("pause");

    return 0;
}