#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <utility>
using namespace std;

constexpr auto PI = 3.14159265358979323846;

// ʹ��������������������
random_device rd;
mt19937 gen(rd());

typedef struct Genetic_Struct {
    vector<int> rna;            //�����Ŵ�����
    int protein;                //����������״(������)
    double fitness;             //������Ӧ��
}VECGS;

class GeneticAlgorithm {
public:
    int rna_size = 16;                   // rna����(�����)
    int population_size = 16;           // ��Ⱥ��ģ
    int generation_total = 1000;        // ��Ⱥ����������
    double mutate_rate = 0.02;          // ���ͻ����
    double selection_rate = 0.5;          //��Ȼѡ���е��Ҵ���

    vector<VECGS> population;           //����һ����Ⱥ

    vector<vector<VECGS>> generation;   //��¼ÿһ������Ⱥ��Ϣ

private:
    //argsort����������һ��fitnessֵ���飬������һ����fitnessֵ��������������
    vector<int> argsort(const vector<double> fitness) {
        // ��fitnessֵ������������һ��pair��
        vector<pair<int, double>> items;
        for (int i = 0; i < fitness.size(); ++i)
            items.emplace_back(i, fitness[i]);

        // ʹ��sort������fitnessֵ�Ӵ�С����
        sort(items.begin(), items.end(), [this](const pair<int, double> &a, const pair<int, double> &b)->bool {return a.second > b.second; });

        // ������������������
        vector<int> result;
        for (const auto &item : items)
            result.push_back(item.first);

        return result;
    }

    //��������Ŵ�����(���������Ŵ���ϢתΪʮ����)
    int translation_rna(vector<int> rna) {
        // ����һ�����ͱ������洢ʮ������
        int _protein = 0;
        // ѭ��������������vector
        for (int i = 0; i < rna.size(); ++i) {
            // ����λ������ʮ������
            _protein += rna[i] * pow(2, rna.size() - i - 1);
        }

        return _protein;
    }

    //���������Ӧ��
    double calculate_fitness(int protein) { return target_function(protein); }

    //�����������RNAƬ��
    vector<int> init_rna() {
        vector<int>rna;

        // ʹ��std::uniform_int_distribution���ɽ���0��1֮����������
        uniform_int_distribution<> dis(0, 1);

        // �������rna����
        for (int i = 0; i < rna_size; ++i) rna.push_back(dis(gen));

        return rna;
    }

    //������ʼ��Ⱥ(��ʼ����Ⱥ)
    vector<VECGS> init_population() {
        //��"population"�Ĵ�С��Ϊ"population_size"ָ���Ĵ�С;
        population.resize(population_size);

        // ����"population"�����е�ÿ��Ԫ�أ���������һ�������RNAƬ�Σ�����Ϊ�����ʣ���������Ӧ��
        for (auto &i : population) {
            i.rna = init_rna();//�����������RNAƬ��
            i.protein = translation_rna(i.rna);//��������RNAΪ������
            i.fitness = calculate_fitness(i.protein);//����������Ӧ��
        }

        return population;
    }

    //��Ȼѡ��
    vector<VECGS> select(vector<VECGS> population) {
        //������ѡ����Ӧ�Ƚϸߵĸ���(��Ӧ�ȴӴ�С����)
        vector<VECGS> population_ng;

        //����һ�����������ָ��"population"�е�����"fitness"
        vector<double> fitness;
        for (auto &i : population) fitness.push_back(i.fitness);

        //����һ����fitnessֵ��������������(��С)
        vector<int> population_index = argsort(fitness);

        //ʹ����Ӧ�ĸ����ŵ��б�ǰ��,Ȼ��ȡǰ"%selection_rate"����Ӧ�����ĸ���
        int num_indices = population_index.size() * selection_rate;
        for (auto i = population_index.begin(); i < population_index.begin() + num_indices; ++i)
            population_ng.push_back(population[*i]);

        return population_ng;
    }

    //��������
    vector<VECGS> crossover(vector<VECGS> population) {
        vector<VECGS> population_ng;

        //Ԥͳ����Ȼѡ�����Ⱥ��ʣ��ĸ������±�[0,population.size() - 1]
        int selected_population_size = population.size() - 1;

        for (int i = 0; i < population_size; ++i) {
            //����������б�;
            uniform_int_distribution<int> dis(0, selected_population_size);

            //�����Ⱥ�еĸ��������ڵ���2, ��ô�ͻ�������ɽ���; �����Ⱥ�н�ʣ��һ������, ��ô���彫��������Է�ֳ;
            if (population.size() >= 2) {
                // ���ѡ��������ͬ�ĸ���
                int j, k; do {
                    j = dis(gen);
                    k = dis(gen);
                } while (j == k);

                // ���ѡ�����������rna����λ��
                uniform_int_distribution<int> disk_rs(0, rna_size - 1);
                int rs = disk_rs(gen);

                // ��������������Ŵ���Ϣ
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

    //����ͻ��
    vector<VECGS> mutate(vector<VECGS> population) {
        //����������б�;
        uniform_int_distribution<int> dis(0, 100);

        // �������м��,ÿ��������и��ʻᷢ��ͻ��
        for (auto &vecgs : population) {
            for (auto &rna : vecgs.rna) {
                if (dis(gen) < (mutate_rate * 100)) rna = rna == 0 ? 1 : 0;
            }
        }

        return population;
    }

public:
    //Ŀ�꺯��
    double target_function(double input) {
        double x = input * (PI / 180);
        double output = -pow(x, 2) + 20 * sin(x) + 100;
        return output;
    }

    //�����ݻ�
    vector<VECGS> iteration(bool Record_generation = false) {
        //������ʼ��Ⱥ
        init_population();

        //��ʼ����
        for (int i = 0; i < generation_total; ++i) {
            //��Ȼѡ��
            population = select(population);

            //��������
            population = crossover(population);

            //����ͻ��
            population = mutate(population);

            //����RNA�ͼ�����Ⱥ��Ӧ��
            for (auto &j : population) {
                j.protein = translation_rna(j.rna);//��������RNAΪ������
                j.fitness = calculate_fitness(j.protein);//����������Ӧ��
            }

            if (Record_generation == true) generation.push_back(population);
        }

        //�������Ž���Ⱥ
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