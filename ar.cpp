#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <functional>

using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution <>floatDist(-1, 0);

double K(int i, int j, int num_kernel) {
    if (num_kernel == 0) {
        return pow(i, 0.98) + pow(j, 0.98);
    }
    if (num_kernel == 1) {
        return i * j;
    }
    if (num_kernel == 2) {
        return 1;
    }
    if (num_kernel == 3) {
        return 2;
    }
    if (num_kernel == 4) {
        return i + j;
    }
    if (num_kernel == 5) {
        return pow(double(i) / j, 0.98) + pow(double(j) / i, 0.98);
    }
    if (num_kernel == 6) {
        return pow(double(i) / j, 0.95) + pow(double(j) / i, 0.95);
    }
}


int main(int argc, char *argv[]) // num_kernel, rank, max_time, N_0, lamda, log_period
{
    int num_kernel = std::stoi(argv[1]);
    double cur_time = 0.0,
            max_time = std::atof(argv[3]);
    long long N_0 = std::atof(argv[4]);
    double lamda = std::atof(argv[5]);
    long long period = std::stoi(argv[6]);
    unsigned int start_time = clock(); // системное время
    
    double V = N_0;
    long long N = N_0;
    vector <int> size(N_0, 1); // размер i частицы - size[i]
    double K_max = K(1, 1, num_kernel);
    double max_size = 1;
    double second_moment = N_0; 
    int counter = 0;
    unsigned long long count_coag = 0;

    
    //////////ЛОГИ/////
    ofstream second_moment_out("second_moments_ar.txt");
    ofstream time_out("time_ar.txt");
    ofstream concentrations_end_out("concentrations_end_ar.txt");
    ofstream max_sizes_out("max_sizes_ar.txt");
    ofstream concentration_time_out("concentration_time_ar.txt");
    ofstream special_information_out("special_information_ar.txt");
    ///////////////////


    while (cur_time < max_time) {
        cur_time += ((2  *  V) / (N * (N - 1) * K_max)) * (1 / (1 + lamda));
        // ВЫБОР ЧАСТИЦ 
        int k = -floatDist(gen) * N;
        int l = -floatDist(gen) * N;
        while (k == l) {
            l = int(-floatDist(gen) * N);
        }
        int i = size[k];
        int j = size[l];
        if (i > j) {
            swap(i, j);
            swap(k, l);
        }


        // ПРИНЯТИЕ - ОТКЛОНЕНИЕ
        double random_val = -floatDist(gen);
        if (K_max * random_val > K(i, j, num_kernel)) {
            continue;
        }

        count_coag++;
        // ФРАГМЕНТАЦИЯ
        if (-floatDist(gen) < lamda / (lamda + 1)) {
            second_moment += (i + j) -  i * i -  j * j;
            size[k] = size[N - 1];
            size[l] = size[N - 2];
            size[N - 1] = 1;
            size[N - 2] = 1;
            N += i + j - 2;
            // ОБНОВЛЕНИЕ МАКСИМУМА
            if (i == max_size || j == max_size) {
                max_size = 1;
                for (int w = 0 ; w < N; w++) {
                    max_size = std::max(max_size, double(size[w]));
                }
            }
            size.resize(N, 1); 
        } else {
            // СЛИПАНИЕ
            max_size = (max_size > i + j) ? max_size : i + j; 
            second_moment += (i + j) * (i + j) - i * i - j * j;
            size[k] += size[l];
            size[l] = size[N - 1];
            size[N - 1] = 1;
            N -= 1;
            
            // УДВОЕНИЕ СИСТЕМЫ
            if (N < N_0 / 2) {
                size.resize(2 * N);
                for (int w = 0; w < N; w++) {
                    size[w + N] = size[w];
                }
                N *= 2;
                V *= 2;
                second_moment *= 2;
            }
        }
        K_max = K(max_size, max_size, 1);
        

        // СБОР ЛОГОВ
        counter += 1;
        counter %= period;
        if (counter == 0) {
            max_sizes_out << max_size / V << " ";
            second_moment_out << second_moment / V << " ";
            time_out << cur_time << " ";
            concentration_time_out << double(N) / V << " ";
            cout << "Масса: " << V << " N: " << N <<  " Концентрация: " << double(N) / V  << " Второй момент: " << second_moment / V << " max_size: " << max_size << " Время: " << cur_time << endl;    
        }
    }


    // дозапись логов
    vector <int> N_store_(max_size + 1);
    for (int w = 1 ; w < N; w++) {
	    N_store_[size[w]] += 1;
    }
    for (int w = 1 ; w < max_size; w++) {
        concentrations_end_out << double(N_store_[w]) / V << " ";
    }

    max_sizes_out.close();
    concentration_time_out.close();
    concentrations_end_out.close();
    second_moment_out.close();
    time_out.close();

    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    special_information_out << double(search_time)/CLOCKS_PER_SEC << " ";
    cout << endl << "ВРЕМЯ РАБОТЫ: " << double(search_time)/CLOCKS_PER_SEC << endl;
    special_information_out.close();
    cout << "КОЛИЧЕСТВО СТОЛКНОВЕНИЙ: " << count_coag << endl;
    return 0;
}

