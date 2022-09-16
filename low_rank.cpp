#include <iostream>
#include <cmath>
#include <vector> 
#include <random>
#include <fstream>
#include <functional>
#include <ctime>


using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution <>floatDist(-1, 0);


double K(int i, bool u_v, int num_rank, int num_kernel) {
    if (num_kernel == 0) { // i ^ 0.98 + j ^ 0.98
    	if (num_rank == 1)
            if (u_v == true) return pow(i, 0.98);
                else return 1;
    	if (num_rank == 2)
            if (u_v == true) return 1;
                else return pow(i, 0.98);
    }

    if (num_kernel == 1) // ij
        if (num_rank == 1) 
            if (u_v == true) return i;
                else return i;
    

    if (num_kernel == 2) // 1
    	if (num_rank == 1)
	        if (u_v == true) return 1;
	            else return 1;
    

    if (num_kernel == 3) // 2
    	if (num_rank == 1)
	        if (u_v == true) return 2;
                else return 1;


    
    if (num_kernel == 4) { // i + j
	    if (num_rank == 1)
		    if (u_v == true) return i;
		        else return 1.0;
		
	    if (num_rank == 2) 
		    if (u_v == true) return 1.0;
                else return i;    
    }


    if (num_kernel == 5) { // (i / j) ** 0.98 + (j / i) ** 0.98
        if (num_rank == 1)
            if (u_v == true) return pow(i, 0.98);
                else return pow(1.0 / double(i), 0.98);

        if (num_rank == 2)
            if (u_v == true) return pow(1.0 / double(i), 0.98);
                else return pow(i, 0.98);
    }


    if (num_kernel == 6) { // (i / j) ** 0.95 + (j / i) ** 0.95
        if (num_rank == 1) 
            if (u_v == true) return pow(i, 0.95);
                else return pow(1.0 / double(i), 0.95);

        if (num_rank == 2)
            if (u_v == true) return pow(1.0 / double(i), 0.95);
                else return pow(i, 0.95);
    }

}


void build (int v, int tl, int tr, vector <double> &t, vector <long long> &N_store, bool u_v, int num_rank, int num_kernel) {
	if (tl == tr) {
	    t[v] = (tl > 0) ? double(N_store[tl]) * K(tl, u_v, num_rank, num_kernel) : 0; 
    } else {
	    int tm = (tl + tr) / 2;
	    build (v*2, tl, tm, t, N_store, u_v, num_rank, num_kernel);
	    build (v*2+1, tm+1, tr, t, N_store, u_v, num_rank, num_kernel);
	    t[v] = t[v*2] + t[v*2+1];
	}
}


void update (int v, int tl, int tr, int pos, double new_val, vector <double> &t) {
	if (tl == tr)
		t[v] = new_val;
	else {
		int tm = (tl + tr) / 2;
		if (pos <= tm)
			update (v*2, tl, tm, pos, new_val, t);
		else
			update (v*2+1, tm+1, tr, pos, new_val, t);
		t[v] = t[v*2] + t[v*2+1];
	}
}


int find (vector <double> &t, int M){
    double random = (-floatDist(gen));
    double rand_val = random * t[1];
    int v = 1;
    while (v < M) {
        v *= 2;
        if (t[v] < rand_val) {
            rand_val -= t[v];
            v += 1;
        }
    }
    return v - M;
}


void updating_all_trees(int pos, int M, int rank, vector <long long> &N_store, vector < vector <double> > &all_trees_u, vector < vector <double> > &all_trees_v, int num_kernel) {
    for (int w = 0 ; w < rank ; w++) {
        update(1, 0, M - 1, pos, K(pos, true, w + 1, num_kernel) * N_store[pos] , all_trees_u[w]); 
        update(1, 0, M - 1, pos, K(pos, false, w + 1, num_kernel) * N_store[pos] , all_trees_v[w]);
    }
}


double special_number(vector < vector <double> > &all_trees_u, vector < vector <double> > &all_trees_v) {
    double sum = 0;
    for (int i = 0 ; i < all_trees_u.size(); i++) {
        sum += all_trees_u[i][1] * all_trees_v[i][1];
    }
    return sum;
}


void coagulation_effect(int i, int j, vector <long long> &N_store, long long &N, double &second_moment) {
    N_store[i] -= 1;
    N_store[j] -= 1;
    N_store[i + j] += 1;
    N -= 1;
    second_moment += (i + j) * (i + j) - i * i - j * j;
}


void shattering_effect(int i, int j, vector <long long> &N_store, long long &N, double &second_moment) {
    N_store[i] -= 1;
    N_store[j] -= 1;
    N_store[1] += i + j;
    N += i + j - 2;
    second_moment += (i + j) -  i * i -  j * j;
}


void doubling_system(long long &N, double &V, vector <long long> &N_store, vector < vector <double> > &all_trees_u, vector < vector <double> > &all_trees_v, double &second_moment) {
	N *= 2;
    V *= 2;
    second_moment *= 2;
    for (int w = 1 ; w < N_store.size(); w++) {
        N_store[w] *= 2;
    }
    for (int w = 0 ; w < all_trees_u.size(); w++) {
        for (int s = 0 ; s < all_trees_u[w].size() ; s++) {
            all_trees_u[w][s] *= 2;
            all_trees_v[w][s] *= 2;
        }
    }
}


void update_max_size(long long &max_size, vector <long long> &N_store) {
    max_size = N_store.size() - 1;
    while (N_store[max_size] == 0) {
        max_size--;
    }
}


int main(int argc, char* argv[]) // num_kernel, rank, max_time, N_0, lamda, period
{
    // ИНИЦИАЛИЗАЦИЯ ПАРАМЕТРОВ МОДЕЛИРУЕМОЙ СИСТЕМЫ
    int num_kernel = std::stoi(argv[1]); // номер ядра в фукции: C()
    int rank = std::stoi(argv[2]); // ранг исследумого ядра
    double cur_time = 0.0, 
            max_time = std::atof(argv[3]); // время моделирования
    long long N_0 = std::atof(argv[4]); // первоначальное число частиц
    double lamda = std::atof(argv[5]); // параметр фрагментации
    long long period = std::stoi(argv[6]); // период сбора логов 
    
    long long  M = 2; // граница размера частиц хранимая текущими структурами данных // debug значение 
    double V = N_0; // объём системы
    long long N = N_0; // текушее число частиц
    double second_moment = N_0; // второй момент системы
    long long max_size = 1; // максимальный размер частицы в текущий момент
    int counter = 0; // счётчик количества столкновений
    unsigned int start_time = clock(); // системное время    
    double random;
    unsigned long long count_coag = 0;


    // ИНИЦИАЛИЗАЦИЯ СТРУКТУР ДАННЫХ
    vector <long long> N_store(M, 0); // N_store[i] - количество частиц размера i
    N_store[1] = N_0;
    vector < vector <double> > all_trees_u(rank, vector <double> (2 * M)); // хранилище деревьев отрезков для векторов - u
    vector < vector <double> > all_trees_v(rank, vector <double> (2 * M)); // для векторв - v
    for (int w = 0 ; w < rank ; w++) {
        build(1, 0, M - 1, all_trees_u[w], N_store, true, w + 1, num_kernel);
        build(1, 0, M - 1, all_trees_v[w], N_store, false, w + 1, num_kernel);
    }

    
    //////////ЛОГИ/////
    ofstream second_moment_out("second_moments_lr.txt");
    ofstream time_out("time_lr.txt");
    ofstream concentrations_end_out("concentrations_end_lr.txt");
    ofstream max_sizes_out("max_sizes_lr.txt");
    ofstream concentration_time_out("concentration_time_lr.txt");
    ofstream special_information_out("special_information_lr.txt");
    ///////////////////
   
 
    while (cur_time < max_time) {
        // УДВОЕНИЕ СИСТЕМЫ В СЛУЧАЕ ДЕГРАДАЦИИ  	
	    if (N <= N_0 / 2) {
            doubling_system(N, V, N_store, all_trees_u, all_trees_v, second_moment);
        }
	
        // ШАГ ПО ВРЕМЕНИ
        double total_activity = special_number(all_trees_u, all_trees_v);
        cur_time += ((2 * V) / total_activity) * (1 / (1 + lamda));
        // ВЫБОР КОМПАНЕНТЫ КАНОНИЧЕСКОГО РАЗЛОЖЕНИЯ
        int target_summand = 0;
        random = (-floatDist(gen));
        double random_val = random * total_activity;
        
        for (target_summand = 0; target_summand < rank - 1; target_summand++) {
            random_val -= all_trees_u[target_summand][1] * all_trees_v[target_summand][1];
            if (random_val <= 0) {
                break;
            }
        }
        
        // ВЫБОР РАЗМЕРОВ ЧАСТИЦ ДЛЯ ВЗАИМОДЕЙСТВИЯ
        int i = find(all_trees_u[target_summand], M);
        int j = find(all_trees_v[target_summand], M);
       
        // ПРОВЕРКА СУЩЕСТВОВАНИЯ ЧАСТИЦ  
        if (N_store[i] == 0 || N_store[j] == 0) {
            continue;
        }

        // ПРОВЕРКА СТОЛКНОВЕНИЯ С САМИМ СОБОЙ
        random = (-floatDist(gen));
        if ((i == j) && (random * N_store[i] <= 1)) {
            continue; 
        }
        // СБОР ЛОГОВ 
        counter++;
        counter %= period;
        if (counter == 0) {
            max_sizes_out << max_size / V << " ";
            second_moment_out << second_moment / V << " ";
            time_out << cur_time << " ";
            concentration_time_out << double(N) / V << " ";
	        cout << "Масса: " << V << " N: " << N <<  " Концентрация: " << double(N) / V  << " Второй момент: " << second_moment / V << " max_size: " << max_size << " Время: " << cur_time << endl;
	    }
        
        count_coag++;    
     
        // СЛУЧАЙ ФРАГМЕНТАЦИИ 
        random = (-floatDist(gen));
        if (random < lamda / (1 + lamda)) {
            shattering_effect(i, j, N_store, N, second_moment);
	        if (i == max_size || j == max_size) {
		        update_max_size(max_size, N_store);
	        }

            updating_all_trees(i, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);
            updating_all_trees(j, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);
            updating_all_trees(1, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);            
            continue;
        }
        
        max_size = (max_size > i + j) ? max_size : i + j;

        // РАСШИРЕНИЕ СТРУКТУР ДАННЫХ
        if (i + j >= M) {
            M *= 2;
            N_store.resize(M, 0);
            coagulation_effect(i, j, N_store, N, second_moment);

            for (int w = 0 ; w < rank ; w++) {
                all_trees_u[w].resize(2 * M);
                all_trees_v[w].resize(2 * M);
        	    build(1, 0, M - 1, all_trees_u[w], N_store, true, w + 1, num_kernel);
        	    build(1, 0, M - 1, all_trees_v[w], N_store, false, w + 1, num_kernel);            
	        }
            continue;
        }
        coagulation_effect(i, j, N_store, N, second_moment);
        updating_all_trees(i, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);
        updating_all_trees(j, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);
        updating_all_trees(i + j, M, rank, N_store, all_trees_u, all_trees_v, num_kernel);
    }
   
 
    // ЗАПИСЬ ФИНАЛЬНЫХ КОНЦЕНТРАЦИЙ
    for (int w = 1 ; w < max_size; w++) {
	    concentrations_end_out << double(N_store[w]) / V << " ";
    }


    max_sizes_out.close();
    concentration_time_out.close();
    concentrations_end_out.close();
    second_moment_out.close();
    time_out.close();
     
    unsigned int end_time = clock();
    unsigned int search_time = end_time - start_time;
    special_information_out << double(search_time)/CLOCKS_PER_SEC << " "; 
    cout << endl <<"ВРЕМЯ РАБОТЫ: " <<  double(search_time)/CLOCKS_PER_SEC << endl;
    cout << "КОЛИЧЕСТВО СТОЛКНОВЕНИЙ: "<< count_coag << endl;
    special_information_out.close();

    return 0;
}
