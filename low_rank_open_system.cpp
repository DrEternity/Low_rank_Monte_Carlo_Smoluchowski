#include <iostream>
#include <cmath>
#include <vector> 
#include <random>
#include <fstream>
#include <functional>
#include <ctime>
#include <map>
#include <cmath>

using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution <>floatDist(-1, 0);


class Modeling
{
    double K(long long i, bool u_v, int num_rank) {
    
        if (num_rank == 1) 
            if (u_v == true) return pow(i, 0.1);
                else return pow(1.0 / double(i), 0.1);

        if (num_rank == 2)
            if (u_v == true) return pow(1.0 / double(i), 0.1);
                else return pow(i, 0.1);
        return 0;
    }
    


    void build (long long v, long long tl, long long tr, vector <double> &t, bool u_v, int num_rank) {
	if (tl == tr) {
	    t[v] = (tl > 0) ? double(N_store[tl]) * K(tl, u_v, num_rank) : 0;
        } else {
	    long long tm = (tl + tr) / 2;
	    build (v*2, tl, tm, t, u_v, num_rank);
	    build (v*2+1, tm+1, tr, t, u_v, num_rank);
	    t[v] = t[v*2] + t[v*2+1];
	}
    }


    void update (long long v, long long tl, long long tr, long long pos, double new_val, vector <double> &t) {
	if (tl == tr) {
	    t[v] = new_val;
        } else {
	    long long tm = (tl + tr) / 2;
	    if (pos <= tm)
		update (v*2, tl, tm, pos, new_val, t);
	    else
		update (v*2+1, tm+1, tr, pos, new_val, t);
	    t[v] = t[v*2] + t[v*2+1];
	}
    }


    long long find (vector <double> &t){
        double random = (-floatDist(gen));
        double rand_val = random * t[1];
        long long v = 1;
        while (v < M) {
            v *= 2;
            if (t[v] < rand_val) {
                rand_val -= t[v];
                v += 1;
            }
        }
        return v - M;
    }


    void updating_all_trees(long long pos) {
        total_activity = 0;
        for (int w = 0 ; w < rank ; w++) {
            update(1, 0, M - 1, pos, K(pos, true, w + 1) * N_store[pos] , all_trees_u[w]);
            update(1, 0, M - 1, pos, K(pos, false, w + 1) * N_store[pos] , all_trees_v[w]);
            total_activity += all_trees_u[w][1] * all_trees_v[w][1];
        }
    }

    bool particle_selection(long long &i, long long &j) {
        // ВЫБОР КОМПАНЕНТЫ КАНОНИЧЕСКОГО РАЗЛОЖЕНИЯ
        int target_summand = 0;
        double random_val = (-floatDist(gen)) * total_activity;

        for (target_summand = 0; target_summand < rank - 1; target_summand++) {
            random_val -= all_trees_u[target_summand][1] * all_trees_v[target_summand][1];
            if (random_val <= 0) {
                break;
            }
        }

        // ВЫБОР РАЗМЕРОВ ЧАСТИЦ ДЛЯ ВЗАИМОДЕЙСТВИЯ
        i = find(all_trees_u[target_summand]);
        j = find(all_trees_v[target_summand]);

        // ПРОВЕРКА СУЩЕСТВОВАНИЯ ЧАСТИЦ  
        if (N_store[i] == 0 || N_store[j] == 0) {
            return false;
        }

        // ПРОВЕРКА СТОЛКНОВЕНИЯ С САМИМ СОБОЙ
        if ((i == j) && ((-floatDist(gen)) * N_store[i] <= 1)) {
            return false;
        }
        return true;
    }

    bool fragmentation(long long i, long long j) {
        double random = (-floatDist(gen));
        if (random < lamda / (1 + lamda)) {
            particle_change(i, -1);
            particle_change(j, -1);
            particle_change(1, i + j);
            return true;
        }
        return false;
    }

    private:
        ///// DATA /////////
        double V, lamda, cur_time;
        vector <long long> N_store;
        map <long long, double> p_k;
        map <long long, double> last_time_add;
        long long limit_size;
        int rank;

        long long mass, second_moment, largest_particle, M, N, N_0;
        double total_activity;
        unsigned long long count_coag;
        vector < vector <double> > all_trees_u;
        vector < vector <double> > all_trees_v;
         

        void particle_change(long long i, long long delta_count) {
            N += delta_count;
            mass += delta_count * i;
            second_moment += i * i * delta_count;
            if (i >= M) {
                M = 1 << (static_cast<long long>(log2(i)) + 1);
                N_store.resize(M, 0);
                for (int w = 0 ; w < rank ; w++) {
                    all_trees_u[w].resize(2 * M);
                    all_trees_v[w].resize(2 * M);
                    build(1, 0, M - 1, all_trees_u[w], true, w + 1);
                    build(1, 0, M - 1, all_trees_v[w], false, w + 1);
                }
            }
            N_store[i] += delta_count;
            if (i == largest_particle && N_store[i] == 0) {
                int j = i;
                for (; j >= 1 && N_store[j] == 0 ; j--); 
                largest_particle = j;
            }
            if (i > largest_particle && N_store[i] > 0) {
                largest_particle = i;
            } 
            updating_all_trees(i);
        }

        void scale_system(unsigned int alpha) {
            V *= alpha;
            mass *= alpha;
            second_moment *= alpha;
            N *= alpha;
            total_activity *= alpha * alpha;
            for (int w = 1; w < N_store.size(); w++) N_store[w] *= alpha;
            for (int w = 0 ; w < all_trees_u.size(); w++) {
                for (int s = 0 ; s < all_trees_u[w].size() ; s++) {
                    all_trees_u[w][s] *= alpha;
                    all_trees_v[w][s] *= alpha;
                }
            }
        }

    public:
        Modeling(double _V, double _lamda, vector <long long> &_N_store, long long _limit_size, map <long long, double> &_p_k) {
            cur_time = mass = second_moment = largest_particle = N = N_0 = count_coag = total_activity = 0;
            M = 1;
            rank = 2;
            all_trees_u = vector < vector <double> > (rank, vector <double> (2 * M));
            all_trees_v = vector < vector <double> > (rank, vector <double> (2 * M));
            lamda = _lamda;
            V = _V;
            p_k = _p_k;
            limit_size = _limit_size;
            for (long long i = 1; i < _N_store.size(); i++) {
                particle_change(i, _N_store[i]); 
            }
            N_0 = N;             
        }

        void simulate_process(double time, double V_max, long long period_logging) {
            double max_time = cur_time + time;
            long long int i, j;
            if (cur_time == 0) {
                long long max_num = -1;
                for (auto w: p_k) {
                    if (max_num == -1)
                        max_num = w.first;
                    if (p_k[max_num] < p_k[w.first])
                        max_num = w.first;
                }
                for (auto w: p_k) {
                    if (p_k[max_num] == p_k[w.first])
                        particle_change(w.first, 1);
                }
                if (max_num != -1) {
                    cur_time += 1 / (p_k[max_num] * V);
                }
            }

            while (cur_time < max_time) {
                if (total_activity > 0) {
                    cur_time += ((2 * V) / (total_activity)) * (1 / (1 + lamda));
                }
                for (auto w: p_k)
                    if ((cur_time - last_time_add[w.first]) * V * p_k[w.first] >= 1) {
                        long long extra = static_cast<long long>((cur_time - last_time_add[w.first]) * V * p_k[w.first]);
                        particle_change(w.first, extra);
                        last_time_add[w.first] += extra / (V * p_k[w.first]); 
                    }

                if (!particle_selection(i, j))
                    continue;

                count_coag++;
                if (count_coag % period_logging == 0) {
                     print_parametres_system();
                }

                if (fragmentation(i, j))
                    continue;

                particle_change(i, -1);
                particle_change(j, -1);
                if (i + j > limit_size) {
                    continue;
                }
                particle_change(i + j, 1);

                if (V < V_max) {
                    scale_system(2);
                }
            }
        }

        void print_parametres_system() {
            cout << " Volume " << V
                 << " Mass: " <<  mass
                 << " Density: " << mass / V 
                 << " Second_Moment: " << second_moment / V
                 << " Largest: " << largest_particle
                 << " N: " << N
                 << " total_activity " << total_activity
                 << " count_coag " << count_coag
                 << " Time: " << cur_time
                 << endl;
        }

        void write_to_file() {
            ofstream concentrations_end_out("concentrations_end_lr.txt");
            //concentrations_end_out << V << " " << largest_particle << endl; 
            for (int w = 1 ; w <= largest_particle; w++) {
                concentrations_end_out << double(N_store[w]) / V << " ";
            }
            concentrations_end_out.close();
        }

        void system_state(vector <long long> &_N_store, double &_V) {
            _N_store = N_store;
            _V = V;
        }
};


int main() {
    map <long long, double> p_k;
    p_k[1] = 1;
    p_k[100] = 0.1;
    vector <long long> N_store;
    Modeling test(1, 0, N_store, 1 << 15, p_k); // V, lamda, N_store, limit_size, p_k
    test.simulate_process(50, 1 << 20, 100000); // time, V_max, period_logging
    
    test.write_to_file();
    return 0;
}
