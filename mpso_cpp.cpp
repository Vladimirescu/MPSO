#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;

const int n_particles = 20;
const int n_dims = 2;

double schaffer_f6(double particle[2]){
    double sum_sq = pow(particle[0],2) + pow(particle[1], 2);
    double numarator = pow(sin(sqrt(sum_sq)), 2) - 0.5;
    double numitor = pow(1 + 0.001*sum_sq, 2);
    return 0.5 + numarator/numitor;
}

double get_norm(double x[n_dims]){
    double sum = 0.0;
    for(int d = 0; d<n_dims; d++) sum += pow(x[d],2);
    return sqrt(sum);
}

int run_mpso(double X_init[n_particles][n_dims], double V_init[n_particles][n_dims],
        double w, int max_iters, double tol, double c1, double c2, int verbose=1, double w_decrease=0){

    double X[n_particles][n_dims];
    double V[n_particles][n_dims];
    double X_best[n_particles][n_dims];
    double fitness_best[n_particles];

    for(int i=0; i<n_particles; i++){
        for(int j=0; j<n_dims; j++){
            X[i][j] = X_init[i][j];
            V[i][j] = V_init[i][j];
            X_best[i][j] = X_init[i][j];
        }
        fitness_best[i] = schaffer_f6(X[i]);
    }

    bool found_min = false;
    int n = 0;
    int best_index;
    double best_fitness;
    double rand1, rand2;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    while(n < max_iters){

        // update valori fitness si gaseste cea mai buna particula pana acum
        best_index = 0;
        best_fitness = schaffer_f6(X[0]);
        for(int i=0; i<n_particles; i++){
            // verifica daca noua pozitie este mai buna decat precedenta
            // daca da, modifica local best
            if(schaffer_f6(X[i]) < fitness_best[i]){
                fitness_best[i] = schaffer_f6(X[i]);
                for(int j=0; j<n_dims; j++){
                    X_best[i][j] = X[i][j];
                }
            }
            if(schaffer_f6(X[i]) < best_fitness){
                best_fitness = schaffer_f6(X[i]);
                best_index = i;
            }
        }

        if(best_fitness < tol && n != 0){
            found_min = true;
            if(verbose){
                cout << "Minimum found in " << n << " iterations" << endl;
                cout << "Best fitness value : " << best_fitness << endl;
                cout << "Best particle : " << endl;
                for(int i=0; i<n_dims; i++){
                    cout << X[best_index][i] << " ";
                }
            }
            break;
        }

        // calcul viteze + update pozitii
        for(int i=0; i<n_particles; i++){
            rand1 = distribution(generator);
            rand2 = distribution(generator);
            for(int j=0; j<n_dims; j++){
                V[i][j] = w * V[i][j] + c1 * rand1 * (X_best[i][j] - X[i][j]) +
                        c2 * rand2 * (X_best[best_index][j] - X[i][j]);
                // prag maxim pentru veolcity
                if(V[i][j] > 2){
                    V[i][j] = 2;
                }

                // update pozitie
                X[i][j] += V[i][j];

                // verificare conditii pentru pozitie
                if(X[i][j] > 100) X[i][j] = 100;
                if(X[i][j] < -100) X[i][j] = -100;
            }
        }

        w -= w_decrease;

        n++;
    }

    if(!found_min){
        if(verbose) cout << "Minimul nu a fost gasit." << endl;
        return -1;
    }else{
        return n;
    }
}


int main() {

    fstream fout;
    fout.open("results.csv", ios::out | ios::app);

    double X_i[n_particles][n_dims];
    double V_i[n_particles][n_dims];

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution_init_X(-100,100);
    std::uniform_real_distribution<double> distribution_init_V(0,2);

    double w[] = {1.8, 1.6, 1.4, 1.2, 1.1, 1.05, 1, 0.95, 0.9, 0.85, 0.8, 0.0};
    int n_failures[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int iter_sums[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 01};
    int max_iters = 4000;
    double tol = 0.001;
    double c1 = 2.0;
    double c2 = 2.0;
    int verbose=1;
    int iters;
    int n_experiments = 30;
    double w_decr = 1.4 / max_iters;
    int iter_sum_decr = 0;
    int n_failures_decr = 0;

    for(int e=0; e<n_experiments; e++) {
        cout << "Experiment " << e+1 << " >>>>> " << endl;
        for (int i = 0; i < n_particles; i++) {
            for (int j = 0; j < n_dims; j++) {
                X_i[i][j] = distribution_init_X(generator);
                V_i[i][j] = distribution_init_V(generator);
            }
        }
        fout << e+1 << ",";
        for(int k = 0; k < sizeof(w) / sizeof(w[0]); k++) {
            iters = run_mpso(X_i, V_i, w[k], max_iters, tol, c1, c2, verbose, 0);
            if (iters == -1) {
                cout << "W = " << w[k] << " : " << "---" << endl;
                n_failures[k] ++;
                fout << " " << ",";
            } else {
                cout << "W = " << w[k] << " : " << iters << endl;
                iter_sums[k] += iters;
                fout << iters << ",";
            }
        }

        // decrementare liniara w
        iters = run_mpso(X_i, V_i, 1.4, max_iters, tol, c1, c2, verbose, w_decr);
        if (iters == -1) {
            cout << "W descrestere liniara : " << "---" << endl;
            n_failures_decr ++;
            fout << " " << "\n";
        } else {
            cout << "W descrestere liniara : " << iters << endl;
            iter_sum_decr += iters;
            fout << iters << "\n";
        }
    }

    cout << "Final results : "<<endl;
    for(int k = 0; k < sizeof(w) / sizeof(w[0]); k++) {
        cout << "w = " << w[k] << " : " << n_failures[k] << " failures ; Avg. iters : "<<
        iter_sums[k] / (n_experiments - n_failures[k]) <<endl;
    }
    cout << "w descrestere liniara : " << n_failures_decr << " failures ; Avg. iters : "<<
    iter_sum_decr / (n_experiments - n_failures_decr);
}

