#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <optional>

#include <future>
#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "molecule.h"
#include "vector3d.h"

template <typename T>
class Environment{
private:
    char* config_path;
    char* coords_path;
    char* distribution_path;
    char* char_path;

    T molecule_mass;
    T eps;
    T sigma;
    T temperatute;
    T env_size;
    unsigned molecules_count;

    unsigned max_step_count;
    unsigned save_period;
    unsigned current_step = 0;
    T dt;
    std::vector<Molecule<T>> molecules;
    std::optional<std::default_random_engine> reng;
    unsigned n_treads;


public:
    Environment(char* config_path, char* coords_path, char* distribution_path, char* char_path): 
        config_path(config_path), coords_path(coords_path), distribution_path(distribution_path), char_path(char_path)
    {
        std::ifstream ifs(config_path);
        json jf = json::parse(ifs);

        std::ofstream out1(coords_path);
        out1.close();
        std::ofstream out3(distribution_path);
        out3.close();
        std::ofstream out2(char_path);
        out3.close();

        molecule_mass = jf["molecule_mass"];
        eps = jf["eps"];
        sigma = jf["sigma"];
        temperatute = jf["temperatute"];
        env_size = jf["env_size"];
        molecules_count = jf["molecules_count"];
        bool grid_init_mode = jf["grid_init_mode"];
        max_step_count = jf["max_step_count"];
        save_period = jf["save_period"];
        dt = jf["dt"];
        *reng = std::default_random_engine(jf["seed"]);
        n_treads = jf["n_treads"];
        
        unsigned mc3 = static_cast<unsigned>(std::pow(molecules_count, 1./3));
        this->molecules_count = mc3 * mc3 * mc3;
        T grid_size = static_cast<T>(env_size) / static_cast<T>(mc3);
        
        std::uniform_real_distribution<T> dstr(0, 1);
        std::normal_distribution<T> dstr_normal(0, std::sqrt(temperatute / molecule_mass));

        for (unsigned i = 0; i < mc3; i++){
            for (unsigned j = 0; j < mc3; j++){
                for (unsigned k = 0; k < mc3; k++){
                    if (grid_init_mode) {
                        Vector3d<T> position(grid_size * (i + 1. / 2), 
                                                    grid_size * (j + 1. / 2), 
                                                    grid_size * (k + 1. / 2));
                        
                        Vector3d<T> velocity_n(dstr_normal(*reng), 
                                            dstr_normal(*reng), 
                                            dstr_normal(*reng));
                        
                        molecules.push_back(Molecule<T>(position, velocity_n, molecule_mass));
                    }
                    else {
                        Vector3d<T> position(grid_size * (i + dstr(*reng) - 1. / 2), 
                                                    grid_size * (j + dstr(*reng) - 1. / 2), 
                                                    grid_size * (k + dstr(*reng) - 1. / 2));
                        
                        Vector3d<T> velocity_n(dstr_normal(*reng), 
                                            dstr_normal(*reng), 
                                            dstr_normal(*reng));
                        
                        molecules.push_back(Molecule<T>(position, velocity_n, molecule_mass));
                    }
            }
        }
        }
    }
    void simulate(){
        update_forces();
        for (current_step = 0; current_step < max_step_count; current_step++){
            if (current_step % save_period == 0){
                save_data();
            }
            step();
        }
    }

private:
    Vector3d<T> calc_force(const Molecule<T> &v1, const Molecule<T> &v2) {
        Vector3d<T> r = v2.position - v1.position % env_size;
        T p = sigma / r.norm();
        return 4 * eps / sigma / sigma * (6 * std::pow(p, 8) - 12 * std::pow(p, 14)) * r;
    }

    T calc_energy(const Molecule<T> &v1, const Molecule<T> &v2) {
        Vector3d<T> r = v2.position - v1.position % env_size;
        T p = std::pow(sigma / r.norm(), 6);
        return 4 * eps * (p * p - p);
    }
    
    void save_data(){
        std::ofstream outfile(coords_path, std::ios::app);
        outfile << molecules_count << std::endl;
        outfile << "Frame " << current_step << std::endl;
        for (unsigned i = 0; i < molecules_count; i++) {
            outfile << "O " << molecules[i].position.x << " " << molecules[i].position.y << " " << molecules[i].position.z << std::endl;
        }
        outfile.close();

        std::ofstream distribution_file(distribution_path, std::ios::app);
        distribution_file << molecules_count << std::endl;
        outfile << "Frame " << current_step << std::endl;
        for (unsigned i = 0; i < molecules_count; i++){
            distribution_file << molecules[i].velocity.x << " " << molecules[i].velocity.y << " " << molecules[i].velocity.z << std::endl;
        }
        distribution_file << std::endl;
        distribution_file.close();

        std::ofstream char_file(char_path, std::ios::app);
        char_file << calc_energy() << ' ' << calc_std_distanse() << ' ' << calc_temperature() << std::endl;
        distribution_file.close();
    }

    T calc_temperature(){
        T temp = 0;
        for (unsigned i = 0; i < molecules_count; ++i){
            T velocity = molecules[i].velocity.norm();
            temp += velocity * velocity;
        }

        return temp / molecules_count / 3;
    }

    T calc_std_distanse(){
        T std = 0;
        for (unsigned i = 0; i < molecules_count; ++i){
            T position = (molecules[i].position - Vector3d<T>(env_size / 2, env_size / 2, env_size / 2)).norm();
            std += position * position;
        }

        return std::sqrt(std / molecules_count);
    }

    T calc_subenergy(unsigned start, unsigned stop){
        T energy = 0;
        for (unsigned i = start; i < stop; ++i){
            energy += molecules[i].mass * std::pow(molecules[i].velocity.norm(), 2) / 2;
            for (unsigned j = 0; j < molecules_count; ++j){
                for (int k1 = -1; k1 < 2; ++k1){
                    for (int k2 = -1; k2 < 2; ++k2){
                        for (int k3 = -1; k3 < 2; ++k3){
                            if (i == j && k1 == 0 && k2 == 0 && k3 == 0) continue;
                             Molecule<T> virtual_m(
                                molecules[j].position % env_size + env_size * Vector3d<T>(k1, k2, k3), 
                                molecules[j].velocity,
                                molecules[j].mass
                            );
                            energy += calc_energy(molecules[i], virtual_m) / 2;
                        }
                    }
                }
            }
        }
        return energy;
    }

    void calc_subforce(unsigned start, unsigned stop){
        for (unsigned i = start; i < stop; ++i){
            molecules[i].reset_force();
            for (unsigned j = 0; j < molecules_count; ++j){
                if (i == j) continue;
                for (int k1 = -1; k1 < 2; ++k1){
                    for (int k2 = -1; k2 < 2; ++k2){
                        for (int k3 = -1; k3 < 2; ++k3){
                             Molecule<T> virtual_m(
                                molecules[j].position % env_size + env_size * Vector3d<T>(k1, k2, k3), 
                                molecules[j].velocity,
                                molecules[j].mass
                            );
                            molecules[i].add_force(calc_force(molecules[i], virtual_m));
                        }
                    }
                }
            }
        }
    }

    T calc_energy(){
        T energy = 0;
        std::future<T>* futures = new std::future<T>[n_treads];
        unsigned a = molecules_count / n_treads;
        unsigned b = molecules_count % n_treads;
        for (unsigned tread_k = 0; tread_k < n_treads; ++tread_k){
            unsigned start = tread_k*a;
            unsigned stop = (tread_k+1)*a;
            if (tread_k == n_treads - 1) stop += b;

            futures[tread_k] = std::async([this, start, stop]() {
                return this->calc_subenergy(start, stop);
                });
        }
        for (unsigned tread_k = 0; tread_k < n_treads; ++tread_k){
            energy += futures[tread_k].get();
            }
        delete[] futures;
        return energy;
    }

    void update_forces(){
        std::thread* futures = new std::thread[n_treads];
        unsigned a = molecules_count / n_treads;
        unsigned b = molecules_count % n_treads;
        for (unsigned tread_k = 0; tread_k < n_treads; ++tread_k){
            unsigned start = tread_k*a;
            unsigned stop = (tread_k+1)*a;
            if (tread_k == n_treads - 1) stop += b;

            futures[tread_k] = std::thread([this, start, stop]() {
                return this->calc_subforce(start, stop);
                });
        }
        for (unsigned tread_k = 0; tread_k < n_treads; ++tread_k){
            futures[tread_k].join();
            }
        delete[] futures;
        }

    void step() {
        for (unsigned i = 0; i < molecules_count; i++){
            molecules[i].position = molecules[i].position + molecules[i].velocity*dt + 0.5 * molecules[i].get_accelerate() * dt * dt;
        }
        update_forces();
        for (unsigned i = 0; i < molecules_count; i++){
            molecules[i].velocity = molecules[i].velocity + 0.5 * dt * (molecules[i].get_accelerate() + molecules[i].get_prev_accelerate());
        }
    }
};

int main(int argc, char* argv[]){
    if (argc != 5) {
        std::cout << "Неправильный ввод" << std::endl;
        return -1;
    }

    Environment<double> env(argv[1], argv[2], argv[3], argv[4]);
    env.simulate();
    return 0;
}
