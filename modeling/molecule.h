#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include "vector3d.h"

template <typename T>
class Molecule{
public:
    Vector3d<T> position;
    Vector3d<T> velocity;
    T mass;

    Molecule(Vector3d<T> position, Vector3d<T> velocity, T mass): 
        position(position), velocity(velocity), mass(mass) {}
    
    void add_force(Vector3d<T> force){
        force_ = force_ + force;
    }

    void reset_force(){
        force_prev_ = force_;
        force_ = Vector3d<T>();
    }

    Vector3d<T> get_accelerate(){
        return force_ / mass;
    }

    Vector3d<T> get_prev_accelerate(){
        return force_prev_ / mass;
    }
private:
    Vector3d<T> force_ = Vector3d<T>();
    Vector3d<T> force_prev_ = Vector3d<T>();
};

#endif