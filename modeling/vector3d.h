#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <string>

template <typename T>
class Vector3d{
public:
    T x = 0;
    T y = 0;
    T z = 0;

    Vector3d(): x(0), y(0), z(0) {}
    Vector3d(T x, T y, T z): x(x), y(y), z(z) {}

    T norm() {
        return std::sqrt(x*x + y*y + z*z);
    }

    std::string represent() {
        std::ostringstream oss;
        oss << x << " " << y << " " << z;
        return oss.str();
    }
};

template <typename T>
Vector3d<T> operator%(const Vector3d<T> &a, const T b){
    return Vector3d<T>((a.x > 0) ? std::fmod(a.x, b) : std::fmod(a.x, b) + b, 
                    (a.y > 0) ? std::fmod(a.y, b) : std::fmod(a.y, b) + b, 
                    (a.z > 0) ? std::fmod(a.z, b) : std::fmod(a.z, b) + b);
}

template <typename T>
Vector3d<T> operator+(const Vector3d<T> &a, const Vector3d<T> &b){
    return Vector3d<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename T>
Vector3d<T> operator-(const Vector3d<T> &a, const Vector3d<T> &b){
    return Vector3d<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename T>
Vector3d<T> operator-(const Vector3d<T> &a){
    return Vector3d<T>(-a.x, -a.y, -a.z);
}

template <typename T>
Vector3d<T> operator*(const T a, const Vector3d<T> &vec){
    return Vector3d<T>(a*vec.x, a*vec.y, a*vec.z);
}

template <typename T>
Vector3d<T> operator*(const Vector3d<T> &vec, const T a){
    return Vector3d<T>(a*vec.x, a*vec.y, a*vec.z);
}

template <typename T>
Vector3d<T> operator/(const Vector3d<T> &vec, const T a){
    return Vector3d<T>(vec.x / a, vec.y / a, vec.z / a);
}

template <typename T>
Vector3d<T> operator*(const Vector3d<T> &a, const Vector3d<T> &b){
    return Vector3d<T>(a.x * b.x, a.y * b.y, a.z * b.z);
}

template <typename T>
Vector3d<T>& operator+=(Vector3d<T>& lhs, const Vector3d<T>& rhs) {
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

template <typename T>
Vector3d<T>& operator-=(Vector3d<T>& lhs, const Vector3d<T>& rhs) {
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
}

#endif