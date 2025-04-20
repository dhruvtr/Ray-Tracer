#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

struct vec3
{
    double x, y, z;

    vec3()
    {
        x = 0;
        y = 0;
        z = 0;
    }
    vec3(double a, double b, double c)
    {
        x = a;
        y = b;
        z = c;
    }

    double size() const { return std::sqrt(x * x + y * y + z * z); }
};

// Defne the operators for vec3
inline vec3 operator+(const vec3 &u, const vec3 &v) { return vec3(u.x + v.x, u.y + v.y, u.z + v.z); }
inline vec3 operator-(const vec3 &u) { return vec3(-u.x, -u.y, -u.z); }
inline vec3 operator-(const vec3 &u, const vec3 &v) { return vec3(u.x - v.x, u.y - v.y, u.z - v.z); }
inline vec3 operator*(const vec3 &u, const vec3 &v) { return vec3(u.x * v.x, u.y * v.y, u.z * v.z); }
inline vec3 operator*(double t, const vec3 &v) { return vec3(t * v.x, t * v.y, t * v.z); }
inline vec3 operator*(const vec3 &v, double t) { return t * v; }
inline vec3 operator/(const vec3 &v, double t) { return (1 / t) * v; }

inline double dot(const vec3 &u, const vec3 &v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline vec3 cross(const vec3 &u, const vec3 &v)
{
    return vec3(u.y * v.z - u.z * v.y,
                u.z * v.x - u.x * v.z,
                u.x * v.y - u.y * v.x);
}

inline vec3 normalise(const vec3 &v)
{
    return v / v.size();
}

#endif