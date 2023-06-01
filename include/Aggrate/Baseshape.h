#ifndef _BASESHAPE_H_
#define _BASESHAPE_H_

#include <cmath>
template <typename T>
struct Point3
{
    T x, y, z;
    Point3(){};
    Point3(T X, T Y, T Z) : x(X), y(Y), z(Z){};

    void normalize()
    {
        T length= sqrt(double(x*x+y*y+z*z));
        x = x/length;
        y = y/length;
        z= z/length;
    }

    Point3<T> operator-(const Point3<T> &a) const
    {
        return Point3<T>(x - a.x, y - a.y, z - a.z);
    }
    Point3<T> operator+(const Point3<T> &a) const
    {
        return Point3<T>(x + a.x, y + a.y, z + a.z);
    }
    friend Point3<T> operator*(T a, const Point3<T> &b)
    {
        return Point3<T>(a * b.x, a * b.y, a * b.z);
    }

    friend Point3<T> operator/( const Point3<T> &b,T a)
    {
        return Point3<T>(b.x / a, b.y / a, b.z / a);
    }

    Point3<float> operator()(const Point3<T> &a)
    {
        return Point3f(a.x, a.y, a.z);
    }

    T operator[](const int i) const
    {
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        return z;
    }

    T &operator[](const int i)
    {
        if (i == 0)
            return x;
        if (i == 1)
            return y;
        return z;
    }
    /* data */
};

template <typename T>
Point3<T> Min(const Point3<T> &p1, const Point3<T> &p2)
{
    return Point3<T>(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
                     std::min(p1.z, p2.z));
}

template <typename T>
Point3<T> Max(const Point3<T> &p1, const Point3<T> &p2)
{
    return Point3<T>(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
                     std::max(p1.z, p2.z));
}

template <typename T1, typename T2>
double Disanct_nosqrt(const Point3<T1> &p1, const Point3<T2> &p2)
{
    return ((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) + (p2.z - p1.z) * (p2.z - p1.z));
}

typedef Point3<double> Point3d;
typedef Point3<float> Point3f;
typedef Point3<int> Point3i;

struct Ray
{
    Point3d o, d;
    Ray(Point3d &O, Point3d &D) : o(O), d(D){};
    float tMax;
};

#endif