#ifndef _BOUNDS_H_
#define _BOUNDS_H_
#include "Aggrate/Common.h"
#include "Aggrate/Baseshape.h"
#include <glog/logging.h>

inline double gamma(int n);
template <typename T>
struct Bounds3
{
public:
    // Bounds3 Public Methods
    Bounds3()
    {
        T minNum = std::numeric_limits<T>::lowest();
        T maxNum = std::numeric_limits<T>::max();
        pMin = Point3<T>(maxNum, maxNum, maxNum);
        pMax = Point3<T>(minNum, minNum, minNum);
    }
    explicit Bounds3(const Point3<T> &p) : pMin(p), pMax(p) {}
    // inline bool Bounds3<T>::IntersectPP(const Ray &ray, const Point3d &invDir) const;
    inline bool IntersectP(const Ray &ray, const Point3d &invDir,
                           const int dirIsNeg[3], double scale = 1 + 2. * gamma(3)) const;
   inline double IntersectPD(const Ray &ray, const Point3d &invDir,
                                     const int dirIsNeg[3],const double& scale) const;
        Bounds3(const Point3<T> &p1, const Point3<T> &p2)
        : pMin(std::min(p1.x, p2.x), std::min(p1.y, p2.y),
               std::min(p1.z, p2.z)),
          pMax(std::max(p1.x, p2.x), std::max(p1.y, p2.y),
               std::max(p1.z, p2.z)) {}
    // const Point3<T> &operator[](int i) const;
    // Point3<T> &operator[](int i);

    inline Point3<T> &operator[](int i)
    {
        return (i == 0 ? pMin : pMax);
    }

    inline Point3<T> operator[](int i) const
    {
        return (i == 0 ? pMin : pMax);
    }
    // inline Bounds3<T> operator=(Bounds3<T> & another)
    // {
    //     return Bounds3<T>(another.pMin,another.pMax);
    // }

    inline Point3<T> center() const
    {
        Point3<T> ret = pMax - pMin;
        ret = ret / 2.;
        return ret;
    }

    bool operator==(const Bounds3<T> &b) const
    {
        return b.pMin == pMin && b.pMax == pMax;
    }
    bool operator!=(const Bounds3<T> &b) const
    {
        return b.pMin != pMin || b.pMax != pMax;
    }

    Point3<T> Corner(int corner) const
    {
        CHECK(corner >= 0 && corner < 8);
        return Point3<T>((*this)[(corner & 1)].x,
                         (*this)[(corner & 2) ? 1 : 0].y,
                         (*this)[(corner & 4) ? 1 : 0].z);
    }
    Point3<T> Diagonal() const { return pMax - pMin; }
    T SurfaceArea() const
    {
        Point3<T> d = Diagonal();
        return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
    }
    T Volume() const
    {
        Point3<T> d = Diagonal();
        return d.x * d.y * d.z;
    }
    int MaximumExtent() const
    {
        Point3<T> d = Diagonal();
        if (d.x > d.y && d.x > d.z)
            return 0;
        else if (d.y > d.z)
            return 1;
        else
            return 2;
    }

    Point3<T> Offset(const Point3<T> &p) const
    {
        Point3<T> o = p - pMin;
        if (pMax.x > pMin.x)
            o.x /= pMax.x - pMin.x;
        if (pMax.y > pMin.y)
            o.y /= pMax.y - pMin.y;
        if (pMax.z > pMin.z)
            o.z /= pMax.z - pMin.z;
        return o;
    }

    friend std::ostream &operator<<(std::ostream &os, const Bounds3<T> &b)
    {
        os << "[ " << b.pMin << " - " << b.pMax << " ]";
        return os;
    }

    // Bounds3 Public Data
    Point3<T> pMin, pMax;
};
typedef Bounds3<float> Bound3f;
typedef Bounds3<double> Bound3d;

template <typename T>
Bounds3<T> Union(const Bounds3<T> &b, const Point3<T> &p)
{
    Bounds3<T> ret;
    ret.pMin = Min(b.pMin, p);
    ret.pMax = Max(b.pMax, p);
    return ret;
}

template <typename T>
Bounds3<T> Union(const Bounds3<T> &b1, const Bounds3<T> &b2)
{
    Bounds3<T> ret;
    ret.pMin = Min(b1.pMin, b2.pMin);
    ret.pMax = Max(b1.pMax, b2.pMax);
    return ret;
}

template <typename T>
inline bool Bounds3<T>::IntersectP(const Ray &ray, const Point3d &invDir,
                                   const int dirIsNeg[3], double scale) const
{
    const Bounds3<T> &bounds = *this;
    // Check for ray intersection against $x$ and $y$ slabs
    double tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
    double tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
    double tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
    double tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

    // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
    // tMax *= 1 + 2.*  gamma(3);
    // tyMax *= 1 + 2. * gamma(3);
    tMax *= scale;
    tyMax *= scale;
    if (tMin > tyMax || tyMin > tMax)
        return false;
    if (tyMin > tMin)
        tMin = tyMin;
    if (tyMax < tMax)
        tMax = tyMax;

    // Check for ray intersection against $z$ slab
    double tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
    double tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

    // Update _tzMax_ to ensure robust bounds intersection
    // tzMax *= 1 + 2. * gamma(3);
    tzMax *= scale;
    if (tMin > tzMax || tzMin > tMax)
        return false;
    if (tzMin > tMin)
        tMin = tzMin;
    if (tzMax < tMax)
        tMax = tzMax;
    return (tMax > 0);
}

template <typename T>
inline double Bounds3<T>::IntersectPD(const Ray &ray, const Point3d &invDir,
                                     const int dirIsNeg[3],const double& scale) const
{
    const Bounds3<T> &bounds = *this;
    // Check for ray intersection against $x$ and $y$ slabs
    double tMin = (bounds[dirIsNeg[0]].x - ray.o.x) * invDir.x;
    double tMax = (bounds[1 - dirIsNeg[0]].x - ray.o.x) * invDir.x;
    double tyMin = (bounds[dirIsNeg[1]].y - ray.o.y) * invDir.y;
    double tyMax = (bounds[1 - dirIsNeg[1]].y - ray.o.y) * invDir.y;

    // Update _tMax_ and _tyMax_ to ensure robust bounds intersection
    // tMax *= 1 + 2.*  gamma(3);
    // tyMax *= 1 + 2. * gamma(3);
    tMax *= scale;
    tyMax *= scale;
    if (tMin > tyMax || tyMin > tMax)
        return false;
    if (tyMin > tMin)
        tMin = tyMin;
    if (tyMax < tMax)
        tMax = tyMax;

    // Check for ray intersection against $z$ slab
    double tzMin = (bounds[dirIsNeg[2]].z - ray.o.z) * invDir.z;
    double tzMax = (bounds[1 - dirIsNeg[2]].z - ray.o.z) * invDir.z;

    // Update _tzMax_ to ensure robust bounds intersection
    // tzMax *= 1 + 2. * gamma(3);
    tzMax *= scale;
    if (tMin > tzMax || tzMin > tMax)
        return false;
    if (tzMin > tMin)
        tMin = tzMin;
    if (tzMax < tMax)
        tMax = tzMax;
    return tMax;
}

// template <typename T>
// inline bool Bounds3<T>::IntersectPP(const Ray &ray, const Point3d &invDir) const
// {
//     double a, b, c, d, e, f;
//     a = pMin.x;
//     b = pMin.y;
//     c = pMin.z;
//     d = pMax.x;
//     e = pMax.y;
//     f = pMax.z;

// }

#endif