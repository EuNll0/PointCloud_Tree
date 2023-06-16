#include <iostream>
#include "Aggrate/BVH_1.h"
#include "Aggrate/Memory.h"
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <glog/logging.h>
#include <chrono>

#define _pclshow

#ifdef _pclshow
template <typename T>
void AddCubeLine(pcl::visualization::PCLVisualizer::Ptr viewer, Bounds3<T> &bound)
{
    static int k = 1;
    static bool color = false;
    using namespace std;
    Point3<T> pMin = bound.pMin;
    Point3<T> pMax = bound.pMax;
    float width = (float)(pMax.x - pMin.x);
    float length = (float)(pMax.y - pMin.y);
    float height = (float)(pMax.z - pMin.z);
    pcl::PointXYZ p0((float)pMin.x, (float)pMin.y, (float)pMin.z);
    pcl::PointXYZ p1((float)pMax.x, (float)pMax.y, (float)pMax.z);
    pcl::PointXYZ p2(p0.x + width, p0.y, p0.z);
    pcl::PointXYZ p3(p0.x, p0.y + length, p0.z);
    pcl::PointXYZ p4(p0.x, p0.y, p0.z + height);
    pcl::PointXYZ p5(p1.x - width, p1.y, p1.z);
    pcl::PointXYZ p6(p1.x, p1.y - length, p1.z);
    pcl::PointXYZ p7(p1.x, p1.y, p1.z - height);
    if (color)
    {
        viewer->addLine(p0, p2, 254, 254, 254, "line" + to_string(k));
        viewer->addLine(p0, p3, 254, 254, 254, "line" + to_string(k + 1));
        viewer->addLine(p0, p4, 254, 254, 254, "line" + to_string(k + 2));
        viewer->addLine(p1, p5, 254, 254, 254, "line" + to_string(k + 3));
        viewer->addLine(p1, p6, 254, 254, 254, "line" + to_string(k + 4));
        viewer->addLine(p1, p7, 254, 254, 254, "line" + to_string(k + 5));
        viewer->addLine(p2, p7, 254, 254, 254, "line" + to_string(k + 6));
        viewer->addLine(p2, p6, 254, 254, 254, "line" + to_string(k + 7));
        viewer->addLine(p3, p5, 254, 254, 254, "line" + to_string(k + 8));
        viewer->addLine(p3, p7, 254, 254, 254, "line" + to_string(k + 9));
        viewer->addLine(p4, p5, 254, 254, 254, "line" + to_string(k + 10));
        viewer->addLine(p4, p6, 254, 254, 254, "line" + to_string(k + 11));
    }
    else
    {
        viewer->addLine(p0, p2, 254, 254, 254, "line" + to_string(k));
        viewer->addLine(p0, p3, 254, 254, 254, "line" + to_string(k + 1));
        viewer->addLine(p0, p4, 254, 254, 254, "line" + to_string(k + 2));
        viewer->addLine(p1, p5, 254, 254, 254, "line" + to_string(k + 3));
        viewer->addLine(p1, p6, 254, 254, 254, "line" + to_string(k + 4));
        viewer->addLine(p1, p7, 254, 254, 254, "line" + to_string(k + 5));
        viewer->addLine(p2, p7, 254, 254, 254, "line" + to_string(k + 6));
        viewer->addLine(p2, p6, 254, 254, 254, "line" + to_string(k + 7));
        viewer->addLine(p3, p5, 254, 254, 254, "line" + to_string(k + 8));
        viewer->addLine(p3, p7, 254, 254, 254, "line" + to_string(k + 9));
        viewer->addLine(p4, p5, 254, 254, 254, "line" + to_string(k + 10));
        viewer->addLine(p4, p6, 254, 254, 254, "line" + to_string(k + 11));
    }
    // viewer->addLine(p0,p1,0,0,0,"line" + to_string(k+6));
    k += 12;
    color = !color;
}
#endif

int main(int argc, char *argv[])
{
    using namespace pcl;
    using namespace std;
    google::InitGoogleLogging(argv[0]);

    FLAGS_alsologtostderr = 1;
    // FLAGS_log_dir = "./log";
    FLAGS_colorlogtostderr = true;

    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::io::loadPCDFile<pcl::PointXYZRGB>("/home/eunll0/Desktop/work/0421/rabbit.pcd", *cloud);
    // pcl::io::loadPCDFile<pcl::PointXYZ>("/home/bay/Desktop/work/map.pcd", *cloud);

#ifdef _pclshow
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("viewer"));
    viewer->addCoordinateSystem(1);
#endif

    auto &points = cloud->points;
    vector<Point3f> pts;
    for (auto &p : points)
    {
        p.x *= 100;
        p.y *= 100;
        p.r = 254;
        p.g = 254;
        p.b = 254;
        p.z *= 100;
        pts.push_back(Point3f(p.x, p.y, p.z));
    }

    LOG(INFO) << "Point sizes : " << points.size();
    auto beforeTime = std::chrono::steady_clock::now();
    // BVH_ACC1 bvh(pts, 0.1, 0.2,BVH_ACC1::SplitMethod::EqualCounts);
    // BVH_ACC1 bvh(pts,0.1,0.2,BVH_ACC1::SplitMethod::Middle,8);
    BVH_ACC1 bvh(pts,0.05,0.09,5,BVH_ACC1::SplitMethod::SAH,8);
    auto afterTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "The time to build tree :" << duration_millsecond;

    vector<Bound3f> ret;

    pcl::PointXYZ o(1, 2, 3);

    Point3f oo(o.x, o.y, o.z), dd(-3, 2, -5);
    dd.normalize();
    Ray ray(oo, dd);
    pcl::PointXYZ d = pcl::PointXYZ(o.x + 5 * dd.x, o.y + 5 * dd.y, o.z + 5 * dd.z);

    beforeTime = std::chrono::steady_clock::now();
    bvh.IntersectPB(ray, ret,-1);
    // bvh.GetBoundTree(2,ret);
    afterTime = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "Search block: " << ret.size();
    LOG(INFO) << "The time to search tree :" << duration_millsecond;

    // bvh.GetBoundTree(2,ret);

#ifdef _pclshow
    vector<uint> rett;
    vector<uint> px;

    beforeTime = std::chrono::steady_clock::now();
    bvh.IntersectP(ray, rett, -1);
    afterTime = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "The time to search tree :" << duration_millsecond<<"\t"<<rett.size();

  
    DLOG(INFO) << rett.size();
    for (uint i = 0; i < rett.size(); i++)
    {
        cloud->points[rett[i]].r = 253;
        cloud->points[rett[i]].g = 2;
        cloud->points[rett[i]].b = 2;
    }

    viewer->addSphere(o, 0.5);
    viewer->addLine(o, d);
    viewer->addPointCloud(cloud);
    for (uint i = 0; i < ret.size(); i++)
        AddCubeLine(viewer, ret[i]);
    viewer->spin();
#endif

    google::ShutdownGoogleLogging();
    return 0;
}
