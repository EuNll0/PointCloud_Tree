#include <iostream>
#include <memory>
#include "Aggrate/Common.h"
#include "Aggrate/Memory.h"
#include "Aggrate/BVH_1.h"
#include <pcl/io/ply_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <glog/logging.h>

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

int main()
{
    using namespace std;
    using namespace pcl;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_1(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("viewer"));
    viewer->addCoordinateSystem(1);
    pcl::io::loadPLYFile("/home/eunll0/Desktop/work/0419/pointcloud.ply", *cloud);
    // pcl::io::loadPLYFile("./test.ply", *cloud_1);
    // {
    //     viewer->addPointCloud(cloud_1);
    //     viewer->spin();
    //     return 0;
    // }

    auto &points = cloud->points;
    vector<Point3f> pts;
    pts.reserve(points.size());
    for (auto &p : points)
    {
        pts.push_back(Point3f(p.x, p.y, p.z));
    }

    pcl::PointCloud<pcl::PointXYZ>().swap(*cloud);
    LOG(INFO) << "Point sizes : " << pts.size();

    auto beforeTime = std::chrono::steady_clock::now();
    // BVH_ACC1 bvh(pts,0.1,0.5,BVH_ACC1::SplitMethod::EqualCounts);
    BVH_ACC1 bvh(pts,0.01,0.02,BVH_ACC1::SplitMethod::EqualCounts);
    auto afterTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "The time to build tree :" << duration_millsecond;

    pcl::PointXYZ o(1, 2, 3);

    Point3d oo(o.x, o.y, o.z), dd(-3, 2, -5);
    dd.normalize();
    Ray ray(oo, dd);
    pcl::PointXYZ d = pcl::PointXYZ(o.x + 5 * dd.x, o.y + 5 * dd.y, o.z + 5 * dd.z);

    vector<Bound3f> ret;
    beforeTime = std::chrono::steady_clock::now();
    bvh.IntersectPB(ray, ret);
    // bvh.GetBoundTree(2,ret);
    afterTime = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "Block found: " << ret.size();
    LOG(INFO) << "The time of IntersectPB :" << duration_millsecond;

    viewer->addSphere(o, 0.1);
    viewer->addLine(o, d);
    for (uint i = 0; i < ret.size(); i++)
        AddCubeLine(viewer, ret[i]);

    vector<uint> rett;
    vector<uint> px;
    beforeTime = std::chrono::steady_clock::now();
    bvh.IntersectP(ray, rett,1);
    afterTime = std::chrono::steady_clock::now();
    duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    LOG(INFO) << "The time of IntersectP :" << duration_millsecond;

    for (int i = 0; i < rett.size(); i++)
    {
        bvh.GetPointsFromNode(rett[i], px);
    }
    DLOG(INFO) << rett.size();
    for (uint i = 0; i < px.size(); i++)
    {
        auto &node = bvh._pointcloud[px[i]];
        cloud_1->points.push_back(pcl::PointXYZRGB(node.x, node.y, node.z, 0, 250, 0));
    }

    viewer->addPointCloud(cloud_1);
    // viewer->addPointCloud(cloud,"kkk");

    viewer->spin();
    pcl::io::savePLYFile("./test.ply", *cloud_1);
    return 0;
}