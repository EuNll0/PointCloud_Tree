#include <pcl/io/ply_io.h>
#include <iostream>
#include <chrono>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>

int main()
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_0(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_1(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_2(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_4(new pcl::PointCloud<pcl::PointXYZRGB>);
    pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("viewer"));
    pcl::io::loadPLYFile("/home/eunll0/Desktop/work/0419/pointcloud_1.ply", *cloud_4);
    pcl::io::loadPLYFile("./test.ply", *cloud_2);
    viewer->addCoordinateSystem(1);
    viewer->initCameraParameters();

    pcl::PointXYZ o(1, 2, 3);
    pcl::PointXYZ d = pcl::PointXYZ(1 + 1 * -3, 2 + 1 * 2, 3 + 1 * -5);

    viewer->addSphere(o,0.3);
    viewer->addLine(o, d);

    viewer->addPointCloud(cloud_4, "1");
    viewer->addPointCloud(cloud_2, "2");
    viewer->spin();
    return 0;

    pcl::io::loadPLYFile("./t2.ply", *cloud_2);
    pcl::io::loadPLYFile("./t1.ply", *cloud_1);
    pcl::io::loadPLYFile("./t0.ply", *cloud_0);
    *cloud_4 += (*cloud_1);
    *cloud_4 += (*cloud_0);
    *cloud_4 += (*cloud_2);

    auto beforeTime = std::chrono::steady_clock::now();
    auto afterTime = std::chrono::steady_clock::now();
    double duration_millsecond = std::chrono::duration<double, std::milli>(afterTime - beforeTime).count();
    std::cout << "The time to build tree :" << duration_millsecond;
    return 0;
    for (size_t i = 0; i < cloud_4->points.size(); i++)
    {
        for (size_t j = 0; j < cloud_2->points.size(); j++)
        {
            if (cloud_4->points[i].x == cloud_2->points[j].x)
            {
                cloud_4->points[i].g = 254;
                cloud_4->points[i].b = 2;
                cloud_4->points[i].r = 0;
            }
        }
    }

    for (size_t i = 0; i < cloud_4->points.size(); i++)
    {
        for (size_t j = 0; j < cloud_1->points.size(); j++)
        {
            if (cloud_4->points[i].x == cloud_1->points[j].x)
            {
                cloud_4->points[i].r = 254;
                cloud_4->points[i].g = 2;
                cloud_4->points[i].b = 0;
            }
        }
    }

    for (size_t i = 0; i < cloud_4->points.size(); i++)
    {
        for (size_t j = 0; j < cloud_0->points.size(); j++)
        {
            if (cloud_4->points[i].x == cloud_0->points[j].x)
            {
                cloud_4->points[i].b = 254;
                cloud_4->points[i].g = 2;
                cloud_4->points[i].r = 0;
            }
        }
    }

    viewer->addPointCloud(cloud_4);
    viewer->addCoordinateSystem(1);

    // pcl::PointXYZ o(1, 2, 3);
    // pcl::PointXYZ d = pcl::PointXYZ(1 + 1 * -3, 2 + 1 * 2, 3 + 1 * -5);

    viewer->addLine(o, d);

    viewer->spin();
    pcl::io::savePLYFile("./final.ply", *cloud_4);
    return 0;
}