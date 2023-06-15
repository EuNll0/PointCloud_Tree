#ifndef _BVH_1_H_
#define _BVH_1_H_
#include "Aggrate/Common.h"
#include <utility>
#include <limits>
#include <glog/logging.h>
#include <sstream>
#include <string>

struct BucketInfo
{
    int count = 0;
    Bound3f bounds;
};

struct BVHBuildNode
{
    bool is_only;
    Bounds3<float> bounds;
    BVHBuildNode *children[2];
    int splitAxis, firstPrimOffset, nPrimitives;

    void InitLeaf(int first, int n, const Bounds3<float> &b)
    {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = std::move(b);
        children[0] = children[1] = nullptr;
    }

    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1)
    {
        children[0] = c0;
        children[1] = c1;
        // DCHECK(c0->bounds != nullptr);
        // DCHECK(c1->bounds != nullptr);
        bounds = std::move(Union(c0->bounds, c1->bounds));
        splitAxis = axis;
        nPrimitives = 0;
    }
};
struct LinearBVHNode
{
    Bounds3<float> bounds;
    union
    {
        uint primitivesOffset;  // leaf
        uint secondChildOffset; // interior
    };
    uint16_t nPrimitives; // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[1];
};

class BVH_ACC1
{
public:
    enum class SplitMethod
    {
        SAH,
        Middle,
        EqualCounts
    };
    const double _voxel_length;
    const double _minBoundLength;
    const int nBuckets = 14;

    const SplitMethod _method;
    LinearBVHNode *nodes = nullptr;
    std::vector<Point3<float>> &_pointcloud;
    std::vector<uint> orderdata;
    uint totalNodes = 0;

    BVH_ACC1(std::vector<Point3<float>> &pointcloud, double voxel_length = 0.5, double minBoundLength = 1,
             SplitMethod method = SplitMethod::Middle, uint MemorySize = 128) : _pointcloud(pointcloud), _method(method),
                                                                                _voxel_length(voxel_length / 2.), _minBoundLength(minBoundLength)
    {
        LOG(INFO) << "Begin to build the tree";
        MemoryArena area(MemorySize * 1024 * 1024);
        BVHBuildNode *root;
        std::vector<uint> pointInfo;
        pointInfo.resize(_pointcloud.size());
        // uint *p = new uint[_pointcloud.size()];

        // #   pragma omp parallel for
        for (uint i = 0; i < pointcloud.size(); i++)
        {
            pointInfo[i] = i;
        }

        orderdata.resize(pointInfo.size());
        uint offset = 0;
        if (method == SplitMethod::SAH)
        {
            root = recursiveBuild_SAH(area, 0, pointInfo.size(), pointInfo, &totalNodes);
        }
        else
        {
            root = recursiveBuild(area, 0, pointInfo.size(), pointInfo, &totalNodes);
        }
        nodes = AllocAligned<LinearBVHNode>(totalNodes);
        flattenBVHTree(root, &offset);
        char output[1024];
        sprintf(output, "BVH created with %u nodes for %lu "
                        "points (%.2f MB), arena allocated %.2f MB",
                totalNodes, pointInfo.size(),
                float(totalNodes * sizeof(LinearBVHNode)) /
                    (1024.f * 1024.f),
                float(area.TotalAllocated()) /
                    (1024.f * 1024.f));
        LOG(INFO) << std::string(output);
        // delete p;
        DCHECK_EQ(totalNodes, offset);
    }

    ~BVH_ACC1()
    {
        FreeAligned(nodes);
    };

    void IntersectP(const Ray &ray, std::vector<uint> &ret_pointIndex, uint depth = 1, double scale = 1 + 2 * gamma(3))
    {
        if (!nodes)
            return;
        Point3d invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
        uint nodesToVisit[1024] = {0};
        nodesToVisit[0] = 0;
        uint toVisitOffset = 1, currentNodeIndex = 0;
        uint leftnode_index = 1, rightnode_index = nodes->secondChildOffset;
        uint currentdepth = 0;
        while (toVisitOffset != 0)
        {
            LinearBVHNode *node = &nodes[currentNodeIndex];
            if (node->nPrimitives > 0)
            {
                uint block = 0;
                DLOG(INFO) << "DCHECK with leave!"
                           << "\t" << int(node->pad[0]);
                if (node->pad[0])
                {
                    for (int i = 0; i < node->nPrimitives; i++)
                    {
                        auto tmp = orderdata[i+node->primitivesOffset];
                        auto &point = _pointcloud[tmp];
                        auto vx = Point3f(_voxel_length, _voxel_length, _voxel_length);
                        if (Bound3f(point - vx, point + vx).IntersectP(ray, invDir, dirIsNeg, scale))
                        {
                            ret_pointIndex.push_back(tmp);
                            block++;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < node->nPrimitives; i++)
                    {
                        ret_pointIndex.push_back(orderdata[i + node->primitivesOffset]);
                    }
                }
                if (block > 0)
                    currentdepth++;
                if (currentdepth == depth)
                {
                    break;
                }
            }
            else
            {
                leftnode_index = currentNodeIndex + 1;
                rightnode_index = node->secondChildOffset;

                const LinearBVHNode *leftnode = &nodes[leftnode_index];
                const LinearBVHNode *rightnode = &nodes[rightnode_index];

                double leftInsect_ = leftnode->bounds.IntersectPD(ray, invDir, dirIsNeg, scale);
                double rightInsect_ = rightnode->bounds.IntersectPD(ray, invDir, dirIsNeg, scale);
                bool leftInsect = leftInsect_ > 0 ? 1 : 0;
                bool rightInsect = rightInsect_ > 0 ? 1 : 0;

                if (leftInsect && !rightInsect)
                {
                    nodesToVisit[++toVisitOffset] = leftnode_index;
                }
                else if (!leftInsect && rightInsect)
                {
                    nodesToVisit[++toVisitOffset] = rightnode_index;
                }

                else if (leftInsect && rightInsect)
                {
                    if (leftInsect_ < rightInsect_)
                    {
                        nodesToVisit[++toVisitOffset] = rightnode_index;
                        nodesToVisit[++toVisitOffset] = leftnode_index;
                    }
                    else
                    {
                        nodesToVisit[++toVisitOffset] = leftnode_index;
                        nodesToVisit[++toVisitOffset] = rightnode_index;
                    }
                    DLOG(INFO) << "DCHECK all intersect:  " << leftnode_index << "\t" << rightnode_index;
                }

                // DLOG(INFO) << "currentNodeIndex is " << currentNodeIndex << "\t the node status:\t" << leftInsect << "\t" << rightInsect;
                // DLOG(INFO) << "To visit offset is  " << toVisitOffset << "\n";
            }
            currentNodeIndex = nodesToVisit[toVisitOffset--];
        }
    }

    bool IntersectPB(const Ray &ray, std::vector<Bounds3<float>> &ret_bounds, uint depth = 1, double scale = 1 + 2 * gamma(3)) const
    {
        if (!nodes)
            return false;
        Point3d invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
        int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
        uint nodesToVisit[1024] = {0};
        nodesToVisit[0] = 0;
        uint toVisitOffset = 1, currentNodeIndex = 0;
        uint leftnode_index = 1, rightnode_index;
        uint currentdepth = 0;
        while (toVisitOffset != 0)
        {
            const LinearBVHNode *node = &nodes[currentNodeIndex];
            if (node->nPrimitives > 0)
            {
                DLOG(INFO) << "DCHECK with leave!"
                           << "\t" << int(node->pad[0]);
                uint block = 0;
                if (node->pad[0])
                {
                    for (int i = 0; i < node->nPrimitives; i++)
                    {
                        auto tmp = orderdata[i+node->primitivesOffset];
                        auto &point = _pointcloud[tmp];
                        auto vx = Point3f(_voxel_length, _voxel_length, _voxel_length);
                        auto voxel = Bound3f(point - vx, point + vx);
                        if (voxel.IntersectP(ray, invDir, dirIsNeg, scale))
                        {
                            ret_bounds.push_back(voxel);
                            block++;
                        }
                    }
                }
                else
                {
                    ret_bounds.push_back(node->bounds);
                }
                if (block > 0)
                    currentdepth++;
                if (currentdepth == depth)
                {
                    break;
                }
            }
            else
            {
                leftnode_index = currentNodeIndex + 1;
                rightnode_index = node->secondChildOffset;

                const LinearBVHNode *leftnode = &nodes[leftnode_index];
                const LinearBVHNode *rightnode = &nodes[rightnode_index];
                // ret_bounds.push_back(leftnode->bounds);
                // ret_bounds.push_back(rightnode->bounds);

                double leftInsect_ = leftnode->bounds.IntersectPD(ray, invDir, dirIsNeg, scale);
                double rightInsect_ = rightnode->bounds.IntersectPD(ray, invDir, dirIsNeg, scale);
                bool leftInsect = leftInsect_ > 0 ? 1 : 0;
                bool rightInsect = rightInsect_ > 0 ? 1 : 0;

                if (leftInsect && !rightInsect)
                {
                    nodesToVisit[++toVisitOffset] = leftnode_index;
                }
                else if (!leftInsect && rightInsect)
                {
                    nodesToVisit[++toVisitOffset] = rightnode_index;
                }

                else if (leftInsect && rightInsect)
                {
                    if (leftInsect_ < rightInsect_)
                    {
                        nodesToVisit[++toVisitOffset] = rightnode_index;
                        nodesToVisit[++toVisitOffset] = leftnode_index;
                    }
                    else
                    {
                        nodesToVisit[++toVisitOffset] = leftnode_index;
                        nodesToVisit[++toVisitOffset] = rightnode_index;
                    }
                    // DLOG(INFO) << "DCHECK all intersect:  " << leftnode_index << "\t" << rightnode_index;
                }

                // DLOG(INFO) << "currentNodeIndex is " << currentNodeIndex << "\t the node status:\t" << leftInsect << "\t" << rightInsect;
                // DLOG(INFO) << "To visit offset is  " << toVisitOffset << "\n";
            }
            currentNodeIndex = nodesToVisit[toVisitOffset--];
        }
        return false;
    }

    void GetBoundTree(const uint depth, std::vector<Bounds3<float>> &ret)
    {
        uint nodesToVisit[1024] = {0};
        uint nodeDepth[1024] = {0};
        nodesToVisit[0] = 0;
        uint toVisitOffset = 0, currentNodeIndex = 0, currentDepth = 0;
        while (true)
        {
            const LinearBVHNode *node = &nodes[currentNodeIndex];

            if (currentDepth == depth || node->nPrimitives > 0)
            {
                ret.push_back(node->bounds);
                if (toVisitOffset == 0)
                    break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
                currentDepth = nodeDepth[toVisitOffset];
            }
            else
            {
                currentDepth++;
                nodeDepth[toVisitOffset++] = currentDepth;
                nodesToVisit[toVisitOffset - 1] = currentNodeIndex + 1;
                currentNodeIndex = node->secondChildOffset;
            }
        }
    }

    void GetPointsFromNode(uint position, std::vector<uint> &ret)
    {
        if (position > totalNodes)
        {
            return;
        }
        LinearBVHNode *node = &nodes[position];
        int primitives_offset = node->primitivesOffset;
        uint16_t nprimitive = node->nPrimitives;
        DLOG(INFO) << "GetPointsFromNode(): primitives_offset   nprimitive" << primitives_offset << "   " << nprimitive;
        if (nprimitive > 0)
        {
            for (uint16_t i = 0; i < nprimitive; i++)
            {
                ret.push_back(orderdata[primitives_offset + i]);
            }
        }
        else
        {
            // Todo
        }
    }

    void SearchPoints(const Ray &ray, std::vector<uint> &input_points, double ridius, std::vector<uint> &output_points)
    {
        double ridius_sqrt = ridius * ridius;
        for (uint i = 0; i < input_points.size(); i++)
        {
            // Todo
        }
    }

private:
    uint flattenBVHTree(BVHBuildNode *node, uint *offset)
    {
        LinearBVHNode *linearNode = &nodes[*offset];
        linearNode->bounds = node->bounds;
        int myOffset = (*offset)++;
        if (node->nPrimitives > 0)
        {
            linearNode->primitivesOffset = node->firstPrimOffset;
            linearNode->nPrimitives = node->nPrimitives;
            linearNode->pad[0] = node->is_only ? 0 : 1;
        }
        else
        {
            linearNode->axis = node->splitAxis;
            linearNode->nPrimitives = 0;
            flattenBVHTree(node->children[0], offset);
            linearNode->secondChildOffset = flattenBVHTree(node->children[1], offset);
        }
        return myOffset;
    }

    BVHBuildNode *recursiveBuild(MemoryArena &area, uint start, uint end, std::vector<uint> &pointInfo,
                                 uint *tatalnodes)
    {
        // DCHECK((*tatalnodes) < alloc_all_memory) << "NOT ENOUGH MEMORY!!!";
        BVHBuildNode *node = area.Alloc<BVHBuildNode>();
        DCHECK(end > start);
        (*tatalnodes)++;
        int dim;

        static float xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = ymin = zmin = std::numeric_limits<float>::max();
        xmax = ymax = zmax = -std::numeric_limits<float>::max();

        // float sumx = 0, sumy = 0, sumz = 0;
        for (uint i = start; i < end; i++)
        {
            Point3<float> &p = _pointcloud[pointInfo[i]];
            xmin = std::min(xmin, p.x);
            ymin = std::min(ymin, p.y);
            zmin = std::min(zmin, p.z);
            xmax = std::max(xmax, p.x);
            ymax = std::max(ymax, p.y);
            zmax = std::max(zmax, p.z);
            // sumx += p.x;
            // sumy += p.y;
            // sumz += p.z;
        }

        // sumx /= (end-start);
        // sumy /= (end-start);
        // sumz /= (end-start);
        // float corx = 0,cory=0,corz=0;
        // for (uint i = start; i < end; i++)
        // {
        //     Point3<float> &p = *pointInfo[i].second;
        //     corx += (p.x - sumx) * (p.x - sumx) ;
        //     cory += (p.y - sumy) * (p.y - sumy) ;
        //     corz += (p.z - sumz) * (p.z - sumz) ;
        // }
        // if(corx >= cory && cory >= corz) {dim = 0;}
        // else if(cory >= corz ) {dim = 1;}
        // else  {dim = 2;}

        Bounds3<float> bounds(Point3<float>(xmin, ymin, zmin), Point3<float>(xmax, ymax, zmax));
        dim = bounds.MaximumExtent();
        uint nPrimitives = end - start;
        int mid = (start + end) / 2;
        auto &&diagonal = (bounds.Diagonal());
        if (diagonal.x < _minBoundLength && diagonal.y < _minBoundLength && diagonal.z < _minBoundLength)
        {
            bounds.pMax.x += _voxel_length;
            bounds.pMax.y += _voxel_length;
            bounds.pMax.z += _voxel_length;
            bounds.pMin.x -= _voxel_length;
            bounds.pMin.y -= _voxel_length;
            bounds.pMin.z -= _voxel_length;
            uint firstposition = orderdata.size();
            for (uint i = start; i < end; i++)
            {
                orderdata.push_back(pointInfo[i]);
            }
            node->InitLeaf(firstposition, nPrimitives, bounds);
            node->is_only = true;
            return node;
        }
        if (nPrimitives <= 10)
        {
            bounds.pMax.x += _voxel_length;
            bounds.pMax.y += _voxel_length;
            bounds.pMax.z += _voxel_length;
            bounds.pMin.x -= _voxel_length;
            bounds.pMin.y -= _voxel_length;
            bounds.pMin.z -= _voxel_length;
            uint firstposition = orderdata.size();
            for (uint i = start; i < end; i++)
            {
                orderdata.push_back(pointInfo[i]);
            }
            node->InitLeaf(firstposition, nPrimitives, bounds);
            node->is_only = false;
            return node;
        }

        auto &pointcloud = _pointcloud;
        switch (_method)
        {
        case SplitMethod::Middle:
        {
            float pmid = (bounds.pMax[dim] + bounds.pMin[dim]) / 2;
            auto p = std::partition(&pointInfo[start], &pointInfo[end - 1] + 1, [dim, pmid, &pointcloud](const auto &pi)
                                    { return (pointcloud[pi])[dim] < pmid; });
            mid = p - &pointInfo[0];
            // LOG(INFO) << "START MID END IS :   " << start << "\t" << mid << "\t" << end << "\t" << dim << "\n";
            break;
        }
        case SplitMethod::EqualCounts:
        {
            std::nth_element(&pointInfo[start], &pointInfo[mid], &pointInfo[end - 1] + 1,
                             [dim, &pointcloud](const auto &a, const auto &b)
                             { return (pointcloud[a])[dim] < (pointcloud[b])[dim]; });
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(area, start, mid, pointInfo, tatalnodes),
                           recursiveBuild(area, mid, end, pointInfo, tatalnodes));
        return node;
    }

    BVHBuildNode *recursiveBuild_SAH(MemoryArena &area, uint start, uint end, std::vector<uint> &info, uint *tatalnodes)
    {
        BVHBuildNode *node = area.Alloc<BVHBuildNode>();
        DCHECK(end > start);
        (*tatalnodes)++;
        static float xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = ymin = zmin = std::numeric_limits<float>::max();
        xmax = ymax = zmax = -std::numeric_limits<float>::max();

        for (uint i = start; i < end; i++)
        {
            Point3<float> &p = _pointcloud[info[i]];
            xmin = std::min(xmin, p.x);
            ymin = std::min(ymin, p.y);
            zmin = std::min(zmin, p.z);
            xmax = std::max(xmax, p.x);
            ymax = std::max(ymax, p.y);
            zmax = std::max(zmax, p.z);
        }
        float min_tmp[3] = {xmin, ymin, zmin};
        float max_tmp[3] = {xmax, ymax, zmax};

        Bounds3<float> bounds(Point3<float>(xmin, ymin, zmin), Point3<float>(xmax, ymax, zmax));
        auto &&diagonal = (bounds.Diagonal());
        if ((diagonal.x < _minBoundLength && diagonal.y < _minBoundLength && diagonal.z < _minBoundLength))
        {
            bounds.pMax.x += _voxel_length;
            bounds.pMax.y += _voxel_length;
            bounds.pMax.z += _voxel_length;
            bounds.pMin.x -= _voxel_length;
            bounds.pMin.y -= _voxel_length;
            bounds.pMin.z -= _voxel_length;
            uint firstposition = orderdata.size();
            for (uint i = start; i < end; i++)
            {
                orderdata.push_back(info[i]);
            }
            node->InitLeaf(firstposition, end - start, bounds);
            node->is_only = true;
            return node;
        }
        if (end - start <= nBuckets)
        {
            bounds.pMax.x += _voxel_length;
            bounds.pMax.y += _voxel_length;
            bounds.pMax.z += _voxel_length;
            bounds.pMin.x -= _voxel_length;
            bounds.pMin.y -= _voxel_length;
            bounds.pMin.z -= _voxel_length;
            uint firstposition = orderdata.size();
            for (uint i = start; i < end; i++)
            {
                orderdata.push_back(info[i]);
            }
            node->InitLeaf(firstposition, end - start, bounds);
            node->is_only = false;
            return node;
        }

        uint dim = bounds.MaximumExtent();

        BucketInfo buckets[nBuckets];

        float tobetween0to1 = 1.f / (max_tmp[dim] - min_tmp[dim]);
        for (uint i = start; i < end; i++)
        {
            int b = int(nBuckets * ((_pointcloud[info[i]][dim] - min_tmp[dim]) * tobetween0to1));
            if (b == nBuckets)
                b--;
            DCHECK_GE(b, 0);
            DCHECK_LT(b, nBuckets);
            buckets[b].count++;
            buckets[b].bounds = Union(buckets[b].bounds, _pointcloud[info[i]]);
        }

        float cost[nBuckets - 1];
        for (int i = 0; i < nBuckets - 1; i++)
        {
            Bound3f b0, b1;
            int count0 = 0, count1 = 0;
            for (int j = 0; j <= i; j++)
            {
                b0 = Union(b0, buckets[j].bounds);
                count0 += buckets[j].count;
            }

            for (int j = i + 1; j < nBuckets; j++)
            {
                b1 = Union(b1, buckets[j].bounds);
                count1 += buckets[j].count;
            }
            cost[i] = (count0 * b0.SurfaceArea()) + (count1 * b1.SurfaceArea());
        }

        float mincost = cost[0];
        int minCostSplitBucket = 0;
        for (int i = 1; i < nBuckets - 1; i++)
        {
            if (cost[i] < mincost)
            {
                mincost = cost[i];
                minCostSplitBucket = i;
            }
        }
        // LOG(INFO)<<minCostSplitBucket;
        auto p = std::partition(&(info[start]), &(info[end - 1]) + 1,
                                [&](const auto &a)
                                {
                                    int b = nBuckets * (_pointcloud[a][dim] - min_tmp[dim]) * tobetween0to1;
                                    if (b == nBuckets)
                                        b--;
                                    return b <= minCostSplitBucket;
                                    // return _pointcloud[a][dim] <pmid;
                                });
        int mid = p - &(info[0]);
        // LOG(INFO) << "START MID END IS " << start << "\t" << mid << "\t" << end << "\n";
        node->InitInterior(dim,
                           recursiveBuild_SAH(area, mid, end, info, tatalnodes),
                           recursiveBuild_SAH(area, start, mid, info, tatalnodes));
        return node;
    }
};

#endif
