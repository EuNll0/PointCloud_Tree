# PointCloud_Tree
ponint cloud library 

```
mkdir bin
cd bin
cmake ..
make -j3
```
最终版本。
构树的输入参数voxel_length代表点云分辨率，minBoundLength代表最小包围盒宽度，与main版本不同，点云分辨率不一定要小于最小包围盒宽度。
在构树时，若点云包成的包围盒宽度小于最小包围盒宽度则生成叶子，若该包围盒内的点云个数小于10时同样生成叶子，通过该种方式可以极大地降低内存占用情况。

使用LinearBVHNode的pad变量标识该叶子是最小包围盒生成的叶子还是最小个数生成的叶子。
在进行光线碰撞时，若碰撞的盒子标识为最小包围盒，则将该叶子所有点云输出，否则则针对其中每一个点云进行光线碰撞。
