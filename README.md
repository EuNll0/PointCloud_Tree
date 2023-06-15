# PointCloud_Tree
ponint cloud library 

```
mkdir bin
cd bin
cmake ..
make -j3
```
针对点云构建包围盒，包围盒参数voxel_length代表点云分辨率，minBoundLength代表最小包围盒宽度。
main分支下，点云分辨率必须小于最小包围盒长度，在构建树时若点云围成的包围盒小于最小包围盒长度则生成叶子节点。
