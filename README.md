# Interval-Based
Interval-Based Continuous Collision Detection 

This is an implementation of continuous collision detection methods described in [Snyder, John, Interval Analysis For Computer Graphics, 1992] and [Redon, Stephane and Kheddar, Abderrahmane and Coquillart, Sabine, Fast Continuous Collision Detection between Rigid Bodies, 2002]

:warning: CAUTION: The Code Is Tested Using GCC, Probably Also Works On MSVC and Clang


---

## Compiling Instruction 

To compile the code, first you need to install 
* CMake (https://cmake.org/), 
* Boost (https://www.boost.org/) 

in your system. 

To build the executable file, you can use CMake
```sh
cd Interval-Based/
mkdir build
cd build
cmake ../  -DCMAKE_BUILD_TYPE=Release
make
```


---
 
## Usage
Include `#include <interval_ccd/interval_ccd.hpp>`

To use the method described in [Snyder, 1992], use the function 

`bool vertexFaceCCD_Interval(const Eigen::Vector3d& vertex_start, const Eigen::Vector3d& face_vertex0_start, const Eigen::Vector3d& face_vertex1_start, const Eigen::Vector3d& face_vertex2_start, const Eigen::Vector3d& vertex_end, const Eigen::Vector3d& face_vertex0_end, const Eigen::Vector3d& face_vertex1_end, const Eigen::Vector3d& face_vertex2_end, double& toi);` 
for vertex-face ccd, 

use the function

`bool edgeEdgeCCD_Interval(const Eigen::Vector3d& edge0_vertex0_start, const Eigen::Vector3d& edge0_vertex1_start, const Eigen::Vector3d& edge1_vertex0_start, const Eigen::Vector3d& edge1_vertex1_start, const Eigen::Vector3d& edge0_vertex0_end, const Eigen::Vector3d& edge0_vertex1_end, const Eigen::Vector3d& edge1_vertex0_end, const Eigen::Vector3d& edge1_vertex1_end, double& toi);`
for edge-edge ccd.

To use the method described in [Redon, 2002], use the function 

`bool vertexFaceCCD_Redon(const Eigen::Vector3d& vertex_start, const Eigen::Vector3d& face_vertex0_start, const Eigen::Vector3d& face_vertex1_start, const Eigen::Vector3d& face_vertex2_start, const Eigen::Vector3d& vertex_end, const Eigen::Vector3d& face_vertex0_end, const Eigen::Vector3d& face_vertex1_end, const Eigen::Vector3d& face_vertex2_end, double& toi);` 
for vertex-face ccd, 

use the function

`bool edgeEdgeCCD_Redon(const Eigen::Vector3d& edge0_vertex0_start, const Eigen::Vector3d& edge0_vertex1_start, const Eigen::Vector3d& edge1_vertex0_start, const Eigen::Vector3d& edge1_vertex1_start, const Eigen::Vector3d& edge0_vertex0_end, const Eigen::Vector3d& edge0_vertex1_end, const Eigen::Vector3d& edge1_vertex0_end, const Eigen::Vector3d& edge1_vertex1_end, double& toi);`
for edge-edge ccd.


For vertex-face ccd, the inputs are the vertex and the triangle vertices in the begining and the end of the time step.

For edge-edge ccd, the inputs are the vertices of the two edges in the begining and the end of the time step.

 `toi` is the output impact time.



