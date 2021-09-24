%%Demonstration Code for Data Visualization Seminar Offered by the Chair of Computer Graphics and Visualization at TUM
%%Seminor Topic: Topology-based Tensor Field Visualization
%%Author: Junpeng Wang (junpeng.wang@tum.de)
%%Date: 2021-09-14
clear
clc

stressfileName = './data/cantilever2D_R500_iLoad5.vtk';

TensorTopologyAnalysis2D(stressfileName);
