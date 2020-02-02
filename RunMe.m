%  MATLAB Source Codes for the book "Cooperative Dedcision and Planning for
%  Connected and Automated Vehicles" published by Mechanical Industry Press
%  in 2020.
% 《智能网联汽车协同决策与规划技术》书籍配套代码
%  Copyright (C) 2020 Bai Li
%  2020.02.02
% ==============================================================================
%  第二章. 混合A*做路径决策，S-T图搜索做速度决策，数值求解最优控制问题做轨迹优化，并呈现结果
% ==============================================================================
%  备注：
%  1. 请读者支持正版AMPL
%  2. 基于该部分代码的研究成果须引用以下参考文献：
%  a) B. Li, and Z. Shao, “A unified motion planning method for
%     parking an autonomous vehicle in the presence of irregularly placed
%     obstacles,” Knowledge-Based Systems, vol. 86, pp. 11–20, 2015.
%  b) B. Li, K. Wang, and Z. Shao, “Time-optimal maneuver planning in
%     automatic parallel parking using a simultaneous dynamic optimization
%     approach,” IEEE Transactions on Intelligent Transportation Systems, vol.
%     17, no. 11, pp. 3263–3274, 2016.
% ==============================================================================
close all
clc

% % 参数设置
global vehicle_geometrics_ % 车辆轮廓几何尺寸
vehicle_geometrics_.vehicle_wheelbase = 2.8;
vehicle_geometrics_.vehicle_front_hang = 0.96;
vehicle_geometrics_.vehicle_rear_hang = 0.929;
vehicle_geometrics_.vehicle_width = 1.942;
vehicle_geometrics_.vehicle_length = vehicle_geometrics_.vehicle_wheelbase + vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_rear_hang;
vehicle_geometrics_.radius = hypot(0.25 * vehicle_geometrics_.vehicle_length, 0.5 * vehicle_geometrics_.vehicle_width);
vehicle_geometrics_.r2x = 0.25 * vehicle_geometrics_.vehicle_length - vehicle_geometrics_.vehicle_rear_hang;
vehicle_geometrics_.f2x = 0.75 * vehicle_geometrics_.vehicle_length - vehicle_geometrics_.vehicle_rear_hang;
global vehicle_kinematics_ % 车辆运动能力参数
vehicle_kinematics_.vehicle_v_max = 2.5;
vehicle_kinematics_.vehicle_a_max = 0.5;
vehicle_kinematics_.vehicle_phy_max = 0.7;
vehicle_kinematics_.vehicle_w_max = 0.5;
vehicle_kinematics_.vehicle_kappa_max = tan(vehicle_kinematics_.vehicle_phy_max) / vehicle_geometrics_.vehicle_wheelbase;
vehicle_kinematics_.vehicle_turning_radius_min = 1 / vehicle_kinematics_.vehicle_kappa_max;
global environment_scale_ % 车辆所在环境范围
environment_scale_.environment_x_min = -20;
environment_scale_.environment_x_max = 20;
environment_scale_.environment_y_min = -20;
environment_scale_.environment_y_max = 20;
environment_scale_.x_scale = environment_scale_.environment_x_max - environment_scale_.environment_x_min;
environment_scale_.y_scale = environment_scale_.environment_y_max - environment_scale_.environment_y_min;

% % 路径决策 + 速度决策过程
global st_graph_search_ % 用于S-T图搜索的A星算法涉及的参数
st_graph_search_.num_nodes_s = 80;
st_graph_search_.num_nodes_t = 40 * 3 + 1;
st_graph_search_.multiplier_H_for_A_star = 2.0;
st_graph_search_.penalty_for_inf_velocity = 4;
global vehicle_TPBV_ obstacle_vertexes_ % 导入既定算例以及静止障碍物分布情况
load TaskSetup.mat
[x, y, theta, path_length, completeness_flag] = ProvideCoarsePathViaHybridAStarSearch(); % 决策路径
st_graph_search_.resolution_s = path_length / st_graph_search_.num_nodes_s;
global dynamic_obs % 布设移动障碍物、指定运动时域并给出移动物体的运动轨迹
st_graph_search_.max_t = round(path_length * 2);
st_graph_search_.max_s = path_length;
st_graph_search_.resolution_t = st_graph_search_.max_t / st_graph_search_.num_nodes_t;
dynamic_obs = GenerateDynamicObstacles();
[t,s] = SearchVelocityInStGraph(x, y, theta); % 速度决策

% % 轨迹规划
[x, y, theta, v, a, phy, w] = FormInitialGuess(x, y, theta, t, s);
WriteInitialGuessForNLP(x, y, theta, v, a, phy, w);
WriteObstacleSetupsForNLP();
WriteBoundaryValues();
!ampl rr.run
% % 轨迹规划结果展示
load opti_flag.txt
if (opti_flag)
    set(0,'DefaultLineLineWidth',1);
    GenerateVideo();
    DrawStaticFigrue();
end