function [t, s] = SearchVelocityInStGraph(x, y, theta)
global st_graph_search_ costmap_ dynamic_obs
costmap_ = zeros(st_graph_search_.num_nodes_t, st_graph_search_.num_nodes_s);
for ii = 1 : st_graph_search_.num_nodes_t
    for jj = 1 : st_graph_search_.num_nodes_s
        cur_x = x(jj);
        cur_y = y(jj);
        cur_theta = theta(jj);
        for kk = 1 : size(dynamic_obs, 2)
            if (IsVehicleCollidingWithMovingObstacle(cur_x, cur_y, cur_theta, dynamic_obs{ii,kk}))
                costmap_(ii, jj) = 1;
                continue;
            end
        end
    end
end
ind_vec = SearchStPathViaAStar();
ind1 = ind_vec(:,1)';
ind2 = ind_vec(:,2)';
ind1n = ind1(1); ind2n = ind2(1);
for ii = 2 : length(ind1)
    if (ind1(ii) ~= ind1(ii-1))
        ind1n = [ind1n, ind1(ii)];
        ind2n = [ind2n, ind2(ii)];
    end
end
t = (ind1n - 1) .* st_graph_search_.resolution_t;
s = (ind2n - 1) .* st_graph_search_.resolution_s;
end

function ind_vec = SearchStPathViaAStar()
global st_graph_search_ costmap_
grid_space_2D_ = cell(st_graph_search_.num_nodes_t, st_graph_search_.num_nodes_s);
init_node = zeros(1,11);
goal_ind = [st_graph_search_.num_nodes_t, st_graph_search_.num_nodes_s];
% Information of each element in each node:
% Dim # | Variable
%  1        null
%  2        null
%  3        f
%  4        g
%  5        h
%  6        is_in_openlist
%  7        is_in_closedlist
%  8-9      index of current node
%  10-11    index of parent node
init_node(4) = 0;
init_node(6) = 1;
init_node(8:9) = [1, 1];
init_node(5) = sum(abs(init_node(8:9) - goal_ind));
init_node(3) = init_node(4) + st_graph_search_.multiplier_H_for_A_star * init_node(5);
init_node(10:11) = [-999,-999];
openlist_ = init_node;
grid_space_2D_{init_node(8), init_node(9)} = init_node;
% % We choose 5 neiborhood expansion, which indicates that we enable that s
% may NOT change monotonously
expansion_pattern = [
    0 1;
    0 -1;
    1 1;
    1 -1;
    1 0];
expansion_length = [
    1 + st_graph_search_.penalty_for_inf_velocity;
    1 + st_graph_search_.penalty_for_inf_velocity;
    1.414;
    1.414;
    1];
completeness_flag = 0;

iter = 0;
while (~isempty(openlist_))
    iter = iter + 1;
    cur_node_order = find(openlist_(:,3) == min(openlist_(:,3))); cur_node_order = cur_node_order(end);
    cur_node = openlist_(cur_node_order, :);
    cur_ind = cur_node(8:9);
    if ((cur_ind(1) == goal_ind(1))&&(cur_ind(2) == goal_ind(2)))
        completeness_flag = 1;
        break;
    end
    cur_g = cur_node(4);
    % Remove cur_node from open list and add it in closed list
    openlist_(cur_node_order, :) = [];
    grid_space_2D_{cur_ind(1), cur_ind(2)}(6) = 0;
    grid_space_2D_{cur_ind(1), cur_ind(2)}(7) = 1;
    for ii = 1 : size(expansion_pattern,1)
        child_node_ind = cur_ind + expansion_pattern(ii,:);
        if ((child_node_ind(1) < 1)||(child_node_ind(1) > st_graph_search_.num_nodes_t)||(child_node_ind(2) < 1)||(child_node_ind(2) > st_graph_search_.num_nodes_s))
            continue;
        end
        % If the child node has been explored ever before, and then if the child has been within the closed list, abandon it and continue.
        if ((~isempty(grid_space_2D_{child_node_ind(1), child_node_ind(2)}))&&(grid_space_2D_{child_node_ind(1), child_node_ind(2)}(7) == 1))
            continue;
        end
        child_g = cur_g + expansion_length(ii);
        child_h = sum(abs(child_node_ind - goal_ind));
        child_f = child_g + st_graph_search_.multiplier_H_for_A_star * child_h;
        child_node_prepare = [0, 0, child_f, child_g, child_h, 1, 0, child_node_ind, cur_ind];
        % If the child node has been explored ever before
        if (~isempty(grid_space_2D_{child_node_ind(1), child_node_ind(2)}))
            % The child must be in the open list now, then check if its
            % recorded parent deserves to be switched as our cur_node.
            if (grid_space_2D_{child_node_ind(1), child_node_ind(2)}(4) > child_g + 0.1)
                child_node_order1 = find(openlist_(:,8) == child_node_ind(1));
                child_node_order2 = find(openlist_(child_node_order1,9) == child_node_ind(2));
                openlist_(child_node_order1(child_node_order2), :) = [];
                grid_space_2D_{child_node_ind(1), child_node_ind(2)} = child_node_prepare;
                openlist_ = [openlist_; child_node_prepare];
            end
        else % Child node has never been explored before
            % If the child node is collison free
            if (costmap_(child_node_ind(1),child_node_ind(2)) == 0)
                openlist_ = [openlist_; child_node_prepare];
                grid_space_2D_{child_node_ind(1), child_node_ind(2)} = child_node_prepare;
            else % If the child node involves collisons
                child_node_prepare(7) = 1;
                child_node_prepare(6) = 0;
                grid_space_2D_{child_node_ind(1), child_node_ind(2)} = child_node_prepare;
            end
        end
    end
end

if (completeness_flag)
    ind_vec = cur_node(8:9);
else
    ind_vec = [];
    return;
end
parent_ind = grid_space_2D_{ind_vec(1), ind_vec(2)}(10:11);
while (parent_ind(1) > -1)
    ind_vec = [parent_ind; ind_vec];
    parent_ind = grid_space_2D_{parent_ind(1), parent_ind(2)}(10:11);
end
end

function is_collided = IsVehicleCollidingWithMovingObstacle(x, y, theta, V)
is_collided = 0;
if (min(hypot(V.x - x, V.y - y)) > 10)
    return;
end
Vcar = CreateVehiclePolygonFull(x, y, theta);
if (any(inpolygon(Vcar.x, Vcar.y, V.x, V.y)))
    is_collided = 1;
    return;
end
if (any(inpolygon(V.x, V.y, Vcar.x, Vcar.y)))
    is_collided = 1;
    return;
end
end

function Vcar = CreateVehiclePolygonFull(x, y, theta)
global vehicle_geometrics_
cos_theta = cos(theta);
sin_theta = sin(theta);
vehicle_half_width = vehicle_geometrics_.vehicle_width * 0.5;
AX = x + (vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_wheelbase) * cos_theta - vehicle_half_width * sin_theta;
BX = x + (vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_wheelbase) * cos_theta + vehicle_half_width * sin_theta;
CX = x - vehicle_geometrics_.vehicle_rear_hang * cos_theta + vehicle_half_width * sin_theta;
DX = x - vehicle_geometrics_.vehicle_rear_hang * cos_theta - vehicle_half_width * sin_theta;
AY = y + (vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_wheelbase) * sin_theta + vehicle_half_width * cos_theta;
BY = y + (vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_wheelbase) * sin_theta - vehicle_half_width * cos_theta;
CY = y - vehicle_geometrics_.vehicle_rear_hang * sin_theta - vehicle_half_width * cos_theta;
DY = y - vehicle_geometrics_.vehicle_rear_hang * sin_theta + vehicle_half_width * cos_theta;
Vcar.x = [AX, BX, CX, DX, AX];
Vcar.y = [AY, BY, CY, DY, AY];
end