function obstacle_cell = GenerateDynamicObstacles()
global xyt_graph_search_
load DynObs
Nobs = size(dynamic_obs,2);
obstacle_cell = cell(xyt_graph_search_.num_nodes_t, Nobs);
for ii = 1 : Nobs
    dx = dynamic_obs{end,ii}.x(1) - dynamic_obs{1,ii}.x(1);
    dy = dynamic_obs{end,ii}.y(1) - dynamic_obs{1,ii}.y(1);
    for jj = 1 : xyt_graph_search_.num_nodes_t
        temp.x = dynamic_obs{1,ii}.x + dx / xyt_graph_search_.num_nodes_t * (jj - 1);
        temp.y = dynamic_obs{1,ii}.y + dy / xyt_graph_search_.num_nodes_t * (jj - 1);
        temp.A = dynamic_obs{1,ii}.A;
        obstacle_cell{jj, ii} = temp;
    end
end
end