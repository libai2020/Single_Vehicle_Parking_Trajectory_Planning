function ind = ConvertConfigToIndex(x, y, time)
global xyt_graph_search_ environment_scale_
ind1 = round((x - environment_scale_.environment_x_min) / xyt_graph_search_.resolution_x) + 1;
ind2 = round((y - environment_scale_.environment_y_min) / xyt_graph_search_.resolution_y) + 1;
ind3 = round(time / xyt_graph_search_.resolution_t) + 1;
if (ind1 < 1)
    ind1 = 1;
elseif (ind1 > xyt_graph_search_.num_nodes_x)
    ind1 = xyt_graph_search_.num_nodes_x;
end
if (ind2 < 1)
    ind2 = 1;
elseif (ind2 > xyt_graph_search_.num_nodes_y)
    ind2 = xyt_graph_search_.num_nodes_y;
end
if (ind3 < 1)
    ind3 = 1;
elseif (ind3 > xyt_graph_search_.num_nodes_t)
    ind3 = xyt_graph_search_.num_nodes_t;
end
ind = [ind1, ind2, ind3];
end