function [x_new, y_new, theta_new, v, a, phy, w] = FormInitialGuess(x, y, theta, t, s)
global st_graph_search_
if (length(t) > st_graph_search_.num_nodes_t)
    tt = t(end);
    ss = s(end);
    for ii = (length(t) - 1) : -1 : 1
        if (t(ii+1) ~= t(ii))
            tt = [t(ii), tt];
            ss = [s(ii), ss];
        end
    end
    tt(1) = 0; ss(1) = 0;
    t = tt; s = ss;
end
if (length(t) ~= st_graph_search_.num_nodes_t)
    error '[FormInitialGuess] Error code 1';
end
[x, y, theta] = ResamplePathWithEqualDistance2(x, y, theta);
ss = zeros(1,length(x));
for ii = 2 : length(x)
    ss(ii) = ss(ii-1) + hypot(x(ii) - x(ii-1), y(ii) - y(ii-1));
end
x_new = [];
y_new = [];
theta_new = [];
for ii = 1 : length(s)
    err = abs(s(ii) - ss);
    ind = find(err == min(err)); ind = ind(end);
    x_new = [x_new, x(ind)];
    y_new = [y_new, y(ind)];
    theta_new = [theta_new, theta(ind)];
end
[v, a, phy, w] = FulfillProfiles(x_new, y_new, theta_new);
end

function [x, y, theta] = ResamplePathWithEqualDistance2(x, y, theta)
for ii = 2 : length(theta)
    while (theta(ii) - theta(ii-1) > pi)
        theta(ii) = theta(ii) - 2 * pi;
    end
    while (theta(ii) - theta(ii-1) < -pi)
        theta(ii) = theta(ii) + 2 * pi;
    end
end
x_extended = [];
y_extended = [];
theta_extended = [];
for ii = 1 : (length(x) - 1)
    distance = hypot(x(ii+1)-x(ii), y(ii+1)-y(ii));
    LARGE_NUM = round(distance * 100);
    temp = linspace(x(ii), x(ii+1), LARGE_NUM);
    temp = temp(1,1:(LARGE_NUM - 1));
    x_extended = [x_extended, temp];
    
    temp = linspace(y(ii), y(ii+1), LARGE_NUM);
    temp = temp(1,1:(LARGE_NUM - 1));
    y_extended = [y_extended, temp];
    
    temp = linspace(theta(ii), theta(ii+1), LARGE_NUM);
    temp = temp(1,1:(LARGE_NUM - 1));
    theta_extended = [theta_extended, temp];
end
x_extended = [x_extended, x(end)];
y_extended = [y_extended, y(end)];
theta_extended = [theta_extended, theta(end)];
index = round(linspace(1, length(x_extended), 1000));
x = x_extended(index);
y = y_extended(index);
theta = theta_extended(index);
end

function [v, a, phy, w] = FulfillProfiles(x, y, theta)
global st_graph_search_
Nfe = st_graph_search_.num_nodes_t;
% Judge velocity direction
vdr = zeros(1, Nfe);
for ii = 2 : (Nfe - 1)
    addtion = (x(ii+1) - x(ii)) * cos(theta(ii)) + (y(ii+1) - y(ii)) * sin(theta(ii));
    if (addtion > 0)
        vdr(ii) = 1;
    else
        vdr(ii) = -1;
    end
end
v = zeros(1, Nfe);
a = zeros(1, Nfe);
dt = st_graph_search_.max_t / Nfe;
for ii = 2 : Nfe
    v(ii) = vdr(ii) * sqrt(((x(ii) - x(ii-1)) / dt)^2 + ((y(ii) - y(ii-1)) / dt)^2);
end
for ii = 2 : Nfe
    a(ii) = (v(ii) - v(ii-1)) / dt;
end
phy = zeros(1, Nfe);
w = zeros(1, Nfe);
global vehicle_kinematics_ vehicle_geometrics_
phy_max = vehicle_kinematics_.vehicle_phy_max;
w_max = vehicle_kinematics_.vehicle_w_max;
for ii = 2 : (Nfe-1)
    phy(ii) = atan((theta(ii+1) - theta(ii)) * vehicle_geometrics_.vehicle_wheelbase / (dt * v(ii)));
    if (phy(ii) > phy_max)
        phy(ii) = phy_max;
    elseif (phy(ii) < -phy_max)
        phy(ii) = -phy_max;
    end
end
for ii = 2 : (Nfe-1)
    w(ii) = (phy(ii+1) - phy(ii)) / dt;
    if (w(ii) > w_max)
        w(ii) = w_max;
    elseif (w(ii) < -w_max)
        w(ii) = -w_max;
    end
end
end