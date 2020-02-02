function WriteInitialGuessForNLP(x, y, theta, v, a, phy, w)
global vehicle_geometrics_
fl = vehicle_geometrics_.vehicle_front_hang + vehicle_geometrics_.vehicle_wheelbase;
rl = vehicle_geometrics_.vehicle_rear_hang;
hw = vehicle_geometrics_.vehicle_width / 2;
x = ConvertSimpleProfileToOCDTForm(x);
y = ConvertSimpleProfileToOCDTForm(y);
theta = ConvertSimpleProfileToOCDTForm(theta);
v = ConvertSimpleProfileToOCDTForm(v);
a = ConvertSimpleProfileToOCDTForm(a);
phy = ConvertSimpleProfileToOCDTForm(phy);
w = ConvertSimpleProfileToOCDTForm(w);
Nfe = size(x,1);
delete('IG.INIVAL');
fid = fopen('IG.INIVAL', 'w');
for ii = 1 : Nfe
    for jj = 0 : 3
        fprintf(fid, 'let x[%g,%g] := %f;\r\n', ii, jj, x(ii,jj+1));
        fprintf(fid, 'let y[%g,%g] := %f;\r\n', ii, jj, y(ii,jj+1));
        fprintf(fid, 'let theta[%g,%g] := %f;\r\n', ii, jj, theta(ii,jj+1));
        fprintf(fid, 'let v[%g,%g] := %f;\r\n', ii, jj, v(ii,jj+1));
        fprintf(fid, 'let phy[%g,%g] := %f;\r\n', ii, jj, phy(ii,jj+1));
        fprintf(fid, 'let egoV[%g, %g, 1, 1] := %f;\r\n', ii, jj, x(ii,jj+1) + fl * cos(theta(ii,jj+1)) - hw * sin(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 2, 1] := %f;\r\n', ii, jj, x(ii,jj+1) + fl * cos(theta(ii,jj+1)) + hw * sin(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 3, 1] := %f;\r\n', ii, jj, x(ii,jj+1) - rl * cos(theta(ii,jj+1)) + hw * sin(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 4, 1] := %f;\r\n', ii, jj, x(ii,jj+1) - rl * cos(theta(ii,jj+1)) - hw * sin(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 1, 2] := %f;\r\n', ii, jj, y(ii,jj+1) + fl * sin(theta(ii,jj+1)) + hw * cos(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 2, 2] := %f;\r\n', ii, jj, y(ii,jj+1) + fl * sin(theta(ii,jj+1)) - hw * cos(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 3, 2] := %f;\r\n', ii, jj, y(ii,jj+1) - rl * sin(theta(ii,jj+1)) - hw * cos(theta(ii,jj+1)));
        fprintf(fid, 'let egoV[%g, %g, 4, 2] := %f;\r\n', ii, jj, y(ii,jj+1) - rl * sin(theta(ii,jj+1)) + hw * cos(theta(ii,jj+1)));
    end
    for jj = 1 : 3
        fprintf(fid, 'let a[%g,%g] := %f;\r\n', ii, jj, a(ii,jj+1));
        fprintf(fid, 'let w[%g,%g] := %f;\r\n', ii, jj, w(ii,jj+1));
    end
end
fclose(fid);
end

function x_new = ConvertSimpleProfileToOCDTForm(x)
Nfe = (length(x) - 1) / 3;
x_new = zeros(Nfe, 4);
x_new(1,1:4) = x(1:4);
x(1:4) = [];
for ii = 2 : Nfe
    x_new(ii,1) = x_new(ii-1,4);
    x_new(ii,2:4) = x(1:3);
    x(1:3) = [];
end
if (length(x) > 0)
    error '[WriteInitialGuessForNLP] code 1'
end
end