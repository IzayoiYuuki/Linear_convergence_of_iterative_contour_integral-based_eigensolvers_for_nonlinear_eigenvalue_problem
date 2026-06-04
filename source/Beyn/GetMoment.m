function [Mom0, Mom1] = GetMoment(z, omega, V, Amu_func, dAmu_func, Sconfig, facA)
% this function solve the linear systems
% used in numerical integral
% with three different methods provided here.

[n, nex] = size(V);
center = (z(1) + z(1 + end/2))/2;

Mom0 = zeros(n, nex);
Mom1 = zeros(n, nex);

if Sconfig.lu == 0
    % slash
    for intp = 1 : length(z)
        Vtmp = Amu_func(z(intp))\(dAmu_func(z(intp))*V);
        Mom0 = Mom0 + omega(intp)*Vtmp;
        Mom1 = Mom1 + (z(intp) - center)*omega(intp)*Vtmp;
    end
else
    % sparse direct with LU decompostion
    for intp = 1 : length(z)
        Vtmp = facA{intp}\(dAmu_func(z(intp))*V);
        Mom0 = Mom0 + omega(intp)*Vtmp;
        Mom1 = Mom1 + (z(intp) - center)*omega(intp)*Vtmp;
    end
end
Mom0 = Mom0/(2*pi*1i);
Mom1 = Mom1/(2*pi*1i);