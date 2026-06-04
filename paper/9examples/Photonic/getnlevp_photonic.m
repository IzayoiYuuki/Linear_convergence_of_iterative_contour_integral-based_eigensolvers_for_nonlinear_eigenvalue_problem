function [COEFFS, FUN, Amu_func] = getnlevp_photonic()
% This function will return COEFFS, FUN, Amu_func just
% like what nlevp will do

file_mat_M17;
file_mat_M16;
file_mat_M15;

COEFFS = cell(3, 1);
COEFFS{1} = mat_M17;
COEFFS{2} = mat_M15;
COEFFS{3} = mat_M16;

p = [1.000e+02, 1.445e+00, 1.670e+02, 0., 0.];
q = [-1.000e+02, 1.111e+00, -1.064e-01];

syms x;
syms f(x) f1(x) f2(x);
lp = length(p); lq = length(q);
f1 = 0;
for i = 1 : lp
    f1 = f1 + p(end - i + 1)*x^(i - 1);
end
f2 = 0;
for i = 1 : lq
    f2 = f2 + q(end - i + 1)*x^(i - 1);
end
f = f1/f2;

FUN3 = cell(40, 1);
FUN3{1} = matlabFunction(f);
for i = 2 : 40
    f = diff(f);
    FUN3{i} = matlabFunction(f);
end

FUN = @(x) allfunc(x, FUN3);

rf = @(x) -(100*x^4 + (289*x^3)/200 + 167*x^2)/(100*x^2 - (1111*x)/1000 + 133/1250);
Amu_func = @(x) COEFFS{1} - x^2*COEFFS{2} + rf(x)*COEFFS{3};

end

function [varargout] = allfunc(lambda, F3)
nout = max(nargout, 1);

[m, n] = size(lambda);
for i = 1 : nout
    varargout{i} = [];
end

for j = 1 : n
    for i = 1 : m
        tmp = lambda(i, j);
        varargout{1} = [varargout{1};[1, -tmp^2, F3{1}(tmp)]];
        if nout > 1
            varargout{2} = [varargout{2};[0, -2*tmp, F3{2}(tmp)]];
        end
        if nout > 2
            varargout{3} = [varargout{3};[0, -2, F3{3}(tmp)]];
        end
        if nout > 3
            for q = 4 : nout
                varargout{q} = [varargout{q};[0, 0, F3{q}(tmp)]];
            end
        end
    end
end

end