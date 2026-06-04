clear;

% Solver settings
tol_res = 1e-12;
Maxiter = 20;

% Problem set
COEFFS{1} = [0  12; -2 14]; COEFFS{1} = blkdiag(COEFFS{1}, 0);
COEFFS{2} = [-1 -6;  2 -9]; COEFFS{2} = blkdiag(COEFFS{2}, 0);
COEFFS{3} = [1   0;  0  1]; COEFFS{3} = blkdiag(COEFFS{3}, 1);
Amu_func = @(x) COEFFS{1} + x*COEFFS{2} + x^2*COEFFS{3};
n = length(COEFFS{1});
dAmu_func = @(lambda) eye(n);

% Accurate solutions
[L0, L1] = CompLinearization(COEFFS);
[U, Lambda] = eig(full(L0), full(L1), 'vector');
U_true = U(1:n, :); Lambda_true = Lambda;

% Contour settings
center = 1.5;
rho = 1;
ls = 1;
N = 512;

% initial guess
nev = 2;
nex = ceil(nev*1);
% rng(1); V = randn(n, nex);
V = eye(n, nex);

% linear solver
Sconfig.lu = 0;

Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 256;

% save file
outname = './figure/same_eigv/Data_same_eigv.txt';
writelines(string(datetime('now')), outname, WriteMode="append");

%% test NLFEAST and Beyn on 3 contours B(1.5, 1) B(2.5, 1) B(3.5, 1)

for i = 1 : 3
[~, Lambda, ~, info] = myNLFEAST_complin(COEFFS, @FUN, Amu_func, V, N, rho, center, ls,...
                                    Maxiter, tol_res, nev, Sconfig);
dlmwrite(outname, Lambda.', '-append', 'delimiter', ',', 'precision', 16);
dlmwrite(outname, info.res_all.', '-append', 'delimiter', ',', 'precision', 4);

[~, Lambda, ~, info] = BeynNEP(Amu_func, dAmu_func, V, N, rho, center, ls, 1, 0, 0, nev, Sconfig);
dlmwrite(outname, Lambda.', '-append', 'delimiter', ',', 'precision', 16);
dlmwrite(outname, info.res_all.', '-append', 'delimiter', ',', 'precision', 4);

center = center + 1;
end

function [varargout] = FUN(lambda)
    nout = max(nargout, 1);
    
    [m, n] = size(lambda);
    for i = 1 : nout
        varargout{i} = [];
    end
    
    for j = 1 : n
        for i = 1 : m
            tmp = lambda(i, j);
            varargout{1} = [varargout{1};[1, tmp, tmp^2]];
            if nout > 1
                varargout{2} = [varargout{2};[0, 1, 2*tmp]];
            end
            if nout > 2
                varargout{3} = [varargout{3};[0, 0, 2]];
            end
            if nout > 3
                for q = 4 : nout
                    varargout{q} = [varargout{q};[0, 0, 0]];
                end
            end
        end
    end

end