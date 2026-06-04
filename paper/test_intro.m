clear;
rng(0);


%% Problem settings
% Problem
% n = 5041;
% [COEFFS, FUN, Amu_func] = nlevp('butterfly', n);
[COEFFS, FUN, Amu_func] = nlevp('butterfly');
n = length(COEFFS{1});
% Contour
center = 0.35+1i*0.25;
rho = 0.1;
ls = 1;
% initial guess
nev = 11;
nex = ceil(nev*1);
rng(0);
V = randn(n, nex) + 1i*randn(n, nex);

%% Solver settings
% Residual
tol_res = 0;
tol_rank = 1e-4;
% dAmu_func(\xi) = T'(\xi), use for Beyn: \int T(\xi)\(T'(\xi)V);
dAmu_func = @(lambda) speye(n);
% store LU decomposition?
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 256;
% Max iteration
Maxiter = 10;

%% file settings
% Output address
outname = './figure/intro/Data_intro.txt';
% mark date
writelines(string(datetime('now')), outname, WriteMode="append");
% save region
dlmwrite(outname, [center, rho, ls], '-append', 'delimiter', ',', 'precision', 16);

%% Test 1: N = 64 - 512, Beyn's method with 10 iterations
for i = 2 : 5
    N = 2^i*16;
    [U, Lambda, res, ~] = BeynNEP(Amu_func, dAmu_func, V, N, rho, center, ls,...
                           Maxiter, tol_res, tol_rank, nev, Sconfig);
    dlmwrite(outname, res, '-append', 'delimiter', ',', 'precision', 4);
end
dlmwrite(outname, Lambda.', '-append', 'delimiter', ',', 'precision', 16);

%% Test 2: N = 64 - 512, NLFEAST with 10 iterations
for i = 2 : 5
    N = 2^i*16;
    [U, Lambda, res, info] = myNLFEAST(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                       Maxiter, tol_res, nev, Sconfig);
    dlmwrite(outname, res, '-append', 'delimiter', ',', 'precision', 4);
end