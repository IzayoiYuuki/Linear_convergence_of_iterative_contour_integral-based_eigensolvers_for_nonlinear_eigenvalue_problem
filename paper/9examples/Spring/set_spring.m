% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
n = 3000;
[COEFFS, FUN, Amu_func] = nlevp('spring', n, 1, 0.6202,  0.4807);
n = length(COEFFS{1});

% Contour settings
alpha = -1.575;
beta = -1.552;
center = (alpha + beta)/2;
rho = (beta - alpha)/2;
ls = 0.5;
N = 8;

% initial guess
nev = 32;
nex = ceil(nev*1);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 32;