% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
n = 20000;
[COEFFS, FUN, Amu_func] = nlevp('loaded_string', n);
n = length(COEFFS{1});

% Contour settings
center = 1000;
rho = 850;
ls = 0.1;
N = 8;

% initial guess
nev = 10;
nex = ceil(nev*1);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 16;