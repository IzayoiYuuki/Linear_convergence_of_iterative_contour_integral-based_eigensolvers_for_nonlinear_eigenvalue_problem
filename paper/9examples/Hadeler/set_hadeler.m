% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
n = 5000;
alpha = 100;
[COEFFS, FUN, Amu_func] = nlevp('hadeler', n, alpha);
n = length(COEFFS{1});

% Contour settings
center = 0.01;
rho = 0.001;
ls = 0.1;
N = 8;

% initial guess
nev = 13;
nex = ceil(nev*1);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 16;