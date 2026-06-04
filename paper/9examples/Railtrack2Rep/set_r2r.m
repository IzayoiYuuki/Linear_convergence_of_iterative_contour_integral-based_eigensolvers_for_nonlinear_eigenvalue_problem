% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
[COEFFS, FUN, Amu_func, XCOEFFS] = nlevp('railtrack2_rep', 705*51, 1000);
n = length(COEFFS{1});

% Contour settings
center = 1;
rho = 0.3;
ls = 1;
N = 16;

% initial guess
nev = 2;
nex = ceil(2);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 16;