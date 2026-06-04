% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
n = 5000;
[COEFFS, FUN, Amu_func] = nlevp('butterfly', n);
n = length(COEFFS{1});

% Contour settings
center = 0.25+1i*0.236;
rho = 0.0035;
ls = 1;
N = 16;

% initial guess
nev = 9;
nex = ceil(nev*1);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 256;