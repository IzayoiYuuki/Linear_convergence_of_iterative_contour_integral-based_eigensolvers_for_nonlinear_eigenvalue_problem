% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
n = 9900;
z = 1;
[COEFFS, FUN, Amu_func] = nlevp('acoustic_wave_2d', n, z);
n = length(COEFFS{1});

% Contour settings
center = 0;
rho = 1.65;
ls = 0.3;
N = 8;

% initial guess
nev = 10;
nex = ceil(nev);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 256;
