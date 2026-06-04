% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
[COEFFS, FUN, Amu_func] = getnlevp_photonic();
n = length(COEFFS{1});
xax = 0.188365;
yax = 0.384419;
center = xax/2 + 1i*yax/2;
rho = xax/2;

% Contour settings
xax = 0.188365;
yax = 0.384419;
center = xax/2 + 1i*yax/2;
rho = xax/2;
ls = yax/xax;
N = 32;

% initial guess
nev = 16;
nex = ceil(16);
rng(0); V = randn(n, nex) + 1i*randn(n, nex);

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 128;