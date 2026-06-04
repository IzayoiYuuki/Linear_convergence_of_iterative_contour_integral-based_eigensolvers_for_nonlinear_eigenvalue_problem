% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
[COEFFS, FUN, Amu_func] = nlevp('pdde_symmetric');
n = length(COEFFS{1});  

% Contour settings
center = 0; 
rho = 1;       
ls = 1;         
N = 8;      

% initial guess
nev = 6;  
nex = ceil(nev*1); 
rng(0); V = randn(n, nex) + 1i*randn(n, nex); 

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 32;

[U, Lambda, res, info] = myNLFEAST(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                   Maxiter, tol_res, nev, Sconfig);
