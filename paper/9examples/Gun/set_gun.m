% Solver settings
tol_res = 1e-12;
Maxiter = 50;

% Problem set
[COEFFS, FUN, Amu_func] = nlevp('gun');
K = COEFFS{1};
M = COEFFS{2};
W1 = COEFFS{3};
W2 = COEFFS{4};
n = length(K);  

% Contour settings
center = 66762; 
rho = 45738;       
ls = 1;         
N = 32;      

% initial guess
nev = 21;  
nex = ceil(nev*1); 
rng(0); V = randn(n, nex) + 1i*randn(n, nex); 

% linear solver
Sconfig.lu = 1;
Sconfig.Beyn.lu = 0;
Sconfig.Beyn.N = 128;