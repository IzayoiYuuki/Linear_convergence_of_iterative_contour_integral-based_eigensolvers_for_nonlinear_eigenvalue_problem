clear;
rng(0);

% Solver settings
tol_res = 1e-12;
Maxiter = 20;

% Problem set
[COEFFS, FUN, Amu_func] = nlevp('butterfly');
n = length(COEFFS{1});

% Contour settings
center = 0.5+1i*0.5;
rho = 0.1;
ls = 1;
N = 64;

% linear solver
Sconfig.lu = 0;

% initial guess
nev = 5;
nex = ceil(nev*1);
V = randn(n, nex) + 1i*randn(n, nex);

tic;

[U, Lambda, res] = myNLFEAST(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                    Maxiter, tol_res, nev, Sconfig);
toc;

figure(1);
poltresult(Lambda, center, rho, ls);

figure(2);
semilogy(1:length(res), res, '-o', 'LineWidth', 3);

function poltresult(Lambda, center, rho, ls)

    Np = 10000;
    qNode = (1:Np)';
    lx = rho;
    sx = ls*rho;
    tt = (2 * pi / Np) * qNode - (pi / Np);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;

    plot(z, 'LineWidth', 1.5);
    hold on;
    plot(Lambda, '*', 'LineWidth', 1, 'MarkerSize', 6);
    xlabel('Real', 'Interpreter', 'latex');
    ylabel('Imag', 'Interpreter', 'latex');
    xlim([real(center) - rho, real(center) + rho]);
    ylim([imag(center) - ls*rho, imag(center) + ls*rho]);
    set(gca, 'color', 'none');
    set(gca, 'TickLabelInterpreter', 'latex');
    hold off;

end



