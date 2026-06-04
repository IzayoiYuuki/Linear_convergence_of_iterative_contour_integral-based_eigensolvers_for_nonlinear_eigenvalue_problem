function [U, Lambda, res, info] = BeynNEP(Amu_func, dAmu_func, V, N, rho, center, ls,...
                                    Maxiter, tol_res, tol_rank, nev,...
                                    Sconfig)
% solving Non-linear Eigenvalue Problems (NEP) with Beyn's method
info = struct();
info.time = zeros(1, 6);
[n, nex] = size(V);

% set integral points and derivatives
[z, omega] = set_contour(rho, center, ls, N);

% prepare LU decomposition if needed
if Sconfig.lu == 1
    lutime = tic;
    [facA] = prep_lu(z, Amu_func);
    info.time(1) = toc(lutime);
else
    facA = 0;
end

for iter = 1 : Maxiter
    
    % Get moment 0 and moment 1
    solvetime = tic;
    [Mom0, Mom1] = GetMoment(z, omega, V, Amu_func, dAmu_func, Sconfig, facA);
    info.time(2) = info.time(2) + toc(solvetime);

    projtime = tic;

    % Estimate the right eigenvectors V
    [V, Sigma, W] = svd(Mom0, 'econ');

    % Numerical rank
    % (for now, we set tol_rank = 1e-14 and do not deflate in any cases)
    ridx = find(diag(Sigma) > tol_rank);
    if length(ridx) ~= nex
        V = V(:, ridx);
        W = W(:, ridx);
        Sigma = Sigma(ridx, ridx);
        nex = length(ridx);
    end
    
    % Solve eig of B, estimate the eigenvalues Lambda
    B = V'*Mom1*W/Sigma;
    [S, Lambda] = eig(B);

    % Approximate eigenvalues and eigenvectors
    Lambda = diag(Lambda) + center;
    U = V*S;

    info.time(3) = info.time(3) + toc(projtime);

    % Check Convergence
    res_vec = zeros(nex, 1);
    for idx = 1 : nex
        res_vec(idx) = norm(Amu_func(Lambda(idx))*U(:, idx))/...
            (normest(Amu_func(Lambda(idx)), 1e-1)*norm(U(:, idx)));
    end
    [tmp, tmpidx] = mink(res_vec, nev);
    U = U(:, tmpidx);
    Lambda = Lambda(tmpidx);
    info.res_all = res_vec;
    res(iter) = tmp(end);
    if res(iter) < tol_res
        break;
    end

end

end

% function [z, omega] = set_contour(rho, center, ls, N)
% % set the integral points z and its derivatives omega
%     lx = rho;
%     sx = ls * rho;
%     qNode = (1:N)';
%     tt = (2 * pi / N) * qNode - (pi / N);
%     z =  lx*cos(tt) + 1i*sx*sin(tt) + center;
%     omega = (-lx*sin(tt) + 1i*sx*cos(tt))*2*pi/N;
% end

function [z, omega] = set_contour(rho, center, ls, N)
% set the integral points z and its derivatives omega
    lx = rho;
    sx = ls * rho;
    qNode = (1:N)';
    tt = (2 * pi / N) * qNode - (pi / N);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;
    omega = (-lx*sin(tt) + 1i*sx*cos(tt))*2*pi/N;
end


