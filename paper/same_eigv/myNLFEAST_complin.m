function [Vout, Lambdaout, res, info] = myNLFEAST_complin(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                    Maxiter, tol_res, nev, Sconfig)
% solving Non-linear Eigenvalue Problems (NEP) with Beyn's method
info = struct();
info.time = zeros(1, 4);
total_time = tic;

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

% prepare norm for checking convergence
for i = 1 : length(COEFFS)
    coeff_norms(i, 1) = normest(COEFFS{i}, 1e-1);
end

% Obtain approximations of U and Lambda from initial guess V
[V, Lambda] = SolveProjectedNEP_complin(rho, center, ls, V, COEFFS, FUN, Sconfig);
res_vec = getres(V, Lambda, coeff_norms, COEFFS, FUN);
[tmp, tmpidx] = mink(res_vec, nev);
Vout = V(:, tmpidx);
Lambdaout = Lambda(tmpidx);
res(1) = tmp(end);
info.total_time(1) = toc(total_time);

for iter = 1 : Maxiter
    
    % Reverse iteration like filter of NLFEAST
    solvetime = tic;
    V = NLFIter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA);
    info.time(2) = info.time(2) + toc(solvetime);

    % Solve projection problem V'*T*V
    projtime = tic;
    [V, Lambda] = SolveProjectedNEP_complin(rho, center, ls, V, COEFFS, FUN, Sconfig);
    info.time(3) = info.time(3) + toc(projtime);

    % Check Convergence
    res_vec = getres(V, Lambda, coeff_norms, COEFFS, FUN);
    [tmp, tmpidx] = mink(res_vec, nev);
    Vout = V(:, tmpidx);
    Lambdaout = Lambda(tmpidx);
    res(iter + 1) = tmp(end);
    info.res_all = res_vec;
    info.total_time(iter + 1) = toc(total_time);
    if res(iter + 1) < tol_res
        break;
    end

end

end

function [z, omega] = set_contour(rho, center, ls, N)
% set the integral points z and its derivatives omega
    lx = rho;
    sx = ls * rho;
    qNode = (1:N)';
    tt = (2 * pi / N) * qNode - (pi / N);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;
    omega = (-lx*sin(tt) + 1i*sx*cos(tt))*2*pi/N;
end

function [Vout] = NLFIter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA)
% A single iteration of filter of NLFEAST
    [n, nex] = size(V);
    Vout = zeros(n, nex);
    
    TV = getTv(V, Lambda, COEFFS, FUN);
    
    if Sconfig.lu == 0
        % slash
        for intp = 1 : length(z)
            Vtmp = (V - Amu_func(z(intp))\TV)/diag(z(intp) - Lambda);
            Vout = Vout + omega(intp)*Vtmp;
        end
    else
        % sparse direct with LU decompostion
        for intp = 1 : length(z)
            Vtmp = (V - facA{intp}\TV)/diag(z(intp) - Lambda);
            Vout = Vout + omega(intp)*Vtmp;
        end
    end
    Vout = Vout/(2*pi*1i);

end

function [Tv] = getTv(V, Xi, COEFFS, FUN)
% compute T(Xi, V) = [T(xi_1)v_1, T(xi_2)v_2, ..., T(xi_k)v_k]
    Tv = zeros(size(V));
    funtmp = FUN(Xi);
    for i = 1 : length(COEFFS)
        Tv = Tv + COEFFS{i}*V*diag(funtmp(:, i));
    end
end

function [res_vec] = getres(V, Lambda, coeff_norms, COEFFS, FUN)
% estimate relevant residuals
    nex = size(V, 2);
    res_vec = zeros(nex, 1);
    Tv = getTv(V, Lambda, COEFFS, FUN);
    for idx = 1 : nex
        res_vec(idx) = norm(Tv(:, idx));
        res_vec(idx) = res_vec(idx)/(abs(FUN(Lambda(idx)))*coeff_norms*norm(V(:, idx)));
    end
end



