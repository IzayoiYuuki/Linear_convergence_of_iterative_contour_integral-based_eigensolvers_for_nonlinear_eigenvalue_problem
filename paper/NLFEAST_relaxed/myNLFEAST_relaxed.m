function [Vout, Lambdaout, res, info] = myNLFEAST_relaxed(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
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
[V, Lambda] = SolveProjectedNEP(rho, center, ls, V, COEFFS, FUN, Sconfig);
res_vec = getres(V, Lambda, coeff_norms, COEFFS, FUN);
[tmp, tmpidx] = mink(res_vec, nev);
Vout = V(:, tmpidx);
Lambdaout = Lambda(tmpidx);
res(1) = tmp(end);
info.total_time(1) = toc(total_time);

for iter = 1 : Maxiter
    
    % Reverse iteration like filter of NLFEAST
    solvetime = tic;
    if iter > 0
        V = NLFAIter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA);
    else
        V = NLFIter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA);
    end
    info.time(2) = info.time(2) + toc(solvetime);

    % Solve projection problem V'*T*V
    projtime = tic;
    [V, ~] = qr(V, 0);
    [V, Lambda] = SolveProjectedNEP(rho, center, ls, V, COEFFS, FUN, Sconfig);
    info.time(3) = info.time(3) + toc(projtime);

    % Check Convergence
    res_vec = getres(V, Lambda, coeff_norms, COEFFS, FUN);
    [tmp, tmpidx] = mink(res_vec, nev);
    Vout = V(:, tmpidx);
    Lambdaout = Lambda(tmpidx);
    res(iter + 1) = tmp(end);
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

function [Vout] = NLFAIter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA)
% A single iteration of filter of NLFEAST with Newton and Arnoldi
% solve T(\lambda)\dT(\lambda)x in S = {x, PT(\lambda)x, (PT(\lambda))^2x, ...}
% that is: solve min ||T(\lambda)Sy - dT(\lambda)x|| 
    
    % dimension of Arnoldi subspace is (mA + 1)
    mA = 1;
    [n, nex] = size(V);
    Vout = zeros(n, nex);
    dTV = getdTv(V, Lambda, COEFFS, FUN);
    
    % S_all = {V, PT(\Lambda, V), PT(\Lambda, T(\Lambda, V)), ...}
    S_all{1} = V;
    for j = 1 : mA
        S_all{j+1} = getPTv(z, omega, S_all{j}, Lambda, COEFFS, FUN, facA); % Checked
    end

    % % Keep this and comment out all the following things
    % % we will get NLFEAST
    % Vout = V - S_all{2};

    % for i = 1 : nex
    %     % Generate search space S
    %     S = V(:, i);
    %     for j = 1 : mA
    %         S = [S, S_all{j+1}(:, i)];
    %     end
    % 
    %     % [S, ~] = qr(S, 0);
    % 
    %     % solve T(\lambda)Sy = dT(\lambda)x (Least square)
    %     TS = getTv(S, ones(1, size(S, 2))*Lambda(i), COEFFS, FUN);
    % 
    %     [~, ~, Vp] = svd(TS);
    %     y = Vp(:, end);
    %     % y(2) = y(2) + 0.1; y = y/norm(y);
    % 
    %     % Approximate solution
    %     Vout(:, i) = S*y;
    % end

    for i = 1 : nex
        % Generate search space S
        S = V(:, i);
        for j = 1 : mA
            S = [S, S_all{j+1}(:, i)];
        end

        % solve T(\lambda)Sy = dT(\lambda)x (Least square)
        TS = getTv(S, ones(1, size(S, 2))*Lambda(i), COEFFS, FUN);
        y = TS\dTV(:, i);

        % Approximate solution
        Vout(:, i) = S*y;
    end

    % % Just Newton
    % for i = 1 : nex
    %     Vout(:, i) = Amu_func(Lambda(i))\dTV(:, i);
    % end

end

function [Tv] = getTv(V, Xi, COEFFS, FUN)
% compute T(Xi, V) = [T(xi_1)v_1, T(xi_2)v_2, ..., T(xi_k)v_k]
    Tv = zeros(size(V));
    funtmp = FUN(Xi);
    for i = 1 : length(COEFFS)
        Tv = Tv + COEFFS{i}*V*diag(funtmp(:, i));
    end
end

function [dTV] = getdTv(V, Xi, COEFFS, FUN)
% compute dT(Xi, V) = [T'(xi_1)v_1, T'(xi_2)v_2, ..., T(xi_k)'v_k]
    dTV = zeros(size(V));
    funtmp = cell(2, 1);
    [funtmp{:}] = FUN(Xi);
    funtmp = funtmp{2};
    for i = 1 : length(COEFFS)
        dTV = dTV + COEFFS{i}*V*diag(funtmp(:, i));
    end
end

function [PTv] = getPTv(z, omega, V, Xi, COEFFS, FUN, facA)
% compute PT(Xi, V) = P[T(xi_1)v_1, T(xi_2)v_2, ..., T(xi_k)v_k]
    PTv = zeros(size(V));    
    Tv = getTv(V, Xi, COEFFS, FUN);
    for intp = 1 : length(z)
        Vtmp = (facA{intp}\Tv)/diag(z(intp) - Xi);
        PTv = PTv + omega(intp)*Vtmp;
    end
    PTv = PTv/(2*pi*1i);
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


%% By Nian
function y = precond(v,facA,z,rho,omega)
y = zeros(size(v));
for intp = 1 : length(z)
    vtmp = (facA{intp}\v)/((z(intp) - rho));
    y = y + omega(intp)*vtmp;
    % y = y + omega(intp)*(facA{intp}\v);
    
end
end

function y = my_gmres(A,q,b,m)
n = size(q,1);
Q = zeros(n,m+1);
H = zeros(m+1,m);
q = q/norm(q);
Q(:,1) = q;
for iter = 1:m
    norm(Q(:,iter))
    q = A(Q(:,iter));
    norm(q)
    if norm(q)<1e-9
        disp('NB')
        iter = iter-1;
        break
    end
    h = Q(:,1:iter)'*q;
    q = q-Q(:,1:iter)*h;
    dh = Q(:,1:iter)'*q;
    q = q-Q(:,1:iter)*dh;
    H(1:iter,iter) = h+dh;
    H(iter+1,iter) = norm(q);
    Q(:,iter+1) = q/H(iter+1,iter);
end
if iter==0
    y = Q(:,1);
else
    y = Q(:,1:iter)*(H(1:iter+1,1:iter)\(Q(:,1:iter+1)'*b));
end
% norm(A(y)-b)
end

function [Vout] = NLFA2Iter(z, omega, V, Lambda, COEFFS, FUN, Amu_func, Sconfig, facA)
% same as NLFAIter but with gmres
    [n, nex] = size(V);
    Vout = zeros(n, nex);
    
    dTV = getdTv(V, Lambda, COEFFS, FUN);

    m = 3;
    for ii = 1:nex
        b = precond(dTV(:,ii),facA,z,Lambda(ii),omega);
        Vout(:,ii) = my_gmres(@(v) precond(getTv(v, Lambda(ii), COEFFS, FUN),facA,z,Lambda(ii),omega),V(:,ii),b,m);
    end
end