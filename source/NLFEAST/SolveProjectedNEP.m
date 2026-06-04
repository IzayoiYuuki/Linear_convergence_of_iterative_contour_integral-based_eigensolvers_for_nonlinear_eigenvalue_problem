function [U, Lambda] = SolveProjectedNEP(rho, center, ls, S, COEFFS, FUN, Sconfig)
% Project T by S'*T*S and solve by Beyn's method

[~, nev] = size(S);

% make a Amu_func
for i = 1 : length(COEFFS)
    COEFFS{i} = S'*COEFFS{i}*S;
end
Amu_func = @(lambda) tmpfunc(lambda, COEFFS, FUN);
dAmu_func = @(lambda) speye(nev);

% random initial guess
V = randn(nev, nev) + 1i*randn(nev, nev);

% Beyn method
[U, Lambda] = BeynNEP(Amu_func, dAmu_func, V, Sconfig.Beyn.N, rho, center, ls, 1, 0, 0, nev, Sconfig.Beyn);
U = S*U;

end

function [Amu] = tmpfunc(mu, COEFFS, FUN)
    fmu = FUN(mu);
    Amu = fmu(1)*COEFFS{1};
    for i = 2 : length(fmu)
        Amu = Amu + fmu(i)*COEFFS{i};
    end
end