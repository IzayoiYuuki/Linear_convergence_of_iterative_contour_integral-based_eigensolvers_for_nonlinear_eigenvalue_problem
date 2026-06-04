function [U, Lambda] = SolveProjectedNEP_complin(rho, center, ls, S, COEFFS, FUN, Sconfig)
% Project T by S'*T*S and solve by Beyn's method

[~, nex] = size(S);

% Project T -> S'TS
for i = 1 : length(COEFFS)
    COEFFS{i} = S'*COEFFS{i}*S;
end

% Solve eig of companion linearized matrix pencil
[L0, L1] = CompLinearization(COEFFS);
[U, Lambda] = eig(full(L0), full(L1), 'vector');

% Companion linearized matrix pencil has p*nex eigs
% pick nex from them which are closest to center of the ellipse contour
idx = inellipse(Lambda, rho, center, ls);
idx = idx(1:nex);
U = U(1:nex, idx); Lambda = Lambda(idx);

% Get projected solution
U = S*U;

end

function [idx] = inellipse(Lambda, rho, center, ls)
    Lambda = Lambda - center;
    % idx = find( (imag(Lambda)/(rho*ls)).^2 + (real(Lambda)/rho).^2 < 1);
    dis = (imag(Lambda)/(rho*ls)).^2 + (real(Lambda)/rho).^2;
    [~, idx] = sort(dis, 'ascend');
end
