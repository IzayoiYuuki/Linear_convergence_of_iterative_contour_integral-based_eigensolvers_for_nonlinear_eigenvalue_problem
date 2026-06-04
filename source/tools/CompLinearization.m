function [L0, L1] = CompLinearization(COEFFS)

n = length(COEFFS{1});
p = length(COEFFS);

L0 = speye(n*(p-1), n*(p-1));
L1 = sparse(n*(p-1), n*(p-1));

L0(1:n, 1:n) = -COEFFS{1};

for i = 1 : p-1
    L1(1:n, n*(i-1)+1:n*i) = COEFFS{i+1};
end
L1(n+1:end, 1:n*(p-2)) = speye(n*(p-2), n*(p-2));






