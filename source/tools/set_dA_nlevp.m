function [dA] = set_dA_nlevp(m, center, Ai, FUN)
% set dA of PEP for infGMRES
% (higher derivatives will be filled with 0)

fl = cell(m + 1, 1);
[fl{:}] = FUN(center);

d = length(Ai);
n = length(Ai{1});

dA = cell(m + 1, 1);
for i = 1 : m + 1
    dA{i} = sparse(n, n);
    for j = 1 : d
        dA{i} = dA{i} + fl{i}(j)*Ai{j};
    end
    dA{i} = (1/factorial(i - 1))*dA{i};
end

end