function [facA] = prep_lu(z, Amu_func)
% prepare lu decomposition for further use

% allocate memory
n = length(Amu_func(z(1)));
facA = cell(length(z), 1);

% for each integral point, do lu decomposition
for intp = 1 : length(z)
    facA{intp} = decomposition(Amu_func(z(intp)));
end