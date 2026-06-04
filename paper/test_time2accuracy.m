clear;
rng(0);
hold off;

TestMatrices = ["aw2", "photonic", "gun"];
SettingFiles = "set_" + TestMatrices;
OutNames = "./figure/time2accuracy/" + TestMatrices + ".txt";

% number of quadrature nodes
N_NLFEAST = 16:16:48;
N_Beyn{1} = 16:8:360; % for aw2
N_Beyn{2} = 16:64:3060; % for gun
N_Beyn{3} = 16:16:720; % for photonics

% warmup
warmup();

for idx = 1 : 3

    % output address
    outname = OutNames(idx);
    writelines(string(datetime('now')), outname, WriteMode="append");
    
    % load problem
    try
        run(SettingFiles(idx));
    catch
        disp(TestMatrices(idx) + ' failed!');
        continue;
    end
    
    % by NLFEAST
    writelines("%%%%%%%%% NLFEAST %%%%%%%%%", outname, WriteMode="append");
    for N = N_NLFEAST
        [~, ~, res, info] = myNLFEAST(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                           150, tol_res, nev, Sconfig);
    
        for iter = 1 : length(res)
            dlmwrite(outname, [N, iter, res(iter), info.total_time(iter)], '-append', 'delimiter', ',', 'precision', 4);
        end 
    end
    
    % by Beyn
    % some additional settings
    tol_rank = 1e-4;
    dAmu_func = @(lambda) speye(n);
    Sconfig.lu = 0;
    Maxiter = 1;
    
    writelines("%%%%%%%%% Beyn's method %%%%%%%%%", outname, WriteMode="append");
    for N = N_Beyn{idx}
        tic;
        [~, ~, res, ~] = BeynNEP(Amu_func, dAmu_func, V, N, rho, center, ls,...
                               Maxiter, tol_res, tol_rank, nev, Sconfig);
        runtime = toc;
    
        dlmwrite(outname, [N, min(res), runtime], '-append', 'delimiter', ',', 'precision', 4);
    end

end