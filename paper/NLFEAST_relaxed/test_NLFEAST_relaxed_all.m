clear;
rng(0);
hold off;

TestMatrices = ["spring", "aw2", "butterfly", "ls", "photonic", "r2r", "hadeler", "gun", "pdde"];

Nlist = [8, 16, 16, 16, 32, 0, 8, 64, 32];

SettingFiles = "set_" + TestMatrices;
OutNames = "./figure/appendix/" + TestMatrices + ".txt";
FigureNames = "./figure/appendix/" + TestMatrices + ".pdf";

% warmup
warmup();

for idx = 1 : 9

    try
        run(SettingFiles(idx));
    catch
        disp(TestMatrices(idx) + ' failed!');
        continue;
    end
    N = Nlist(idx);
    Maxiter = 250;
    
    % solve with NLFEAST
    tic;
    [U, Lambda, res, info] = myNLFEAST_relaxed(COEFFS, FUN, Amu_func, V, N, rho, center, ls,...
                                    Maxiter, tol_res, nev, Sconfig);
    info.time(4) = toc;
    
    %% Save data
    writelines(string(datetime('now')), OutNames(idx), WriteMode="append");
    dlmwrite(OutNames(idx), info.time, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(OutNames(idx), res, '-append', 'delimiter', ',', 'precision', 4);
    dlmwrite(OutNames(idx), Lambda.', '-append', 'delimiter', ',', 'precision', 14);
    dlmwrite(OutNames(idx), [center, rho, ls], '-append', 'delimiter', ',', 'precision', 14);

    %% Display results
    disp("Now " + TestMatrices(idx));
    disp("LU decomposition time: " + info.time(1));
    disp("Filter time: " + info.time(2));
    disp("Solving projected problem time: " + info.time(3));
    disp("Total time: " + info.time(4));

    figure(1);
    poltresult(Lambda, center, rho, ls);

    figure(2);
    semilogy(1:length(res), res, '-o', 'LineWidth', 3);

end

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