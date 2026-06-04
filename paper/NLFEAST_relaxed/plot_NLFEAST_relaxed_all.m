clear;
rng(0);
hold off;

% beautiful colors and markers
marklist = ["-x", "-o", "-<", "-*", "-+", "->"];
colorlist = {[65,105,225]/255, [255,165,0]/255, [255,105,180]/255, [148,0,211]/255, [60,179,113]/255, [128,128,0]/255};

% data
TestMatrices = ["spring", "aw2", "butterfly", "ls", "photonic", "r2r", "hadeler", "gun", "pdde"];
FileNames = "./figure/appendix/" + TestMatrices + ".txt";
FullName = {'{\tt{spring}} (\(8\))', '{\tt{acoustic\_wave\_2d}} (\(16\))', '{\tt{butterfly}} (\(16\))',...
              '{\tt{loaded\_string}} (\(16\))', '{\tt{photonics}} (\(32\))', '{\tt{railtrack2\_rep}} (\(64\))',...
               '{\tt{hadeler}} (\(8\))', '{\tt{gun}} (\(64\))', '{\tt{pdde\_symmetric}} (\(32\))'};

%% Convergence History
gca_conv = figure('Units', 'centimeters', 'Position', [0 0 30 30]);

for fileNo = 1 : 9

    % read from files
    FullFile = readlines(FileNames(fileNo));
    res = str2num(FullFile(3));

    % plot a 3 x 3 figures for 9 examples
    subplot(3, 3, fileNo);
    semilogy(0:length(res)-1, res, marklist(2), 'linewidth', 2, 'MarkerSize', 3, 'Color', colorlist{1});
    yline(10^(-12), '--', 'linewidth', 1.5);
    hold off;
    yticks([10^(-12), 10^(-8), 10^(-4), 1])
    yticklabels({'\(10^{-12}\)', '\(10^{-8}\)', '\(10^{-4}\)', '\(10^{0}\)'})
    ylim([10^(-13), 1]);
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca, 'FontSize', 20);
    title(FullName(fileNo), 'FontSize', 20, 'Interpreter', 'latex')
    ylabel('Residual', 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'color', 'none');

end

exportgraphics(gca_conv, './figure/appendix/appendix.pdf', 'BackgroundColor','none');

function plot_region(Lambda, center, rho, ls)

    Np = 10000;
    qNode = (1:Np)';
    lx = rho;
    sx = ls*rho;
    tt = (2 * pi / Np) * qNode - (pi / Np);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;

    plot(z, 'LineWidth', 1.5, 'Color', [65,105,225]/255);
    hold on;
    plot(Lambda, '*', 'LineWidth', 1.2, 'MarkerSize', 8, 'Color', [255,165,0]/255);
    xlabel('Real', 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('Imag', 'Interpreter', 'latex', 'FontSize', 18);
    xlim([real(center) - rho, real(center) + rho]);
    ylim([imag(center) - ls*rho, imag(center) + ls*rho]);
    set(gca, 'color', 'none');
    set(gca, 'TickLabelInterpreter', 'latex');
    hold off;

end