clear;
rng(0);
hold off;

% plot parameters
marklist = ["-x", "-o", "-<", "-*", "-+", "->"];
colorlist = {[65,105,225]/255, [255,165,0]/255, [255,105,180]/255, [148,0,211]/255, [60,179,113]/255, [128,128,0]/255};

TestMatrices = ["aw2", "photonic", "gun"];
SettingFiles = "set_" + TestMatrices;
OutNames = "./figure/time2accuracy/" + TestMatrices + ".txt";

FullName = {'{\tt{acoustic\_wave\_2d}}', '{\tt{photonics}}',  '{\tt{gun}}'};

gca1 = figure('Units', 'centimeters', 'Position', [0 0 42 25]);
for fileNo = 1 : 3

    %% Load
    FullFile = readlines(OutNames(fileNo));
    
    %% Read NLFEAST
    data_nf = [];
    
    idx = 3;
    line_now = str2num(FullFile(idx));
    while ~isempty(line_now)
        
        data_nf = [data_nf; line_now];
    
        idx = idx + 1;
        line_now = str2num(FullFile(idx));
    end
    
    %% Read Beyn
    data_by = [];
    
    idx = idx + 1;
    line_now = str2num(FullFile(idx));
    while ~isempty(line_now)
        
        data_by = [data_by; line_now];
    
        idx = idx + 1;
        line_now = str2num(FullFile(idx));
    end
    
    %% plot
    subplot(2, 3, fileNo);
    
    % 3 curves of NLFEAST
    idx = find(data_nf(:, 1) == 16);
    semilogy(data_nf(idx, 4), data_nf(idx, 3), marklist(1), 'MarkerIndices', ceil(linspace(1,length(idx), 25)), 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{1}); hold on;
    idx = find(data_nf(:, 1) == 32);
    semilogy(data_nf(idx, 4), data_nf(idx, 3), marklist(2), 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{2}); hold on;
    idx = find(data_nf(:, 1) == 48);
    semilogy(data_nf(idx, 4), data_nf(idx, 3), marklist(3), 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{3}); hold on;
    
    % 1 curve for Beyn
    semilogy(data_by(:, 3), data_by(:, 2), marklist(4), 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{4});
    yline(10^(-12), '--', 'linewidth', 1.5); hold off;
    
    title(FullName{fileNo}, 'Interpreter', 'latex');
    yticks([10^(-12), 10^(-9), 10^(-6), 10^(-3), 1])
    yticklabels({'\(10^{-12}\)', '\(10^{-9}\)', '\(10^{-6}\)', '\(10^{-3}\)', '\(10^{0}\)'})
    endidx = find(data_by(:, 2) > 1e-13);
    xlim([0, max(data_by(endidx(end), 3))]);
    ylim([10^(-13), 1]);
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca, 'FontSize', 20);
    ylabel('Residual', 'Interpreter', 'latex', 'FontSize', 20);
    xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 20);
    set(gca, 'color', 'none');
    

end
% legend('Interpreter', 'latex', 'box', 'off');
    
lgd = legend('NLFEAST with $N=16$', 'NLFEAST with $N=32$', 'NLFEAST with $N=48$', "Beyn's method with \(N=16\)--\(3060\)", "Tolerance",...
            'FontSize', 20, 'Position', [0.27, 0.40, 0.5, 0.05],'Interpreter','latex',...
            'NumColumnsMode', 'manual', 'NumColumns', 2, 'Color','white');

exportgraphics(gca1, './figure/time2accuracy/time2accuracy.pdf', 'BackgroundColor','none');
