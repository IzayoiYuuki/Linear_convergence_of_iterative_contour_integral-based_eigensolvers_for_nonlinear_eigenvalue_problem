clear;
rng(0);
hold off;

% plot parameters
marklist = ["-x", "-o", "-<", "-*", "-+", "->"];
colorlist = {[65,105,225]/255, [255,165,0]/255, [255,105,180]/255, [148,0,211]/255, [60,179,113]/255, [128,128,0]/255};

%% Load
Lambda = [0, 1, 2, 3, 4];
Lambda_same = [3, 4];

%% Read NLFEAST
figure('Units', 'centimeters', 'Position', [0 0 19 14]);

plot_region(1.5, 1, 1, [30,144,255]/255); hold on;
plot_region(2.5, 1, 1, [34,139,34]/255);
plot_region(3.5, 1, 1, [205,92,92]/255);
plot(real(Lambda), imag(Lambda), '*', 'LineWidth', 1, 'MarkerSize', 10, 'Color', colorlist{2});

text(0-0.1, 0.2, '\(\lambda_0\)', 'Interpreter', 'latex', 'FontSize', 20);
text(1-0.1, 0.2, '\(\lambda_1\)', 'Interpreter', 'latex', 'FontSize', 20);
text(2-0.1, 0.2, '\(\lambda_2\)', 'Interpreter', 'latex', 'FontSize', 20);
text(3-0.1, 0.2, '\(\lambda_3\)', 'Interpreter', 'latex', 'FontSize', 20);
text(4-0.1, 0.2, '\(\lambda_4\)', 'Interpreter', 'latex', 'FontSize', 20);

title('Test contours', 'Interpreter', 'latex')
xlabel('Real', 'Interpreter', 'latex');
ylabel('Imag', 'Interpreter', 'latex');
xlim([-0.5, 4.5])
ylim([-2.5, 1.5])
xticks([0, 1, 2, 3, 4, 5]);
yticks([-2, -1, 0 ,1]);
set(gca, 'color', 'none');
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);
legend('Contour \(\mathcal{C}_{1}\)', 'Contour \(\mathcal{C}_{2}\)', 'Contour \(\mathcal{C}_{3}\)', 'Eigenvalues', 'Interpreter', 'latex', 'Location','south', 'box', 'off', 'NumColumnsMode', 'manual', 'NumColumns', 2);
hold off;

exportgraphics(gca, './figure/same_eigv/same_eigv.pdf', 'BackgroundColor','none');

function plot_region(center, rho, ls, mycolor)
% draw ellipse region of interest
    Np = 10000;
    qNode = (1:Np)';
    lx = rho;
    sx = ls*rho;
    tt = (2 * pi / Np) * qNode - (pi / Np);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;
    plot(z, 'LineWidth', 2, 'Color', mycolor); hold on;
end