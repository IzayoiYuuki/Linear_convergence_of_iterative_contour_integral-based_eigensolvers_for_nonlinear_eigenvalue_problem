clear;
rng(0);
hold off;

% plot parameters
marklist = ["-x", "-o", "-<", "-*", "-+", "->"];
colorlist = {[65,105,225]/255, [255,165,0]/255, [255,105,180]/255, [148,0,211]/255, [60,179,113]/255, [128,128,0]/255};

%% Load
FullFile = readlines('./figure/intro/Data_intro.txt');

%% figure 1: Beyn's method, 10 iter, 64, 128, 256, 512 qNodes
figure('Units', 'centimeters', 'Position', [0 0 19 14]);
for i = 1 : 4
    res = str2num(FullFile(i + 2));
    semilogy(0:length(res), [1, res], marklist{i}, 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{i}); hold on;
end
hold off;

title("Beyn's method", 'Interpreter', 'latex')
legend('$N=64$', '$N=128$', '$N=256$', '$N=512$', 'Interpreter', 'latex', 'box', 'off');
xlabel('Iterations', 'Interpreter', 'latex');
xlim([0, 10]);
ylabel('Residual', 'Interpreter', 'latex');
yticks([10^(-16), 10^(-12), 10^(-8), 10^(-4), 1])
yticklabels({'\(10^{-16}\)', '\(10^{-12}\)', '\(10^{-8}\)', '\(10^{-4}\)', '\(10^{0}\)'})
ylim([10^(-16), 1]);
set(gca, 'color', 'none');
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);

% slightly move legend up
h = legend;
p = h.Position;
p(2) = p(2) + 0.02;
h.Position = p;

exportgraphics(gca, './figure/intro/intro1-1.pdf', 'BackgroundColor','none');

%% figure 2: region of interest and eigenvalues to be computed
Lambda_true = [
0.269116796917073 + 0.236990802383966i
0.304852019949294 + 0.220448968829496i
0.284829383301611 + 0.255205421896188i
0.322139826088162 + 0.240048282456614i
0.364150108908551 + 0.188363837242100i
0.306735530841174 + 0.285466354068229i
0.385239732817062 + 0.211448611690893i
0.346363345200310 + 0.272800818668618i
0.330110365868871 + 0.326877051559231i
0.460677802232896 + 0.126426476945846i
0.414337427372999 + 0.250379039094508i
0.372405194621209 + 0.317743029889994i
0.350625139549668 + 0.377179177618922i
0.490854516544653 + 0.162957609479054i
0.445077822903038 + 0.303605915468083i
0.395637467268562 + 0.372346046766981i
0.365202259406344 + 0.432142315147135i
0.525615028746292 + 0.221195148371234i
0.412720094633443 + 0.432008118258471i
0.472414292856460 + 0.367772850997250i
0.372778266450175 + 0.484565226949958i
0.558411920025264 + 0.293560342939760i
0.375005942907254 + 0.523831272772491i
0.422238258332365 + 0.488878601898202i
0.492747300357317 + 0.437711328232297i
0.425581765009511 + 0.531357920644182i
0.586838033752989 + 0.376181125229503i
0.689173106766718 + 0.134763075535075i
0.685616268417965 + 0.175362345453074i
0.504161647567241 + 0.504468476873191i
0.691914824016340 + 0.241545259034048i
0.508021566585071 + 0.554368392875770i
0.607319473462600 + 0.464835005916696i
0.708452708786816 + 0.326781061408727i
0.617186814982521 + 0.549454958734074i
0.804864028131497 + 0.220828026737011i
0.727592288862459 + 0.427243632816257i
0.810178271703325 + 0.255622608575701i
0.618460715650399 + 0.613056604489285i
0.820755730726838 + 0.316437316393255i
0.741163194143420 + 0.537745550182735i
0.836379010253795 + 0.404087655723378i
0.885260989817921 + 0.282262475155174i
0.894552041898822 + 0.319324937764334i
0.744283786331639 + 0.646541111813964i
0.853220341698646 + 0.517099955869046i
0.738844823886403 + 0.731658534006025i
0.864617980453660 + 0.651815654480523i];
Lambda = str2num(FullFile(7));
Region = str2num(FullFile(2));

figure('Units', 'centimeters', 'Position', [0 0 19 14]);
plot(Lambda_true, '*', 'LineWidth', 1, 'MarkerSize', 10, 'Color', colorlist{2}); hold on;
plot(Lambda, 'o', 'LineWidth', 1.2, 'MarkerSize', 10, 'Color', colorlist{4});
plot_region(Region(1), Region(2), Region(3));

title('\texttt{butterfly}', 'Interpreter', 'latex')
xlabel('Real', 'Interpreter', 'latex');
ylabel('Imag', 'Interpreter', 'latex');
xlim([0, 0.8])
ylim([0, 0.6])
set(gca, 'color', 'none');
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);
legend('Eigenvalues', 'Eigenvalues of interest', 'Contour', 'Interpreter', 'latex', 'Location','southwest', 'box', 'off');
hold off;

exportgraphics(gca, './figure/intro/intro1-2.pdf', 'BackgroundColor','none');

%% figure 3: NLFEAST 10 iter, 64, 128, 256, 512 qNodes
figure('Units', 'centimeters', 'Position', [0 0 19 14]);
for i = 1 : 4
    res = str2num(FullFile(i + 7));
    semilogy(0:length(res), [1, res], marklist{i}, 'LineWidth', 3, 'MarkerSize', 8, 'Color', colorlist{i}); hold on;
end
hold off;

title('NLFEAST', 'Interpreter', 'latex')
legend('$N=64$', '$N=128$', '$N=256$', '$N=512$', 'Interpreter', 'latex', 'box', 'off');
xlabel('Iterations', 'Interpreter', 'latex');
xlim([0, 10]);
ylabel('Residual', 'Interpreter', 'latex');
yticks([10^(-16), 10^(-12), 10^(-8), 10^(-4), 1])
yticklabels({'\(10^{-16}\)', '\(10^{-12}\)', '\(10^{-8}\)', '\(10^{-4}\)', '\(10^{0}\)'})
ylim([10^(-16), 1]);
set(gca, 'color', 'none');
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 20);

% slightly move legend up
h = legend;
p = h.Position;
p(2) = p(2) + 0.02;
h.Position = p;

exportgraphics(gca, './figure/intro/intro2.pdf', 'BackgroundColor','none');


%% functions
function plot_region(center, rho, ls)
% draw ellipse region of interest
    Np = 10000;
    qNode = (1:Np)';
    lx = rho;
    sx = ls*rho;
    tt = (2 * pi / Np) * qNode - (pi / Np);
    z =  lx*cos(tt) + 1i*sx*sin(tt) + center;
    plot(z, 'LineWidth', 2, 'Color', [65,105,225]/255); hold on;
end