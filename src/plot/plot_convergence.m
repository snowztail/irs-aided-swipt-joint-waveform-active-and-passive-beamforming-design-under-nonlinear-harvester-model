clear; clc; close all; config_convergence;

%% * Retrieve data
% * Load data and select the WPT point
load('../data/convergence.mat');
sampleIndex = nSamples;

% * Algorithm 1
scaIter = aoSolution{sampleIndex}.scaIter{1}(2, :);
scaIter(end : nIters) = scaIter(end);

% * Algorithm 2
gpIter = aoSolution{sampleIndex}.gpIter{1}(2, :);
gpIter(end : nIters) = gpIter(end);

% * Algorithm 3
mScaIter = lcSolution{sampleIndex}.mScaIter{1}(2, :);
mScaIter(end : nIters) = mScaIter(end);

% * Algorithm 4
bcdIter = aoSolution{sampleIndex}.bcdIter(2, :);
bcdIter(end : nIters) = bcdIter(end);

% * Algorithm 5
lcBcdIter = lcSolution{sampleIndex}.lcBcdIter(2, :);
lcBcdIter(end : nIters) = lcBcdIter(end);

%% * Convergence trends of SCA-based algorithms
figure('name', 'Convergence trends of SCA-based algorithms', 'position', [0, 0, 500, 400]);
scaPlot = tiledlayout(2, 1, 'tilespacing', 'loose');

% * Algorithm 1
nexttile;
plot(1e6 * scaIter, 'color', '#0072BD', 'linestyle', '-', 'marker', 'o');
grid on;
legend('SCA', 'location', 'se');
ylabel('DC current [$\mu$A]');
xlim([1, nIters]);
xticks([1, 5 : 5 : nIters]);
box on;

% * Algorithm 3
nexttile;
plot(1e6 * mScaIter, 'color', '#D95319', 'linestyle', '--', 'marker', '+');
grid on;
legend('M-SCA', 'location', 'se');
xlabel('Iteration index');
ylabel('DC current [$\mu$A]');
xlim([1, nIters]);
xticks([1, 5 : 5 : nIters]);
box on;

savefig('../figures/convergence_sca.fig');
matlab2tikz('../../assets/convergence_sca.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
close;

%% * Convergence trends of GP-based algorithm
figure('name', 'Convergence trends of GP-based algorithm', 'position', [0, 0, 500, 400]);

% * Algorithm 2
plot(1e6 * gpIter, 'color', '#0072BD', 'linestyle', '-', 'marker', 'o');
grid on;
legend('GP', 'location', 'se');
xlabel('Iteration index');
ylabel('DC current [$\mu$A]');
xlim([1, nIters]);
xticks([1, 5 : 5 : nIters]);
box on;

savefig('../figures/convergence_gp.fig');
matlab2tikz('../../assets/convergence_gp.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
close;

%% * Convergence trends of BCD-based algorithms
figure('name', 'Convergence trends of BCD-based algorithms', 'position', [0, 0, 500, 400]);
bcdPlot = tiledlayout(2, 1, 'tilespacing', 'loose');

% * Algorithm 4
nexttile;
plot(1e6 * bcdIter, 'color', '#0072BD', 'linestyle', '-', 'marker', 'o');
grid on;
legend('BCD', 'location', 'se');
ylabel('DC current [$\mu$A]');
xlim([1, nIters]);
xticks([1, 5 : 5 : nIters]);
box on;

% * Algorithm 5
nexttile;
plot(1e6 * lcBcdIter, 'color', '#D95319', 'linestyle', '--', 'marker', '+');
grid on;
legend('LC-BCD', 'location', 'se');
xlabel('Iteration index');
ylabel('DC current [$\mu$A]');
xlim([1, nIters]);
xticks([1, 5 : 5 : nIters]);
box on;

savefig('../figures/convergence_bcd.fig');
matlab2tikz('../../assets/convergence_bcd.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
