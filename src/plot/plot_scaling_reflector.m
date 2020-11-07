clear; clc; close all; config_scaling_reflector;

%% * Load batch data
indexSet = 1 : nBatches;
rateSet = zeros(nBatches, length(Variable.nReflectors));
currentLinearSet = zeros(nBatches, length(Variable.nReflectors));
currentNonlinearSet = zeros(nBatches, length(Variable.nReflectors));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/scaling_reflector/scaling_reflector_%d.mat', iBatch), 'rateInstance', 'currentLinearInstance', 'currentNonlinearInstance');
        rateSet(iBatch, :) = rateInstance;
        currentLinearSet(iBatch, :) = currentLinearInstance;
        currentNonlinearSet(iBatch, :) = currentNonlinearInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
rate = mean(rateSet(indexSet, :), 1);
currentLinear = mean(currentLinearSet(indexSet, :), 1);
currentNonlinear = mean(currentNonlinearSet(indexSet, :), 1);
snrDb = pow2db(2 .^ (rate / nSubbands));
save('../data/scaling_reflector.mat');

%% * SNR and power plots
figure('name', 'Average SNR and harvested DC power vs number of reflectors');
pathlossPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

% * SNR plot
nexttile;
plotHandle = plot(Variable.nReflectors, snrDb);
grid on;
legend('WIT', 'location', 'nw');
xlabel('Number of reflectors');
ylabel('Average subband SNR [dB]');
xlim([Variable.nReflectors(1), Variable.nReflectors(end)]);
xticks(Variable.nReflectors(1 : 2 : end));
box on;
apply_style(plotHandle);

% * Power plot
nexttile;
plotHandle = gobjects(1, 2);
hold all;
plotHandle(1) = plot(Variable.nReflectors, mag2db(currentLinear));
plotHandle(2) = plot(Variable.nReflectors, mag2db(currentNonlinear));
hold off;
grid on;
legend('Linear WPT', 'Nonlinear WPT', 'location', 'nw');
xlabel('Number of reflectors');
ylabel('Average output DC current [dBA]');
xlim([Variable.nReflectors(1), Variable.nReflectors(end)]);
xticks(Variable.nReflectors(1 : 2 : end));
yticks(-100 : 20 : 0);
box on;
apply_style(plotHandle);

savefig('../figures/scaling_reflector.fig');
matlab2tikz('../../assets/scaling_reflector.tex');
