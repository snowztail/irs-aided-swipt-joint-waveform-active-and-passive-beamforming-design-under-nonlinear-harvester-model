clear; clc; close all; config_scaling_tx;

%% * Load batch data
indexSet = 1 : nBatches;
rateSet = zeros(nBatches, length(Variable.nTxs));
currentLinearSet = zeros(nBatches, length(Variable.nTxs));
currentNonlinearSet = zeros(nBatches, length(Variable.nTxs));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/scaling_tx/scaling_tx_%d.mat', iBatch), 'rateInstance', 'currentLinearInstance', 'currentNonlinearInstance');
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
save('../data/scaling_tx.mat');

%% * Rate and current plots
figure('name', 'Average SNR and harvested DC power vs number of transmit antennas');
pathlossPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

% * SNR plot
nexttile;
plotHandle = plot(Variable.nTxs, snrDb);
grid on;
legend('WIT', 'location', 'nw');
xlabel('Number of transmit antennas');
ylabel('Average subband SNR [dB]');
xlim([Variable.nTxs(1), Variable.nTxs(end)]);
xticks(Variable.nTxs([1, 2 : 2 : end]));

apply_style(plotHandle);

% * Power plot
nexttile;
plotHandle = gobjects(1, 2);
hold all;
plotHandle(1) = plot(Variable.nTxs, mag2db(currentLinear));
plotHandle(2) = plot(Variable.nTxs, mag2db(currentNonlinear));
hold off;
grid on;
legend('Linear WPT', 'Nonlinear WPT', 'location', 'nw');
xlabel('Number of transmit antennas');
ylabel('Average output DC current [dBA]');
xlim([Variable.nTxs(1), Variable.nTxs(end)]);
xticks(Variable.nTxs([1, 2 : 2 : end]));
yticks(-100 : 20 : 0)

apply_style(plotHandle);
savefig('../figures/scaling_tx.fig');
matlab2tikz('../../assets/scaling_tx.tex');
