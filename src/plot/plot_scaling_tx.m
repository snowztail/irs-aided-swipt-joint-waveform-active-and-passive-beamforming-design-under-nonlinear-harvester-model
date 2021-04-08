clear; clc; close all; config_scaling_tx;

%% * Load batch data
indexSet = 1 : nBatches;
rateWfSet = zeros(nBatches, length(Variable.nTxs));
currentAssSet = zeros(nBatches, length(Variable.nTxs));
currentSmfSet = zeros(nBatches, length(Variable.nTxs));
currentSdrSet = zeros(nBatches, length(Variable.nTxs));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/scaling_tx/scaling_tx_%d.mat', iBatch), 'rateWf', 'currentAss', 'currentSmf', 'currentSdr');
        rateWfSet(iBatch, :) = rateWf;
        currentAssSet(iBatch, :) = currentAss;
        currentSmfSet(iBatch, :) = currentSmf;
        currentSdrSet(iBatch, :) = currentSdr;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
rateWf = mean(rateWfSet(indexSet, :), 1);
currentAss = mean(currentAssSet(indexSet, :), 1);
currentSmf = mean(currentSmfSet(indexSet, :), 1);
currentSdr = mean(currentSdrSet(indexSet, :), 1);
snrDb = pow2db(2 .^ (rateWf / nSubbands));
save('../data/scaling_tx.mat');

%% * Rate and current plots
figure('name', 'SNR and harvested DC power vs number of transmit antennas', 'position', [0, 0, 500, 400]);
pathlossPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

% * SNR plot
nexttile;
plotHandle = plot(Variable.nTxs, snrDb);
grid on;
legend('WF', 'location', 'se');
xlabel('Number of transmit antennas');
ylabel('SNR [dB]');
xlim([Variable.nTxs(1), Variable.nTxs(end)]);
xticks(Variable.nTxs([1, 2 : 2 : end]));
box on;
apply_style(plotHandle);

% * Power plot
nexttile;
plotHandle = gobjects(1, 3);
hold all;
plotHandle(1) = plot(Variable.nTxs, mag2db(currentAss));
plotHandle(2) = plot(Variable.nTxs, mag2db(currentSmf));
plotHandle(3) = plot(Variable.nTxs, mag2db(currentSdr));
hold off;
grid on;
legend('ASS', 'SMF', 'GP', 'location', 'se');
xlabel('Number of transmit antennas');
ylabel('Output DC current [dBA]');
xlim([Variable.nTxs(1), Variable.nTxs(end)]);
xticks(Variable.nTxs([1, 2 : 2 : end]));
yticks(-100 : 20 : 0);
box on;
apply_style(plotHandle);

savefig('../figures/scaling_tx.fig');
matlab2tikz('../../assets/scaling_tx.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
