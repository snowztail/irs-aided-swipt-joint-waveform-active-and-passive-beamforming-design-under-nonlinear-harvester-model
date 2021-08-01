clear; clc; close all; config_scaling_reflector;

%% * Load batch data
indexSet = 1 : nBatches;
rateWfSet = zeros(nBatches, length(Variable.nReflectors));
currentAssSet = zeros(nBatches, length(Variable.nReflectors));
currentSmfSet = zeros(nBatches, length(Variable.nReflectors));
currentSdrSet = zeros(nBatches, length(Variable.nReflectors));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/scaling_reflector/scaling_reflector_%d.mat', iBatch), 'rateWf', 'currentAss', 'currentSmf', 'currentSdr');
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
save('../data/scaling_reflector.mat');

%% * SNR and power plots
figure('name', 'SNR and harvested DC power vs number of reflectors', 'position', [0, 0, 500, 400]);
pathlossPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

% * SNR plot
nexttile;
plotHandle = plot(Variable.nReflectors, snrDb);
grid on;
legend('WF', 'location', 'se');
ylabel('SNR [dB]');
xlim([Variable.nReflectors(1), Variable.nReflectors(end)]);
xticks(Variable.nReflectors(1 : 2 : end));
box on;
apply_style(plotHandle);

% * Power plot
nexttile;
plotHandle = gobjects(1, 3);
hold all;
plotHandle(1) = plot(Variable.nReflectors, mag2db(currentAss));
plotHandle(2) = plot(Variable.nReflectors, mag2db(currentSmf));
plotHandle(3) = plot(Variable.nReflectors, mag2db(currentSdr));
hold off;
grid on;
legend('LEH', 'SMF', 'GP', 'location', 'se', 'orientation', 'horizontal');
xlabel('Number of reflectors');
ylabel('DC current [dBA]');
xlim([Variable.nReflectors(1), Variable.nReflectors(end)]);
xticks(Variable.nReflectors(1 : 2 : end));
yticks(-100 : 20 : 0);
box on;
apply_style(plotHandle);

savefig('../figures/scaling_reflector.fig');
matlab2tikz('../../assets/scaling_reflector.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
