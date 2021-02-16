clear; clc; close all;

%% * Max eigenvalue over sum eigenvalue vs number of subbands
config_re_subband;

indexSet = 1 : nBatches;
eigRatioSubband = cell(nBatches, length(Variable.nSubbands), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_subband/re_subband_%d.mat', iBatch), 'reSolution');
		for iSubband = 1 : length(Variable.nSubbands)
			for iSample = 1 : nSamples
				eigRatioSubband{iBatch, iSubband, iSample} = reSolution{iSubband}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
clearvars -except eigRatioSubband;

%% * Max eigenvalue over sum eigenvalue vs number of transmit antennas
config_re_tx;

indexSet = 1 : nBatches;
eigRatioTx = cell(nBatches, length(Variable.nTxs), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_tx/re_tx_%d.mat', iBatch), 'reSolution');
		for iTx = 1 : length(Variable.nTxs)
			for iSample = 1 : nSamples
				eigRatioTx{iBatch, iTx, iSample} = reSolution{iTx}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
clearvars -except eigRatioSubband eigRatioTx;

%% * Max eigenvalue over sum eigenvalue vs number of IRS elements
config_re_reflector;

indexSet = 1 : nBatches;
eigRatioReflector = cell(nBatches, length(Variable.nReflectors), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_reflector/re_reflector_%d.mat', iBatch), 'reSolution');
		for iReflector = 1 : length(Variable.nReflectors)
			for iSample = 1 : nSamples
				eigRatioReflector{iBatch, iReflector, iSample} = reSolution{iReflector}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
clearvars -except eigRatioSubband eigRatioTx eigRatioReflector;

%% * Eigenratio plots
figure('name', 'Max eigenvalue of the relaxed solution over sum eigenvalue of the relaxed solution under different configurations');
eigenRatioPlot = tiledlayout(3, 1, 'tilespacing', 'compact');

config_re_subband;
nexttile;
plotHandle = gobjects(1, length(Variable.nSubbands));
legendString = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
	plotHandle(iSubband) = cdfplot([eigRatioSubband{:, iSubband, :}]);
	legendString{iSubband} = sprintf('$N = %d$', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid on;
box on;
title('');
legend(legendString);
apply_style(plotHandle);

config_re_tx;
nexttile;
plotHandle = gobjects(1, length(Variable.nTxs));
legendString = cell(1, length(Variable.nTxs));
for iTx = 1 : length(Variable.nTxs)
	plotHandle(iTx) = cdfplot([eigRatioTx{:, iTx, :}]);
	legendString{iTx} = sprintf('$M = %d$', Variable.nTxs(iTx));
    hold on;
end
hold off;
grid on;
box on;
title('');
legend(legendString);
xlabel('Max eigenvalue over sum eigenvalue');
apply_style(plotHandle);

config_re_reflector;
nexttile;
plotHandle = gobjects(1, length(Variable.nReflectors));
legendString = cell(1, length(Variable.nReflectors));
for iReflector = 1 : length(Variable.nReflectors)
	plotHandle(iReflector) = cdfplot([eigRatioReflector{:, iReflector, :}]);
	legendString{iReflector} = sprintf('$L = %d$', Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid on;
box on;
title('');
legend(legendString);
ylabel('Cumulative probability');
apply_style(plotHandle);

savefig('../figures/cdf_eigenratio.fig');
matlab2tikz('../../assets/cdf_eigenratio.tex');
