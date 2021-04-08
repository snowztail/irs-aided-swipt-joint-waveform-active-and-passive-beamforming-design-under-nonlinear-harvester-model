clear; clc; close all; config_re_noise;

%% * Load batch data
indexSet = 1 : nBatches;
reAoSet = cell(nBatches, length(Variable.noisePower));
reLcSet = cell(nBatches, length(Variable.noisePower));
aoPowerRatio = zeros(nBatches, length(Variable.noisePower), nSamples);
lcPowerRatio = zeros(nBatches, length(Variable.noisePower), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_noise/re_noise_%d.mat', iBatch), 'reAoInstance', 'reLcInstance', 'reAoSolution', 'reLcSolution');
        reAoSet(iBatch, :) = reAoInstance;
        reLcSet(iBatch, :) = reLcInstance;
		for iNoise = 1 : length(Variable.noisePower)
			for iSample = 1 : nSamples
				aoPowerRatio(iBatch, iNoise, iSample) = reAoSolution{iNoise}{iSample}.powerRatio;
				lcPowerRatio(iBatch, iNoise, iSample) = reLcSolution{iNoise}{iSample}.powerRatio;
			end
		end
	catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAoNoise = cell(1, length(Variable.noisePower));
reLcNoise = cell(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    reAoNoise{iNoise} = mean(cat(3, reAoSet{indexSet, iNoise}), 3);
    reLcNoise{iNoise} = mean(cat(3, reLcSet{indexSet, iNoise}), 3);
end
aoPowerRatio = permute(mean(aoPowerRatio(indexSet, :, :), 1), [2 3 1]);
lcPowerRatio = permute(mean(lcPowerRatio(indexSet, :, :), 1), [2 3 1]);
save('../data/re_noise.mat');

%% * R-E and splitting ratio plots
figure('name', 'R-E region vs average noise power', 'position', [0, 0, 500, 400]);
legendString = cell(2, length(Variable.noisePower));
plotHandle = gobjects(2, length(Variable.noisePower));
hold all;
for iNoise = 1 : length(Variable.noisePower)
    plotHandle(1, iNoise) = plot(reAoNoise{iNoise}(1, :) / nSubbands, 1e6 * reAoNoise{iNoise}(2, :));
    plotHandle(2, iNoise) = plot(reLcNoise{iNoise}(1, :) / nSubbands, 1e6 * reLcNoise{iNoise}(2, :));
    legendString{1, iNoise} = sprintf('BCD: $\\sigma_n^2 = %d$ dBm', pow2db(Variable.noisePower(iNoise)) + 30);
    legendString{2, iNoise} = sprintf('LC-BCD: $\\sigma_n^2 = %d$ dBm', pow2db(Variable.noisePower(iNoise)) + 30);
end
hold off;
grid on;
legend(legendString(:));
xlabel('Per-subband rate [bps/Hz]');
ylabel('Output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_group_style(plotHandle(:), 2);
savefig('../figures/re_noise.fig');
matlab2tikz('../../assets/re_noise.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
close;

% * Power splitting ratio
figure('name', 'Splitting ratio vs average noise power', 'position', [0, 0, 500, 400]);
plotHandle = gobjects(2, length(Variable.noisePower));
hold all;
for iNoise = 1 : length(Variable.noisePower)
    plotHandle(1, iNoise) = plot(reAoNoise{iNoise}(1, :) / nSubbands, aoPowerRatio(iNoise, :));
    plotHandle(2, iNoise) = plot(reLcNoise{iNoise}(1, :) / nSubbands, lcPowerRatio(iNoise, :));
end
hold off;
grid on;
xlabel('Per-subband rate [bps/Hz]');
ylabel('Splitting ratio');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_group_style(plotHandle(:), 2);

savefig('../figures/splitting_ratio_noise.fig');
matlab2tikz('../../assets/splitting_ratio_noise.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
