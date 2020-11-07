clear; clc; close all; config_re_noise;

%% * Load batch data
indexSet = 1 : nBatches;
reSet = cell(nBatches, length(Variable.noisePower));
powerRatio = zeros(nBatches, length(Variable.noisePower), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_noise/re_noise_%d.mat', iBatch), 'reInstance', 'reSolution');
        reSet(iBatch, :) = reInstance;
		for iNoise = 1 : length(Variable.noisePower)
			for iSample = 1 : nSamples
				powerRatio(iBatch, iNoise, iSample) = reSolution{iNoise}{iSample}.powerRatio;
			end
		end
	catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reNoise = cell(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    reNoise{iNoise} = mean(cat(3, reSet{indexSet, iNoise}), 3);
end
powerRatio = permute(mean(powerRatio(indexSet, :, :), 1), [2 3 1]);
save('../data/re_noise.mat');

%% * R-E and splitting ratio plots
figure('name', 'Average R-E region and power splitting ratio vs average noise power');
tiledlayout(2, 1, 'tilespacing', 'compact');

% * R-E region
nexttile;
legendString = cell(1, length(Variable.noisePower));
plotHandle = gobjects(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    plotHandle(iNoise) = plot(reNoise{iNoise}(1, :) / nSubbands, 1e6 * reNoise{iNoise}(2, :));
    legendString{iNoise} = sprintf('$\\sigma_n = %d$ dBm', pow2db(Variable.noisePower(iNoise)) + 30);
    hold on;
end
hold off;
grid on;
legend(legendString);
xlabel('Average subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

% * Power splitting ratio
nexttile;
plotHandle = gobjects(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    plotHandle(iNoise) = plot(reNoise{iNoise}(1, :) / nSubbands, powerRatio(iNoise, :));
    hold on;
end
hold off;
grid on;
xlabel('Average subband rate [bps/Hz]');
ylabel('Average power splitting ratio');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

savefig('../figures/re_noise.fig');
matlab2tikz('../../assets/re_noise.tex');
