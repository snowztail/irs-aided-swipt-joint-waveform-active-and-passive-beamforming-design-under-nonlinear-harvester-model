clear; clc; close all; config_noise;

%% * Load batch data
reSet = cell(nBatches, length(Variable.noisePower));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_noise/re_noise_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reNoise = cell(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    reNoise{iNoise} = mean(cat(3, reSet{:, iNoise}), 3);
end
save('../data/re_noise.mat');

%% * R-E plots
figure('name', 'R-E region vs average noise power');
legendString = cell(1, length(Variable.noisePower));
plotHandle = gobjects(1, length(Variable.noisePower));
for iNoise = 1 : length(Variable.noisePower)
    plotHandle(iNoise) = plot(reNoise{iNoise}(1, :) / nSubbands, 1e6 * reNoise{iNoise}(2, :));
    legendString{iNoise} = sprintf('$\\sigma_n = %d$ dB', pow2db(Variable.noisePower(iNoise)));
    hold on;
end
hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);

apply_style(plotHandle);
savefig('../figures/re_noise.fig');
matlab2tikz('../../assets/re_noise.tex');
