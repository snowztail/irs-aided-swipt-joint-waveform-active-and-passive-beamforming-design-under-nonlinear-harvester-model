clear; clc; close all; config_los;

%% * Load batch data
indexSet = 1 : nSamples;
reNlosSet = cell(nBatches, 1);
reLosSet = cell(nBatches, 1);
rateNlosSet = zeros(nBatches, 1);
rateLosSet = zeros(nBatches, 1);
currentNlosSet = zeros(nBatches, 1);
currentLosSet = zeros(nBatches, 1);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_los/re_los_%d.mat', iBatch), 'reNlosInstance', 'reLosInstance');
        reNlosSet{iBatch} = reNlosInstance;
		reLosSet{iBatch} = reLosInstance;
		rateNlosSet(iBatch) = reNlosInstance(1, 1);
		rateLosSet(iBatch) = reLosInstance(1, 1);
		currentNlosSet(iBatch) = reNlosInstance(2, end);
		currentLosSet(iBatch) = reLosInstance(2, end);
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reNlos = mean(cat(3, reNlosSet{indexSet}), 3);
reLos = mean(cat(3, reLosSet{indexSet}), 3);
save('../data/re_los.mat');

%% * R-E plots
figure('name', 'R-E region for IRS-aided NLoS and LoS channels');
plotHandle = gobjects(1, nCases);
hold all;
plotHandle(1) = plot(reNlos(1, :) / nSubbands, 1e6 * reNlos(2, :));
plotHandle(2) = plot(reLos(1, :) / nSubbands, 1e6 * reLos(2, :));
hold off;
grid on;
legend('NLoS', 'LoS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);

apply_style(plotHandle);
savefig('../figures/re_los.fig');
matlab2tikz('../../assets/re_los.tex');
close;

%% * CDF plots
figure('name', 'WIT and WPT CDF for IRS-aided NLoS and LoS channels');
cdfPlot = tiledlayout(2, 1, 'tilespacing', 'compact');

nexttile;
plotHandle = gobjects(1, nCases);
hold all;
plotHandle(1) = cdfplot(rateNlosSet);
plotHandle(2) = cdfplot(rateLosSet);
hold off;
grid on;
legend('NLoS', 'LoS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Cumulative probability');
xlim([0 inf]);
ylim([0 inf]);
apply_style(plotHandle);

nexttile;
plotHandle = gobjects(1, nCases);
hold all;
plotHandle(1) = cdfplot(currentNlosSet);
plotHandle(2) = cdfplot(currentLosSet);
hold off;
grid on;
legend('NLoS', 'LoS');
xlabel('Average output DC current [$\mu$A]');
ylabel('Cumulative probability');
xlim([0 inf]);
ylim([0 inf]);
apply_style(plotHandle);

savefig('../figures/cdf_los.fig');
matlab2tikz('../../assets/cdf_los.tex');
