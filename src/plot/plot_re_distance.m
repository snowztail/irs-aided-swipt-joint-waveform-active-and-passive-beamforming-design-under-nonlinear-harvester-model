clear; clc; close all; config_re_distance;

%% * Load batch data
indexSet = 1 : nBatches;
reSet = cell(nBatches, length(Variable.horizontalDistance));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_distance/re_distance_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reDistance = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    reDistance{iDistance} = mean(cat(3, reSet{indexSet, iDistance}), 3);
end
save('../data/re_distance.mat');

%% * R-E plots
figure('name', 'R-E region vs AP-IRS horizontal distance');
legendString = cell(1, length(Variable.horizontalDistance));
plotHandle = gobjects(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    plotHandle(iDistance) = plot(reDistance{iDistance}(1, :) / nSubbands, 1e6 * reDistance{iDistance}(2, :));
    legendString{iDistance} = sprintf('$d_H = %s$', num2str(Variable.horizontalDistance(iDistance)));
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
savefig('../figures/re_distance.fig');
matlab2tikz('../../assets/re_distance.tex');
