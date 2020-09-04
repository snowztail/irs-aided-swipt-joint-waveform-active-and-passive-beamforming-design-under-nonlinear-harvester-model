clear; clc; config_distance;

%% * Load batch data
reSet = cell(nBatches, length(Variable.horizontalDistance));
for iBatch = 1 : nBatches
    load(sprintf('../data/re_distance_%d.mat', iBatch), 'reInstance');
    reSet(iBatch, :) = reInstance;
end

%% * Average over batches
reDistance = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    reDistance{iDistance} = mean(cat(3, reSet{:, iDistance}), 3);
end
save('../data/re_distance.mat');

%% * R-E plots
figure('name', 'R-E region vs AP-IRS horizontal distance');
legendString = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    plot(reDistance{iDistance}(1, :) / nSubbands, 1e6 * reDistance{iDistance}(2, :));
    legendString{iDistance} = sprintf('d_H = %d', Variable.horizontalDistance(iDistance));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('../figures/re_distance.fig');
