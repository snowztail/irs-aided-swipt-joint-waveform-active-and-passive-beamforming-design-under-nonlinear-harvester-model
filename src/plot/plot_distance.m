clear; clc; load('../data/re_distance.mat');

%% * Average over batches
reDistance = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    reDistance{iDistance} = mean(cat(3, reSet{:, iDistance}), 3);
end
save('../data/re_distance.mat', 'reDistance', '-append');

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
