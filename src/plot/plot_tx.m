clear; clc; load('../data/re_tx.mat');

%% * Average over batches
reTx = cell(1, length(Variable.nTxs));
for iSubband = 1 : length(Variable.nTxs)
    reTx{iSubband} = mean(cat(3, reSet{:, iSubband}), 3);
end
save('../data/re_tx.mat', 'reTx', '-append');

%% * R-E plots
figure('name', 'R-E region vs number of transmit antennas');
legendString = cell(1, length(Variable.nTxs));
for iTx = 1 : length(Variable.nTxs)
    plot(reTx{iTx}(1, :) / nSubbands, 1e6 * reTx{iTx}(2, :));
    legendString{iTx} = sprintf('M = %d', Variable.nTxs(iTx));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('../figures/re_tx.fig');
