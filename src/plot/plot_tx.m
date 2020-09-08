clear; clc; config_tx;

%% * Load batch data
reSet = cell(nBatches, length(Variable.nTxs));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_tx_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reTx = cell(1, length(Variable.nTxs));
for iSubband = 1 : length(Variable.nTxs)
    reTx{iSubband} = mean(cat(3, reSet{:, iSubband}), 3);
end
save('../data/re_tx.mat');

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
