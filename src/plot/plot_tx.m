clear; clc; close all; config_tx;

%% * Load batch data
indexSet = 1 : nSamples;
reSet = cell(nBatches, length(Variable.nTxs));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_tx/re_tx_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reTx = cell(1, length(Variable.nTxs));
for iSubband = 1 : length(Variable.nTxs)
    reTx{iSubband} = mean(cat(3, reSet{indexSet, iSubband}), 3);
end
save('../data/re_tx.mat');

%% * R-E plots
figure('name', 'R-E region vs number of transmit antennas');
legendString = cell(1, length(Variable.nTxs));
plotHandle = gobjects(1, length(Variable.nTxs));
for iTx = 1 : length(Variable.nTxs)
    plotHandle(iTx) = plot(reTx{iTx}(1, :) / nSubbands, 1e6 * reTx{iTx}(2, :));
    legendString{iTx} = sprintf('$M = %d, L = %d$', Variable.nTxs(iTx), nReflectors);
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
savefig('../figures/re_tx.fig');
matlab2tikz('../../assets/re_tx.tex');
