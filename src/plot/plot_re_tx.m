clear; clc; close all; config_re_tx;

%% * Load batch data
indexSet = 1 : nBatches;
reAoSet = cell(nBatches, length(Variable.nTxs));
reLcSet = cell(nBatches, length(Variable.nTxs));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_tx/re_tx_%d.mat', iBatch), 'reAoInstance', 'reLcInstance');
        reAoSet(iBatch, :) = reAoInstance;
        reLcSet(iBatch, :) = reLcInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAoTx = cell(1, length(Variable.nTxs));
reLcTx = cell(1, length(Variable.nTxs));
for iSubband = 1 : length(Variable.nTxs)
    reAoTx{iSubband} = mean(cat(3, reAoSet{indexSet, iSubband}), 3);
    reLcTx{iSubband} = mean(cat(3, reLcSet{indexSet, iSubband}), 3);
end
save('../data/re_tx.mat');

%% * R-E plots
figure('name', 'R-E region vs number of transmit antennas', 'position', [0, 0, 500, 400]);
legendString = cell(2, length(Variable.nTxs));
plotHandle = gobjects(2, length(Variable.nTxs));
hold all;
for iTx = 1 : length(Variable.nTxs)
    plotHandle(1, iTx) = plot(reAoTx{iTx}(1, :) / nSubbands, 1e6 * reAoTx{iTx}(2, :));
    plotHandle(2, iTx) = plot(reLcTx{iTx}(1, :) / nSubbands, 1e6 * reLcTx{iTx}(2, :));
	legendString{1, iTx} = sprintf('BCD: $M = %d$', Variable.nTxs(iTx));
	legendString{2, iTx} = sprintf('LC-BCD: $M = %d$', Variable.nTxs(iTx));
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

savefig('../figures/re_tx.fig');
matlab2tikz('../../assets/re_tx.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
