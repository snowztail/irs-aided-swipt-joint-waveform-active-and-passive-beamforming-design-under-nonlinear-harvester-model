clear; clc; close all; config_re_tx;

%% * Load batch data
indexSet = 1 : nBatches;
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
figure('name', 'Average R-E region vs number of transmit antennas', 'position', [0, 0, 500, 400]);
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
xlabel('Average subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

savefig('../figures/re_tx.fig');
matlab2tikz('../../assets/re_tx.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
