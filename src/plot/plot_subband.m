clear; clc; config_subband;

%% * Load batch data
reSet = cell(nBatches, length(Variable.nSubbands));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_subband_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reSubband = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
    reSubband{iSubband} = mean(cat(3, reSet{:, iSubband}), 3);
end
save('../data/re_subband.mat');

%% * R-E plots
figure('name', 'R-E region vs number of subbands');
legendString = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
    plot(reSubband{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reSubband{iSubband}(2, :), 'linewidth', 2);
    legendString{iSubband} = sprintf('N = %d', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
xlim([0 inf]);
ylim([0 inf]);
savefig('../figures/re_subband.fig');
matlab2tikz('../figures/re_subband.tex');
