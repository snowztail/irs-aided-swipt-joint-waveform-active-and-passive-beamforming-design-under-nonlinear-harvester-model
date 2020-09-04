clear; clc; config_snr;

%% * Load batch data
reSet = cell(nBatches, length(Variable.snr));
for iBatch = 1 : nBatches
    load(sprintf('../data/re_snr_%d.mat', iBatch), 'reInstance');
    reSet(iBatch, :) = reInstance;
end

%% * Average over batches
reSnr = cell(1, length(Variable.snr));
for iSnr = 1 : length(Variable.snr)
    reSnr{iSnr} = mean(cat(3, reSet{:, iSnr}), 3);
end
save('../data/re_snr.mat');

%% * R-E plots
figure('name', 'R-E region vs large-scale SNR');
legendString = cell(1, length(Variable.snr));
for iSnr = 1 : length(Variable.snr)
    plot(reSnr{iSnr}(1, :) / nSubbands, 1e6 * reSnr{iSnr}(2, :));
    legendString{iSnr} = sprintf('SNR = %d dB', pow2db(Variable.snr(iSnr)));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
ylim([0 inf]);
savefig('../figures/re_snr.fig');
