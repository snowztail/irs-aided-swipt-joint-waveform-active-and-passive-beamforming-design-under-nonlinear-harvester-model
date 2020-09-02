clear; clc; load('../data/re_snr.mat');

%% * Average over batches
reSnr = cell(1, length(Variable.snr));
for iSnr = 1 : length(Variable.snr)
    reSnr{iSnr} = mean(cat(3, reSet{:, iSnr}), 3);
end
save('../data/re_snr.mat', 'reSnr', '-append');

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
