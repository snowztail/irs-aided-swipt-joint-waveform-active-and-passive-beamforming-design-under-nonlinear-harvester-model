clear; clc; config_irs;

%% * Load batch data
reSet = cell(nBatches, nCases);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_irs_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reIrs = cell(1, nCases);
for iCase = 1 : nCases
    reIrs{iCase} = mean(cat(3, reSet{:, iCase}), 3);
end
save('../data/re_irs.mat');

%% * R-E plots
figure('name', 'R-E region for adaptive, fixed and no IRS');
for iCase = 1 : nCases
    plot(reIrs{iCase}(1, :) / nSubbands, 1e6 * reIrs{iCase}(2, :), 'linewidth', 2);
    hold on;
end
hold off;
grid on;
legend('Adaptive IRS', 'Ideal FS IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [\muA]');
xlim([0 inf]);
ylim([0 inf]);
savefig('../figures/re_irs.fig');
matlab2tikz('../figures/re_irs.tex');
