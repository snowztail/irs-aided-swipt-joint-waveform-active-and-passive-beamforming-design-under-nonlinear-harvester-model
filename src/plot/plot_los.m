clear; clc; close all; config_los;

%% * Load batch data
reSet = cell(nBatches, nCases);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_los_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reLos = cell(1, nCases);
for iCase = 1 : nCases
    reLos{iCase} = mean(cat(3, reSet{:, iCase}), 3);
end
save('../data/re_los.mat');

%% * R-E plots
figure('name', 'R-E region for IRS-aided NLoS and LoS channels');
for iCase = 1 : nCases
    plot(reLos{iCase}(1, :) / nSubbands, 1e6 * reLos{iCase}(2, :));
    hold on;
end
hold off;
grid on;
legend('NLoS', 'LoS');
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
savefig('../figures/re_los.fig');
matlab2tikz('../../assets/re_los.tex');
