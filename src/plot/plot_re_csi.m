clear; clc; close all; config_re_csi;

%% * Load batch data
indexSet = 1 : nBatches;
reNoIrsSet = cell(nBatches, 1);
reErrorSet = cell(nBatches, length(Variable.cascadedErrorVariance));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_csi/re_csi_%d.mat', iBatch), 'reNoIrsInstance', 'reErrorInstance');
        reNoIrsSet(iBatch) = reNoIrsInstance;
        reErrorSet(iBatch, :) = reErrorInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reNoIrsCsi = mean(cat(3, reNoIrsSet{indexSet}), 3);
reErrorCsi = cell(1, length(Variable.cascadedErrorVariance));
for iError = 1 : length(Variable.cascadedErrorVariance)
    reErrorCsi{iError} = mean(cat(3, reErrorSet{indexSet, iError}), 3);
end
save('../data/re_csi.mat');

%% * R-E plots
figure('name', 'R-E region vs cascaded channel estimation error', 'position', [0, 0, 500, 400]);
legendString = cell(1, length(Variable.cascadedErrorVariance) + 1);
plotHandle = gobjects(1, length(Variable.cascadedErrorVariance) + 1);
hold all;
for iError = 1 : length(Variable.cascadedErrorVariance)
    plotHandle(iError) = plot(reErrorCsi{iError}(1, :) / nSubbands, 1e6 * reErrorCsi{iError}(2, :));
	legendString{iError} = sprintf('$\\epsilon_n^2 = %s \Lambda_I\Lambda_R$', num2str(Variable.cascadedErrorVarianceRatio(iError)));
end
plotHandle(end) = plot(reNoIrsCsi(1, :) / nSubbands, 1e6 * reNoIrsCsi(2, :));
legendString{end} = 'no IRS';
hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

savefig('../figures/re_csi.fig');
matlab2tikz('../../assets/re_csi.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
