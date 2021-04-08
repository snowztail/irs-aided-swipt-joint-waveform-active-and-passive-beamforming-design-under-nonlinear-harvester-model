clear; clc; close all; config_re_csi;

%% * Load batch data
indexSet = 1 : nBatches;
reRandomSet = cell(nBatches, length(Variable.nSubbands));
reErrorSet = cell(nBatches, length(Variable.nSubbands), length(Variable.cascadedErrorVariance));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_csi/re_csi_%d.mat', iBatch), 'reRandomInstance', 'reErrorInstance');
        reRandomSet(iBatch, :) = reRandomInstance;
        reErrorSet(iBatch, :, :) = reErrorInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reRandomCsi = cell(length(Variable.nSubbands), 1);
reErrorCsi = cell(length(Variable.nSubbands), length(Variable.cascadedErrorVariance));
for iSubband = 1 : length(Variable.nSubbands)
	reRandomCsi{iSubband} = mean(cat(3, reRandomSet{indexSet, iSubband}), 3);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		reErrorCsi{iSubband, iError} = mean(cat(4, reErrorSet{indexSet, iSubband, iError}), 4);
	end
end
save('../data/re_csi.mat');

%% * R-E plots
figure('name', 'R-E region vs cascaded channel estimation error', 'position', [0, 0, 500, 400]);
legendString = cell(length(Variable.nSubbands), length(Variable.cascadedErrorVariance) + 1);
plotHandle = gobjects(length(Variable.nSubbands), length(Variable.cascadedErrorVariance) + 1);
hold all;
for iSubband = 1 : length(Variable.nSubbands)
	nSubbands = Variable.nSubbands(iSubband);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		plotHandle(iSubband, iError) = plot(reErrorCsi{iSubband, iError}(1, :) / nSubbands, 1e6 * reErrorCsi{iSubband, iError}(2, :));
		legendString{iSubband, iError} = sprintf('$\\epsilon_n^2 = %s \\Lambda_I\\Lambda_R$ $(N = %s)$', num2str(Variable.cascadedErrorVarianceRatio(iError)), num2str(Variable.nSubbands(iSubband)));
	end
	plotHandle(iSubband, iError + 1) = plot(reRandomCsi{iSubband}(1, :) / nSubbands, 1e6 * reRandomCsi{iSubband}(2, :));
	legendString{iSubband, iError + 1} = sprintf('Random IRS $(N = %s)$', num2str(Variable.nSubbands(iSubband)));
end
hold off;
grid on;
legend(vec(transpose(legendString)));
xlabel('Per-subband rate [bps/Hz]');
ylabel('Output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_group_style(vec(transpose(plotHandle)), length(Variable.cascadedErrorVariance) + 1);

savefig('../figures/re_csi.fig');
matlab2tikz('../../assets/re_csi.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
