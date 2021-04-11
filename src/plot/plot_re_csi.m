clear; clc; close all; config_re_csi;

%% * Load batch data
indexSet = 1 : nBatches;
reRandomSet = cell(nBatches, length(Variable.nReflectors));
reErrorSet = cell(nBatches, length(Variable.nReflectors), length(Variable.cascadedErrorVariance));
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
reRandomCsi = cell(length(Variable.nReflectors), 1);
reErrorCsi = cell(length(Variable.nReflectors), length(Variable.cascadedErrorVariance));
for iReflector = 1 : length(Variable.nReflectors)
	reRandomCsi{iReflector} = mean(cat(3, reRandomSet{indexSet, iReflector}), 3);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		reErrorCsi{iReflector, iError} = mean(cat(4, reErrorSet{indexSet, iReflector, iError}), 4);
	end
end
save('../data/re_csi.mat');

%% * R-E plots
figure('name', 'R-E region vs cascaded channel estimation error', 'position', [0, 0, 500, 400]);
legendString = cell(length(Variable.nReflectors), length(Variable.cascadedErrorVariance) + 1);
plotHandle = gobjects(length(Variable.nReflectors), length(Variable.cascadedErrorVariance) + 1);
hold all;
for iReflector = 1 : length(Variable.nReflectors)
	nSubbands = Variable.nReflectors(iReflector);
	for iError = 1 : length(Variable.cascadedErrorVariance)
		plotHandle(iReflector, iError) = plot(reErrorCsi{iReflector, iError}(1, :) / nSubbands, 1e6 * reErrorCsi{iReflector, iError}(2, :));
		legendString{iReflector, iError} = sprintf('$\\epsilon_n^2 = %s \\Lambda_I\\Lambda_R$ $(L = %s)$', num2str(Variable.cascadedErrorVarianceRatio(iError)), num2str(Variable.nReflectors(iReflector)));
	end
	plotHandle(iReflector, iError + 1) = plot(reRandomCsi{iReflector}(1, :) / nSubbands, 1e6 * reRandomCsi{iReflector}(2, :));
	legendString{iReflector, iError + 1} = sprintf('Random IRS $(L = %s)$', num2str(Variable.nReflectors(iReflector)));
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
