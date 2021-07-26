clear; clc; close all; config_re_lc;

%% * Load batch data
indexSet = 1 : nBatches;
reAoSet = cell(nBatches, length(Variable.nSubbands));
reLcSet = cell(nBatches, length(Variable.nSubbands), length(Variable.alpha));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_lc/re_lc_%d.mat', iBatch), 'reAoInstance', 'reLcInstance');
		reAoSet(iBatch, :) = reAoInstance;
		reLcSet(iBatch, :, :) = reLcInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAoSubband = cell(length(Variable.nSubbands), 1);
reLcSubband = cell(length(Variable.nSubbands), length(Variable.alpha));
for iSubband = 1 : length(Variable.nSubbands)
	reAoSubband{iSubband} = mean(cat(3, reAoSet{indexSet, iSubband}), 3);
	for iAlpha = 1 : length(Variable.alpha)
		reLcSubband{iSubband, iAlpha} = mean(cat(4, reLcSet{indexSet, iSubband, iAlpha}), 4);
	end
end
save('../data/re_lc.mat');

%% * R-E plots
for iSubband = 1 : length(Variable.nSubbands)
	figure('name', sprintf('R-E region vs SMF ratio for $N = %s$', num2str(Variable.nSubbands(iSubband))), 'position', [0, 0, 500, 400]);
	legendString = cell(2, length(Variable.alpha) + 1);
	plotHandle = gobjects(2, length(Variable.alpha) + 1);
	hold all;
	plotHandle(1, 1) = plot(reAoSubband{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reAoSubband{iSubband}(2, :));
	legendString{1, 1} = 'BCD';
	plotHandle(2, 1) = plot([reAoSubband{iSubband}(1, end) / Variable.nSubbands(iSubband), reAoSubband{iSubband}(1, 1) / Variable.nSubbands(iSubband)], [1e6 * reAoSubband{iSubband}(2, end), 1e6 * reAoSubband{iSubband}(2, 1)]);
	legendString{2, 1} = 'AO (TS)';
	for iAlpha = 1 : length(Variable.alpha)
		plotHandle(1, iAlpha + 1) = plot(reLcSubband{iSubband, iAlpha}(1, :) / Variable.nSubbands(iSubband), 1e6 * reLcSubband{iSubband, iAlpha}(2, :));
		legendString{1, iAlpha + 1} = sprintf('LC-BCD: $\\alpha = %s$', num2str(Variable.alpha(iAlpha)));
		plotHandle(2, iAlpha + 1) = plot([reLcSubband{iSubband, iAlpha}(1, end) / Variable.nSubbands(iSubband), reLcSubband{iSubband, iAlpha}(1, 1) / Variable.nSubbands(iSubband)], [1e6 * reLcSubband{iSubband, iAlpha}(2, end), 1e6 * reLcSubband{iSubband, iAlpha}(2, 1)]);
		legendString{2, iAlpha + 1} = sprintf('LC (TS): $\\alpha = %s$', num2str(Variable.alpha(iAlpha)));
	end
	hold off;
	grid on;
	legend(legendString(:));
	xlabel('Per-subband rate [bps/Hz]');
	ylabel('DC current [$\mu$A]');
	xlim([0 inf]);
	ylim([0 inf]);
	box on;
	apply_group_style(plotHandle(:), 2);

	savefig(sprintf('../figures/re_lc_%ssubbands.fig', num2str(Variable.nSubbands(iSubband))));
	matlab2tikz(sprintf('../assets/re_lc_%ssubbands.tex', num2str(Variable.nSubbands(iSubband))), 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
	close;
end
