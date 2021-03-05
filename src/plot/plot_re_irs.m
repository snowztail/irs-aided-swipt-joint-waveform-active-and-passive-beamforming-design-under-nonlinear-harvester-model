clear; clc; close all; config_re_irs;

%% * Load batch data
indexSet = 1 : nBatches;
reIdealIrsSet = cell(nBatches, length(Variable.bandwidth));
reAdaptiveIrsSet = cell(nBatches, length(Variable.bandwidth));
reWitIrsSet = cell(nBatches, length(Variable.bandwidth));
reWptIrsSet = cell(nBatches, length(Variable.bandwidth));
reNoIrsSet = cell(nBatches, length(Variable.bandwidth));
for iBatch = 1 : nBatches
    try
		load(sprintf('../data/re_irs/re_irs_%d.mat', iBatch), 'reIdealIrsInstance', 'reAdaptiveIrsInstance', 'reWitIrsInstance', 'reWptIrsInstance', 'reNoIrsInstance');
        reIdealIrsSet(iBatch, :) = reIdealIrsInstance;
        reAdaptiveIrsSet(iBatch, :) = reAdaptiveIrsInstance;
        reWitIrsSet(iBatch, :) = reWitIrsInstance;
        reWptIrsSet(iBatch, :) = reWptIrsInstance;
        reNoIrsSet(iBatch, :) = reNoIrsInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reIdealIrs = cell(1, length(Variable.bandwidth));
reAdaptiveIrs = cell(1, length(Variable.bandwidth));
reWitIrs = cell(1, length(Variable.bandwidth));
reWptIrs = cell(1, length(Variable.bandwidth));
reNoIrs = cell(1, length(Variable.bandwidth));
for iBandwidth = 1 : length(Variable.bandwidth)
	reIdealIrs{iBandwidth} = mean(cat(3, reIdealIrsSet{indexSet, iBandwidth}), 3);
	reAdaptiveIrs{iBandwidth} = mean(cat(3, reAdaptiveIrsSet{indexSet, iBandwidth}), 3);
	reWitIrs{iBandwidth} = mean(cat(3, reWitIrsSet{indexSet, iBandwidth}), 3);
	reWptIrs{iBandwidth} = mean(cat(3, reWptIrsSet{indexSet, iBandwidth}), 3);
	reNoIrs{iBandwidth} = mean(cat(3, reNoIrsSet{indexSet, iBandwidth}), 3);
end
save('../data/re_irs.mat');

%% * R-E plots
for iBandwidth = 1 : length(Variable.bandwidth)
	figure('name', sprintf('Average R-E region for ideal, adaptive, nonadaptive and no IRS for $B = %d$ MHz', Variable.bandwidth(iBandwidth) / 1e6));
	plotHandle = gobjects(1, length(Variable.bandwidth));
	hold all;
	plotHandle(1) = plot(reIdealIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reIdealIrs{iBandwidth}(2, :));
	plotHandle(2) = plot(reAdaptiveIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reAdaptiveIrs{iBandwidth}(2, :));
	plotHandle(3) = plot(reWitIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reWitIrs{iBandwidth}(2, :));
	plotHandle(4) = plot(reWptIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reWptIrs{iBandwidth}(2, :));
	plotHandle(5) = plot(reNoIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reNoIrs{iBandwidth}(2, :));
	hold off;
	grid on;
	legend('Ideal FS IRS', 'Adaptive IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
	xlabel('Average subband rate [bps/Hz]');
	ylabel('Average output DC current [$\mu$A]');
	xlim([0 inf]);
	ylim([0 inf]);
    box on;
	apply_style(plotHandle);

	savefig(sprintf('../figures/re_irs_%dmhz.fig', Variable.bandwidth(iBandwidth) / 1e6));
	matlab2tikz(sprintf('../../assets/re_irs_%dmhz.tex', Variable.bandwidth(iBandwidth) / 1e6));
	close;
end
