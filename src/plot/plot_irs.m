clear; clc; close all; config_irs;

%% * Load batch data
reAdaptiveIrsSet = cell(nBatches, length(Variable.bandwidth));
reFsIrsSet = cell(nBatches, length(Variable.bandwidth));
reWitIrsSet = cell(nBatches, length(Variable.bandwidth));
reWptIrsSet = cell(nBatches, length(Variable.bandwidth));
reNoIrsSet = cell(nBatches, length(Variable.bandwidth));
for iBatch = 1 : nBatches
    try
		load(sprintf('../data/re_irs/re_irs_%d.mat', bandwidth / 1e6, iBatch), 'reAdaptiveIrsInstance', 'reFsIrsInstance', 'reWitIrsInstance', 'reWptIrsInstance', 'reNoIrsInstance');
        reAdaptiveIrsSet(iBatch, :) = reAdaptiveIrsInstance;
        reFsIrsSet(iBatch, :) = reFsIrsInstance;
        reWitIrsSet(iBatch, :) = reWitIrsInstance;
        reWptIrsSet(iBatch, :) = reWptIrsInstance;
        reNoIrsSet(iBatch, :) = reNoIrsInstance;
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reAdaptiveIrs = cell(1, length(Variable.bandwidth));
reFsIrs = cell(1, length(Variable.bandwidth));
reWitIrs = cell(1, length(Variable.bandwidth));
reWptIrs = cell(1, length(Variable.bandwidth));
reNoIrs = cell(1, length(Variable.bandwidth));
for iBandwidth = 1 : length(Variable.bandwidth)
	reAdaptiveIrs{iBandwidth} = mean(cat(3, reAdaptiveIrsSet{:, iBandwidth}), 3);
	reFsIrs{iBandwidth} = mean(cat(3, reFsIrsSet{:, iBandwidth}), 3);
	reWitIrs{iBandwidth} = mean(cat(3, reWitIrsSet{:, iBandwidth}), 3);
	reWptIrs{iBandwidth} = mean(cat(3, reWptIrsSet{:, iBandwidth}), 3);
	reNoIrs{iBandwidth} = mean(cat(3, reNoIrsSet{:, iBandwidth}), 3);
end
save('../data/re_irs.mat');

%% * R-E plots
figure('name', 'R-E region for adaptive, fixed and no IRS');
rePlot = tiledlayout(length(Variable.bandwidth), 1, 'tilespacing', 'compact');
for iBandwidth = 1 : length(Variable.bandwidth)
	nexttile;
	plotHandle = gobjects(1, length(Variable.bandwidth));
	hold all;
	plotHandle(1) = plot(reAdaptiveIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reAdaptiveIrs{iBandwidth}(2, :));
	plotHandle(2) = plot(reFsIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reFsIrs{iBandwidth}(2, :));
	plotHandle(3) = plot(reWitIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reWitIrs{iBandwidth}(2, :));
	plotHandle(4) = plot(reWptIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reWptIrs{iBandwidth}(2, :));
	plotHandle(5) = plot(reNoIrs{iBandwidth}(1, :) / nSubbands, 1e6 * reNoIrs{iBandwidth}(2, :));
	hold off;
	grid on;
	xlim([0 inf]);
	ylim([0 inf]);
	if iBandwidth == 1
		legend('Adaptive IRS', 'Ideal FS IRS', 'WIT-optimized IRS', 'WPT-optimized IRS', 'No IRS');
	elseif iBandwidth == 2
		ylabel('Average output DC current [$\mu$A]');
	end
	title(sprintf('$B = %d$ MHz', Variable.bandwidth(iBandwidth) / 1e6));
	apply_style(plotHandle);
end
xlabel('Per-subband rate [bps/Hz]');

savefig('../figures/re_irs.fig');
matlab2tikz('../../assets/re_irs.tex');
