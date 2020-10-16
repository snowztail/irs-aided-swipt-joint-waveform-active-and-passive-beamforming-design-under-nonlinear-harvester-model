clear; clc; close all; config_re_subband;

%% * Load batch data
indexSet = 1 : nBatches;
reAdaptiveIrsSet = cell(nBatches, length(Variable.nSubbands));
reFsIrsSet = cell(nBatches, length(Variable.nSubbands));
infoAmplitudeSet = cell(nBatches, length(Variable.nSubbands));
powerAmplitudeSet = cell(nBatches, length(Variable.nSubbands));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_subband/re_subband_%d.mat', iBatch), 'reAdaptiveIrsInstance', 'reFsIrsInstance', 'reAdaptiveIrsSolution', 'reFsIrsSolution');
		reAdaptiveIrsSet(iBatch, :) = reAdaptiveIrsInstance;
		reFsIrsSet(iBatch, :) = reFsIrsInstance;
		for iSubband = 1 : length(Variable.nSubbands)
			infoAmplitudeSet{iBatch, iSubband} = sort(reAdaptiveIrsSolution{iSubband}{end}.infoAmplitude);
			powerAmplitudeSet{iBatch, iSubband} = sort(reAdaptiveIrsSolution{iSubband}{end}.powerAmplitude);
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAdaptiveIrs = cell(1, length(Variable.nSubbands));
reFsIrs = cell(1, length(Variable.nSubbands));
infoAmplitude = cell(1, length(Variable.nSubbands));
powerAmplitude = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
	reAdaptiveIrs{iSubband} = mean(cat(3, reAdaptiveIrsSet{indexSet, iSubband}), 3);
	reFsIrs{iSubband} = mean(cat(3, reFsIrsSet{indexSet, iSubband}), 3);
	infoAmplitude{iSubband} = mean(cat(3, infoAmplitudeSet{indexSet, iSubband}), 3);
	powerAmplitude{iSubband} = mean(cat(3, powerAmplitudeSet{indexSet, iSubband}), 3);
end
save('../data/re_subband.mat');

%% * R-E plots for frequnecy-flat IRS
figure('name', 'R-E region vs number of subbands for frequency-flat IRS');
legendString = cell(1, length(Variable.nSubbands) + 2);
plotHandle = gobjects(1, length(Variable.nSubbands) + 2);
for iSubband = 1 : length(Variable.nSubbands)
    plotHandle(iSubband) = plot(reAdaptiveIrs{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reAdaptiveIrs{iSubband}(2, :));
    legendString{iSubband} = sprintf('PS: $N = %d$', Variable.nSubbands(iSubband));
	hold on;
end

% * Optimal strategy for medium number of subbands (TS + PS)
subbandIndex = 4;
optIndex = convhull(transpose([0, reAdaptiveIrs{subbandIndex}(1, :) / Variable.nSubbands(subbandIndex); 0, 1e6 * reAdaptiveIrs{subbandIndex}(2, :)])) - 1;
optIndex = optIndex(2 : end - 1);
plotHandle(iSubband + 1) = plot(reAdaptiveIrs{subbandIndex}(1, optIndex) / Variable.nSubbands(subbandIndex), 1e6 * reAdaptiveIrs{subbandIndex}(2, optIndex), 'r');
legendString{iSubband + 1} = sprintf('TS + PS: $N = %d$', Variable.nSubbands(subbandIndex));
hold on;

% * Optimal strategy for large number of subbands (TS)
subbandIndex = 5;
plotHandle(iSubband + 2) = plot([reAdaptiveIrs{subbandIndex}(1, end) / Variable.nSubbands(subbandIndex), reAdaptiveIrs{subbandIndex}(1, 1) / Variable.nSubbands(subbandIndex)], [1e6 * reAdaptiveIrs{iSubband}(2, end), 1e6 * reAdaptiveIrs{iSubband}(2, 1)], 'k');
legendString{iSubband + 2} = sprintf('TS: $N = %d$', Variable.nSubbands(subbandIndex));

hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);

apply_style(plotHandle);
savefig('../figures/re_subband_ff_irs.fig');
matlab2tikz('../../assets/re_subband_ff_irs.tex');
close;

%% * R-E plots for frequency-selective IRS
figure('name', 'R-E region vs number of subbands for frequency-selective IRS');
legendString = cell(1, length(Variable.nSubbands) + 2);
plotHandle = gobjects(1, length(Variable.nSubbands) + 2);
for iSubband = 1 : length(Variable.nSubbands)
    plotHandle(iSubband) = plot(reFsIrs{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reFsIrs{iSubband}(2, :));
    legendString{iSubband} = sprintf('PS: $N = %d$', Variable.nSubbands(iSubband));
	hold on;
end

% * Optimal strategy for medium number of subbands (TS + PS)
subbandIndex = 4;
optIndex = convhull(transpose([0, reFsIrs{subbandIndex}(1, :) / Variable.nSubbands(subbandIndex); 0, 1e6 * reFsIrs{subbandIndex}(2, :)])) - 1;
optIndex = optIndex(2 : end - 1);
plotHandle(iSubband + 1) = plot(reFsIrs{subbandIndex}(1, optIndex) / Variable.nSubbands(subbandIndex), 1e6 * reFsIrs{subbandIndex}(2, optIndex), 'r');
legendString{iSubband + 1} = sprintf('TS + PS: $N = %d$', Variable.nSubbands(subbandIndex));
hold on;

% * Optimal strategy for large number of subbands (TS)
subbandIndex = 5;
plotHandle(iSubband + 2) = plot([reFsIrs{subbandIndex}(1, end) / Variable.nSubbands(subbandIndex), reFsIrs{subbandIndex}(1, 1) / Variable.nSubbands(subbandIndex)], [1e6 * reFsIrs{iSubband}(2, end), 1e6 * reFsIrs{iSubband}(2, 1)], 'k');
legendString{iSubband + 2} = sprintf('TS: $N = %d$', Variable.nSubbands(subbandIndex));

hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);

apply_style(plotHandle);
savefig('../figures/re_subband_fs_irs.fig');
matlab2tikz('../../assets/re_subband_fs_irs.tex');
close;

%% * Waveform amplitude
figure('name', 'Sorted waveform amplitude vs number of subbands');
waveformPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	stem(1 : Variable.nSubbands(iSubband), infoAmplitude{iSubband}, 'marker', 'o');
    hold on;
    stem(1 : Variable.nSubbands(iSubband), powerAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 inf]);
    xticks(1 : Variable.nSubbands(iSubband));
	hold off;
	grid on;
	if iSubband == 1
		legend('$s_I$', '$s_P$');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
	end
end
xlabel('Sorted subband index');

savefig('../figures/waveform_subband.fig');
matlab2tikz('../../assets/waveform_subband.tex');
