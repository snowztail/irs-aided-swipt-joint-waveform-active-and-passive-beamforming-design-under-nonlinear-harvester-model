clear; clc; close all; config_subband;

%% * Load batch data
reSet = cell(nBatches, length(Variable.nSubbands));
infoAmplitudeSet = cell(nBatches, length(Variable.nSubbands));
powerAmplitudeSet = cell(nBatches, length(Variable.nSubbands));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_subband_%d.mat', iBatch), 'reInstance', 'reSolution');
		reSet(iBatch, :) = reInstance;
		for iSubband = 1 : length(Variable.nSubbands)
			infoAmplitudeSet{iBatch, iSubband} = sort(reSolution{iSubband}{end}.infoAmplitude);
			powerAmplitudeSet{iBatch, iSubband} = sort(reSolution{iSubband}{end}.powerAmplitude);
		end
    catch
        disp(iBatch);
    end
end

%% * Average over batches
reSubband = cell(1, length(Variable.nSubbands));
infoAmplitude = cell(1, length(Variable.nSubbands));
powerAmplitude = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
	reSubband{iSubband} = mean(cat(3, reSet{:, iSubband}), 3);
	infoAmplitude{iSubband} = mean(cat(3, infoAmplitudeSet{:, iSubband}), 3);
	powerAmplitude{iSubband} = mean(cat(3, powerAmplitudeSet{:, iSubband}), 3);
end
save('../data/re_subband.mat');

%% * R-E plots
figure('name', 'R-E region vs number of subbands');
legendString = cell(1, length(Variable.nSubbands));
for iSubband = 1 : length(Variable.nSubbands)
    plot(reSubband{iSubband}(1, :) / Variable.nSubbands(iSubband), 1e6 * reSubband{iSubband}(2, :));
    legendString{iSubband} = sprintf('$N = %d$', Variable.nSubbands(iSubband));
    hold on;
end
hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
savefig('../figures/re_subband.fig');
matlab2tikz('../../assets/re_subband.tex');
close;

%% * Waveform amplitude
figure('name', 'Sorted waveform amplitude vs number of subbands');
amplitudePlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
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
	end
	if iSubband == 3
		ylabel('Waveform amplitude');
	end
end
xlabel('Sorted subband index');
savefig('../figures/waveform_subband.fig');
matlab2tikz('../../assets/waveform_subband.tex');
