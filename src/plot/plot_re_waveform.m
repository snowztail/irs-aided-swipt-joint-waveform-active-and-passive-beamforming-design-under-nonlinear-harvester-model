clear; clc; close all; config_re_waveform;

%% * Load data
load('../data/re_waveform.mat', 'reIrsSolution', 'reNoIrsSolution');

witIrsInfoAmplitude = cell(1, length(Variable.nSubbands));
witNoIrsInfoAmplitude = cell(1, length(Variable.nSubbands));
witIrsPowerAmplitude = cell(1, length(Variable.nSubbands));
witNoIrsPowerAmplitude = cell(1, length(Variable.nSubbands));

wptIrsInfoAmplitude = cell(1, length(Variable.nSubbands));
wptNoIrsInfoAmplitude = cell(1, length(Variable.nSubbands));
wptIrsPowerAmplitude = cell(1, length(Variable.nSubbands));
wptNoIrsPowerAmplitude = cell(1, length(Variable.nSubbands));

noIrsChannelAmplitude = cell(1, length(Variable.nSubbands));
witIrsChannelAmplitude = cell(1, length(Variable.nSubbands));
wptIrsChannelAmplitude = cell(1, length(Variable.nSubbands));

for iSubband = 1 : length(Variable.nSubbands)
	witIrsInfoAmplitude{iSubband} = sort(reIrsSolution{iSubband}{1}.infoAmplitude);
	witNoIrsInfoAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{1}.infoAmplitude);
	witIrsPowerAmplitude{iSubband} = sort(reIrsSolution{iSubband}{1}.powerAmplitude);
	witNoIrsPowerAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{1}.powerAmplitude);

	wptIrsInfoAmplitude{iSubband} = sort(reIrsSolution{iSubband}{2}.infoAmplitude);
	wptNoIrsInfoAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{2}.infoAmplitude);
	wptIrsPowerAmplitude{iSubband} = sort(reIrsSolution{iSubband}{2}.powerAmplitude);
	wptNoIrsPowerAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{2}.powerAmplitude);

	noIrsChannelAmplitude{iSubband} = sort(abs(reNoIrsSolution{iSubband}{2}.channel));
	witIrsChannelAmplitude{iSubband} = sort(abs(reIrsSolution{iSubband}{1}.compositeChannel));
	wptIrsChannelAmplitude{iSubband} = sort(abs(reIrsSolution{iSubband}{2}.compositeChannel));
end

%% * Channel amplitude
figure('name', 'Sorted channel amplitude with and without IRS vs number of subbands', 'position', [0, 0, 500, 400]);
channelPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	bar(1 : Variable.nSubbands(iSubband), [noIrsChannelAmplitude{iSubband}, witIrsChannelAmplitude{iSubband}, wptIrsChannelAmplitude{iSubband}]);
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 1e-2]);
    xticks(1 : Variable.nSubbands(iSubband));
	yticks([0 1e-2]);
	hold off;
	grid on;
    if iSubband == 1
		legend('No IRS', 'IRS: WIT', 'IRS: WPT', 'location', 'northoutside', 'orientation', 'horizontal');
	elseif iSubband == 3
		ylabel('Channel amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/channel_amplitude.fig');
matlab2tikz('../../assets/channel_amplitude.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}, ' 'scaled y ticks=false, ' 'yticklabel=\pgfkeys{/pgf/number format/.cd,fixed,precision=2}\pgfmathprintnumber{\tick}']);
close;

%% * Waveform amplitude
% * WIT modulated
figure('name', 'Sorted WIT modulated amplitude with and without IRS vs number of subbands', 'position', [0, 0, 500, 400]);
witInfoPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'loose');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), witIrsInfoAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), witNoIrsInfoAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	yticks([0 2 4]);
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'location', 'northoutside', 'orientation', 'horizontal');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wit_modulated.fig');
matlab2tikz('../../assets/waveform_wit_modulated.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
close;

% * WPT modulated
figure('name', 'Sorted WPT modulated amplitude with and without IRS vs number of subbands', 'position', [0, 0, 500, 400]);
wptInfoPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), wptIrsInfoAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), wptNoIrsInfoAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	yticks([0 2 4]);
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'location', 'northoutside', 'orientation', 'horizontal');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wpt_modulated.fig');
matlab2tikz('../../assets/waveform_wpt_modulated.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
close;

% * WPT multisine
figure('name', 'Sorted WPT multisine amplitude with and without IRS vs number of subbands', 'position', [0, 0, 500, 400]);
wptPowerPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), wptIrsPowerAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), wptNoIrsPowerAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	yticks([0 2 4]);
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{P}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{P}}$', 'location', 'northoutside', 'orientation', 'horizontal');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wpt_multisine.fig');
matlab2tikz('../../assets/waveform_wpt_multisine.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
