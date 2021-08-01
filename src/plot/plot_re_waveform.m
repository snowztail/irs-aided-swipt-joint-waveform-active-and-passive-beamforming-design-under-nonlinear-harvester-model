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

for iSubband = 1 : length(Variable.nSubbands)
	witIrsInfoAmplitude{iSubband} = sort(reIrsSolution{iSubband}{1}.infoAmplitude);
	witNoIrsInfoAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{1}.infoAmplitude);
	witIrsPowerAmplitude{iSubband} = sort(reIrsSolution{iSubband}{1}.powerAmplitude);
	witNoIrsPowerAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{1}.powerAmplitude);

	wptIrsInfoAmplitude{iSubband} = sort(reIrsSolution{iSubband}{2}.infoAmplitude);
	wptNoIrsInfoAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{2}.infoAmplitude);
	wptIrsPowerAmplitude{iSubband} = sort(reIrsSolution{iSubband}{2}.powerAmplitude);
	wptNoIrsPowerAmplitude{iSubband} = sort(reNoIrsSolution{iSubband}{2}.powerAmplitude);
end

%% * Waveform amplitude
figure('name', 'Sorted WIT/WPT waveform amplitude with and without IRS vs number of subbands', 'position', [0, 0, 500, 400]);
% * WIT modulated
witInfoPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), witIrsInfoAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), witNoIrsInfoAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{I}}$');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wit_modulated.fig');
matlab2tikz('../../assets/waveform_wit_modulated.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);

% * WPT modulated
wptInfoPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), wptIrsInfoAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), wptNoIrsInfoAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{I}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{I}}$');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wpt_modulated.fig');
matlab2tikz('../../assets/waveform_wpt_modulated.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);

% * WPT multisine
wptPowerPlot = tiledlayout(length(Variable.nSubbands), 1, 'tilespacing', 'compact');
for iSubband = 1 : length(Variable.nSubbands)
	nexttile;
	hold all;
	stem(1 : Variable.nSubbands(iSubband), wptIrsPowerAmplitude{iSubband}, 'marker', 'o');
    stem(1 : Variable.nSubbands(iSubband), wptNoIrsPowerAmplitude{iSubband}, 'marker', 'x');
	xlim([0 Variable.nSubbands(iSubband) + 1]);
	ylim([0 sqrt(2 * txPower)]);
    xticks(1 : Variable.nSubbands(iSubband));
	hold off;
	grid on;
    if iSubband == 1
		legend('IRS: {\boldmath${s}$}$_{\mathrm{P}}$', 'No IRS: {\boldmath${s}$}$_{\mathrm{P}}$');
	elseif iSubband == 3
		ylabel('Waveform amplitude');
    end
    box on;
end
xlabel('Sorted subband index');

savefig('../figures/waveform_wpt_multisine.fig');
matlab2tikz('../../assets/waveform_wpt_multisine.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
