clear; clc; close all; config_re_reflector;

%% * Load batch data
indexSet = 1 : nBatches;
reAoSet = cell(nBatches, length(Variable.nReflectors));
reLcSet = cell(nBatches, length(Variable.nReflectors));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_reflector/re_reflector_%d.mat', iBatch), 'reAoInstance', 'reLcInstance');
        reAoSet(iBatch, :) = reAoInstance;
        reLcSet(iBatch, :) = reLcInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAoReflector = cell(1, length(Variable.nReflectors));
reLcReflector = cell(1, length(Variable.nReflectors));
for iReflector = 1 : length(Variable.nReflectors)
    reAoReflector{iReflector} = mean(cat(3, reAoSet{indexSet, iReflector}), 3);
    reLcReflector{iReflector} = mean(cat(3, reLcSet{indexSet, iReflector}), 3);
end
save('../data/re_reflector.mat');

%% * R-E plots
figure('name', 'R-E region vs number of reflectors', 'position', [0, 0, 500, 400]);
legendString = cell(2, length(Variable.nReflectors));
plotHandle = gobjects(2, length(Variable.nReflectors));
hold all;
for iReflector = 1 : length(Variable.nReflectors)
    plotHandle(1, iReflector) = plot(reAoReflector{iReflector}(1, :) / nSubbands, 1e6 * reAoReflector{iReflector}(2, :));
    plotHandle(2, iReflector) = plot(reLcReflector{iReflector}(1, :) / nSubbands, 1e6 * reLcReflector{iReflector}(2, :));
	legendString{1, iReflector} = sprintf('BCD: $L = %d$', Variable.nReflectors(iReflector));
	legendString{2, iReflector} = sprintf('LC-BCD: $L = %d$', Variable.nReflectors(iReflector));
end
hold off;
grid on;
legend(legendString(:));
xlabel('Per-subband rate [bps/Hz]');
ylabel('Output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_group_style(plotHandle(:), 2);

savefig('../figures/re_reflector.fig');
matlab2tikz('../../assets/re_reflector.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
