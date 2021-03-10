clear; clc; close all; config_re_reflector;

%% * Load batch data
indexSet = 1 : nBatches;
reSet = cell(nBatches, length(Variable.nReflectors));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_reflector/re_reflector_%d.mat', iBatch), 'reInstance');
        reSet(iBatch, :) = reInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reReflector = cell(1, length(Variable.nReflectors));
for iReflector = 1 : length(Variable.nReflectors)
    reReflector{iReflector} = mean(cat(3, reSet{indexSet, iReflector}), 3);
end
save('../data/re_reflector.mat');

%% * R-E plots
figure('name', 'Average R-E region vs number of reflectors', 'position', [0, 0, 500, 400]);
legendString = cell(1, length(Variable.nReflectors));
plotHandle = gobjects(1, length(Variable.nReflectors));
for iReflector = 1 : length(Variable.nReflectors)
    plotHandle(iReflector) = plot(reReflector{iReflector}(1, :) / nSubbands, 1e6 * reReflector{iReflector}(2, :));
	legendString{iReflector} = sprintf('$M = %d, L = %d$', nTxs, Variable.nReflectors(iReflector));
    hold on;
end
hold off;
grid on;
legend(legendString);
xlabel('Average subband rate [bps/Hz]');
ylabel('Average output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

savefig('../figures/re_reflector.fig');
matlab2tikz('../../assets/re_reflector.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
