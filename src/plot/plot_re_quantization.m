clear; clc; close all; config_re_quantization;

%% * Load batch data
indexSet = 1 : nBatches;
reNoIrsSet = cell(nBatches, 1);
reIrsSet = cell(nBatches, 1);
reQuantizedSet = cell(nBatches, length(Variable.nQuantizeBits));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_quantization/re_quantization_%d.mat', iBatch), 'reNoIrsInstance', 'reIrsInstance', 'reQuantizedInstance');
        reNoIrsSet{iBatch} = reNoIrsInstance;
        reIrsSet{iBatch} = reIrsInstance;
        reQuantizedSet(iBatch, :) = reQuantizedInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reNoIrs = mean(cat(3, reNoIrsSet{indexSet}), 3);
reIrs = mean(cat(3, reIrsSet{indexSet}), 3);
reQuantized = cell(1, length(Variable.nQuantizeBits));
for iBit = 1 : length(Variable.nQuantizeBits)
    reQuantized{iBit} = mean(cat(3, reQuantizedSet{indexSet, iBit}), 3);
end
save('../data/re_quantization.mat');

%% * R-E plots
figure('name', 'R-E region vs IRS quantization bits', 'position', [0, 0, 500, 400]);
legendString = cell(1, length(Variable.nQuantizeBits) + 2);
plotHandle = gobjects(1, length(Variable.nQuantizeBits) + 2);
hold all;
plotHandle(1) = plot(reNoIrs(1, :) / nSubbands, 1e6 * reNoIrs(2, :));
legendString{1} = 'No IRS';
for iBit = 1 : length(Variable.nQuantizeBits)
    plotHandle(iBit + 1) = plot(reQuantized{iBit}(1, :) / nSubbands, 1e6 * reQuantized{iBit}(2, :));
	legendString{iBit + 1} = sprintf('$b = %s$', num2str(Variable.nQuantizeBits(iBit)));
end
plotHandle(end) = plot(reIrs(1, :) / nSubbands, 1e6 * reIrs(2, :));
legendString{end} = 'Unquantized IRS';
hold off;
grid on;
legend(legendString);
xlabel('Per-subband rate [bps/Hz]');
ylabel('Output DC current [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_style(plotHandle);

savefig('../figures/re_quantization.fig');
matlab2tikz('../../assets/re_quantization.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
