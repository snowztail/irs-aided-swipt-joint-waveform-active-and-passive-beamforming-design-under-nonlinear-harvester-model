clear; clc; close all; config_re_distance;

%% * Load batch data
indexSet = 1 : nBatches;
reAoSet = cell(nBatches, length(Variable.horizontalDistance));
reLcSet = cell(nBatches, length(Variable.horizontalDistance));
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_distance/re_distance_%d.mat', iBatch), 'reAoInstance', 'reLcInstance');
        reAoSet(iBatch, :) = reAoInstance;
        reLcSet(iBatch, :) = reLcInstance;
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

%% * Average over batches
reAoDistance = cell(1, length(Variable.horizontalDistance));
reLcDistance = cell(1, length(Variable.horizontalDistance));
for iDistance = 1 : length(Variable.horizontalDistance)
    reAoDistance{iDistance} = mean(cat(3, reAoSet{indexSet, iDistance}), 3);
    reLcDistance{iDistance} = mean(cat(3, reLcSet{indexSet, iDistance}), 3);
end
save('../data/re_distance.mat');

%% * R-E plots
figure('name', 'R-E region vs AP-IRS horizontal distance', 'position', [0, 0, 500, 400]);
legendString = cell(2, length(Variable.horizontalDistance));
plotHandle = gobjects(2, length(Variable.horizontalDistance));
hold all;
for iDistance = 1 : length(Variable.horizontalDistance)
    plotHandle(1, iDistance) = plot(reAoDistance{iDistance}(1, :) / nSubbands, 1e6 * reAoDistance{iDistance}(2, :));
    plotHandle(2, iDistance) = plot(reLcDistance{iDistance}(1, :) / nSubbands, 1e6 * reLcDistance{iDistance}(2, :));
    legendString{1, iDistance} = sprintf('BCD: $d_H = %s$ m', num2str(Variable.horizontalDistance(iDistance)));
    legendString{2, iDistance} = sprintf('LC-BCD: $d_H = %s$ m', num2str(Variable.horizontalDistance(iDistance)));
end
hold off;
grid on;
% legend(legendString(:));
xlabel('Per-subband rate [bps/Hz]');
ylabel('DC [$\mu$A]');
xlim([0 inf]);
ylim([0 inf]);
box on;
apply_group_style(plotHandle(:), 2);

dummyHandle = gobjects(1, 5);
hold all;
dummyHandle(1) = plot(-1, -1, 'color', '#000000', 'linestyle', '-');
dummyHandle(2) = plot(-1, -1, 'color', '#000000', 'linestyle', '--');
dummyHandle(3) = plot(-1, -1, 'color', '#0072BD', 'linestyle', 'none', 'marker', 'o');
dummyHandle(4) = plot(-1, -1, 'color', '#D95319', 'linestyle', 'none', 'marker', '+');
dummyHandle(5) = plot(-1, -1, 'color', '#EDB120', 'linestyle', 'none', 'marker', 's');
hold off;
legend(dummyHandle, 'BCD', 'LC-BCD', sprintf('$d_{\\mathrm{H}} = %s$ m', num2str(Variable.horizontalDistance(1))), sprintf('$d_{\\mathrm{H}} = %s$ m', num2str(Variable.horizontalDistance(2))), sprintf('$d_{\\mathrm{H}} = %s$ m', num2str(Variable.horizontalDistance(3))));

savefig('../figures/re_distance.fig');
matlab2tikz('../../assets/re_distance.tex', 'extraaxisoptions', ['title style={font=\huge}, ' 'label style={font=\huge}, ' 'ticklabel style={font=\LARGE}, ' 'legend style={font=\LARGE}']);
