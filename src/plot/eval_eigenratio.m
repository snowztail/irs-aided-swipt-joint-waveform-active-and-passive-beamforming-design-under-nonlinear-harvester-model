%% * Max eigenvalue over sum eigenvalue vs number of transmit antennas
clear; clc; close all;
config_re_tx;
indexSet = 1 : nBatches;
eigRatioAoTx = cell(nBatches, length(Variable.nTxs), nSamples);
eigRatioLcTx = cell(nBatches, length(Variable.nTxs), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_tx/re_tx_%d.mat', iBatch), 'reAoSolution', 'reLcSolution');
		for iTx = 1 : length(Variable.nTxs)
			for iSample = 1 : nSamples
				eigRatioAoTx{iBatch, iTx, iSample} = reAoSolution{iTx}{iSample}.eigRatio;
				eigRatioLcTx{iBatch, iTx, iSample} = reLcSolution{iTx}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
eigRatioTx = zeros(length(Variable.nTxs), 2);
for iTx = 1 : length(Variable.nTxs)
	eigRatioTx(iTx, 1) = min([eigRatioAoTx{:, iTx, :}]);
	eigRatioTx(iTx, 2) = min([eigRatioLcTx{:, iTx, :}]);
end
1 - eigRatioTx

%% * Max eigenvalue over sum eigenvalue vs number of subbands
clear; clc; close all;
config_re_subband;
indexSet = 1 : nBatches;
eigRatioAoSubband = cell(nBatches, length(Variable.nSubbands), nSamples);
eigRatioLcSubband = cell(nBatches, length(Variable.nSubbands), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_subband/re_subband_%d.mat', iBatch), 'reAoSolution', 'reLcSolution');
		for iSubband = 1 : length(Variable.nSubbands)
			for iSample = 1 : nSamples
				eigRatioAoSubband{iBatch, iSubband, iSample} = reAoSolution{iSubband}{iSample}.eigRatio;
				eigRatioLcSubband{iBatch, iSubband, iSample} = reLcSolution{iSubband}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end

eigRatioSubband = zeros(length(Variable.nSubbands), 2);
for iSubband = 1 : length(Variable.nSubbands)
	eigRatioSubband(iSubband, 1) = min([eigRatioAoSubband{:, iSubband, :}]);
	eigRatioSubband(iSubband, 2) = min([eigRatioLcSubband{:, iSubband, :}]);
end
1 - eigRatioSubband

%% * Max eigenvalue over sum eigenvalue vs number of IRS elements
clear; clc; close all;
config_re_reflector;
indexSet = 1 : nBatches;
eigRatioAoReflector = cell(nBatches, length(Variable.nReflectors), nSamples);
eigRatioLcReflector = cell(nBatches, length(Variable.nReflectors), nSamples);
for iBatch = 1 : nBatches
    try
        load(sprintf('../data/re_reflector/re_reflector_%d.mat', iBatch), 'reAoSolution', 'reLcSolution');
		for iReflector = 1 : length(Variable.nReflectors)
			for iSample = 1 : nSamples
				eigRatioAoReflector{iBatch, iReflector, iSample} = reAoSolution{iReflector}{iSample}.eigRatio;
				eigRatioLcReflector{iBatch, iReflector, iSample} = reLcSolution{iReflector}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
eigRatioReflector = zeros(length(Variable.nReflectors), 2);
for iReflector = 1 : length(Variable.nReflectors)
	eigRatioReflector(iReflector, 1) = min([eigRatioAoReflector{:, iReflector, :}]);
	eigRatioReflector(iReflector, 2) = min([eigRatioLcReflector{:, iReflector, :}]);
end
1 - eigRatioReflector

%% * Max eigenvalue over sum eigenvalue vs number of IRS elements
clear; clc; close all;
config_re_irs;
indexSet = 1 : nBatches;
eigRatioAoBandwidth = cell(nBatches, length(Variable.bandwidth), nSamples);
eigRatioLcBandwidth = cell(nBatches, length(Variable.bandwidth), nSamples);
for iBatch = 2 : nBatches
    try
        load(sprintf('../data/re_irs/re_irs_%d.mat', iBatch), 'reAdaptiveIrsSolution', 'reLcSolution');
		for iBandwidth = 1 : length(Variable.bandwidth)
			for iSample = 1 : nSamples
				eigRatioAoBandwidth{iBatch, iBandwidth, iSample} = reAdaptiveIrsSolution{iBandwidth}{iSample}.eigRatio;
				eigRatioLcBandwidth{iBatch, iBandwidth, iSample} = reLcSolution{iBandwidth}{iSample}.eigRatio;
			end
		end
    catch
		indexSet(indexSet == iBatch) = [];
        disp(iBatch);
    end
end
eigRatioBandwidth = zeros(length(Variable.bandwidth), 2);
for iBandwidth = 1 : length(Variable.bandwidth)
	eigRatioBandwidth(iBandwidth, 1) = min([eigRatioAoBandwidth{:, iBandwidth, :}]);
	eigRatioBandwidth(iBandwidth, 2) = min([eigRatioLcBandwidth{:, iBandwidth, :}]);
end
1 - eigRatioBandwidth
