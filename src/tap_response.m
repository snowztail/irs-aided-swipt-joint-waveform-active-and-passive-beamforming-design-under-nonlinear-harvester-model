clear; clc; setup;

%% * Dataset parameters
% max number of transmit antennas
nTxs = 10;
% number of receive antenna
nRxs = 1;
% max number of reflecting elements in IRS
nReflectors = 100;
% number of clusters in channel model
nClusters = 4;
% number of taps in channel model
nTaps = 18;

%% * Generate CSCG variables as uncorrelated tap response
directVariable = zeros(nClusters, nTaps, nRxs, nTxs);
incidentVariable = zeros(nClusters, nTaps, nReflectors, nTxs);
reflectiveVariable = zeros(nClusters, nTaps, nRxs, nReflectors);
for iCluster = 1 : nClusters
    for iTap = 1 : nTaps
        directVariable(iCluster, iTap, :, :) = sqrt(1 / 2) * (randn(nRxs, nTxs) + 1i * randn(nRxs, nTxs));
        incidentVariable(iCluster, iTap, :, :) = sqrt(1 / 2) * (randn(nReflectors, nTxs) + 1i * randn(nReflectors, nTxs));
        reflectiveVariable(iCluster, iTap, :, :) = sqrt(1 / 2) * (randn(nRxs, nReflectors) + 1i * randn(nRxs, nReflectors));
    end
end

%% * Generate LOS matrix
directLosMatrix = exp(1i * 2 * pi * rand(nRxs, nTxs));
incidentLosMatrix = exp(1i * 2 * pi * rand(nReflectors, nTxs));
reflectiveLosMatrix = exp(1i * 2 * pi * rand(nRxs, nReflectors));

save('data/variable.mat', 'directVariable', 'incidentVariable', 'reflectiveVariable', 'directLosMatrix', 'incidentLosMatrix', 'reflectiveLosMatrix');
