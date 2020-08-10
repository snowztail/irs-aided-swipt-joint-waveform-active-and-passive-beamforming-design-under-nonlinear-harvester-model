%% * Transceiver
% diode k-parameter
k2 = 0.0034;
k4 = 0.3829;
% antenna resistance
resistance = 50;
% coefficients on current terms
beta2 = k2 * resistance;
beta4 = k4 * resistance ^ 2;
% number of transmit and receive antennas
nTxs = 1;
nRxs = 1;
% number of users
nUsers = 1;
% average transmit power
txPower = db2pow(6);
% average noise power
noisePower = db2pow(-70);

%% * Channel
% AP-user distance
directDistance = 15;
% vertical distance from the IRS to the AP-user path
verticalDistance = 2;
% projection of AP-IRS distance to the AP-user path
horizontalDistance = 2;
% AP-IRS and IRS-user distance
[incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance);
% center frequency
centerFrequency = 5.18e9;
% bandwidth
bandwidth = 1e6;
% number of frequency bands
nSubbands = 4;
% channel fading mode ('flat' or 'selective')
fadingMode = 'selective';
% carrier frequency
[subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);

%% * Algorithm
% minimum gain per iteration
tolerance = 1e-8;
% number of CSCG random vectors to generate
nCandidates = 1e3;
% number of samples in R-E curves
nSamples = 40;

%% * Variable
% number of reflecting elements in IRS
Variable.nReflectors = 0 : 10 : 50;
% tap gains and delays
Variable.directTapGain = cell(length(Variable.nReflectors), 1);
Variable.incidentTapGain = cell(length(Variable.nReflectors), 1);
Variable.reflectiveTapGain = cell(length(Variable.nReflectors), 1);

%% * Taps
load('data/variable.mat');
for iReflector = 1 : length(Variable.nReflectors)
    nReflectors = Variable.nReflectors(iReflector);
    % no spatial correlation
    corTx = eye(nTxs);
    corRx = eye(nRxs);
    corIrs = eye(nReflectors);
    % tap gains and delays
    [Variable.directTapGain{iReflector}, directTapDelay] = tap_tgn(corTx, corRx, directLosMatrix(1 : nRxs, 1 : nTxs), directVariable(:, :, 1 : nRxs, 1 : nTxs), 'nlos');
    [Variable.incidentTapGain{iReflector}, incidentTapDelay] = tap_tgn(corTx, corIrs, incidentLosMatrix(1 : nReflectors, 1 : nTxs), incidentVariable(:, :, 1 : nReflectors, 1 : nTxs), 'nlos');
    [Variable.reflectiveTapGain{iReflector}, reflectiveTapDelay] = tap_tgn(corIrs, corRx, reflectiveLosMatrix(1 : nRxs, 1 : nReflectors), reflectiveVariable(:, :, 1 : nRxs, 1 : nReflectors), 'nlos');
end
clearvars nReflectors directVariable incidentVariable reflectiveVariable directLosMatrix incidentLosMatrix reflectiveLosMatrix;
