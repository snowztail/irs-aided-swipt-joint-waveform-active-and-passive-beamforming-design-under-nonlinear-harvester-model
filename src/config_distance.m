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
nTxs = 2;
nRxs = 1;
% number of users
nUsers = 1;
% average transmit power
txPower = 1;
% average noise power
noisePower = db2pow(-60);

%% * Channel
% AP-user distance
directDistance = 10;
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
% number of reflecting elements in IRS
nReflectors = 10;

%% * Algorithm
% minimum gain ratio per iteration
tolerance = 1e-3;
%  number of CSCG random vectors to generate
nCandidates = 1e4;
% number of samples in R-E curves
nSamples = 40;

%% * Variable
% AP-IRS distance
Variable.incidentDistance = 1 : 5;
