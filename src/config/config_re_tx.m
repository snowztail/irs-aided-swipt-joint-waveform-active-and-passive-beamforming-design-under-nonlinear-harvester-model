%% * Transceiver
% diode k-parameter
k2 = 0.0034;
k4 = 0.3829;
% antenna resistance
resistance = 50;
% scale ratio of SMF
alpha = 2;
% coefficients on current terms
beta2 = k2 * resistance;
beta4 = k4 * resistance ^ 2;
% number of receive antennas
nRxs = 1;
% number of users
nUsers = 1;
% average transmit power
txPower = db2pow(10);
% average noise power
noisePower = db2pow(-70);
% IRS element gain
irsGain = 1;
% receive antenna gain
rxGain = db2pow(3);

%% * Channel
% AP-user distance
directDistance = 12;
% vertical distance from the IRS to the AP-user path
verticalDistance = 0.2;
% projection of AP-IRS distance to the AP-user path
horizontalDistance = 0.2;
% AP-IRS and IRS-user distance
[incidentDistance, reflectiveDistance] = coordinate(directDistance, verticalDistance, horizontalDistance);
% center frequency
centerFrequency = 2.4e9;
% bandwidth
bandwidth = 1e6;
% number of frequency bands
nSubbands = 16;
% channel fading mode ('flat' or 'selective')
fadingMode = 'selective';
% carrier frequency
[subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
% number of reflecting elements in IRS
nReflectors = 20;
% spatial correlation
corRx = eye(nRxs);
corIrs = eye(nReflectors);

%% * Algorithm
% minimum gain per iteration
tolerance = 1e-8;
% number of CSCG random vectors to generate
nCandidates = 1e3;
% number of samples in R-E curves
nSamples = 20;
% number of channel realizations
nChannels = 1;

%% * Variable
% number of transmit antennas
Variable.nTxs = [1, 2, 4];

%% * PBS
% number of individual jobs
nBatches = 3e2;
