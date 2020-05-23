%% * Transceiver
% diode k-parameter
k2 = 0.0034;
k4 = 0.3829;
% antenna resistance
resistance = 50;
% number of transmit and receive antennas
nTxs = 1;
nRxs = 1;
% number of users
nUsers = 1;
% average transmit power
txPower = 1;
% average noise power
noisePower = db2pow(- 30);

%% * Channel
% AP-user distance
directDistance = 10;
incidentDistance = 1;
reflectiveDistance = directDistance - incidentDistance;
% center frequency
centerFrequency = 5.18e9;
% bandwidth
bandwidth = 1e6;
% number of frequency bands
nSubbands = 4;
% channel fading mode ("flat" or "selective")
fadingMode = "selective";
% carrier frequency
[subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
% gain on each reflecting element
irsGain = db2pow(3);
% number of reflecting elements in IRS
nReflectors = 5;

%% * Algorithm
% output DC current constraint
currentConstraint = 0;
% minimum rate increase per iteration
tolerance = 1e-8;

%% * Variables
% number of reflecting elements in IRS
Variable.nReflectors = 5 : 5 : 25;
