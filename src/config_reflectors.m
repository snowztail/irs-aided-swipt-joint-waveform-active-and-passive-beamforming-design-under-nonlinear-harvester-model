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
% average transmit and receive power
txPower = 10;
% average noise power
noisePower = db2pow(- 40 - 30);

%% * Channel
% AP-user distance
directDistance = 10;
incidentDistance = 1;
reflectiveDistance = directDistance - incidentDistance;
% pathlosses
[directPathloss] = large_scale_fading("direct", directDistance);
[incidentPathloss] = large_scale_fading("incident", incidentDistance);
[reflectivePathloss] = large_scale_fading("reflective", reflectiveDistance);
% center frequency
centerFrequency = 5.18e9;
% bandwidth
bandwidth = 1e6;
% number of frequency bands
nSubbands = 4;
% channel fading type ('flat' or 'selective')
fadingType = 'selective';
% carrier frequency
[carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
% IRS gain on each reflecting element
irsGain = db2pow(3);

%% * Algorithm
% rate constraint per subband
rateConstraint = 0: 10;
% minimum current gain per iteration
tolerance = 1e-8;

%% * Variables
% number of reflecting elements in IRS
Variable.nReflectors = 5 : 5 : 25;
