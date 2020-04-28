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
txPower = db2pow(- 20 - 30);
rxPower = db2pow(- 20 - 30);
% SNR
snrDb = 20;
% average noise power
noisePower = rxPower / db2pow(snrDb);
% distances
directDistance = 10;
incidentDistance = 9;
reflectiveDistance = directDistance - incidentDistance;
% pathlosses
[directPathloss] = large_scale_fading(directDistance);
[incidentPathloss] = large_scale_fading(incidentDistance);
[reflectivePathloss] = large_scale_fading(reflectiveDistance);

%% * Channel
% center frequency
centerFrequency = 5.18e9;
% bandwidth
bandwidth = 1e6;
% number of frequency bands
nSubbands = 8;
% channel fading type ('flat' or 'selective')
fadingType = 'selective';
% carrier frequency
[carrierFrequency] = carrier_frequency(centerFrequency, bandwidth, nSubbands);
% number of reflecting elements in IRS
nReflectors = 10;

%% * Algorithm
% rate constraint per subband
rateConstraint = 0: 0.5: 10;
% minimum current gain per iteration
tolerance = 5e-8;

%% * Variables
