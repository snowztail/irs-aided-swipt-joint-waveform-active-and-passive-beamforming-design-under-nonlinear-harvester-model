clear; clc; setup;

%% * Transceiver
% number of transmit and receive antennas
nTxs = 10;
nRxs = 1;
% max number of reflecting elements in IRS
nReflectors = 100;

%% * NLOS tap response
[directTapGain, directTapDelay] = tap_tgn(nTxs, nRxs, 'nlos');
[incidentTapGain, incidentTapDelay] = tap_tgn(nTxs, nReflectors, 'nlos');
[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(nReflectors, nRxs, 'nlos');
save('data/tap_nlos.mat', 'directTapGain', 'directTapDelay', 'incidentTapGain', 'incidentTapDelay', 'reflectiveTapGain', 'reflectiveTapDelay');

%% * LOS tap response
[incidentTapGain, incidentTapDelay] = tap_tgn(nTxs, nReflectors, 'los');
[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(nReflectors, nRxs, 'los');
save('data/tap_los.mat', 'directTapGain', 'directTapDelay', 'incidentTapGain', 'incidentTapDelay', 'reflectiveTapGain', 'reflectiveTapDelay');
