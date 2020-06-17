clear; clc; setup;

%% * Transceiver
% number of transmit and receive antennas
nTxs = 1;
nRxs = 1;
% max number of reflecting elements in IRS
nReflectors = 25;

%% * Tap response
[directTapGain, directTapDelay] = tap_tgn(nTxs, nRxs);
[incidentTapGain, incidentTapDelay] = tap_tgn(nTxs, nReflectors);
[reflectiveTapGain, reflectiveTapDelay] = tap_tgn(nReflectors, nRxs);

clearvars nTxs nRxs nReflectors;
save('data/tap.mat');
