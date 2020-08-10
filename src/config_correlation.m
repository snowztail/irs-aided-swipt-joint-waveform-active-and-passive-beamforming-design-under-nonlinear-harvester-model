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
% channel distances and angles
[incidentDistance, reflectiveDistance, directAzimuth, incidentAzimuth, reflectiveAzimuth] = coordinate(directDistance, verticalDistance, horizontalDistance);
% center frequency
centerFrequency = 5.18e9;
% bandwidth
bandwidth = 1e6;
% wavelength
wavelength = 3e8 / centerFrequency;
% number of frequency bands
nSubbands = 4;
% channel fading mode ('flat' or 'selective')
fadingMode = 'selective';
% carrier frequency
[subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
% number of reflecting elements in IRS
nReflectors = 20;

%% * Taps
load('data/variable.mat');
% extract uncorrelated tap gains
directVariable = directVariable(:, :, 1 : nRxs, 1 : nTxs);
incidentVariable = incidentVariable(:, :, 1 : nReflectors, 1 : nTxs);
reflectiveVariable = reflectiveVariable(:, :, 1 : nRxs, 1 : nReflectors);
% extract LOS matrices
directLosMatrix = directLosMatrix(1 : nRxs, 1 : nTxs);
incidentLosMatrix = incidentLosMatrix(1 : nReflectors, 1 : nTxs);
reflectiveLosMatrix = reflectiveLosMatrix(1 : nRxs, 1 : nReflectors);

% * No spatial correlation
corTx = eye(nTxs);
corRx = eye(nRxs);
corIrs = eye(nReflectors);

[directUncorrelatedTapGain, directTapDelay] = tap_tgn(corTx, corRx, directLosMatrix, directVariable, 'nlos');
[incidentUncorrelatedTapGain, incidentTapDelay] = tap_tgn(corTx, corIrs, incidentLosMatrix, incidentVariable, 'nlos');
[reflectiveUncorrelatedTapGain, reflectiveTapDelay] = tap_tgn(corIrs, corRx, reflectiveLosMatrix, reflectiveVariable, 'nlos');

% * Spatial correlation
[directCorTx, directCorRx] = spatial_correlation(wavelength, nTxs, nRxs, directAzimuth, 'transmitter', 'receiver');
[incidentCorTx, incidentCorIrs] = spatial_correlation(wavelength, nTxs, nReflectors, incidentAzimuth, 'transmitter', 'irs');
[reflectiveCorIrs, reflectiveCorRx] = spatial_correlation(wavelength, nReflectors, nRxs, reflectiveAzimuth, 'irs', 'receiver');

[directCorrelatedTapGain, ~] = tap_tgn(directCorTx, directCorRx, directLosMatrix, directVariable, 'nlos');
[incidentCorrelatedTapGain, ~] = tap_tgn(incidentCorTx, incidentCorIrs, incidentLosMatrix, incidentVariable, 'nlos');
[reflectiveCorrelatedTapGain, ~] = tap_tgn(reflectiveCorIrs, reflectiveCorRx, reflectiveLosMatrix, reflectiveVariable, 'nlos');

%% * Algorithm
% minimum gain per iteration
tolerance = 1e-8;
% number of CSCG random vectors to generate
nCandidates = 1e3;
% number of samples in R-E curves
nSamples = 40;
