%% * Initialize script
clear; close all; clc; config_reflectors;

% K * N
[directChannel] = channel_tgn_e(directPathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
% L * N
[incidentChannel] = channel_tgn_e(incidentPathloss, nReflectors, nSubbands, nUsers, carrierFrequency, fadingType);
% K * NL
reflectiveChannel = zeros(nUsers, nSubbands * nReflectors);
for iReflector = 1 : nReflectors
    [reflectiveChannel(iReflector : nReflectors : end)] = channel_tgn_e(reflectivePathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
end
% N * 1
[irs] = irs_selective(directChannel, incidentChannel, reflectiveChannel);
% K * N
[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
