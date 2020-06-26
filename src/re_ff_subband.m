clear; clc; setup; config_subband; load('data/tap.mat');

% * R-E region vs number of subbands
niSample = cell(length(Variable.nSubbands), 1);
ffSample = cell(length(Variable.nSubbands), 1);
for iSubband = 1 : length(Variable.nSubbands)
    % * Update channels
    nSubbands = Variable.nSubbands(iSubband);
    [subbandFrequency] = subband_frequency(centerFrequency, bandwidth, nSubbands);
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

    ni_sdr;
    niSample{iSubband} = niSdrSample;

    ff_sdr;
    ffSample{iSubband} = ffSdrSample;
end
