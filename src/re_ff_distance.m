clear; clc; setup; config_distance; load('data/tap.mat');

%% ! IRS: R-E region vs AP-IRS distance
ffSample = cell(length(Variable.incidentDistance), 1);
for iDistance = 1 : length(Variable.incidentDistance)
    % * Update channels
    incidentDistance = Variable.incidentDistance(iDistance);
    reflectiveDistance = directDistance - incidentDistance;
    [directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
    [incidentChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
    [reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");
    ff_sdr;
    ffSample{iDistance} = ffSdrSample;
end
