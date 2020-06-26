clear; clc; setup; config_reflector; load('data/tap.mat');

% * Construct channels
[directChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, directTapGain, directTapDelay, "direct");
[incidentChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, incidentDistance, incidentTapGain, incidentTapDelay, "incident");
[reflectiveChannel_] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, reflectiveDistance, reflectiveTapGain, reflectiveTapDelay, "reflective");

%% ! IRS: R-E region vs number of IRS elements
ffSample = cell(length(Variable.nReflectors), 1);
for iReflector = 1 : length(Variable.nReflectors)
    % * Update channels
    nReflectors = Variable.nReflectors(iReflector);
    incidentChannel = incidentChannel_(:, :, 1 : nReflectors);
    reflectiveChannel = reflectiveChannel_(:, 1 : nReflectors, :);
    ff_sdr;
    ffSample{iReflector} = ffSdrSample;
end
