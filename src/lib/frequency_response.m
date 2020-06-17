function [directChannel, incidentChannel, reflectiveChannel] = frequency_response(nSubbands, subbandFrequency, fadingMode, nReflectors, directDistance, incidentDistance, reflectiveDistance)
    % Function:
    %   - get frequency response of direct, incident and reflective channels
    %
    % Input:
    %   - nSubbands (N): number of frequency bands
    %   - subbandFrequency (f_n) [nSubbands]: the center frequency of subbands
    %   - fadingMode: fading mode "flat" or "selective"
    %   - nReflectors: number of reflecting elements in IRS
    %   - directDistance: AP-user distance
    %   - incidentDistance: AP-IRS distance
    %   - reflectiveDistance: IRS-user distance
    %
    % Output:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: AP-IRS channel
    %   - directChannel (h_R) [nSubbands * nReflectors * nRxs]: IRS-user channel
    %
    % Comment:
    %   - based on generated tap data
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Jun 20



    load('data/tap.mat');

    % * Direct link
    [directFading] = fading_tgn(directTapGain, directTapDelay, nSubbands, subbandFrequency, fadingMode);
    [directPathloss] = path_loss(directDistance, "direct");
    directChannel = directFading / sqrt(directPathloss);

    % * Incident link
    incidentTapGain = incidentTapGain(:, :, 1 : nReflectors);
    [incidentFading] = fading_tgn(incidentTapGain, incidentTapDelay, nSubbands, subbandFrequency, fadingMode);
    [incidentPathloss] = path_loss(incidentDistance, "incident");
    incidentChannel = incidentFading / sqrt(incidentPathloss);

    % * Reflective link
    reflectiveTapGain = reflectiveTapGain(:, 1 : nReflectors, :);
    [reflectiveFading] = fading_tgn(reflectiveTapGain, reflectiveTapDelay, nSubbands, subbandFrequency, fadingMode);
    [reflectivePathloss] = path_loss(reflectiveDistance, "reflective");
    reflectiveChannel = reflectiveFading / sqrt(reflectivePathloss);

end
