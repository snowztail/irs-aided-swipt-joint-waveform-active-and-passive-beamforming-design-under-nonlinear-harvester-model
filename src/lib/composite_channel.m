function [compositeChannel, concatChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs)
    % Function:
    %   - obtain the composition of direct channel and IRS-aided channel
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - irs (v) [nReflectors]: the IRS gain vector
    %
    % Output:
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: the composite channel
    %   - concatChannel (\Phi) []
    %
    % Comment:
    %   - for frequency-flat IRS
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    [nSubbands, ~, ~] = size(directChannel);

    concatChannel = cell(nSubbands, 1);
    compositeChannel = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        concatChannel{iSubband} = reflectiveChannel(iSubband, :, :) * diag(squeeze(incidentChannel(iSubband, :, :)));
        compositeChannel{iSubband} = directChannel(iSubband, :, :) + concatChannel{iSubband} * irs;
    end

end
