function [compositeChannel, concatVector, concatMatrix] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs)
    % Function:
    %   - obtain the composition of direct channel and IRS-aided channel
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - irs (\phi) [nReflectors]: the IRS gain vector
    %
    % Output:
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: the composite channel
    %   - concatVector (M) [(nReflectors + 1) * nSubbands]: concatenated channel vector
    %   - concatMatrix (R_n) [(nReflectors + 1) * (nReflectors + 1)]: rate SDR matrix
    %
    % Comment:
    %   - for frequency-flat IRS
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    [nSubbands, ~, nReflectors] = size(incidentChannel);
    compositeChannel = zeros(size(directChannel));
    concatMatrix = cell(nSubbands, 1);
    concatVector = zeros(nReflectors + 1, nSubbands);
    for iSubband = 1 : nSubbands
        concatChannel = ctranspose(reflectiveChannel(iSubband, :, :) * diag(squeeze(incidentChannel(iSubband, :, :))));
        compositeChannel(iSubband, :, :) = directChannel(iSubband, :, :) + concatChannel' * irs;
        concatMatrix{iSubband} = [concatChannel * concatChannel', concatChannel * directChannel(iSubband, :, :); directChannel(iSubband, :, :)' * concatChannel', directChannel(iSubband, :, :)' * directChannel(iSubband, :, :)];
        concatVector(1 : end - 1, iSubband) = concatChannel;
    end
    concatVector(end, :) = directChannel';

end
