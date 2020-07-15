function [compositeChannel, concatChannel, concatSubchannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs)
    % Function:
    %   - obtain the composition of direct channel and IRS-aided channel
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - irs (\phi) [nReflectors]: IRS reflection coefficient
    %
    % Output:
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: superposition of direct and extra channels
    %   - concatChannel (V) [nReflectors * (nTxs * nSubbands)]: AP-IRS-user concatenated channel
    %   - concatSubchannel (V_n) {nSubbands}[nReflectors * nTxs]: AP-IRS-user concatenated subchannel
    %
    % Comment:
    %   - \boldsymbol{h}_{D,n}^H = directChannel(iSubband, :)
    %   - \boldsymbol{H}_{I,n} = permute(incidentChannel(iSubband, :, :), [2 3 1])'
    %   - \boldsymbol{h}_{R,n}^H = reflectiveChannel(iSubband, :)
    %   - \boldsymbol{\Phi} = diag(irs')
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    % * Get data
    [nSubbands, nTxs, ~] = size(incidentChannel);

    % * Construct IRS-aided extra channel and AP-IRS-user concatenated channel
    extraChannel = zeros(nSubbands, nTxs);
    concatSubchannel = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        concatSubchannel{iSubband} = diag(reflectiveChannel(iSubband, :)) * permute(incidentChannel(iSubband, :, :), [2 3 1])';
        extraChannel(iSubband, :) = irs' * concatSubchannel{iSubband};
    end
    concatChannel = cat(2, concatSubchannel{:});

    % * Combine for composite channel
    compositeChannel = directChannel + extraChannel;

end
