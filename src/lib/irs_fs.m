function [irs, compositeChannel] = irs_fs(directChannel, incidentChannel, reflectiveChannel)
    % Function:
    %   - adjust frequency-selective IRS to maximize composite channel strength
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %
    % Output:
    %   - irs (\phi) [nReflectors, nSubbands]: IRS reflection coefficients per subband
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: superposition of direct and extra channels
    %
    % Comment:
    %   - align the extra channel with the direct channel at each subband
    %   - for single-user SISO system
    %   - \boldsymbol{h}_{D,n}^H = directChannel(iSubband, :)
    %   - \boldsymbol{H}_{I,n} = permute(incidentChannel(iSubband, :, :), [2 3 1])'
    %   - \boldsymbol{h}_{R,n}^H = reflectiveChannel(iSubband, :)
    %   - \boldsymbol{\Phi}_n = diag(irs(:, iSubband)')
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Apr 20



    % * Get data
    [nSubbands, ~, nReflectors] = size(incidentChannel);

    % * Each IRS element aligns AP-IRS and IRS-user channel per subband
    irs = zeros(nReflectors, nSubbands);
    for iReflector = 1 : nReflectors
        for iSubband = 1 : nSubbands
            irs(iReflector, iSubband) = exp(1i * (angle(conj(directChannel(iSubband))) - angle(incidentChannel(iSubband, :, iReflector)) - angle(conj(reflectiveChannel(iSubband, iReflector)))));
        end
    end

    % * Obtain composite channel
    compositeChannel = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        compositeChannel(iSubband) = directChannel(iSubband) + reflectiveChannel(iSubband, :) * diag(irs(:, iSubband)') * permute(incidentChannel(iSubband, :, :), [2 3 1])';
    end

end
