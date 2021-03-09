function [cascadedChannel] = cascaded_channel(incidentChannel, reflectiveChannel)
    % Function:
    %   - cascade the AP-IRS and IRS-UE channels
    %
    % Input:
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %
    % Output:
    %   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %
    % Comment:
    %   - \boldsymbol{h}_{D,n}^H = directChannel(iSubband, :)
    %   - \boldsymbol{H}_{I,n} = permute(incidentChannel(iSubband, :, :), [2 3 1])'
    %   - \boldsymbol{h}_{R,n}^H = reflectiveChannel(iSubband, :)
    %   - \boldsymbol{\Phi} = diag(irs')
    %   - \boldsymbol{h}_n^H = compositeChannel(iSubband, :)
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20


    % * Get data
    [nSubbands, nTxs, nReflectors] = size(incidentChannel);

    % * Construct AP-IRS-UE cascaded channel
	cascadedChannel = zeros(nReflectors, nTxs, nSubbands);
	for iSubband = 1 : nSubbands
		cascadedChannel(:, :, iSubband) = diag(reflectiveChannel(iSubband, :)) * permute(incidentChannel(iSubband, :, :), [2 3 1])';
	end

end
