function [compositeChannel] = composite_channel(directChannel, cascadedChannel, irs)
    % Function:
    %   - obtain the composition of direct channel and IRS-aided channel
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
	%   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %   - irs (\phi) [nReflectors * 1]: IRS reflection coefficient
    %
    % Output:
    %   - compositeChannel (h) [nSubbands * nTxs]: equivalent composite channel
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
    [~, nTxs, nSubbands] = size(cascadedChannel);

    % * IRS-aided auxiliary channel
	auxiliaryChannel = zeros(nSubbands, nTxs);
	for iSubband = 1 : nSubbands
		auxiliaryChannel(iSubband, :) = irs' * cascadedChannel(:, :, iSubband);
	end

    % * Composite channel
    compositeChannel = directChannel + auxiliaryChannel;

end
