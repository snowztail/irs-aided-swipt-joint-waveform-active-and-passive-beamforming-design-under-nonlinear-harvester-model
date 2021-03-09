function [idealCompositeChannel] = composite_channel_ideal(directChannel, cascadedChannel)
    % Function:
	%	- obtain the ideal composite channel by frequency-selective IRS
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs]: the AP-user channel
	%   - cascadedChannel (V) [nReflectors * nTxs * nSubbands]: AP-IRS-user concatenated channel
    %
    % Output:
    %   - idealCompositeChannel (h) [nSubbands * nTxs]: ideal equivalent composite channel
    %
	% Comment:
	%   - assume each IRS element has independent reflection coefficient at each subband (DoF = NL)
	%	- for SISO case, the ideal frequency-selective IRS should align the direct and auxiliary channels at all subbands
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Sep 20


	% * Get data
    [~, nTxs, nSubbands] = size(cascadedChannel);

    % * IRS-aided auxiliary channel
	auxiliaryChannel = zeros(nSubbands, nTxs);
	for iSubband = 1 : nSubbands
		auxiliaryChannel(iSubband, :) = sum(abs(cascadedChannel(:, :, iSubband))) .* exp(1i * angle(directChannel(iSubband, :)));
	end

    % * Composite channel
    idealCompositeChannel = directChannel + auxiliaryChannel;

end
