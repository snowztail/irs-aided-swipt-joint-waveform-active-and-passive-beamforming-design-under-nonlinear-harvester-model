function [fsIrsCompositeChannel] = fs_irs_composite_channel(directChannel, incidentChannel, reflectiveChannel)
    % Function:
	%	- obtain the optimal composite SISO channel by frequency-selective IRS to maximize the R-E region
    %
    % Input:
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %
    % Output:
    %   - fsIrsCompositeChannel (h) [nSubbands * nTxs * nRxs]: superposition of direct and extra channels
    %
	% Comment:
	%	- the closed form IRS reflection coefficients is obtained for SISO
	%   - assume each IRS element has independent reflection coefficient at each subband (DoF = NL)
	%	- IRS reflection coefficient aligns direct channel and extra channel to maximize composite channel strength
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Sep 20


	% * Get data
	[nSubbands, nTxs, nReflectors] = size(incidentChannel);
	if nTxs ~= 1
		error(sprintf('Sorry, this function is only for SISO.'));
	end

	% * Each IRS element aligns AP-IRS-user channel with AP-user channel at each subband
	irs = zeros(nReflectors, nSubbands);
	for iReflector = 1 : nReflectors
		for iSubband = 1 : nSubbands
			irs(iReflector, iSubband) = exp(1i * (angle(conj(directChannel(iSubband))) - angle(incidentChannel(iSubband, :, iReflector)) - angle(conj(reflectiveChannel(iSubband, iReflector)))));
		end
	end

	% * Obtain AP-IRS-user channel
	extraChannel = zeros(nSubbands, nTxs);
	for iSubband = 1 : nSubbands
		extraChannel(iSubband, :) = irs(:, iSubband)' * diag(reflectiveChannel(iSubband, :)) * permute(incidentChannel(iSubband, :, :), [2 3 1])';
	end

    % * Combine for composite channel
    fsIrsCompositeChannel = directChannel + extraChannel;

end
