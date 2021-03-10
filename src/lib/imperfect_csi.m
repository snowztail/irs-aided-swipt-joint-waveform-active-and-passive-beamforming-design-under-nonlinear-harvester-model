function [imperfectChannel] = imperfect_csi(channel, errorVariance)
    % Function:
    %   - emulate imperfect channel with estimation error
    %
    % Input:
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel frequency response
	%	- errorVariance (\epsilon ^ 2): variance of estimation error
    %
    % Output:
    %   - imperfectChannel (\hat{h}) [nSubbands * nTxs * nRxs]: channel estimation
    %
    % Comment:
    %   - assume the estimation error is independent to the channel and the entries follow zero-mean CSCG distribution
    %
    % Author & Date: Yang (i@snowztail.com) - 10 Mar 21


	imperfectChannel = channel + sqrt(errorVariance / 2) * (randn(size(channel)) + 1i * randn(size(channel)));

end
