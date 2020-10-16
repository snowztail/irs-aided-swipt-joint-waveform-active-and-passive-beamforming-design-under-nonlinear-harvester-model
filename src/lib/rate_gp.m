function [rate, rateMonomial, rateExponent] = rate_gp(channelAmplitude, infoAmplitude, infoRatio, noisePower)
    % Function:
    %   - formulate rate as a function of waveform amplitudes
    %   - decompose rate as sum of monomials
    %
    % Input:
    %   - channelAmplitude (a) [nSubbands * 1]: channel amplitude response
    %   - infoAmplitude (s_I) [1 * nSubbands]: amplitude on information carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %	- rate (R): achievable sum rate of all subbands
    %   - rateMonomial (g_I): GM monomials
    %   - rateExponent (\gamma_I): GM exponents
    %
    % Comment:
    %   - there is a constant term (i.e. 1) in each posynomial
    %   - for the single-user case, the optimum phases are given by matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 04 Jun 19


    % * Get data
    nSubbands = size(channelAmplitude, 1);
    nTerms = 2;

    % * Type of variables
    isKnown = isa(infoAmplitude, 'double');

    % * Constant term 1 in rate expression
    if isKnown
        rateMonomial = ones(nSubbands, nTerms);
    else
        rateMonomial = cvx(ones(nSubbands, nTerms));
    end

    % * SNR term in rate expression
    for iSubband = 1 : nSubbands
        rateMonomial(iSubband, 2) = infoRatio * (infoAmplitude(iSubband) ^ 2 * channelAmplitude(iSubband) ^ 2) / noisePower;
    end

    if isKnown
        ratePosynomial = sum(rateMonomial, 2);
        rateExponent = rateMonomial ./ repmat(ratePosynomial, [1 nTerms]);
        rate = log2(prod(ratePosynomial));
    else
        rateExponent = NaN;
        rate = NaN;
    end

end
