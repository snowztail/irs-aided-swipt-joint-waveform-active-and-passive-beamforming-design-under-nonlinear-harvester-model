function [rate, rateMonomial, rateExponent] = rate_gp(channelAmplitude, infoWaveform, infoRatio, noisePower)
    % Function:
    %   - formulate rate as a function of waveform amplitudes
    %   - decompose rate as sum of monomials
    %
    % Input:
    %   - channelAmplitude (a) [nSubbands * nTxs * nRxs]: channel amplitude response
    %   - infoWaveform (s_I) [nSubbands]: amplitude on information carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - noisePower (\sigma_n^2): average noise power
    %
    % Output:
    %   - rate (R): maximum achievable rate
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
    isKnown = isa(infoWaveform, 'double');

    % * Constant term 1 in rate expression
    if isKnown
        rateMonomial = ones(nSubbands, nTerms);
    else
        rateMonomial = cvx(ones(nSubbands, nTerms));
    end

    % * SNR term in rate expression
    for iSubband = 1 : nSubbands
        rateMonomial(iSubband, 2) = infoRatio * (infoWaveform(iSubband) ^ 2 * norm(channelAmplitude(iSubband, :)) ^ 2) / noisePower;
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
