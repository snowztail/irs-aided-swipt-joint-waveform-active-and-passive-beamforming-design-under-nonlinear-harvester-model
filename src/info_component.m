function [rate, rateMonomial, rateExponent] = info_component(noisePower, channelAmplitude, infoAmplitude, infoRatio)
    % Function:
    %   - formulate rate as a function of amplitudes of multicarrier modulated information waveform
    %
    % InputArg(s):
    %   - noisePower: average noise power
    %   - channelAmplitude: amplitude of channel impulse response
    %   - infoAmplitude: optimum amplitude assigned to information waveform
    %   - infoRatio: information splitting ratio
    %
    % OutputArg(s):
    %   - rate: maximum achievable mutual information
    %   - rateMonomial: monomial components of posynomials that contribute to mutual information
    %   - rateExponent: exponent of the mutual information in the geometric mean
    %
    % Comments:
    %   - returns per-subband rate that decreases as the number of subbands increases
    %   - there is a constant term 1 in each posynomial
    %
    % Author & Date: Yang (i@snowztail.com) - 04 Jun 19



    nSubbands = size(channelAmplitude, 1);

    % number of terms (Kn) in the result posynomials (a constant and a monomial)
    nTerms = 2;

    % type of variables
    isKnown = isa(infoAmplitude, 'double');

    % initialize (a constant term 1 exists in each posynomial)
    if isKnown
        % placeholder for actual values (doubles)
        rateMonomial = ones(nSubbands, nTerms);
    else
        % placeholder for CVX variables (expressions)
        rateMonomial = cvx(ones(nSubbands, nTerms));
    end

    for iSubband = 1: nSubbands
        rateMonomial(iSubband, 2) = infoRatio / noisePower * (infoAmplitude(iSubband) ^ 2 * norm(channelAmplitude(iSubband, :)) ^ 2);
    end

    if isKnown
        ratePosynomial = sum(rateMonomial, 2);
        rateExponent = rateMonomial ./ repmat(ratePosynomial, [1 nTerms]);
        rate = log(prod(ratePosynomial)) / log(2);
    else
        rateExponent = NaN;
        rate = NaN;
    end
    % return per-subband rate
    rate = rate / nSubbands;

end
