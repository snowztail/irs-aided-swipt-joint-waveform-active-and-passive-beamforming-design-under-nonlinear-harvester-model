function [current, currentMonomial, currentExponent] = power_component(k2, k4, resistance, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio)
    % Function:
    %   - formulate output current as a function of amplitudes of multicarrier unmodulated (multisine) power waveform and modulated information waveform
    %
    % InputArg(s):
    %   - k2: diode k-parameter
    %   - k4: diode k-parameter
    %   - resistance: antenna resistance
    %   - channelAmplitude: amplitude of channel impulse response
    %   - infoAmplitude: optimum amplitude assigned to information waveform
    %   - powerAmplitude: optimum amplitude assigned to power waveform
    %   - powerRatio: power splitting ratio
    %
    % OutputArg(s):
    %   - current: target posynomial (zDc) that proportional to the output current
    %   - currentMonomial: monomials (g) as components of the target posynomial
    %   - currentExponent: exponent of the target function in the geometric mean
    %
    % Comments:
    %   - only consider the most fundamental nonlinear model (i.e. truncate at the fourth order)
    %
    % Author & Date: Yang (i@snowztail.com) - 04 Jun 19



    nSubbands = size(channelAmplitude, 1);
    % number of terms in each expression
    nTermsI2 = nSubbands;
    nTermsI4 = nSubbands ^ 2;
    nTermsP2 = nSubbands;
    nTermsP4 = nSubbands * (2 * nSubbands ^ 2 + 1) / 3;
    nTermsI2P2 = nSubbands ^ 2;

    % type of variables
    isKnown = isa(infoAmplitude, 'double');

    % initialize
    if isKnown
        % placeholder for actual values (doubles)
        currentMonomialI2 = zeros(1, nTermsI2);
        currentMonomialI4 = zeros(1, nTermsI4);
        currentMonomialP2 = zeros(1, nTermsP2);
        currentMonomialP4 = zeros(1, nTermsP4);
        currentMonomialI2P2 = zeros(1, nTermsI2P2);
    else
        % placeholder for CVX variables (expressions)
        currentMonomialI2 = cvx(zeros(1, nTermsI2));
        currentMonomialI4 = cvx(zeros(1, nTermsI4));
        currentMonomialP2 = cvx(zeros(1, nTermsP2));
        currentMonomialP4 = cvx(zeros(1, nTermsP4));
        currentMonomialI2P2 = cvx(zeros(1, nTermsI2P2));
    end

    % monomials related to the expectation of time-average of information (modulated) signal to the second order
    iTermI2 = 0;
    for iSubband = 1: nSubbands
        iTermI2 = iTermI2 + 1;
        currentMonomialI2(iTermI2) = norm(channelAmplitude(iSubband, :)) ^ 2 * infoAmplitude(iSubband) ^ 2;
    end
    clearvars iSubband;

    % monomials related to the expectation of time-average of information (modulated) signal to the fourth order
    iTermI4 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            iTermI4 = iTermI4 + 1;
            currentMonomialI4(iTermI4) = ...
                (infoAmplitude(iSubband1) ^ 2 * norm(channelAmplitude(iSubband1, :)) ^ 2) * ...
                (infoAmplitude(iSubband2) ^ 2 * norm(channelAmplitude(iSubband2, :)) ^ 2);
        end
    end
    clearvars iSubband1 iSubband2;

    % monomials related to the time-average of power (multisine) signal to the second order
    iTermP2 = 0;
    for iSubband = 1: nSubbands
        iTermP2 = iTermP2 + 1;
        currentMonomialP2(iTermP2) = norm(channelAmplitude(iSubband, :)) ^ 2 * powerAmplitude(iSubband) ^ 2;
    end
    clearvars iSubband;

    % monomials related to the time-average of power (multisine) signal to the fourth order
    iTermP4 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            for iSubband3 = 1: nSubbands
                iSubband4 = iSubband1 + iSubband2 - iSubband3;
                isValid = iSubband4 >= 1 && iSubband4 <= nSubbands;
                if isValid
                    iTermP4 = iTermP4 + 1;
                    currentMonomialP4(iTermP4) = ...
                        (powerAmplitude(iSubband1) * norm(channelAmplitude(iSubband1, :))) * ...
                        (powerAmplitude(iSubband2) * norm(channelAmplitude(iSubband2, :))) * ...
                        (powerAmplitude(iSubband3) * norm(channelAmplitude(iSubband3, :))) * ...
                        (powerAmplitude(iSubband4) * norm(channelAmplitude(iSubband4, :)));
                else
                    continue;
                end
            end
        end
    end
    clearvars iSubband1 iSubband2 iSubband3 iSubband4;

    % monomials related to the combined power-info terms
    iTermI2P2 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            iTermI2P2 = iTermI2P2 + 1;
            currentMonomialI2P2(iTermI2P2) = ...
                (norm(channelAmplitude(iSubband1, :)) ^ 2 * powerAmplitude(iSubband1) ^ 2) * ...
                (norm(channelAmplitude(iSubband2, :)) ^ 2 * infoAmplitude(iSubband2) ^ 2);
        end
    end
    clearvars iSubband1 iSubband2;

    currentMonomialI2 = currentMonomialI2 * 0.5 * k2 * powerRatio * resistance;
    currentMonomialI4 = currentMonomialI4 * 0.75 * k4 * powerRatio ^ 2 * resistance ^ 2;
    currentMonomialP2 = currentMonomialP2 * 0.5 * k2 * powerRatio * resistance;
    currentMonomialP4 = currentMonomialP4 * 0.375 * k4 * powerRatio ^ 2 * resistance ^ 2;
    currentMonomialI2P2 = currentMonomialI2P2 * 1.5 * k4 * powerRatio ^ 2 * resistance ^ 2;

    % group monomials
    currentMonomial = [currentMonomialI2 currentMonomialI4 currentMonomialP2 currentMonomialP4 currentMonomialI2P2];

    if isKnown
        current = sum(currentMonomial);
        currentExponent = currentMonomial / current;
    else
        current = NaN;
        currentExponent = NaN;
    end

end
