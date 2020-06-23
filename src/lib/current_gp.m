function [current, currentMonomial, currentExponent] = current_gp(beta2, beta4, channel, infoWaveform, powerWaveform, powerRatio)
    % Function:
    %   - formulate output DC current as a function of waveform amplitudes
    %   - decompose current as sum of monomials
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - channel (h) [nSubbands * nTxs * nRxs]: channel amplitude response
    %   - infoWaveform (w_I) [nSubbands]: amplitude on information carriers
    %   - powerWaveform (w_P) [nSubbands]: amplitude on power carriers
    %   - powerRatio (\rho): power splitting ratio
    %
    % Output:
    %   - current (z): posynomial proportional to current
    %   - currentMonomial (g_P): GM monomials
    %   - currentExponent (\gamma_P): GM exponents
    %
    % Comment:
    %   - only consider the most fundamental nonlinear model (i.e. truncate to the fourth order)
    %   - for the single-user case, the optimum phases are given by matched filter
    %
    % Author & Date: Yang (i@snowztail.com) - 04 Jun 19


    % * Initialize algorithm
    nSubbands = size(infoWaveform, 1);
    nTermsP2 = nSubbands;
    nTermsP4 = nSubbands * (2 * nSubbands ^ 2 + 1) / 3;
    nTermsI2 = nSubbands;
    nTermsI4 = nSubbands ^ 2;
    nTermsP2I2 = nSubbands ^ 2;

    % * Type of variables
    isKnown = isa(infoWaveform, 'double');
    if isKnown
        currentMonomialP2 = zeros(1, nTermsP2);
        currentMonomialP4 = zeros(1, nTermsP4);
        currentMonomialI2 = zeros(1, nTermsI2);
        currentMonomialI4 = zeros(1, nTermsI4);
        currentMonomialP2I2 = zeros(1, nTermsP2I2);
    else
        currentMonomialP2 = cvx(zeros(1, nTermsP2));
        currentMonomialP4 = cvx(zeros(1, nTermsP4));
        currentMonomialI2 = cvx(zeros(1, nTermsI2));
        currentMonomialI4 = cvx(zeros(1, nTermsI4));
        currentMonomialP2I2 = cvx(zeros(1, nTermsP2I2));
    end

    % * Monomials
    iTermP2 = 0;
    for iSubband = 1: nSubbands
        iTermP2 = iTermP2 + 1;
        currentMonomialP2(iTermP2) = norm(channel(iSubband, :)) ^ 2 * powerWaveform(iSubband) ^ 2;
    end
    clearvars iSubband;

    iTermP4 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            for iSubband3 = 1: nSubbands
                for iSubband4 = 1: nSubbands
                    if iSubband1 + iSubband2 == iSubband3 + iSubband4
                        iTermP4 = iTermP4 + 1;
                        currentMonomialP4(iTermP4) = ...
                            (powerWaveform(iSubband1) * norm(channel(iSubband1, :))) * ...
                            (powerWaveform(iSubband2) * norm(channel(iSubband2, :))) * ...
                            (powerWaveform(iSubband3) * norm(channel(iSubband3, :))) * ...
                            (powerWaveform(iSubband4) * norm(channel(iSubband4, :)));
                    end
                end
            end
        end
    end
    clearvars iSubband1 iSubband2 iSubband3 iSubband4;

    iTermI2 = 0;
    for iSubband = 1: nSubbands
        iTermI2 = iTermI2 + 1;
        currentMonomialI2(iTermI2) = norm(channel(iSubband, :)) ^ 2 * infoWaveform(iSubband) ^ 2;
    end
    clearvars iSubband;

    iTermI4 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            iTermI4 = iTermI4 + 1;
            currentMonomialI4(iTermI4) = ...
                (infoWaveform(iSubband1) ^ 2 * norm(channel(iSubband1, :)) ^ 2) * ...
                (infoWaveform(iSubband2) ^ 2 * norm(channel(iSubband2, :)) ^ 2);
        end
    end
    clearvars iSubband1 iSubband2;

    iTermP2I2 = 0;
    for iSubband1 = 1: nSubbands
        for iSubband2 = 1: nSubbands
            iTermP2I2 = iTermP2I2 + 1;
            currentMonomialP2I2(iTermP2I2) = ...
                (norm(channel(iSubband1, :)) ^ 2 * powerWaveform(iSubband1) ^ 2) * ...
                (norm(channel(iSubband2, :)) ^ 2 * infoWaveform(iSubband2) ^ 2);
        end
    end
    clearvars iSubband1 iSubband2;

    currentMonomialP2 = (1 / 2) * beta2 * powerRatio * currentMonomialP2;
    currentMonomialP4 = (3 / 8) * beta4 * powerRatio ^ 2 * currentMonomialP4;
    currentMonomialI2 = (1 / 2) * beta2 * powerRatio * currentMonomialI2;
    currentMonomialI4 = (3 / 4) * beta4 * powerRatio ^ 2 * currentMonomialI4;
    currentMonomialP2I2 = (3 / 2) * beta4 * powerRatio ^ 2 *  currentMonomialP2I2;

    currentMonomial = [currentMonomialP2 currentMonomialP4 currentMonomialI2 currentMonomialI4 currentMonomialP2I2];

    % * Posynomial and exponents
    if isKnown
        current = sum(currentMonomial);
        currentExponent = currentMonomial / current;
    else
        current = NaN;
        currentExponent = NaN;
    end

end
