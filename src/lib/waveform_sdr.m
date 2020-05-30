function [infoWaveform, powerWaveform, infoRatio, powerRatio, rate, current] = waveform_sdr(k2, k4, resistance, txPower, noisePower, rateConstraint, tolerance, compositeChannel, infoWaveform, powerWaveform, infoRatio, powerRatio)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %   - compute the output DC current and user rate
    %
    % Input:
    %   - k2: diode k-parameters
    %   - k4: diode k-parameters
    %   - resistance (R_ant): antenna resistance
    %   - txPower (P): average transmit power
    %   - noisePower (\sigma_n^2): average noise power
    %   - rateConstraint (\bar{R}): average output DC current constraint
    %   - tolerance (\epsilon): minimum rate gain ratio per iteration
    %   - compositeChannel (h) [nSubbands * nTxs * nRxs]: total composite channel
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers (previous solution)
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers (previous solution)
    %   - infoRatio (\bar{\rho}): information splitting ratio (previous solution)
    %   - powerRatio (\rho): power splitting ratio (previous solution)
    %
    % Output:
    %   - infoWaveform (w_I) [nSubbands]: weight on information carriers
    %   - powerWaveform (w_P) [nSubbands]: weight on power carriers
    %   - infoRatio (\bar{\rho}): information splitting ratio
    %   - powerRatio (\rho): power splitting ratio
    %   - rate: user rate
    %   - current: average output DC current
    %
    % Comment:
        % TODO
    %
    % Author & Date: Yang (i@snowztail.com) - 21 May 20



    [nSubbands, ~, ~] = size(compositeChannel);


















    isConverged = false;
    current = NaN;
    rate = NaN;
    [channelAmplitude, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize(txPower, compositeChannel);
    [~, ~, rateExponent] = info_component(noisePower, channelAmplitude, infoAmplitude, infoRatio);
    [~, ~, currentExponent] = power_component(k2, k4, resistance, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

    counter = 0;
    while ~isConverged
        counter = counter + 1;
        cvx_begin gp
            cvx_solver mosek
            variable t
            variable infoAmplitude(nSubbands, 1)
            variable powerAmplitude(nSubbands, 1)
            variable infoRatio
            variable powerRatio
            % formulate expressions
            [~, rateMonomial, ~] = info_component(noisePower, channelAmplitude, infoAmplitude, infoRatio);
            [~, currentMonomial, ~] = power_component(k2, k4, resistance, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

            minimize (1 / t)
            subject to
                0.5 * (infoAmplitude' * infoAmplitude + powerAmplitude' * powerAmplitude) <= txPower;
                2 ^ (nSubbands * rateConstraint) * prod(prod((rateMonomial ./ rateExponent) .^ (- rateExponent))) <= 1;
                t * prod((currentMonomial ./ currentExponent) .^ (- currentExponent)) <= 1;
                infoRatio + powerRatio <= 1;
        cvx_end

        if cvx_status == "Solved"
            [rate_, ~, rateExponent] = info_component(noisePower, channelAmplitude, infoAmplitude, infoRatio);
            [current_, ~, currentExponent] = power_component(k2, k4, resistance, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);
        else
            break;
        end

        isConverged = (current_ - current) <= tolerance || counter >= 1e2;
        rate = rate_;
        current = current_;
    end
    infoWaveform = infoAmplitude .* exp(- 1i * angle(compositeChannel));
    powerWaveform = powerAmplitude .* exp(- 1i * angle(compositeChannel));

end
