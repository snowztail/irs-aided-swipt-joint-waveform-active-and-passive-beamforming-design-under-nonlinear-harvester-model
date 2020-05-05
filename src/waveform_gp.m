function [infoWaveform, powerWaveform, rate, current] = waveform_gp(k2, k4, resistance, txPower, noisePower, rateConstraint, tolerance, compositeChannel)
    % Function:
    %   - optimize the information and power waveform to maximize the R-E region
    %   - compute the output DC current
    %
    % InputArg(s):
    %   - k2: diode k-parameters
    %   - k4: diode k-parameters
    %   - resistance: antenna resistance
    %   - txPower: average transmit power
    %   - noisePower: average noise power
    %   - rateConstraint: rate constraint per subband
    %   - tolerance: minimum current gain per iteration
    %   - compositeChannel [\boldsymbol{H}] (nSubbands * 1): total composite channel
    %
    % OutputArg(s):
    %   - infoWaveform: [\boldsymbol{w}_I] (nSubbands * 1): complex information weight
    %   - powerWaveform: [\boldsymbol{w}_P] (nSubbands * 1): complex power weight
    %   - rate: maximum rate per subband
    %   - current: average output current
    %
    % Comment(s):
    %   - the phase is matched to the composite channel
    %
    % Author & Date: Yang (i@snowztail.com) - 29 Apr 20



    nSubbands = size(compositeChannel, 1);
    isConverged = false;
    isSolvable = true;
    current = NaN;
    rate = NaN;
    [channelAmplitude, infoAmplitude, powerAmplitude, infoRatio, powerRatio] = initialize(txPower, compositeChannel);
    [~, ~, rateExponent] = info_component(noisePower, channelAmplitude, infoAmplitude, infoRatio);
    [~, ~, currentExponent] = power_component(k2, k4, resistance, channelAmplitude, infoAmplitude, powerAmplitude, powerRatio);

    while ~isConverged && isSolvable
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

        isConverged = (current_ - current) <= tolerance;
        rate = rate_;
        current = current_;
    end
    infoWaveform = infoAmplitude .* exp(- 1i * angle(compositeChannel));
    powerWaveform = powerAmplitude .* exp(- 1i * angle(compositeChannel));

end
