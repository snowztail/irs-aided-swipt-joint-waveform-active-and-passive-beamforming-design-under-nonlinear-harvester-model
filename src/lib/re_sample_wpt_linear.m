function [sample, solution] = re_sample_wpt_linear(beta2, beta4, directChannel, incidentChannel, reflectiveChannel, txPower, noisePower, nCandidates, tolerance)
    % Function:
    %   - optimize the waveform and IRS reflection coefficients based on linear harvetser model to maximize average output DC current
    %
    % Input:
    %   - beta2: coefficients on second-order current terms
    %   - beta4: coefficients on fourth-order current terms
    %   - directChannel (h_D) [nSubbands * nTxs * nRxs]: the AP-user channel
    %   - incidentChannel (h_I) [nSubbands * nTxs * nReflectors]: the AP-IRS channel
    %   - reflectiveChannel (h_R) [nSubbands * nReflectors * nRxs]: the IRS-user channel
    %   - txPower (P): average transmit power budget
    %   - noisePower (\sigma_n^2): average noise power
    %   - nCandidates (Q): number of CSCG random vectors to generate
    %   - tolerance (\epsilon): minimum current gain per iteration
    %
    % Output:
    %   - sample [2 * nSamples]: rate-energy sample
    %   - solution: IRS reflection coefficient, composite channel, waveform, splitting ratio and eigenvalue ratio
    %
    % Comment:
    %   - based on linear harvester model
    %
    % Author & Date: Yang (i@snowztail.com) - 10 Oct 20


    % * Get data
	[nSubbands, ~, nReflectors] = size(incidentChannel);

    % * Initialize IRS and composite channel
    irs = exp(1i * 2 * pi * rand(nReflectors, 1));
    [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);

    % * Initialize waveform and splitting ratio
	infoAmplitude = zeros(1, nSubbands) + eps;
    powerAmplitude = sqrt(2 * txPower / nSubbands) * ones(1, nSubbands);
	[infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);
    infoRatio = eps;
    powerRatio = 1 - infoRatio;

    % * AO
    isConverged = false;
    current_ = 0;
	rateConstraint = 0;
	eigRatio = [];
    while ~isConverged
		[irs, eigRatio(end + 1)] = irs_linear(beta2, directChannel, incidentChannel, reflectiveChannel, irs, infoWaveform, powerWaveform, infoRatio, powerRatio, noisePower, rateConstraint, nCandidates);
		[compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs);
		[infoAmplitude, powerAmplitude, current] = waveform_linear(beta2, beta4, compositeChannel, txPower);
		[infoWaveform, powerWaveform] = precoder_mrt(compositeChannel, infoAmplitude, powerAmplitude);
        isConverged = abs(current - current_) <= tolerance;
        current_ = current;
    end

	sample = [eps; current];
	solution = variables2struct(irs, compositeChannel, infoAmplitude, powerAmplitude, infoRatio, powerRatio, eigRatio);

end
