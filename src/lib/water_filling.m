function [streamPower] = water_filling(streamStrength, txPower, noisePower)
    % Function:
    %   - iterative stream power allocation based on water-filling algorithm
    %
    % Input:
    %   - streamStrength: eigenvalues of the channel matrix corresponding to valid non-zero eigenvectors
    %   - txPower: average transmit power
    %   - noisePower: average noise power
    %
    % Output:
    %   - streamPower: power allocated to each stream
    %
    % Restraint:
    %   - this function is based on iterative method which can be inefficient for a large number of streams
    %
    % Comment:
    %   - non-zero eigenvalues correspond to feasible eigenvectors and streams
    %   - favour multiple eigenmode transmission at high SNR and dominant eigenmode transmission at low SNR
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Feb 19



    % calculate SNR
    snr = txPower / noisePower;
    % initialise power allocation
    streamPower_ = zeros(length(streamStrength), 1);
    % obtain number of available streams
    nStreams = length(streamStrength(abs(streamStrength) > eps));
    % sort in descending order for water filling
    [streamStrength, index] = sort(streamStrength, 'descend');
    % calculate the stream base level (status)
    baseLevel = 1 ./ (snr * streamStrength);
    % begin iteration
    for iStream = 1 : nStreams
        % update quasi water level
        waterLevel = 1 / (nStreams - iStream + 1) * (1 + sum(baseLevel(1 : (nStreams - iStream + 1))));
        % try to allocate power with this new water level
        streamPower_(1 : iStream) = waterLevel - baseLevel(1 : iStream);
        % negative allocation is invalid
        isValid = all(streamPower_ >= 0);
        if isValid
            % update power allocation
            streamPower = streamPower_;
        else
            % invalid power allocation, return the latest valid solution
            break;
        end
    end
    % match corresponding eigenvalues
    streamPower = streamPower(index);
    % match transmit power constraint
    streamPower = txPower * streamPower ./ sum(streamPower);

end
