function [corTx, corRx] = spatial_correlation(wavelength, nTxs, nRxs, azimuthAngle, txType, rxType)
    % Function:
    %   - obtain spatial correlation matrices of transmitter, receiver and IRS
    %   - obtain LoS matrices of direct, incident and reflective links
    %
    % Input:
    %   - wavelength: wavelength
    %   - nTxs: number of transmit antennas
    %   - nRxs: number of receive antennas
    %   - azimuthAngle [2 * 1]: azimuth AoD, AoA of channel
    %   - txType: transmitter type ('transmitter', 'irs')
    %   - rxType: receiver type ('irs', 'receiver')
    %
    % Output:
    %   - corTx (R_t) [nTxs * nTxs]: transmit spatial correlation matrix
    %   - corRx (R_r) [nRxs * nRxs]: receive spatial correlation matrix
    %
    % Comment:
    %   - assume half wavelength spacing for most adjacent elements
    %   - assume transmit antennas in a row, receive antennas in a row
    %   - assume IRS reflectors in a plane, each row has 5 elements
    %
    % Reference:
    %   - Q.-U.-A. Nadeem, A. Kammoun, M. Debbah, and M.-S. Alouini, “A Generalized Spatial Correlation Model for 3D MIMO Channels Based on the Fourier Coefficients of Power Spectrums,” IEEE Transactions on Signal Processing, vol. 63, no. 14, pp. 3671–3686, 2015
    %
    % Author & Date: Yang (i@snowztail.com) - 9 Aug 20


    azimuthAod = azimuthAngle(1);
    azimuthAoa = azimuthAngle(2);

    % * Relative coordinates
    switch txType
    case 'transmitter'
        txCoordinate = wavelength / 2 * (0 : nTxs - 1);
    case 'irs'
        txCoordinate = zeros(2, nTxs);
        txCoordinate(1, :) = wavelength / 2 * mod((0 : nTxs - 1), 5);
        txCoordinate(2, :) = wavelength / 2 * floor((0 : nTxs - 1) / 5);
    end

    switch rxType
    case 'irs'
        rxCoordinate = zeros(2, nRxs);
        rxCoordinate(1, :) = wavelength / 2 * mod((0 : nRxs - 1), 5);
        rxCoordinate(2, :) = wavelength / 2 * floor((0 : nRxs - 1) / 5);
    case 'receiver'
        rxCoordinate = wavelength / 2 * (0 : nRxs - 1);
    end

    % * Correlation matrices
    corTx = zeros(nTxs);
    for iTx = 1 : nTxs
        for jTx = 1 : nTxs
            corTx(iTx, jTx) = exp(1i * 2 * pi * norm(txCoordinate(:, iTx) - txCoordinate(:, jTx)) / wavelength * sin(azimuthAod));
        end
    end

    corRx = zeros(nRxs);
    for iRx = 1 : nRxs
        for jRx = 1 : nRxs
            corRx(iRx, jRx) = exp(1i * 2 * pi * norm(rxCoordinate(:, iRx) - rxCoordinate(:, jRx)) / wavelength * sin(azimuthAoa));
        end
    end

end
