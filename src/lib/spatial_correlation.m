function [corMatrix] = spatial_correlation(nElements, azimuthAngle, elevationAngle, wavelength, elementType)
    % Function:
    %   - obtain spatial correlation matrices of transmitter, receiver and IRS
    %   - obtain LOS matrices of direct, incident and reflective links
    %
    % Input:
    %   - nElements: number of antennas or reflectors
    %   - azimuthAngle: horizontal angle between Tx-Rx and y-axis
    %   - elevationAngle: vertical angle between Tx-Rx and z-axis
    %   - wavelength: wavelength
    %   - elementType: element type 'transmitter', 'receiver', or 'irs'
    %
    % Output:
    %   - corMatrix (R) [nElements * nElements]: spatial correlation matrix
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


    corMatrix = zeros(nElements);

    % * Set coordinate
    switch elementType
    case {'transmitter', 'receiver'}
        coordinate = wavelength / 2 * (0 : nElements - 1);
    case 'irs'
        coordinate = zeros(2, nElements);
        coordinate(1, :) = wavelength / 2 * mod((1 : nElements), 5);
        coordinate(2, :) = wavelength / 2 * floor((1 : nElements) / 5);
    end

    % * Obtain correlation matrix
    for iElement = 1 : nElements
        for jElement = 1 : nElements
            corMatrix(iElement, jElement) = exp(1i * 2 * pi * norm(coordinate(:, iElement) - coordinate(:, jElement)) * sin(azimuthAngle) * sin(elevationAngle));
        end
    end

end
