function [pathloss] = path_loss(distance)
    % Function:
    %   - calculate the large-scale signal attenuation based on IEEE TGn channel model D
    %
    % Input:
    %   - distance (d): distance between the transmitter and the receiver
    %
    % Output:
    %   - pathloss (\Lambda): large-scale channel strength reduction
    %
    % Comment:
    %   - consists of the free space loss (exponent = 2) up to 10 m and typical urban loss (exponent = 3.5) onwards
    %   - set reference path loss at 1 m as -30 dB
    %
    % Reference:
    %   - V. Erceg et al., "TGn channel models," in Version 4. IEEE 802.11–03/940r4, May 2004.
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20


    if distance <= 10
        pathloss = db2pow(-30) * distance ^ (-2);
    else
        pathloss = db2pow(-30) * 10 ^ (-2) * (distance / 10) ^ (-3.5);
    end

end
