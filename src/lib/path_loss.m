function [pathloss] = path_loss(distance, linkMode)
    % Function:
    %   - calculate the large-scale fading
    %
    % Input:
    %   - distance (d): distance between the transmitter and the receiver
    %   - linkMode: link mode 'direct', 'incident', or 'reflective'
    %
    % Output:
    %   - pathloss (\Lambda): large-scale channel strength reduction
    %
    % Comment:
    %   - assume direct link is blocked so that its exponent is larger than indirect links
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20


    switch linkMode
    case 'direct'
        exponent = 3.8;
    case {'incident', 'reflective'}
        exponent = 2.2;
    end
    pathloss = db2pow(30 + 10 * exponent * log10(distance / 1));

end
