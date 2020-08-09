function [pathloss] = path_loss(distance, linkMode)
    % Function:
    %   - calculate the large-scale signal attenuation
    %
    % Input:
    %   - distance (d): distance between the transmitter and the receiver
    %   - linkMode: link mode 'direct', 'incident', or 'reflective'
    %
    % Output:
    %   - pathloss (\Lambda): large-scale channel strength reduction
    %
    % Comment:
    %   - assume direct link is blocked while incident and reflective links are not blocked
    %   - assume 3 dBi gain at each IRS element
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20


    switch linkMode
    case 'direct'
        pathloss = db2pow(-30) * distance ^ (-3.2);
    case 'incident'
        pathloss = db2pow(-30) * distance ^ (-2.2) * db2pow(3);
    case 'reflective'
        pathloss = db2pow(-30) * distance ^ (-2.6) * db2pow(3);
    end

end
