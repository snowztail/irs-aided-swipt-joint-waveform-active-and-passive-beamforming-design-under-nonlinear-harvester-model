function [sumPathloss] = sum_pathloss(directDistance, incidentDistance, reflectiveDistance)
    % Function:
    %   - calculate the sum pathloss of all links
    %
    % Input:
    %   - sumPathloss (\Lambda): the sum pathloss of all links
    %
    % Output:
    %   - directDistance (d_D): AP-user distance
    %   - incidentDistance (d_I): AP-IRS distance
    %   - reflectiveDistance (d_R): IRS-user distance
    %
    % Comment:
    %   - pathloss can be interpreted as fraction of power received at the receiver
    %   - with IRS, the equivalent pathloss is the pathloss of AP-user plus the pathloss of AP-IRS-user path
    %   - don't multiply pathloss of incident and reflective links as pathloss is not a linear function of distance
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Aug 20


    sumPathloss = path_loss(directDistance) + path_loss(incidentDistance + reflectiveDistance);

end
