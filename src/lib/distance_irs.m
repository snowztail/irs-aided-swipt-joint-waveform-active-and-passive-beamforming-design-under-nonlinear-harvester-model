function [incidentDistance, reflectiveDistance] = distance_irs(directDistance, verticalDistance, horizontalDistance)
    % Function:
    %   - obtain the AP-IRS (incident) and IRS-user (reflective) distance
    %
    % Input:
    %   - directDistance (d_D): AP-user distance
    %   - verticalDistance (d_V): vertical distance from the IRS to the AP-user path
    %   - horizontalDistance (d_H): projection of AP-IRS distance to the AP-user path
    %
    % Output:
    %   - directDistance (d_I): AP-IRS distance
    %   - reflectiveDistance (d_R): IRS-user distance
    %
    % Comment:
    %   - fix AP and user, move IRS
    %
    % Author & Date: Yang (i@snowztail.com) - 7 Aug 20


    incidentDistance = sqrt(verticalDistance ^ 2 + horizontalDistance ^ 2);
    reflectiveDistance = sqrt(verticalDistance ^ 2 + (directDistance - horizontalDistance) ^ 2);

end
