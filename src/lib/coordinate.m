function [incidentDistance, reflectiveDistance, directAzimuth, incidentAzimuth, reflectiveAzimuth] = coordinate(directDistance, verticalDistance, horizontalDistance)
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
    %   - directAzimuth [2 * 1]: azimuth AoD, AoA of direct channel
    %   - incidentAzimuth [2 * 1]: azimuth AoD, AoA of incident channel
    %   - reflectiveAzimuth [2 * 1]: azimuth AoD, AoA of reflective channel
    %
    % Comment:
    %   - fix AP and user, move IRS
    %   - set AP-user path as reference
    %   - note direction of vectors
    %
    % Author & Date: Yang (i@snowztail.com) - 7 Aug 20


    % * Distances
    incidentDistance = sqrt(verticalDistance ^ 2 + horizontalDistance ^ 2);
    reflectiveDistance = sqrt(verticalDistance ^ 2 + (directDistance - horizontalDistance) ^ 2);

    % * AoD
    directAzimuth(1) = 0;
    incidentAzimuth(1) = atan(verticalDistance / horizontalDistance);
    reflectiveAzimuth(1) = - atan(verticalDistance / (directDistance - horizontalDistance));

    % * AoA
    directAzimuth(2) = directAzimuth(1) + pi;
    incidentAzimuth(2) = incidentAzimuth(1) + pi;
    reflectiveAzimuth(2) = reflectiveAzimuth(1) + pi;

end
