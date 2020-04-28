function [largeScaleFading] = large_scale_fading(distance)
    % Function:
    %   - calculate the large-scale fading
    %
    % InputArg(s):
    %   - distance: separation between the transmitter and the receiver
    %
    % OutputArg(s):
    %   - largeScaleFading [\boldsymbol{\Lambda}]: large-scale channel strength reduction
    %
    % Comment(s):
    %   - large scale fading includes path loss and shadowing
    %   - shadowing is not considered in this project
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20



    % * pathloss
    pathlossExponent = 2;
    pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10));

    % * shadowing
    % shadowingSd = 3;
    % shadowing = db2pow(shadowingSd * randn);
    shadowing = 1;

    % * large-scale fading
    largeScaleFading = pathloss * shadowing;

end
