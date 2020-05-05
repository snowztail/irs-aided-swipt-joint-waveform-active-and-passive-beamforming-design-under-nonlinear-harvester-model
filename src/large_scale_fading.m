function [largeScaleFading] = large_scale_fading(mode, distance)
    % Function:
    %   - calculate the large-scale fading
    %
    % InputArg(s):
    %   - mode: channel mode (direct, incident and reflective)
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
    if mode == "direct"
        pathlossExponent = 3.8;
    elseif mode == "incident" || mode == "reflective"
        pathlossExponent = 2.2;
    end
    pathloss = db2pow(30 + 10 * pathlossExponent * log10(distance / 1));

    % * shadowing
    % shadowingSd = 3;
    % shadowing = db2pow(shadowingSd * randn);
    shadowing = 1;

    % * large-scale fading
    largeScaleFading = pathloss * shadowing;

end
