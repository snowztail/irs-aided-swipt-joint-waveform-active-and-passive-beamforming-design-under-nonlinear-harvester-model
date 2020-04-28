function [compositeChannel] = composite_channel(directChannel, incidentChannel, reflectiveChannel, irs)
    % Function:
    %   - compute the composite channel per subband
    %
    % InputArg(s):
    %   - directChannel [\boldsymbol{h}_D] (nUsers * nSubbands): direct AP-user channel
    %   - incidentChannel [\boldsymbol{h}_I] (nReflectors * nSubbands): incident AP-IRS channel
    %   - reflectiveChannel [\boldsymbol{h}_R] (nUsers * (nSubbands * nReflectors)): reflective IRS-user channel
    %   - irs [\boldsymbol{\v}] ((nSubbands * nReflectors) * 1): IRS reflection coefficients
    %
    % OutputArg(s):
    %   - compositeChannel [\boldsymbol{H}] (nUsers * nSubbands): total composite channel
    %
    % Comment(s):
    %   - the direct channel overlaps with the IRS-aided channel
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Apr 20



    [nUsers, nSubbands] = size(directChannel);
    nReflectors = size(incidentChannel, 1);
    compositeChannel = zeros(nUsers, nSubbands);
    for iSubband = 1 : nSubbands
        phi = reflectiveChannel(:, (iSubband - 1) * nReflectors + 1 : iSubband * nReflectors) * diag(incidentChannel(:, iSubband));
        compositeChannel(:, iSubband) = directChannel(:, iSubband) + phi * irs((iSubband - 1) * nReflectors + 1 : iSubband * nReflectors);
    end

end
