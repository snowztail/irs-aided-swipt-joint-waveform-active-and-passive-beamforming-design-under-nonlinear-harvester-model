function [irs] = irs_selective(directChannel, incidentChannel, reflectiveChannel, irsGain)
    % Function:
    %   - adjust intelligent reflecting surface to maximize the composite channel strength
    %   - align the extra channel with the direct channel at each subband
    %
    % InputArg(s):
    %   - directChannel [\boldsymbol{h}_D] (nUsers * nSubbands): direct AP-user channel
    %   - incidentChannel [\boldsymbol{h}_I] (nReflectors * nSubbands): incident AP-IRS channel
    %   - reflectiveChannel [\boldsymbol{h}_R] (nUsers * (nSubbands * nReflectors)): reflective IRS-user channel
    %   - irsGain: IRS gain on each reflecting element
    %
    % OutputArg(s):
    %   - irs [\boldsymbol{\v}] ((nSubbands * nReflectors) * 1): IRS reflection coefficients
    %
    % Comment(s):
    %   - assume frequency selective IRS (impractical)
    %   - IRS matrix is diagonal with every nReflector elements paired (nSubbands pairs in total)
    %   - pairs are different for selective IRS while same for flat IRS
    %
    % Author & Date: Yang (i@snowztail.com) - 28 Apr 20



    [nReflectors, nSubbands] = size(incidentChannel);
    irsPhase = zeros(nReflectors, nSubbands);
    for iReflector = 1 : nReflectors
        for iSubband = 1 : nSubbands
            irsPhase(iReflector, iSubband) = angle(directChannel(iSubband)) - angle(incidentChannel(iReflector, iSubband)) - angle(reflectiveChannel(iReflector + (iSubband - 1) * nReflectors));
        end
    end
    irs = irsGain * exp(1i * vec(irsPhase));

end
