function [quantizedIrs] = quantized_irs(irs, nQuantizeBits)
    % Function:
    %   - quantize each reflection coefficient by b bits
    %
    % Input:
	%   - irs (\phi) [nReflectors * 1]: IRS reflection coefficients
	%	- nQuantizeBits (b): number of bits used for quantization
    %
    % Output:
    %   - quantizedIrs [nReflectors * 1]: quantized IRS reflection coefficients
    %
    % Comment:
    %   - precision 1 / 2^b, require 2^b impedances
    %
    % Author & Date: Yang (i@snowztail.com) - 10 Mar 21


	quantizedIrs = exp(1i * floor(angle(irs) * 2 ^ nQuantizeBits) / 2 ^ nQuantizeBits);

end
