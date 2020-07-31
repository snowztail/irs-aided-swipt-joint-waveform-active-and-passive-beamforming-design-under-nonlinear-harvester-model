function [blkDiagComponents] = block_diagonal(matrix, nTxs, nSubbands)
    % Function:
    %   - decompose matrix into block-diagonal components
    %
    % Input:
    %   - matrix [(nTxs * nSubbands) * (nTxs * nSubbands)]: target square matrix
    %   - nTxs (M): number of transmit antennas
    %   - nSubbands (N): number of subbands
    %
    % Output:
    %   - blkDiagComponents {2 * nSubbands - 1}[(nTxs * nSubbands) * (nTxs * nSubbands)]: block-diagonal components with null off-diagonal entries
    %
    % Comment:
    %   - for general use
    %
    % Author & Date: Yang (i@snowztail.com) - 16 Jul 20


    blkDiagComponents = cell(2 * nSubbands - 1, 1);
    for blkDiagIndex = - nSubbands + 1 : nSubbands - 1
        blkDiagComponents{-blkDiagIndex + nSubbands} = zeros(nTxs * nSubbands);
        for xIndex = 1 : nSubbands
            yIndex = xIndex - blkDiagIndex;
            isValid = (1 <= yIndex) && (yIndex <= nSubbands);
            if isValid
                blkDiagComponents{-blkDiagIndex + nSubbands}((xIndex - 1) * nTxs + 1 : xIndex * nTxs, (yIndex - 1) * nTxs + 1 : yIndex * nTxs) = matrix((xIndex - 1) * nTxs + 1 : xIndex * nTxs, (yIndex - 1) * nTxs + 1 : yIndex * nTxs);
            end
        end
    end

end
