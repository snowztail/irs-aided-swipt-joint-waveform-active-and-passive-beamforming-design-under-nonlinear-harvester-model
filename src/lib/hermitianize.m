function [hermitianTerm] = hermitianize(term)
    % Function:
    %   - make an object strictly hermitian
    %
    % Input:
    %   - term: something not strictly hermitian
    %
    % Output:
    %   - hermitianTerm: hermitianized term
    %
    % Comment:
    %   - use to solve precision issue
    %
    % Author & Date: Yang (i@snowztail.com) - 17 Jun 20


    hermitianTerm = (1 / 2) * (term + term');

end
