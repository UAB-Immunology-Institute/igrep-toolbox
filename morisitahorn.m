function o = morisitahorn(c)

% MORISITAHORN Morisita-Horn overlap index between two communities
%
%    MORISITAHORN(c) returns a square matrix containing the Morisita-
%    Horn overlap indices for all possible pairs of repertoire of
%    community comparisons (columns of input matrix).  Each column of the
%    input matrix represents a repertoire and each row represents a species
%    (e.g. a grouped lineage).  The values of the matrix represent the
%    number of individuals (e.g. sequences) for the corresponding species
%    and repertoire.  For a pairwise repertoire comparison, the index
%    varies from 0 (no overlap) to 1 (identical species and corresponding
%    frequencies).

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.
        
    % compute overlap between all pairs of columns
    n = size(c, 2);
    X = sum(c);
    sx = repmat(sum(c .^ 2) ./ (X .^ 2), n, 1);
    num = 2 * (c' * c);
    denom = (sx + sx') .* (X' * X);
    o = num ./ denom;
    
return


