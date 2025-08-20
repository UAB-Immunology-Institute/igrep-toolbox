function m = nt2mat(b)

% NT2MAT Convert cell array of nucleotide strings to numeric matrix
%
%    NT2MAT(seq) takes a cell array of strings of nucleotide sequences and
%    returns a numeric matrix, with each number corresponding to a base as
%    follows:
%        A = 1
%        C = 2
%        G = 3
%        T = 4
%    Spaces and "." map to 0.
%
%    Usage:
%
%        a = {...
%                '...AACTG';...
%                '..AAACCGT'};
%
%        b = NT2MAT(a);
%
%        b = [...
%                0 0 0 1 1 2 4 3 0;...
%                0 0 1 1 1 2 2 3 4]

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % replace nucleotides with numbers
    b = regexprep(b, 'a', '1');
    b = regexprep(b, 'c', '2');
    b = regexprep(b, 'g', '3');
    b = regexprep(b, 't', '4');
    b = regexprep(b, '\.', '0');
    b = regexprep(b, 'n', '0');

    % convert cell array of strings to character matrix
    c = char(b);

    % deal with trailing spaces due to sequences of different lengths
    for i = 1:size(c, 1), c(i, :) = regexprep(c(i, :), ' ', '0'); end

    % vectorize
    [ns, nb] = size(c);
    c = c(:);

    % convert to numbers
    d = str2num(c);

    % reshape back into matrix
    m = reshape(d, ns, nb);

return

