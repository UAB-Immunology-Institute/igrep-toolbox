function validateseq(t, s, varargin)

% VALIDATESEQ Test a variable to see if it is a valid sequence
%
%    VALIDATESEQ(t, s) throws an error if the input argument is not a valid
%    sequence.  The first argument is either "nt" for nucleotide or "aa"
%    for amino acid sequence.  Valid options are char, string of lists of
%    each.  Examples of valid nucleotide sequences are:
%
%        'aacccgctcag'
%        "aacgcgcta"
%        ["accgcgtca", "cctgca"]
%        {'ccgtca', 'ttgcgact'}
%
%    This function is a wrapper for isseq and suitable for use in
%    addParameter for error checking routines.

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % test for number of expected sequences optional input
    na = NaN;
    if nargin > 2
        na = varargin{1};
        if ~(isnumeric(na) && isscalar(na) && na > 0)
            error('second optional argument must be numeric scalar > 0');
        end
    end 

    % test for valid sequence
    if isequal(t, 'nt'), x = 'nucleotide'; else, x = 'amino acid'; end
    if ~isseq(t, s), error(['input is not a valid ' x ' sequence']); end
    
    % test for specified number of sequences
    if nargin > 2
        if (isstring(s) || iscellstr(s)) && length(s) ~= na
            error(['expecting ' num2str(na) ' sequences']);
        end
    end

return