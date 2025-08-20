function r = isseq(t, seq)

% ISSEQ Test whether input is a valid nucleotide or amino acid sequence
%
%    ISSEQ(t, seq) returns a 1 if the input sequence is a valid nucleotide
%    or amino acid sequence.  The first argument is one of "nt" for
%    nucleotide" or "aa" for amino acid.  The second argument is the
%    sequence (or list of sequences).  If the second rgument is list of
%    sequences, then the output is a vector of logicals.  Valid inputs are:
%        character array
%        string
%        vector of strings
%        cell vector of chars
%    This function is case-insensitive.  Ignored characters include: "."
%    (a placeholder for alignments), and "%", "_" (wildcards for SQLite
%    pattern matching).
%
%    Usage:
%
%        r = ISSEQ('nt', 'aacccgctcag')
%        r = ISSEQ('nt', "aacgcgcta")
%        r = ISSEQ('nt', ["accgcgtca", "cctgca"])
%        r = ISSEQ('nt', {'ccgtca', 'ttgcgact'})
%        r = ISSEQ('aa', 'ARRDVY')

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % second arg
    if ~ismember(t, {'nt', 'aa'})
        error('sequence type must be "nt" or "aa"');
    end
        
    % valid values
    if isequal(t, 'nt')        
        valid = 'acgt._%';
    else
        valid = 'arndbceqzghilkmfpstwyv._%';
    end
    
    % test sequence input
    if ischar(seq)
        n = 1;
    elseif isstring(seq)
        n = length(seq);
        if n == 1
            seq = char(seq);
        else
            if ~isvector(seq)
                error('for string list input, must be vector');
            end
            seq = cellstr(seq);
        end
    elseif iscellstr(seq)
        if ~isvector(seq)
            error('for cell array of string input, must be vector');
        end
        n = length(seq);
    else
        error('input must be string, string array, char, cell array of chars');
    end

    % test validity
    if ischar(seq)
        r = isempty(setdiff(lower(unique(char(seq))), valid));
    else   
        r = nan(n, 1);
        for i = 1:n
            r(i) = isempty(setdiff(lower(unique(char(seq{i}))), valid));
        end
    end
    
return


