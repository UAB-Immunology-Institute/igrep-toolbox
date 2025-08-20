function r = cdr3len(obj, pop, seq, varargin)
    % IGREP/CDR3LEN CDR3 lengths for sequences of specified pop
    %
    % r = CDR3LEN(obj, pop, f) returns CDR3 nucleotide or amino acid
    % lengths for sequences in the specified population.  The second
    % argument must be 'nt' for nucleotide or 'aa' for amino acid.  The
    % return value is a Matlab table (for n lineages) including
    % lineageIDs, isotypes, vgenes, jgenes and CDR3 lengths as well as the
    % number of sequences per lineage.
    %
    % r = CDR3LEN(obj, pop, f, filt) returns CDR3 lengths for the specified
    % population under the specified filter conditions.  The optional
    % argument is a three-element cell array of strings: the first element
    % is the isotype, the second is the vgene and the third is the jgene -
    % values to limit the sequences by.  For the vgene and jgene, if
    % immediately preceded by a "!", the filter will return all sequences
    % but those with the specified value.  For the isotype, if you want to
    % match an isotype prefix, add a wildcard ("%") to the end of an
    % isotype - for example, "G" will match G1, G1A, G1B if you want to
    % collapse subtypes into a supertype.   To avoid filtering by one of
    % the parameters, use an empty character array.
    %
    % Dot notation usage:
    %     r = D.CDR3LEN('pop', 'nt')
    %     r = D.CDR3LEN('pop', 'nt', {'A', '', ''})
    %     r = D.CDR3LEN('pop', 'aa', {'A', 'IGHV4-34', ''})
    %     r = D.CDR3LEN('pop', 'aa', {'G%', '', 'IGHJ3'})
    %     r = D.CDR3LEN('pop', 'nt', {'', '!IGHV1-2', ''})    
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin > 3
        filt = makefilt(obj, varargin{1});
    else
        filt = '';
    end            
    validatepop(obj, pop);
    validatetextarg(obj, seq);
    if ~ismember(seq, {'nt', 'aa'})
        error('second argument must be "nt" or "aa"');
    end
    if isequal(char(seq), 'nt')
        f = 'nt_cdr3';
    else
        f = 'aa_cdr3';
    end
    r = query(obj, [...
        'SELECT lineageID AS lineage, isotype AS isotype, ',...
        'vgene AS v, jgene AS j, ',...
        'length(' f ') AS cdr3len_' seq ', ',...
        'count(*) AS n ',...
        'FROM Sequence WHERE population = "' char(pop) '" ' filt ' ',...
        'GROUP BY lineageID']);
end
