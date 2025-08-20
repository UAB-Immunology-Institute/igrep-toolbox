function cdrh = cdr3lenhist(obj, pop, seq, method, varargin)
    % IGREP/CDR3LENHIST cdr3 length histogram
    %
    % [cdr3, s] = CDR3LENHIST(obj, pop, seq, method) returns a histogram of 
    % cdr3 lengths for a particular population (first argument).  The 
    % second argument specifies whether length is for nucleotides (nt) or 
    % amino acids (aa).  The third argument specifies whether to compute
    % the histogram over sequences ('sequence') or lineages ('lineage'). 
    % The first output is a Matlab table where the first three columns are
    % the bin (start, end and center), and the fourth column contains the 
    % binned cdr3 lengths.  
    %
    % [cdr3, s] = CDR3LENHIST(obj, pop, seq, method, filt) only considers
    % sequences meeting the filter constraints as described in the
    % igrep/cdr3len method.
    %
    % Dot notation usage:
    %     h = D.CDR3LENHIST('pop', 'aa', 'sequence')
    %     h = D.CDR3LENHIST('pop', 'nt', 'lineage', {'G', '', ''})
    %     h = D.CDR3LENHIST('pop', 'aa', 'lineage', {'', '!IGHV1-2', ''}) 
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.

    validatepop(obj, pop);  
    if nargin > 3
        c = cdr3len(obj, pop, seq, varargin{1});
    else
        c = cdr3len(obj, pop, seq);
    end
    if strcmp(seq, 'nt')
        e = 1.5:3:97.5;
        bctr = e(2:end) - 1.5;
    else
        e = .5:1:32.5;
        bctr = e(2:end) - .5;
    end  
    val = c{:, 5};
    if isequal(method, 'sequence')
        val = repelem(val, c.n);
    end
    h = c.Properties.VariableNames{5};
    head = [{'start', 'end', 'center'}, h];  
    cdrh = [e(1:(end - 1))', e(2:end)', bctr', histcounts(val, e)'];  
    cdrh = array2table(cdrh, 'VariableNames', head);

end      