function r = hist(obj, pop, varargin)
    % IGREP/HIST Lineage size histogram
    %
    % r = HIST(obj, pop) returns a Matlab table where the first column is
    % the size of a lineage (in sequences) and the second column is the
    % number of occurrences of lineage of the corresponding sizes for the
    % specified population.  The sum of the second column is the number of
    % lineages in that population.  The sum of the row products of the
    % result is the number of sequences for that population.  The rows are
    % sorted by the first column, from smallest to largest.
    %
    % r = HIST(obj, pop, iso) returns the same matrix except only
    % considering sequences of a specified isotype.
    %
    % Dot notation usage:
    %     r = D.HIST('mypop')
    %     r = D.HIST('mypop', 'M')      
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin > 2
        a = ranklin(obj, pop, varargin{1});            
    else
        a = ranklin(obj, pop);
    end
    t = tabulate(a.size);
    r = array2table(t(t(:, 2) > 0, 1:2),...
        'VariableNames', {'size', 'occurrences'});
end      
