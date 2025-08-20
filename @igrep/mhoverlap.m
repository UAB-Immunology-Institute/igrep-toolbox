function r = mhoverlap(obj, varargin)
    % IGREP/MHOVERLAP All-by-all Morisita-Horn overlap matrix
    %
    % r = MHOVERLAP(obj) returns square Matlab table containing the values
    % for an all-by-all, pairwise Morisita-Horn overlap matrix.
    %
    % r = MHOVERLAP(obj, iso) restricts the comparisons to sequences of a
    % specified isotype.
    %
    % Dot notaqtion usage:
    %     r = D.MHOVERLAP
    %     r = D.MHOVERLAP('G')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin == 1
        t = lintable(obj);
        tab = t{:, 2:end};
        pop = t.Properties.VariableNames(2:end);
    else                
        validateiso(obj, varargin{1});     
        t = lintable(obj, 1);
        alltab = t{:, 2:end};
        allpop = t.Properties.VariableNames(2:end);        
        jkeep = ~cellfun('isempty', regexp(allpop, ['_' varargin{1} '$']));
        pop = allpop(jkeep);
        tab = alltab(:, jkeep);
    end
    o = morisitahorn(tab); 
    r = array2table(o, 'VariableNames', pop, 'RowNames', pop);
end
