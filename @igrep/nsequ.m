function r = nsequ(obj, pop)
    % IGREP/NSEQU Number of unique VDJ sequences
    %
    % r = NSEQU(obj, pop) returns the number of unique VDJ sequences for a
    % specified population. Note that this may be an overestimate because
    % of variation in sequence starting positions that otherwise might be
    % the same.
    %
    % Dot notation usage:
    %    r = D.NSEQU('mypop')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatepop(obj, pop);    
    a = query(obj, [...
        'SELECT nt_VDJ FROM Sequence ',...
        'WHERE population = "' char(pop) '"']);
    r = length(unique(a.NT_VDJ));

end  
