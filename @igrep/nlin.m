function r = nlin(obj, varargin)
    % IGREP/NLIN Number of lineages for populations
    %
    % r = NLIN(obj) returns and nx2 Matlab table where the first column
    % contains the populations represented in the database and the second
    % column contains the corresponding numbers of lineages.
    %
    % r = NLIN(obj, pop) returns the number of lineages for the population
    % specified in the second argument.
    %
    % Dot notation usage:
    %    r = D.NLIN
    %    r = D.NLIN('mypop')   
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin == 1
        r = obj.lineages;
    else
        validatepop(obj, varargin{1});            
        r = obj.lineages.n(obj.lineages.population == varargin{1});
    end
end
