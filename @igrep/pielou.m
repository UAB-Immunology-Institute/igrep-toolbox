function r = pielou(obj, pop, varargin)
    % IGREP/PIELOU Pielou evenness value for population
    %
    % r = PIELOU(obj, pop) returns the pielou evenness value for the
    % specified population.  Values range from 0 (very oligoclonal) to 1
    % (very polyclonal).
    %
    % r = PIELOU(obj, pop, iso) returns the evenness value only considering
    % sequences of the specified isotype.
    %
    % Dot notation usage:
    %     r = D.PIELOU('mypop')
    %     r = D.PIELOU('mypop', 'A')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin > 2
        lin = ranklin(obj, pop, varargin{1});
    else
        lin = ranklin(obj, pop);
    end            
    p = lin.size ./ sum(lin.size);      
    r = (-1 * sum(p .* log(p))) / log(length(p));
end
