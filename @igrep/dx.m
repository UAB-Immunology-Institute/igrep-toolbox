function r = dx(obj, x, pop, varargin)
    % IGREP/DX Measure of clonal expandedness for a population
    %
    % r = DX(obj, x, pop) returns D(x) value(s) for each desired "x" and
    % the specified population.  For instance, the D20 (x = 20) is the
    % number of lineages that span the top 20% of sequences; i.e. those
    % from 80 to 100% of all sequences when lineages are ordered by size.
    % A population with a few large, expanded clones will have a small D20
    % and polyclonal populations (e.g. naive cells) will have a large D20.
    % The second argument can be a vector, the the return value will have
    % the same number of elements for the results.
    %
    % r = DX(obj, x, pop, iso) returns the D(x) values only considering
    % sequences of the specified isotype.
    %
    % Dot notation usage:
    %     r = D.DX([20 50], 'mypop')
    %     r = D.DX([20 50], 'mypop', 'A')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin > 3
        lin = ranklin(obj, pop, varargin{1});
    else
        lin = ranklin(obj, pop);
    end
    
    if(isempty(lin.size))
        r = 0;
    else
        r = nan(1, length(x));
        if ~isequal(lin, 'No Data')
            p = 100 * cumsum(lin.size) / sum(lin.size);
            for i = 1:length(x)
                r(i) = length(p) - find(p <= (100 - x(i)), 1, 'last');
            end
        end
    end
end 
