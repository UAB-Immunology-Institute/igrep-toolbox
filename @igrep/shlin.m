function [p, nc] = shlin(obj, pop1, pop2, varargin)
    % IGREP/SHLIN numbers, percentages shared lineages between pops
    %
    % [prc, ns] = SHLIN(obj, pop1, pop2) returns percentages of shared
    % lineages of each of a pair of specified populations.  The first
    % return value is a 2-element vector: percentage of shared lineages
    % relative to number of lineages in populations one and two.  The
    % second return value is the number of shared lineages between the two
    % populations.
    %
    % [prc, ns] = SHLIN(obj, pop1, pop2, n) returns values as above except
    % ignores shared lineages where the number of sequences in either
    % population is n or below.
    %
    % Dot notation usage:
    %     [prc, nc] = D.SHLIN('pop1', 'pop2')  
    %     [prc, nc] = D.SHLIN('pop1', 'pop2', 1)
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatetextarg(obj, pop1);
    validatetextarg(obj, pop2);
    t = lintable(obj, cellstr({pop1, pop2}));
    m = t{:, 2:3};
    n1 = length(find(m(:, 1) > 0));
    n2 = length(find(m(:, 2) > 0));
    if ~isempty(varargin)
        validateattributes(varargin{1}, {'numeric'},...
            {'nonempty', 'scalar', '>=', 0});
        m(m <= varargin{1}) = 0;
    end
    nc = length(find(prod(m, 2) > 0));
    p = 100 * nc * [1 / n1, 1 / n2];  
end
