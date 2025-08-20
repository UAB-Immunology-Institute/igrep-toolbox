function r = ranklin(obj, pop, varargin)
    % IGREP/RANKLIN Lineages and their sizes ordered by size
    %
    % r = RANKLIN(obj, pop) returns an nx2 table where the first column
    % contains the lineage ID and the second column contains the number of
    % sequences in that lineage for the specified population.  Rows are 
    % sorted by lineage size (second column)
    %
    % r = RANKLIN(obj, pop, iso) returns the same output except that the
    % number of sequences for each lineage is restricted to a particular
    % isotype.
    %
    % r = RANKLIN(obj, pop, m) returns an nx2 table as above, except that
    % it based on a random sub-sample of sequences from the specified
    % population.
    %
    % r = RANKLIN(obj, pop, iso, m) returns an nx2 table based on a random
    % sub-sample of seqeunces of a particular isotype.
    %
    % Dot notation usage:
    %     r = D.RANKLIN('mypop')
    %     r = D.RANKLIN('mypop', 'A')
    %     r = D.RANKLIN('mypop', 10000)
    %     r = D.RANKLIN('mypop', 'G', 15000)
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatepop(obj, pop);   
    if nargin == 3
        if isnumeric(varargin{1})                    
            validateattributes(varargin{1}, {'numeric'},...
                {'nonempty', 'scalar', '>', 0});
            n = varargin{1};
            subsamp = 1;
            iso = '';
        else
            subsamp = 0;
            validateiso(obj, varargin{1});
            iso = [' AND isotype = "' char(varargin{1}) '"'];
        end
    elseif nargin == 4
            subsamp = 1;
            validateiso(obj, varargin{1});
            validateattributes(varargin{2}, {'numeric'},...
                {'nonempty', 'scalar', '>', 0});
            n = varargin{2};
            iso = [' AND isotype = "' char(varargin{1}) '"'];
    elseif nargin == 2
            subsamp = 0;
            iso = '';
    end
    if subsamp == 0
        q = [...
            'SELECT lineageID AS lineage, count(*) AS size ',...
            'FROM Sequence WHERE ',...
            'population = "' char(pop) '"' iso,...
            ' GROUP BY lineageID ORDER BY size'];
    else
        q = [...
            'SELECT lineageID AS lineage, count(*) AS size ',...
            'FROM Sequence WHERE ',...
            'sequenceID IN ',...
                '(SELECT sequenceID FROM Sequence WHERE ',...
                'population = "' char(pop) '"' iso,...
                ' ORDER BY random() LIMIT ' num2str(n),...
            ') GROUP BY lineageID ORDER BY size'];
    end        
    r = query(obj, q);
end
