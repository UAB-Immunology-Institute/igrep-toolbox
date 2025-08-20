function r = lintable(obj, varargin)
    % IGREP/LINTABLE Lineage table for specified populations
    %
    % t = LINTABLE(obj) returns a lineage table (as a Matlab table) for all
    % populations in a database.  The output is a matlab table where the
    % rows are lineages and the columns are populations.  The first column
    % of the table is the lineageID.
    %
    % t = LINTABLE(obj, poplist) returns a lineage table for just the
    % populations specified by the second argument - a cell array of
    % strings.
    %
    % t = LINTABLE(obj, 1) returns a lineage table for all populations in
    % the database broken down by isotype.  The column headings are a
    % concatenation of unique population/isotype combination.
    %
    % t = LINTABLE(obj, poplist, 1) returns a lineage table for the
    % specified populations broken down by isotype.
    %
    % t = LINTABLE(obj, 2) returns a lineage table for all populations in
    % the database.  Two additional columns in the table (2nd and 3rd) are
    % vgene annotations and jgene annotations.
    %
    % t = LINTABLE(obj, poplist, 2) returns a lineage table for the
    % specified populations.  Two additional columns in the table are vgene
    % and jgene annotations.
    %
    % t = LINTABLE(obj, 3) returns a lineage table for all populations in
    % the database broken down by isotype.  Two additional columns in the
    % table are vgene annotations and jgene annotations.
    %
    % t = LINTABLE(obj, poplist, 3) returns a lineage table for the
    % specified populations in the database broken down by isotype.  Two
    % additional columns are vgene annotations and jgene annotations.
    %
    % Summary of flags (optional last input argument):
    %     1 = break down table by isotype
    %     2 = add v and j gene annotations as return values
    %     3 = break down by isotype and return v and j genes
    %
    % Dot notation usage:
    %     t = D.LINTABLE
    %     t = D.LINTABLE({'pop1', 'pop2'})
    %     t = D.LINTABLE(1)
    %     t = D.LINTABLE({'pop1', 'pop2'}, 1)
    %     t = D.LINTABLE(2)
    %     t = D.LINTABLE({'pop1', 'pop2'}, 2)

    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.

    isoflag = ...
        (nargin == 3 && ...
        (isequal(varargin{2}, 1) || isequal(varargin{2}, 3))) ...
        || ...
        (nargin == 2 && ...
        (isequal(varargin{1}, 1) || isequal(varargin{1}, 3)));                        
    vjflag = ...
        (nargin == 3 && ...
        (isequal(varargin{2}, 2) || isequal(varargin{2}, 3))) ...
        || ...
        (nargin == 2 && ...
        (isequal(varargin{1}, 2) || isequal(varargin{1}, 3)));                   
    popflag = ...
        (nargin == 3 && (iscell(varargin{1}) || isstring(varargin{1}))) || ...
        (nargin == 2 && (iscell(varargin{1}) || isstring(varargin{1}))); 
    if popflag
        pop = cellstr(varargin{1}(:)');
        validatepoplist(obj, pop);
        w = [' WHERE population IN ("',...
            strjoin(pop, '", "') '")'];  
    else
        w = '';
        pop = cellstr(obj.upop');
    end              
    n = length(pop);    
    if vjflag
        a = query(obj,...
            ['SELECT distinct(lineageID) AS lineage, ',...
            'vgene AS v, jgene AS j FROM Sequence' w]);  
        [ulin, linsort] = sort(a.lineage);
        vgene = a.v(linsort);
        jgene = a.j(linsort);
    else 
        a = query(obj,...
            ['SELECT distinct(lineageID) AS lineage FROM Sequence' w]);
        ulin = sort(a.lineage);
    end    
    if isoflag                                            
        iso = obj.uiso;
        ni = length(iso);
        mat = zeros(length(ulin), n * ni);
        popiso = cell(1, n * ni);
        for i = 1:n
            for j = 1:ni
                p1 = sortrows(ranklin(obj, pop{i}, iso{j}), 'lineage');
                [~, ~, idx] = intersect(p1{:, 'lineage'}, ulin);
                k = ((i - 1) * ni) + j;
                mat(idx, k) = p1{:, 'size'};
                popiso(k) = cellstr([pop{i} '_' iso{j}]);
            end
        end
        pop = popiso;
    else
        mat = zeros(length(ulin), n);
        for i = 1:n
            p1 = sortrows(ranklin(obj, pop{i}), 'lineage');
            [~, ~, idx] = intersect(p1{:, 'lineage'}, ulin);
            mat(idx, i) = p1{:, 'size'};
        end
    end                        
    [~, isort] = sort(sum(mat, 2), 'descend');
    r = table;
    r.lineage = ulin(isort);
    if vjflag
        r.v = vgene(isort);
        r.j = jgene(isort);
    end
    r = [r, array2table(mat(isort, :), 'VariableNames', pop)];
end
