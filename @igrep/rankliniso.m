function r = rankliniso(obj, pop)
    % IGREP/RANKLINISO Lineages, sizes and isotype compositions
    %
    % r = RANKLINISO(obj, pop) returns a Matlab table where the first
    % column is lineage ID, the second column is number of sequences in the
    % corresponding lineages for the specified population, and the
    % remaining columns are the numbers of sequences in the corresponding
    % lineages for the specified population for each isotype.  Rows are
    % ordered by lineage size (second column).
    %
    % Dot notation usage:
    %    r = D.RANKLINISO('mypop')
    
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
        'SELECT lineageID, isotype, count(lineageID) AS n ',...
        'FROM Sequence WHERE population = "' char(pop) '" ',...
        'GROUP BY lineageID, isotype']);  
    ulin = unique(a.lineageID);
    ui = sort(obj.uiso);
    mat = zeros(length(ulin), length(ui));
    [~, idxlin] = ismember(a.lineageID, ulin);
    [~, idxiso] = ismember(a.isotype, ui);
    mat(sub2ind(size(mat), idxlin, idxiso)) = a.n;
    r0 = [ulin sum(mat, 2) mat];
    [~, isort] = sort(r0(:, 2));
    r = array2table(r0(isort, :),...
        'VariableNames', [{'lineage', 'total'}, ui']);
    
end
