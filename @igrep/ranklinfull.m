function r = ranklinfull(obj, pop) 
    % IGREP/RANKLINFULL Lineages and attributes ordered by size
    %
    % r = RANKLINFULL(obj, pop) returns Matlab table containing information
    % about each lineage in a populationp specified by the first input
    % argument.  The output columns are:
    %     lineage ID
    %     V-gene
    %     J-gene
    %     number of sequences
    %     average total mutation frequency
    %     average non-silent mutation frequency
    %     average silent mutation frequency
    %     numbers of sequences per each isotype
    %
    % Dot notation usage:
    %     r = D.RANKLINFULL('mypop')    
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatepop(obj, pop);
    q1 = [...
        'SELECT lineageID AS lineage, isotype, count(*) AS n ',...
        'FROM Sequence WHERE population = "' pop,...
        '" GROUP BY lineage, isotype'];   
    q2 = [...
        'SELECT lineageID AS lineage, vgene AS v, jgene AS j, ',...
        'avg(100.0 * V_Mutations / V_Nucleotides) AS avg, ',...
        'avg(100.0 * V_NonsilentMutations / V_Nucleotides) AS avgNS, ',...
        'avg(100.0 * V_SilentMutations / V_Nucleotides) AS avgS ',...
        'FROM Sequence WHERE population = "' pop,...
        '" GROUP BY lineage ORDER BY lineage'];
    a1 = query(obj, q1);
    a2 = query(obj, q2);
    ui = unique(a1.isotype);
    ni = length(ui);
    ul = sort(unique(a1.lineage));   % can count on order??   
    imat = zeros(length(ul), ni);
    for i = 1:ni
        temp = a1(a1.isotype == ui(i), :);
        [~, idx, ~] = intersect(ul, temp.lineage);
        imat(idx, i) = temp.n;
    end
    ntot = sum(imat, 2);
    r = table;
    r.lineage = ul;
    r.v = a2.v;
    r.j = a2.j;
    r.sequences = ntot;
    r.avgMut = a2.avg;
    r.avgNonSilentMut = a2.avgNS;
    r.avgSilentMut = a2.avgS;
    r = sortrows([r array2table(imat, 'VariableNames', ui')],...
        'sequences', 'descend');
end
