function r = linreport(obj, lin)
    % IGREP/LINREPORT Brief summary report for a lineage
    %
    % r = LINREPORT(obj, lin) returns a data structure with information
    % about the specified lineage, including database name, v-gene and
    % j-gene calls.  It also includes a Matlab table with summary
    % information: number of sequences, numbers of unique V(D)J and CDR3
    % nucleotide and amino acid sequences, as well as average total,
    % non-silent and silent mutations. These summaries are also broken
    % down by isotype and/or population.
    %
    % Dot notation usage:
    %     r = D.LINREPORT(12345)
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    v = query(obj, ['SELECT distinct(VGene) FROM Sequence ',...
        'WHERE lineageID = ' num2str(lin)]);
    j = query(obj, ['SELECT distinct(JGene) FROM Sequence ',...
        'WHERE lineageID = ' num2str(lin)]);
    v = v{1, 1};
    j = j{1, 1};
    qm = [...
        'count(*) AS sequences, ',...
        'count(distinct(NT_VDJ)) AS uniqueVDJnt, ',...
        'count(distinct(NT_CDR3)) AS uniqueCDR3nt, ',...
        'count(distinct(AA_VDJ)) AS uniqueVDJaa, ',...
        'count(distinct(AA_CDR3)) AS uniqueCDR3aa, ',...
        'avg(100.0 * V_Mutations / V_Nucleotides) AS avgMutation, ',...
        'avg(100.0 * V_SilentMutations / V_Nucleotides) AS avgSilentMutation, ',...
        'avg(100.0 * V_NonsilentMutations / V_Nucleotides) AS avgNonsilentMutation ',...
        'FROM Sequence ',...
        'WHERE lineageID = ' num2str(lin)];                
    at = [array2table(["any", "any"], 'VariableNames', {'population', 'isotype'}), query(obj, ['SELECT ' qm])];
    ap = query(obj, [...
        'SELECT population, ' qm ' GROUP BY population ',...
        'ORDER BY population']);
    ap.isotype = repmat("any", size(ap, 1), 1);
    ap = movevars(ap, 'isotype', 'After', 'population');
    ai = query(obj, [...
        'SELECT isotype, ' qm ' GROUP BY isotype ',...
        'ORDER BY isotype']);
    ai.population = repmat("any", size(ai, 1), 1);
    ai = movevars(ai, 'population', 'Before', 'isotype');
    aa = query(obj, [...
        'SELECT population, isotype, ' qm ' GROUP BY ',...
        'population, isotype ORDER BY population, isotype']);    
    a = [at; ap; ai; aa];
    r.datasource = obj.datasource;
    r.lineage = lin;
    r.vgene = v;
    r.jgene = j;
    r.summary = a;
end    
