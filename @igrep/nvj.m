function r = nvj(obj, pop, method, varargin)
    % IGREP/NVJ V-gene and J-gene usage for a population
    %
    % r = NVJ(obj, pop, method) returns a table containing the combined
    % V- and J-gene usage for a given population (first argument). The
    % second argument specifies how the gene usage should be tallied -
    % valid values are "sequence" (where the number corresponding to the
    % V-gene and J-gene for a population is the number of sequences with
    % that combination) or "lineage" (where the number is the number of
    % lineages).  The output is a table where the rows are V-genes and the
    % columns are J-genes (and the V- and J-gene names are the RowNames and
    % VariableNames, respectively of the table).  The J-genes shown are the
    % unique V- or J- genes for the V- and entire database (all
    % populations), so it is possible for a specific row or column
    % (specific V- or J-gene) to have all "zeros".
    %
    % r = NVJ(obj, pop, method, flatten) returns distribution of v-j
    % combinations as a three column table: V-gene, J-gene and counts.  The
    % value of the third argument should be either "v" or "j" specifying
    % how to sort the v-j combinations.
    %
    % Dot notation usage:
    %     r = D.NVJ('mypop', 'sequence')  
    %     r = D.NVJ('mypop', 'lineage')
    %     r = D.NVJ('mypop', 'lineage', 'j')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatetextarg(obj, method);
    validatepop(obj, pop);
    if ~isempty(varargin)
        validatetextarg(obj, varargin{1});
        if ~ismember(varargin{1}, {'v', 'j'})
            error('optional arg must be "v" or "j"');
        end
    end
    switch method
        case 'sequence', s = ' count(*) AS n';
        case 'lineage',  s = ' count(distinct(lineageID)) as n';
        otherwise,       error('must be "sequence" or "lineage"');
    end
    a = query(obj, [...
        'SELECT VGene AS v, JGene AS j,' s, ' FROM Sequence WHERE ',...
        'population = "' char(pop) '"GROUP BY VGene, JGene']);
    v = sort(obj.uv);
    j = sort(obj.uj);
    nv = length(v);
    nj = length(j);
    mat = zeros(nv, nj); 
    [~, idxv] = ismember(a.v, v);
    [~, idxj] = ismember(a.j, j);         
    mat(sub2ind(size(mat), idxv, idxj)) = a.n;
    if ~isempty(varargin)
        if isequal(varargin{1}, 'v')
            vcol = reshape(repmat(v, 1, nj)', nv * nj, 1);
            jcol = repmat(j, nv, 1);
            rc = reshape(mat', nv * nj, 1);                    
        else                
            vcol = repmat(v, nj, 1);
            jcol = reshape(repmat(j, 1, nv)', nv * nj, 1);
            rc = mat(:);                
        end
        r = table;
        r.v = vcol;
        r.j = jcol;
        r.count = rc;
    else
        r = array2table(mat);
        r.Properties.VariableNames = j;
        r.Properties.RowNames = v;
    end
end
