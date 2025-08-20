function r = nvx(obj, poplist, method, fam, varargin)
    % IGREP/NVX V-gene (or family) usage for a list of populations
    %
    % r = NVX(obj, poplist, method, family) returns a table containing the
    % counts of sequences or lineages for each V-gene family (or each
    % V-gene within a V-gene family).  The first argument is a list of
    % populations in the database to consider.  The second argument
    % specifies how the gene usage should be tallied - valid values are
    % "sequence" or "lineage".  The third argument is the V-gene family of
    % V-genes to consider (valid values are the numbers 1-7), or "0" to
    % tally across V-gene families.  Columns are the populations (in the
    % same order as specified in the first input argument).
    %
    % r = NVX(obj, poplist, method, family, 1) returns a table of 
    % frequencies (based on the total repertoires) rather than counts.
    %     
    % Dot notation usage:
    %     r = D.NVX({'mypop1', 'mypop2'}, 'sequence', 3)
    %     r = D.NVX({'mypop1', 'mypop2'}, 'lineage', 0, 1);
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatetextarg(obj, method);
    validatepoplist(obj, poplist);
    if ~isnumeric(fam) || ~isscalar(fam) || ~ismember(fam, 0:7)
        error('family must be scalar integer between 0 and 7');
    end     
    switch method
        case 'sequence', s = 'count(*) as n';
        case 'lineage',  s = 'count(distinct(lineageID)) as n';
        otherwise,       error('method must be "sequence" or "lineage"');
    end      
    nrm = NaN;
    if ~isempty(varargin)
        if ~isscalar(varargin{1}) || ~isnumeric(varargin{1}) || varargin{1} ~= 1
            error('optional arg must 1');
        end        
        if isequal(method, 'lineage')
            [~, ip] = ismember(poplist, obj.nlin.population);
            nrm = obj.nlin.n(ip)';
        else
            [~, ip] = ismember(poplist, obj.nseq.population);
            nrm = obj.nseq.total(ip)';
        end            
    end
    p = ['("' strjoin(poplist, '", "') '")'];
    if fam > 0
        v = ['"IGHV' num2str(fam) '%"'];
        a = query(obj, [...
            'SELECT population, vgene AS v, ' s, ' FROM Sequence WHERE ',...
            'population IN ' p ' AND v LIKE ' v ' '...
            'GROUP BY v, population']);     
        uv = unique(a.v);        
        r = zeros(length(uv), length(poplist));
        [~, idxv] = ismember(a.v, uv);
        [~, idxp] = ismember(a.population, poplist);
        r(sub2ind(size(r), idxv, idxp)) = a.n;   
        if ~isnan(nrm), r = 100 * r ./ repmat(nrm, size(r, 1), 1); end
        r = array2table(r, 'VariableNames', poplist, 'RowNames', uv);
    else
        r = zeros(7, length(poplist));        
        for i = 1:7
            a = query(obj, [...
                'SELECT population, ' s, ' FROM Sequence WHERE ',...
                'population IN ' p ' AND ',...
                'vgene LIKE "IGHV' num2str(i) '%" ',...
                'GROUP BY population ORDER BY population']);    
            [~, idxp] = ismember(a.population, poplist);            
            r(i, idxp) = a.n';
        end
        if ~isnan(nrm), r = 100 * r ./ repmat(nrm, size(r, 1), 1); end
        r = array2table(r, 'VariableNames', poplist, 'RowNames',...
            cellstr([repmat('IGHV', 7, 1) num2str((1:7)')]));
    end
end
