function r = mut(obj, pop, region, t, varargin)
    % IGREP/MUT Mutation frequencies for specified population
    %
    % m = MUT(obj, pop, r, t) returns mutation frequencies for the
    % specified population.  The first argument is the desired population.
    % The second argument is one of "v", "fr", "cdr" or "all", indicating
    % the region(s) to report mutation frequencies - the entire v-region,
    % the FR's (FR1 + FR2 + FR3), the CDR's (CDR1 + CDR2), or all three.
    % The third argument is one of "total", "nonsilent", "silent" or "all",
    % indicating the type of mutation to report.  The return value is a
    % Matlab table.  Rows are sequences and columns are lineageID, isotype,
    % vgene, jgene and then the mutation statistics - the number of columns
    % depends on second and third input arguments.
    %
    % m = MUT(obj, pop, r, t, filt) returns mutation frequencies for the
    % specified population under the specified filter specs. The optional
    % argument is a three-element cell array of strings: the first element
    % is the isotype, the second is the vgene and the third is the jgene -
    % values to limit the sequences by.  For the vgene and jgene,  if
    % immediately preceded by a "!", the filter will return all sequences
    % but those with the specified value.  For the isotype, if you want to
    % match an isotype prefix, add a wildcard ("%") to the end of an
    % isotype - for example, "G" will match G1, G1A, G1B if you want to
    % collapse subtypes into a supertype.   To avoid filtering by one of
    % the parameters, use an empty character array.
    %
    % Dot notation usage:
    %     m = D.MUT('pop', 'all', 'all')
    %     m = D.MUT('pop', 'all', 'total', {'A', '', ''})
    %     m = D.MUT('pop', 'cdr', 'all', {'A', 'IGHV4-34', ''})
    %     m = D.MUT('pop', 'cdr', 'all', {'G%', '', 'IGHJ3'})
    %     m = D.MUT('pop', 'fr', 'silent', {'', '!IGHV1-2', ''})
    %     m = D.MUT('pop', 'all', 'all', {'M', 'IGHV3-7', 'IGHJ4'})
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validatetextarg(obj, region);
    v1 = {'all', 'v', 'fr', 'cdr'};
    if ~ismember(region, v1)
        error(['arg 2 must be one of: ' strjoin(v1, ' | ')]);
    end
    validatetextarg(obj, t);
    v2 = {'all', 'total', 'nonsilent', 'silent'};
    if ~ismember(t, v2)
        error(['arg 3 must be one of: ' strjoin(v2, ' | ')]);
    end
    filt = '';
    if nargin > 4, filt = makefilt(obj, varargin{1}); end             
    validatepop(obj, pop);
    key = table;
    key.region = {...
        'v';   'v';   'v';...
        'fr';  'fr';  'fr';...
        'cdr'; 'cdr'; 'cdr'};
    key.type = repmat({'total', 'nonsilent', 'silent'}', 3, 1);
    key.numerator = {...
        'V_Mutations';...
        'V_NonsilentMutations';...
        'V_SilentMutations';...
        'FR1_Mutations + FR2_Mutations + FR3_Mutations';...
        'FR1_NonsilentMutations + FR2_NonsilentMutations + FR3_NonsilentMutations';...
        'FR1_SilentMutations + FR2_SilentMutations + FR3_SilentMutations';...
        'CDR1_Mutations + CDR2_Mutations';...
        'CDR1_NonsilentMutations + CDR2_NonsilentMutations';...
        'CDR1_SilentMutations + CDR2_SilentMutations'};
    key.header = {...
        'V_total';   'V_nonsilent';   'V_silent';...
        'FR_total';  'FR_nonsilent';  'FR_silent';...
        'CDR_total'; 'CDR_nonsilent'; 'CDR_silent'};
    rg = cellstr(region); if isequal(region, 'all'), rg = v1(2:end)'; end            
    ty = cellstr(t); if isequal(t, 'all'), ty = v2(2:end)'; end          
    myregions = reshape(repmat(rg', length(ty), 1),...
        length(rg) * length(ty), 1);
    mytypes = repmat(ty, length(rg), 1);
    s = cell(length(myregions), 1);
    h = cell(length(myregions), 1);
    for i = 1:length(myregions)
        k = intersect(...
            find(strcmp(myregions{i}, key.region)),...
            find(strcmp(mytypes{i}, key.type)));
        s(i) = key.numerator(k);
        h(i) = key.header(k);
    end                                   
    if ~isempty(find(strncmp(s, 'V',   1), 1))
        s = [s; 'V_Nucleotides'];
        h = [h; 'V_nucleotides'];
    end
    if ~isempty(find(strncmp(s, 'FR',  1), 1))
        s = [s; 'FR1_Nucleotides + FR2_Nucleotides + FR3_Nucleotides'];
        h = [h; 'FR_nucleotides'];                    
    end
    if ~isempty(find(strncmp(s, 'CDR', 1), 1))
        s = [s; 'CDR1_Nucleotides + CDR2_Nucleotides'];
        h = [h; 'CDR_nucleotides'];
    end                            
    s0 = {'lineageID AS lineage', 'isotype', 'vgene AS v', 'jgene AS j'}';
    a = query(obj, [...
            'SELECT ' strjoin([s0; s], ', ') ' FROM Sequence ',...
            'WHERE population = "' char(pop) '" ' filt]);
    x = a{:, 5:end};
    inum = find(~cellfun('isempty', regexp(s, 'Mutations')));
    iden = find(~cellfun('isempty', regexp(s, 'Nucleotides')));            
    iv   = find(strncmp(s, 'V', 1));
    ifr  = find(strncmp(s, 'F', 1));
    icdr = find(strncmp(s, 'C', 1));            
    inumv   = intersect(inum, iv);
    inumfr  = intersect(inum, ifr);
    inumcdr = intersect(inum, icdr);            
    idenv   = intersect(iden, iv);
    idenfr  = intersect(iden, ifr);
    idencdr = intersect(iden, icdr);            
    frq = [];
    header = {};            
    if ~isempty(iv)
        for i = 1:length(inumv)
            frq = [frq, 100 * x(:, inumv(i)) ./ x(:, idenv)];
            header = [header, h(inumv(i))];
        end
        frq = [frq, x(:, idenv)];
        header = [header, h(idenv)];                
    end            
    if ~isempty(ifr)
        for i = 1:length(inumfr)
            frq = [frq, 100 * x(:, inumfr(i)) ./ x(:, idenfr)];
            header = [header, h(inumfr(i))];
        end
        frq = [frq, x(:, idenfr)];
        header = [header, h(idenfr)];                
    end                
    if ~isempty(icdr)
        for i = 1:length(inumcdr)
            frq = [frq, 100 * x(:, inumcdr(i)) ./ x(:, idencdr)];
            header = [header, h(inumcdr(i))];
        end
        frq = [frq, x(:, idencdr)];
        header = [header, h(idencdr)];                
    end    
    r = [a(:, 1:4) array2table(frq, 'VariableNames', header)];
end
