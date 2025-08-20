classdef igrep < handle & matlab.mixin.CustomDisplay
    % IGREP is a handle class that encapsulates and provides an interface
    % for an open igseq sqlite database connection.  Public properties are
    % read-only.  Database connectivity is established upon instantiation
    % and is maintained as long as the object exists.
    %
    % IGREP Properties (read-only):
    %
    %   datasource       - name of database (with full path)
    %   dbname           - name of database (without full path)
    %   upop             - cell array of strings of unique population names
    %   uiso             - cell array of strings of unique isotypes
    %   uv               - cell array of strings of unique v-genes
    %   uj               - cell array of strings of unique j-genes
    %   nrow             - number of rows in database
    %   sequences        - number of sequences for each population, isotype
    %   lineages         - number of lineages for each population
    % 
    % IGREP Methods:
    %   
    %   Utilities
    %
    %     igrep          - constructor
    %     delete         - destructor
    %     sql            - run a custom query
    %
    %   Numbers of Things
    %
    %     nseq           - get # of sequences for populations in database
    %     nsequ          - get # of unique VDJ sequences for a population
    %     nlin           - get # of lineages for populations in database
    %     nvj            - return V- and J-gene usage matrix for population
    %     nvx            - return V-family usage for selected populations
    %
    %   Lineage Utilities
    %
    %     ranklin        - get lineages and order by their sizes
    %     rankliniso     - get lineages, sizes and isotype compositions
    %     ranklinfull    - get lineages, v-, j-genes, counts, mutations
    %     lintable       - create a lineage table for specified populations
    %     shlin          - percent shared lineages between two populations
    %     hist           - create a lineage size occurrence histogram
    %     linreport      - return a brief summary for a particular lineage
    %     exportinfo     - export info for selected lineage to file
    %     exportfasta    - export fasta file for selected lineage
    %     exportcircos   - export files used by Circos
    %
    %   Other Characteristics
    %
    %     dx             - return D(x) value(s) for a population
    %     pielou         - return pielou evenness value for a population
    %     cdr3len        - return CDR3 lengths (nt or aa) of population
    %     cdr3lenhist    - return CDR3 length histogram of population
    %     mut            - return mutation freqs of populations, isotypes
    %     muthist        - return histogram of mutation frequencies
    %     mhoverlap      - return all-by-all Morisita-Horn overlap indices
    %     vdj5trim       - return vdj of lineage/pop with 5' end trimmed
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
        
    properties(GetAccess = 'public', SetAccess = 'private')
        
        datasource;
        dbname;
        upop;
        uiso;
        uv;
        uj;
        nrow;
        sequences;
        lineages;
        
    end
    
    properties(GetAccess = 'private', SetAccess = 'private')
        
        handle;
        
    end
    
    methods
        
        function obj = igrep(d)
            % IGREP/IGREP Constructor for IGREP object 
            %
            % obj = IGREP(datasource) creates a IGREP object which is a
            % container for an igseq analysis unit stored in an sqlite
            % database.  The argument (a char array or a string scalar) is
            % an igseq sqlite database (including full path).
            %
            % Usage:
            %
            %    D = igrep('/path/to/my.db')
            validatetextarg(obj, d);
            obj.datasource = char(d);
            try
                h = sqlite(d);   % use native function instead of jdbc 7/29/25
            catch
                error('could not establish database connection');
            end
            obj.handle = h;
            disp(['     open connection to: ' d]);
            a = strsplit(obj.datasource, '/');
            obj.dbname = a{end};
            qp = 'SELECT distinct(population) FROM Sequence';
            qi = 'SELECT distinct(isotype) FROM Sequence';
            qv = 'SELECT distinct(vgene) FROM Sequence ORDER BY vgene';
            qj = 'SELECT distinct(jgene) FROM Sequence ORDER BY jgene';
            qs = 'SELECT count(*) AS n FROM Sequence';
            obj.upop = query(obj, qp).population;
            obj.uiso = query(obj, qi).isotype; 
            obj.uv   = query(obj, qv).Vgene;
            obj.uj   = query(obj, qj).JGene;
            obj.nrow = query(obj, qs).n;
            b = query(obj, [...
                'SELECT population, isotype, count(*) AS n ',...
                'FROM Sequence GROUP BY population, isotype']);
            up = unique(b.population);
            ui = unique(b.isotype);            
            r = zeros(length(up), length(ui));
            for k = 1:size(b, 1)
                irow = strcmp(b.population(k), up);
                icol = strcmp(b.isotype(k),    ui);
                r(irow, icol) = b.n(k);
            end
            t = table;
            t.population = string(up);
            t.total = sum(r, 2);
            obj.sequences = [t array2table(r, 'VariableNames', ui)];
            obj.lineages = query(obj, [...                    
                'SELECT population, count(distinct(lineageID)) AS n ',...
                'FROM Sequence GROUP BY population ',...
                'ORDER BY population']);
        end
            
        function delete(obj)
            % IGREP/DELETE Destructor
            %
            % DELETE(obj) closes the database connection and destroys the
            % object.  This is inherited from the abstract Handle class.
            % This method runs when called explicitly or, more typically,
            % when the variable is cleared or overwritten.
            disp(['    close connection to: ' obj.datasource]);
            close(obj.handle);
        end
                
        % methods in individual files in @igrep directory
        r = sql(obj, q)        
        r = nseq(obj, varargin)        
        r = nsequ(obj, pop)
        r = nlin(obj, varargin)
        r = nvj(obj, pop, method, varargin)
        r = nvx(obj, poplist, method, fam, varargin)
        r = ranklin(obj, pop, varargin)
        r = rankliniso(obj, pop)
        r = ranklinfull(obj, pop) 
        r = hist(obj, pop, varargin)
        r = linreport(obj, lin)
        r = dx(obj, x, pop, varargin)
        r = pielou(obj, pop, varargin)
        r = mut(obj, pop, region, t, varargin)
        r = cdr3len(obj, pop, flag, varargin)        
        r = mhoverlap(obj, varargin)
        r = vdj5trim(obj, offset, filterval)
        [m, lin, pop, varargout] = lintable(obj, varargin)
        [p, nc] = shlin(obj, pop1, pop2, varargin)
        [muth, stat] = muthist(obj, pop, r, t, bin, maxf, varargin)
        [cdrh, stat] = cdr3lenhist(obj, pop, seq, varargin)
        exportinfo(obj, lin)
        exportfasta(obj, lin)
        exportcircos(obj, config)
                                           
    end
    
    methods(Access = 'private')
        
        function a = makefilt(obj, f)
            % construct a "where" clause for filtering on isotype, v, j
            if length(f) ~= 3
                error('optional argument must have three values');
            end
            if ~isempty(f{1})
                
                if ~isequal(f{1}(end), '%')
                    validateiso(obj, f{1});                
                elseif ~strncmp(f{1}(1:(end - 1)), obj.uiso,...
                        length(f{1}) - 1)                
                    error(['isotype ' char(f{1}) ': no database match']); 
                end
                isoclause = [' AND isotype LIKE "' f{1} '"'];
            else
                isoclause = '';
            end                
            if ~isempty(f{2})
                validatetextarg(obj, f{2});
                temp = f{2};
                if strncmp(temp, '!', 1)
                    vclause = [' AND vgene != "' temp(2:end) '"']; 
                else
                    vclause = [' AND vgene = "' temp '"'];
                end
            else
                vclause = '';
            end                
            if ~isempty(f{3})
                validatetextarg(obj, f{3});
                temp = f{3};
                if strncmp(temp, '!', 1)
                    jclause = [' AND jgene != "' temp(2:end) '"']; 
                else
                    jclause = [' AND jgene = "' temp '"'];
                end
            else
                jclause = '';
            end
            a = [isoclause, vclause, jclause];
        end
        
        function a = query(obj, q)
            % IGREP/QUERY Run an SQL query 
                a = fetch(obj.handle, q);
                % convert INT64 columns to DOUBLES for native interface
                idx = varfun(@(x) isa(x, 'int64'), a,...
                    'OutputFormat', 'uniform');
                a = convertvars(a, idx, @double);
                % convert STRING columns to CELLARRAY for native interface
                % jdx = varfun(@(x) isa(x, 'string'), a,...
                %     'OutputFormat', 'uniform');
                % a = convertvars(a, jdx, @cellstr);
        end
        
        function validatetextarg(~, x)
            % test if string scalar or character array
            if ~((isstring(x) && isscalar(x)) || ischar(x))
                error('input must be char array or string scalar');
            end
        end
        
        function validatepop(obj, p)
            % test if population in database - can be a list of pops
            validatetextarg(obj, p);
            if ~ismember(p, obj.upop)
                error(['population ' char(p) ' not in database']);
            end
        end
        
        function validatepoplist(obj, ps)
            % test whether all items in poplist are in database
            if ~(iscellstr(ps) && min(size(ps)) == 1)
                error('poplist is not cell vector of character arrays')
            end
            if ~isempty(setdiff(ps, obj.upop))
                error('one or more populations not in database');
            end
        end
        
        function validateiso(obj, i)
            % test if isotype in database
            validatetextarg(obj, i);
            if ~ismember(i, obj.uiso)
                error(['isotype ' char(i) ' not in database']);
            end
        end
        
    end
    
end
    
