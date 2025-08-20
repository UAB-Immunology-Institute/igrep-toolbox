function result = cdr3search(AU, t, cdr3, seqsim, varargin)

% CDR3SEARCH find sequences in database that match a cdr3 sequence
%
%    CDR3SEARCH(datasource, type, sequence, similarity) queries a igseq
%    database to find sequences that have a specified cdr3 nucleotide or
%    amino acid sequence within a specified sequence similarity tolerance.
%    The first argument is an IGREP object, the second argument is the type
%    of sequence (either "nt" for nucleotide or "aa" for amino acid), the
%    third argument is the CDR3 sequence to search for, and the fourth
%    argument is the percent sequence similarity threshold for the search.
%    The output is a table of hits, where the fields are sequenceID,
%    population, lineageID, isotype, v-gene, j-gene, percent similarity,
%    cdr3 sequence (nt or aa), and vdj sequence (nt or aa).
%
%    CDR3SEARCH(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'v'            restrict search based on particular v-gene
%        'j'            restrict search based on particular j-gene
%        'population'   restrict search to specified population
%
%    Usage:
%
%        CDR3SEARCH(D, 'aa', 'AKDQANYYDSSGYYKRKDYFDY', 85)
%
%        CDR3SEARCH(D, 'aa', 'AKDQANYYDSSGYYKRKDYFDY', 85,...
%            'population', 'PB')
%
%        CDR3SEARCH(D, 'aa', 'AKDQANYYDSSGYYKRKDYFDY', 85,...
%            'v', 'IGHV3-23', 'j', 'IGHJ4', 'population', 'PB')
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % argument validation and attribute-value pair option handling
    pa = inputParser;
    addRequired(pa, 'AU',...
        @(x) validateattributes(x, {'igrep'}, {'nonempty'}));    
    addRequired(pa, 't',...
        @(x) ~isempty(validatestring(x, {'nt', 'aa'})));
    addRequired(pa, 'cdr3',...
        @(x) validateseq(t, x, 1));
    addRequired(pa, 'seqsim',...
        @(x) validateattributes(x, {'numeric'},...
        {'scalar','>', 0, '<=', 100}));
    addParameter(pa, 'v',   '', @(x) validatetext(x, 1));
    addParameter(pa, 'j',   '', @(x) validatetext(x, 1));    
    addParameter(pa, 'population', '', @(x) validatetext(x, 1));
    parse(pa, AU, t, cdr3, seqsim, varargin{:});    
    v = char(pa.Results.v);
    j = char(pa.Results.j);
    p = pa.Results.population;
    
    % test for valid v/j genes, populations
    msg = ' not in database';
    if ~isempty(v) && ~ismember(v, AU.uv),   error([v ' vgene' msg]); end
    if ~isempty(j) && ~ismember(j, AU.uj),   error([j ' jgene' msg]); end
    if ~isempty(p) && ~ismember(p, AU.upop), error([p msg]);          end
    
    % where clauses
    w = '';
    if ~isempty(v), w = [w ' AND vgene = "' v '"'];      end
    if ~isempty(j), w = [w ' AND jgene = "' j '"'];      end
    if ~isempty(p), w = [w ' AND population = "' p '"']; end
    
    % field to length filter on
    f = [t '_cdr3'];
    f1 = [t '_vdj'];
    
    % run query
    r = AU.sql([...
        'SELECT sequenceID, population, lineageID, isotype, ',...
        'vgene, jgene, ' f ', ' f1 ' FROM Sequence ',...
        'WHERE length(' f ') = ' num2str(length(cdr3)) w]);
    d = eq(char(r{:, 7}), repmat(cdr3, size(r, 1), 1));
    ss = 100 * sum(d, 2) / size(d, 2);
    r = [r num2cell(ss)];
    
    % convert output to table
    result = r(ss >= seqsim, :);
    result.Properties.VariableNames([5 6 9]) = {'v', 'j', 'similarity'};
    result = result(:, [1 2 3 4 5 6 9 7 8]);
    
return
    
    
    
    

