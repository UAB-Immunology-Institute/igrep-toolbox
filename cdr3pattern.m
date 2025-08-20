function result = cdr3pattern(AU, pattern)

% CDR3PATTERN find sequences in database that match a cdr3 sequence
%
%    CDR3PATTERN(datasource, pattern) searches for a specific amino acid
%    pattern in the CDR3 heavy chain region.  Acceptable wildcards are "%"
%    for "one or more residues" or "_" for "one residue".  Return value is 
%    a table of hits indicating lineageID, population, isotype, v-gene,
%    j-gene and the CDR3 amino acid sequence.
%
%    Usage:
%
%        CDR3PATTERN(D, '%GATA%F_F%')
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
    addRequired(pa, 'pattern',...
        @(x) validateseq('aa', x, 1));
    parse(pa, AU, pattern);    
        
    % run query
    result = AU.sql([...
        'SELECT lineageID, population, isotype, vgene, jgene, AA_CDR3 ',...
        'FROM Sequence WHERE AA_CDR3 LIKE "' pattern '"']);
      
return
    
    
    
    

