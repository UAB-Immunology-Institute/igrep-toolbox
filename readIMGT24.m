function T = readIMGT24(filename)

% READIMGT24 read in IMGT High/V-Quest 2 or 4 alignment file
%
%    READIMGT24(filename) reads in an IMGT alignment file (a "2_" file for
%    nucleotides, or a "4_" file for amino acids).  The output is a Matlab
%    table of the file contents.  Two additional columns are added: vgene
%    and jgene which are derived from the vGeneAndAllele and jGeneAndAllele
%    respectively.  These use the first annotation hit if there are more
%    than one and strips off the allele information.
%
%    Usage:
%
%        T = READIMGT24('4_IMGT-gapped-AA-sequences.txt')

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % read in IMGT 2 or 4 file
    fmt = ['%f' repmat('%s', 1, 17)];
    fid = fopen(filename);
    h = strsplit(fgetl(fid), char(9));
    C = textscan(fid, fmt, 'delimiter', char(9));
    fclose(fid);

    % remap header names to be suitable for table field names
    r = {...
        'Sequence number',        'seqNum';...
        'Sequence ID',            'seqID';...
        'V-DOMAIN Functionality', 'functionality';...
        'V-GENE and allele',      'vGeneAndAllele';...
        'J-GENE and allele',      'jGeneAndAllele';...
        'D-GENE and allele',      'dGeneAndAllele';...
        'V-D-J-REGION',           'vdj';...
        'V-J-REGION',             'vj';...
        'V-REGION',               'v';...
        'FR1-IMGT',               'fr1';...
        'CDR1-IMGT',              'cdr1';...
        'FR2-IMGT',               'fr2';...
        'CDR2-IMGT',              'cdr2';...
        'FR3-IMGT',               'fr3';...
        'CDR3-IMGT',              'cdr3';...
        'JUNCTION',               'jct';...
        'J-REGION',               'j';...
        'FR4-IMGT',               'fr4'};

    % check file fields
    if ~isequal(sort(h(:)), sort(r(:, 1)))
        error('does not appear to be an IMGT 2- or 4-file: check fields');
    end
    
    % remap headers
    hnew = cell(length(h), 1);
    for i = 1:length(h), hnew(i) = r(strcmp(r(:, 1), h{i}), 2); end

    % create table of IMGT data
    T = table; for i = 1:length(h), T.(hnew{i}) = C{i}; end

    % clean up v- and j-genes, make new columns
    T.vGene = regexprep(T.vGeneAndAllele, '(^Homsap |\*.*$)', '');
    T.jGene = regexprep(T.jGeneAndAllele, '(^Homsap |\*.*$)', '');
        
    % reorder table
    T = T(:, [1:4 19 5 20 6:18]);

return


    
    
    