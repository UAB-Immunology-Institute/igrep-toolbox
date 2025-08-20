function r = vdj5trim(obj, start, filterval)
    % IGREP/VDJ5TRIM Number of unique VDJ sequences
    %
    % r = VDJ5TRIM(obj, start, filterval) returns a table of vdj nucleotide
    % sequences for either a specified lineage or population.  The vdj
    % sequence is broken up into two columns based on the specified 5' 
    % start position.  The output table also returns the lineageID,
    % population, v gene, j gene and isotype.  This is useful because for
    % some datasets the 5' end of the v gene contains random hexamers.
    %
    % t = VDJ5TRIM(obj, 5, 'mypop') returns a table of trimmed vdj
    % sequences for the population 'mypop' starting at base 5.
    %
    % t = VDJ5TRIM(obj, 7, 13023) returns a table of trimmed vdj sequences
    % for the lineage 13023 starting at base 7.
    %
    % Dot notation usage:
    %    r = D.VDJ5TRIM(5, 'mypop')
    %    r = D.VDJ5TRIM(7, 13023)
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.

    if isnumeric(filterval)
        if isscalar(filterval)
            w = ['lineageID = ' num2str(filterval)];
        else
            error('only on lineageID can be specified');
        end
    else
        validatetext(filterval, 1);
        w = ['population = "' char(filterval) '"'];
    end
    a = query(obj, [...
        'SELECT population, lineageID, isotype, vgene, jgene, nt_VDJ ',...
        'FROM Sequence WHERE ' w]);
    seq = char(a.NT_VDJ);
    r = a(:, 1:5);
    r.Properties.VariableNames = {'population', 'lineage', 'isotype', 'v', 'j'};
    r.vdjDiscard = string(seq(:, 1:(start - 1)));
    r.vdjTrimmed = strip(string(seq(:, start:end)));

end
