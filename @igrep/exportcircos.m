function exportcircos(obj, config)
    % IGREP/EXPORTCIRCOS Export files for generating circos plot
    %
    % EXPORTCIRCOS(obj, f) generates configuration and data files suitable
    % for use with Circos.  Input argument is a configuration file. Creates
    % an output folder called "dbCircosOutput" in the current wirking
    % directory.
    %
    % Dot notation usage:
    %     D.EXPORTCIRCOS('/path/to/my/igseq/config/circos.inf')
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    perl('dbCircos.pl', obj.datasource, config)            
end  
