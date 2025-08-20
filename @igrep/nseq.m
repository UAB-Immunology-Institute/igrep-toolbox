function r = nseq(obj, varargin)
    % IGREP/NSEQ Number of sequences for populations
    %
    % r = NSEQ(obj) returns an nx2 Matlab table where the first column
    % contains the populations represented in the database and the second
    % column contains the corresponding numbers of sequences.
    %
    % r = NSEQ(obj, pop) returns the number of sequences for the
    % population specified in the second argument.
    %
    % Dot notation usage:
    %     r = D.NSEQ
    %     r = D.NSEQ('mypop')   
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    if nargin == 1
        r = obj.sequences(:, [1 2]);
    else
        validatepop(obj, varargin{1});
        r = obj.sequences.total(obj.sequences.population == varargin{1});
    end   
end
