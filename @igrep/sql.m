function r = sql(obj, q)
    % IGREP/SQL Run a custom query
    %
    % r = SQL(obj, query) returns the result of the user-specified query.
    % The input should be a char (bound by single quotes). This allows the
    % use of double quotes within the query string. The return value is a
    % cell array.
    %
    % Dot notation usage:
    %     a = 'SELECT count(*) FROM Sequence WHERE isotype = "M"'
    %     r = D.SQL(a)
            
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
            
    if ~ischar(q), error('input must be a char array'); end
    r = query(obj, q);
    
end
