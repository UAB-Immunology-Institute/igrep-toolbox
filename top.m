function r = top(t, varargin)

% TOP show top n rows of table
%
%    r = TOP(t) returns the first 8 rows of an array, cell array or table.
%
%    r = TOP(t, n) returns the first n rows of an array, cell array or
%    table.
%
%    Usage:
%
%        TOP(t)
%        TOP(t, 20)

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % rows of input
    n = size(t, 1);
    
    % if specifying number of rows
    if nargin == 1
        nrow = 8;
    else
        if ~isscalar(varargin{1}) || ~isnumeric(varargin{1})
            error('first argument must be scalar numeric');
        end
        nrow = varargin{1};
    end
    
    % if asked for more rows than there are
    if nrow > n, nrow = n; end
    
    % return value
    r = t(1:nrow, :);
    
return
    
    
