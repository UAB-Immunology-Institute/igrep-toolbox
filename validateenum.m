function validateenum(arg, list)

% VALIDATEENUM Throw error if not a valid text item in a pick list
%
%    VALIDATEENUM(arg, list) throws an error if the input item (first
%    argument - a string or char array) is not in the list (second
%    argument - a string array or cell array of chars).

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    if ~(isstring(list) || iscellstr(list))
        error('list must be string array for cell vector of chars');
    end
    if numel(list) < 2
        error('list must contain at least two items');
    end
    if min(size(list)) > 1
        error('list must be vector');
    end
    validatetext(arg, 1);
    list2 = cellstr(list);
    err = strjoin(list2, ' | ');
    if ~ismember(char(arg), list2)
        error(['must be one of: ' err]);
    end

return

