function validatetext(arg, varargin)

% VALIDATETEXT Throw error if not a valid text item or list of items
%
%    VALIDATETEXT(arg) throws an error if the input argument is not a valid
%    text item: either a char array (text enclosed by single quotes), a
%    string (text enclosed by double quotes), a cell array of char arrays
%    or an array of strings.
%
%    VALIDATETEXT(arg, n) throws an error if there are not n items.  If n =
%    1, the argument can be a char array.

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    n = NaN;
    if nargin > 1
        n = varargin{1};
        if ~(isnumeric(n) && isscalar(n) && n > 0)
            error('second optional argument must be numeric scalar > 0');
        end
    end 
            
    a = ischar(arg);
    b = iscellstr(arg);
    c = isstring(arg);
    m = numel(arg);
    s = min(size(arg));
    if ~(a || ((b || c) && (s == 1)))
        error(['input must be (char array) | (string vector) | ',...
            '(cell vector of char arrays)']);
    else
        if n == 1
            if ~(a || ((b || c) && m == 1))
                error('expecting one text item');
            end
        elseif n > 1
            if ~((b || c) && m == n)
                error(['expecting ' num2str(n) ' text items']);
            end
        end
    end
    

return