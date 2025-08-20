function validatecolor(arg, varargin)

% VALIDATECOLOR Test a variable to see if it is a valid color designation
%
%    VALIDATECOLOR(c) throws an error if the input argument is not a valid
%    color designation.  Valid options are a three-element vector (where
%    each value is between 0 and 1), or a single character (r, g, b, c, y,
%    m, k, w).  Input can be a cell vector of colors.  Examples of valid
%    colors are:
%
%        'r'
%
%        [.2 .5 0]
%
%        {'r', 'b', [0 .7 0]}
%
%    VALIDATECOLOR(c, n) throws an error if there are not n items.  If n =
%    1 then the input can be a cell array, a character or a three element
%    vector.

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    na = NaN;
    if nargin > 1
        na = varargin{1};
        if ~(isnumeric(na) && isscalar(na) && na > 0)
            error('second optional argument must be numeric scalar > 0');
        end
    end 

    if iscell(arg)        
        if min(size(arg)) > 1
            error('cell array must be a vector');
        else
            n = length(arg);
            if isfinite(na) && na ~= n
                error(['expecting ' num2str(na) ' colors']);
            end
            for i = 1:n
                test(arg{i})
            end
        end
    else
        if isfinite(na) && na > 1
            error(['expecting a cell vector of ' num2str(na) ' colors']);
        end
        test(arg);
    end
    
return
        
    function test(c)
        
        if isnumeric(c)        
            if ~((numel(c) == 3) && (min(c) >= 0) && (max(c) <= 1))
                error(['numeric color designation should be 3-element ',...
                    'vector with values between 0 and 1']);
            end        
        elseif ischar(c)        
            v = {'r', 'g', 'b', 'c', 'y', 'm', 'k', 'w'};
            if ~ismember(c, v)
                error(['single character color should be one of: ',...
                    strjoin(v, ', ')]);
            end        
        else        
            error('color should be a 3-element vector or character');        
        end
        
    return
        
       