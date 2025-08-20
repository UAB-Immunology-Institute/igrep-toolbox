function [c, b] = ntcolor(varargin)

% NTCOLOR Colormap for nucleotides
%
%    NTCOLOR returns two values - a 4x3 or 5x3 colormap and a corresponding
%    row labels.  The colormap is the same as that used by Jalview.  The
%    optional input argument is either 4 (default) or 5 and specifies
%    whether to include the placeholder color (gray) as the first row.
%
%    A = green
%    C = orange
%    G = red
%    T = blue
%    placeholder = gray
%
%    Usage:
%
%        [cmap, bases] = ntcolor
%
%        [cmap, bases] = ntcolor(5)

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % input argument error check
    n = 4;
    if nargin > 0
        n = varargin{1};
        if ~isnumeric(n) || ~isscalar(n) || ~ismember(n, [4 5])
            error('argument must be 4 or 5');
        end
    end
            
    % jalview colormap
    c = [...
        200 200 200
        100 247  63;...
        255 179  64;...
        235  65  60;...
        60  136 238] / 255;
    b = {'.', 'a', 'c', 'g', 't'}';
    
    % exclude gray
    if n == 4
        c = c(2:end, :);
        b = b(2:end);
    end

return


