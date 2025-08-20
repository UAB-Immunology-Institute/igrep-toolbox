function [x, y] = hbez(x1, y1a, y1b, x2, y2a, y2b, spread)

% HBEZ Coords for filled bezier useful for horizontal alluvial plots
%
%    [xb, yb] = HBEZ(x1, y1bot, y1top, x2, y2bot, y2top, spread) returns x
%    and y coordinates for filled bezier ribbon between two x-axis 
%    positions.  The first three arguments define the vertical extent of
%    the ribbon for the first x-axis position, and the next three arguments
%    define the vertical extent for the second x-axis position.  The spread 
%    controls how curved the ribbon appears: a value of 0 creates a
%    straight ribbon and a value of 1 creates a sharp elbow.  This function 
%    uses a cubic bezier (two end points, two control points).
%
%    This code is based in part from Mathworks blog entry (posted by Mike
%    Garrity on October 13, 2014) on beziers which can be found at:
%    https://blogs.mathworks.com/graphics/2014/10/13/bezier-curves/

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % coords for bottom edge of bezier ribbon
    p1 = bez([x1; y1a], [x2; y2a], spread); 
    
    % coords for top edge of bezier ribbon
    p2 = bez([x1; y1b], [x2; y2b], spread);

    % coords for entire ribbon suitable for use with "fill" command
    x = [p1(1, :) fliplr(p2(1, :)) p1(1, 1)];
    y = [p1(2, :) fliplr(p2(2, :)) p1(2, 1)];

return

    % helper function for bezier curve: a and d are endpoints
    function pts = bez(a, d, s)
       
        % b and c are control points
        b = [a(1) + (s * (d(1) - a(1))); a(2)];
        c = [d(1) - (s * (d(1) - a(1))); d(2)];
        
        % points along bezier (100 points; change change if necessary)
        t = linspace(0,1,101);
        
        % compute cubic bezier (see reference above)
        pts = ...
            kron((1 - t) .^ 3, a) + ...
            kron(3 * (1 - t) .^ 2 .* t, b) + ...
            kron(3 * (1 - t) .* t .^2, c) + ...
            kron(t .^ 3, d);

    return

