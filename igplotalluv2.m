function igplotalluv2(AU, pop, varargin)

% IGPLOTALLUV2 Create a 2-population alluvial clonality plot
%
%    IGPLOTALLUV2(datasource, pop) creates two clonality plots where common
%    lineages in the two populations are connected by bezier ribbons.  The
%    first argument is an IGREP object.  The second argument is a 2-element
%    cell vector of population names as appears in the specified data
%    source.
%
%    IGPLOTALLUV2(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        alternate names to display on plot (instead of
%                       population names) - a two element cell vector of
%                       char arrays
%        'method'       specifies how to sort clonality bars:
%                           'separate' (default) - sort each clonality bar
%                               separately
%                           'bar1all' - sort by connected lineages first,
%                               sorted by bar 1, then non-connected
%                               lineages separately
%                           'bar2all' - sort by connected lineages first,
%                               sorted by bar 2, then non-connected
%                               lineages separately
%                           'bar1common' - only show connected lineages for
%                               both populations and sort by bar 1
%                           'bar2common' - only show connected lineages for
%                               both populations and sort by bar 2
%        'noconn'       specify max number of sequences to not count as a
%                       connection.  For example, "1" means that the
%                       clonotype in both populations needs to be more than
%                       1 sequence for a ribbon to be drawn.  Default is 0.
%
%    Usage:
%
%        IGPLOTALLUV2(D, {'p7', 'p11'})
%
%        IGPLOTALLUV2(D, {'p7', 'p11'}, 'label', {'mem', 'PC'},...
%            'method', 'bar2all')
%
%    Requires Matlab Database Toolbox

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    % allowed sort methods
    mlist = {'separate', 'bar1all', 'bar2all', 'bar1common', 'bar2common'};

    % argument validation and attribute-value pair option handling
    pa = inputParser;
    addRequired(pa, 'AU',...
        @(x) validateattributes(x, {'igrep'}, {'nonempty'}));
    addRequired(pa, 'pop',...
        @(x) validatetext(x, 2));
    addParameter(pa, 'label', pop,...
        @(x) validatetext(x, 2));
    addParameter(pa, 'method', 'separate',...
        @(x) validateenum(x, mlist));
    addParameter(pa, 'noconn', 0,...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0}));
    parse(pa, AU, pop, varargin{:});    
    pop = cellstr(pop);
    label = cellstr(pa.Results.label);
    method = char(pa.Results.method);
    t = pa.Results.noconn;

    % get data
    Z = table2array(AU.lintable({pop{1} pop{2}}));
    
    % to prevent artifact by deafult sort for singletons
    Z = Z(randperm(size(Z, 1)), :);
        
    % case for connected connected lineages only
    if ismember(method, {'bar1common', 'bar2common'})
        ylab = 'Cumulative Percentage of (Connected) Sequences';
        Z = Z(prod(Z(:, 2:3), 2) > 0, :);
    else
        ylab = 'Cumulative Percentage of Sequences';
    end

    % lineage IDs of common lineages
    linconn = Z(min(Z(:, [2 3]), [], 2) > t, 1);
    
    % express as percent total sequences
    Z(:, 2) = 100 * Z(:, 2) / sum(Z(:, 2));
    Z(:, 3) = 100 * Z(:, 3) / sum(Z(:, 3));
    
    % set up data for clonality bars (only non-zero in each population)
    C1 = Z(Z(:, 2) > 0, [1 2]);
    C2 = Z(Z(:, 3) > 0, [1 3]);
    
    % add third column indicating whether connected or not
    C1(:, 3) = zeros(size(C1, 1), 1);
    C2(:, 3) = zeros(size(C2, 1), 1);
    [~, ii1] = intersect(C1(:, 1), linconn);
    [~, ii2] = intersect(C2(:, 1), linconn);
    C1(ii1, 3) = 1;
    C2(ii2, 3) = 1; 
    
    % sort each of the clonality datasets according to fourth argument
    if isequal(method, 'separate')
        [~, isort1] = sort(C1(:, 2));
        [~, isort2] = sort(C2(:, 2));
        S1 = C1(isort1, :);
        S2 = C2(isort2, :);
    elseif ismember(method, {'bar1common', 'bar2common'})
        if isequal(method, 'bar1common')
            [~, isort] = sort(C1(:, 2));
        else
            [~, isort] = sort(C2(:, 2));
        end
        S1 = C1(isort, :);
        S2 = C2(isort, :);
    elseif ismember(method, {'bar1all', 'bar2all'})
        MP1 = C1(C1(:, 3) == 1, :);
        MU1 = C1(C1(:, 3) == 0, :);
        MP2 = C2(C2(:, 3) == 1, :);
        MU2 = C2(C2(:, 3) == 0, :);
        if isequal(method, 'bar1all')
            [~, isort1] = sort(MP1(:, 2));            
            isort2 = zeros(size(MP1, 1), 1);
            for i = 1:length(isort2)
                isort2(i) = find(MP1(isort1(i), 1) == MP2(:, 1));
            end
        else
            [~, isort2] = sort(MP2(:, 2));
            isort1 = zeros(size(MP2, 1), 1);
            for i = 1:length(isort1)
                isort1(i) = find(MP2(isort2(i), 1) == MP2(:, 1));
            end
        end
        [~, isortu1] = sort(MU1(:, 2));
        [~, isortu2] = sort(MU2(:, 2));
        S1 = [MU1(isortu1, :); MP1(isort1, :)];
        S2 = [MU2(isortu2, :); MP2(isort2, :)];
    end
    
    % data for clonality bars
    CLON1 = [0; cumsum(S1(:, 2))];
    CLON2 = [0; cumsum(S2(:, 2))];         
    n1 = length(CLON1);
    n2 = length(CLON2);
    x1 = reshape([zeros(n1, 1) ones(n1, 1) zeros(n1, 1)]', n1 * 3, 1);
    x2 = reshape([zeros(n2, 1) ones(n2, 1) zeros(n2, 1)]', n2 * 3, 1);
    y1 = reshape((CLON1 * [1 1 1])', n1 * 3, 1);
    y2 = reshape((CLON2 * [1 1 1])', n2 * 3, 1);    
           
    % figure constants
    cw = 60;
    bw = 400;
    left = 100;
    right = 100;
    bot = 50;
    top = 50;
    ph = 840;
    fw = left + cw + bw + cw + right;
    fh = bot + ph + top;
    posname = [0, 0, 1, bot / fh];   
    pos1t = [left             / fw, (bot + ph) / fh, cw / fw, top / fh];
    pos2t = [(left + cw + bw) / fw, (bot + ph) / fh, cw / fw, top / fh];
    pos1  = [left             / fw, bot / fh, cw / fw, ph / fh];    
    pos2  = [(left + cw + bw) / fw, bot / fh, cw / fw, ph / fh];
    posc  = [(left + cw)      / fw, bot / fh, bw / fw, ph / fh];
    pref = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'linewidth', 1,...
        'box', 'on',...
        'xlim', [0 1],...
        'ylim', [0 100]);
    textpref = struct(...
        'fontname', 'arial',...
        'fontsize', 20,...
        'fontweight', 'bold',...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    titpref = struct(...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    maxcol = 300;
    cmap = winter(maxcol);
        
    % make figure
    figure('position', [5 5 fw fh]);

    % database name
    subplot('position', posname);
    text(.5, 50, AU.dbname,...
        'fontname', 'courier',...
        'fontsize', 24,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca, pref, titpref);

    % titles
    subplot('position', pos1t);
    text(.5, 70, label{1}, textpref);
    set(gca, pref, titpref);
    subplot('position', pos2t);
    text(.5, 70, label{2}, textpref);
    set(gca, pref, titpref);

    % make arcs
    subplot('position', posc, 'nextplot', 'add');
    cmap = cmap(randperm(size(cmap, 1)), :);
    for i = 1:length(linconn)
        i1 = find(S1(:, 1) == linconn(i));
        i2 = find(S2(:, 1) == linconn(i));
        [x12, y12] = hbez(...
            0, CLON1(i1), CLON1(i1 + 1),...
            1, CLON2(i2), CLON2(i2 + 1), .5);
        icol = mod(i, maxcol) + 1;
        hf = fill(x12, y12, cmap(icol, :),...
            'edgecolor', cmap(icol, :),...
            'linewidth', .25);
        alpha(hf, .5);
    end
    set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');
    uistack(gca, 'bottom');    

    % make clonality plot
    yt = 0:5:100;
    ytl = cell(length(yt), 1);
    for i = 1:length(yt), ytl(i) = cellstr(num2str(yt(i))); end
    h1 = subplot('position', pos1, 'nextplot', 'add');
    plot(x1, y1, 'k-', 'linewidth', .25);
    set(gca, pref);
    subplot('position', pos2, 'nextplot', 'add');
    plot(x2, y2, 'k-', 'linewidth', .25);
    set(gca, pref);
    set(h1, 'ytick', yt, 'yticklabel', ytl,...
        'fontname', 'arial','fontsize', 14);
    axes(h1);
    hy = ylabel(ylab, 'fontname', 'arial', 'fontsize', 18);
    set(hy, 'position', [-.8 50 1]);

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');
    
return
