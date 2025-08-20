function igplotalluv3(AU, pop, varargin)

% IGPLOTALLUV3 Create a 3-population alluvial clonality plot
%
%    IGPLOTALLUV3(datasource, pop) creates three clonality plots where
%    common lineages between pairs of populations are connected by bezier
%    ribbons.  The first argument is an IGREP object.  The second argument
%    is a 3-element cell vector of population names as appears in 
%    the specified data source.  Ribbons that connect the first and middle
%    populations only, or the middle and third populations only are
%    rendered in shades of blue-green.  Ribbons that connect all three are
%    rendered in shades of orange-yellow.
%
%    NOTE: this plot does not depict connections between the first and
%    third populations only.  Circular graph layout would be more
%    appropriate for illustrating all pairwise relationships.
%
%    IGPLOTALLUV3(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        alternate names to display on plot (instead of
%                       population names) - a two element cell vector of
%                       char arrays
%        'noconn'       specify max number of sequences to not count as a
%                       connection.  For example, "1" means that the
%                       clonotype in both populations needs to be more than
%                       1 sequence for a ribbon to be drawn.  Default is 0.
%  
%    Usage:
%
%        IGPLOTALLUV3(D, {'p1', 'p2', 'p3'})
%
%        IGPLOTALLUV3(D, {'p1', 'p2', 'p3'}, 'label', {'Mem', 'PC', 'PB'})
%
%        IGPLOTALLUV3(D, {'p1', 'p2', 'p3'}, 'noconn', 2);
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

    % argument validation and attribute-value pair option handling
    pa = inputParser;
    addRequired(pa, 'AU',...
        @(x) validateattributes(x, {'igrep'}, {'nonempty'}));
    addRequired(pa, 'pop',...
        @(x) validatetext(x, 3));
    addParameter(pa, 'label', pop,...
        @(x) validatetext(x, 3));
    addParameter(pa, 'noconn', 0,...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0}));
    parse(pa, AU, pop, varargin{:});    
    pop = cellstr(pop);
    label = cellstr(pa.Results.label);
    t = pa.Results.noconn;

    % get lineage table for all three populations
    lintab = AU.lintable(pop);
    lin = lintab.lineage;
    m = lintab{:, 2:4};
    Z = [lin, 100 * m ./ repmat(sum(m), length(lin), 1), m];

    % to prevent visual artifact by deafult sort for singletons
    Z = Z(randperm(size(Z, 1)), :);

    % make structure array
    L = struct;
    for i = 1:3
        k = find(Z(:, i + 1) > 0);
        L(i).lin = Z(k, 1);
        L(i).sz = Z(k, i + 1);
        L(i).count = Z(k, i + 4);
    end

    % figure constants
    cw = 60;
    bw = 400;
    left = 100;
    right = 100;
    bot = 50;
    top = 50;
    ph = 840;
    fw = left + cw + bw + cw + bw + cw + right;
    fh = bot + ph + top;
    pref = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'linewidth', 1,...
        'box', 'on',...
        'xlim', [0 1],...
        'ylim', [0 100]);
    maxcol = 300;
    cmap = winter(maxcol);
    cmap2 = autumn(maxcol);

    % sort each of the clonality datasets
    S = struct;
    for i = 1:3
        [~, isort] = sort(L(i).sz);
        S(i).lin   = L(i).lin(isort);
        S(i).sz    = L(i).sz(isort);
        S(i).count = L(i).count(isort);
        S(i).y     = [0; cumsum(S(i).sz)];
        n = length(S(i).y);
        z = ones(n, 1);
        S(i).xb    = reshape([z * 0, z, z * 0]', n * 3, 1);
        S(i).yb    = reshape((S(i).y * [1 1 1])', n * 3, 1);
        S(i).pos   = [...
            (left + ((i - 1) * (cw + bw))) / fw,...
            bot / fh, cw / fw, ph / fh];
        S(i).post  = [...
            (left + ((i - 1) * (cw + bw))) / fw,...
            (bot + ph) / fh, cw / fw, top / fh];
    end
    
    % make figure
    figure('position', [5 5 fw fh]);

    % database name
    subplot('position', [0, 0, 1, bot / fh]);
    text(.5, 50, AU.dbname,...
        'fontname', 'courier',...
        'fontsize', 24,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');

    % titles
    for i = 1:3
        subplot('position', S(i).post);
        text(.5, 70, label{i},...
            'fontname', 'arial',...
            'fontsize', 20,...
            'fontweight', 'bold',...
            'horizontalalignment', 'center',...
            'interpreter', 'none');            
        set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');
    end
    
    % different arc categories
%     clin12  = setdiff(intersect(S(1).lin, S(2).lin), S(3).lin);
%     clin23  = setdiff(intersect(S(2).lin, S(3).lin), S(1).lin);
%     clin123 = intersect(intersect(S(3).lin, S(2).lin), S(1).lin);
     clin12  = setdiff(intersect(S(1).lin(S(1).count > t), S(2).lin(S(2).count > t)), S(3).lin(S(3).count > t));
     clin23  = setdiff(intersect(S(2).lin(S(2).count > t), S(3).lin(S(3).count > t)), S(1).lin(S(1).count > t));
     clin123 = intersect(intersect(S(1).lin(S(1).count > t), S(2).lin(S(2).count > t)), S(3).lin(S(3).count > t));     

    % axis positions
    pos12 = [(left + cw) / fw, bot / fh, bw / fw, ph / fh];
    pos23 = [(left + cw + bw + cw) / fw, bot / fh, bw / fw, ph / fh];
    
    % make left to middle arcs
    cmap = cmap(randperm(size(cmap, 1)), :);
    hc12 = subplot('position', pos12, 'nextplot', 'add');
    ribbons(S(1), S(2), clin12, cmap);
    set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');

    % make middle to right arcs
    cmap = cmap(randperm(size(cmap, 1)), :);
    hc23 = subplot('position', pos23, 'nextplot', 'add');
    ribbons(S(2), S(3), clin23, cmap);
    set(gca, pref, 'xcolor', 'w', 'ycolor', 'w', 'box', 'off');

    % make all three arcs
    cmap2 = cmap2(randperm(size(cmap, 1)), :);
    axes(hc12); ribbons(S(1), S(2), clin123, cmap2);
    axes(hc23); ribbons(S(2), S(3), clin123, cmap2);

    % make clonality plot
    yt = 0:5:100;
    ytl = cell(length(yt), 1);
    for i = 1:length(yt), ytl(i) = cellstr(num2str(yt(i))); end
    for i = 1:3
        subplot('position', S(i).pos, 'nextplot', 'add');
        plot(S(i).xb, S(i).yb, 'k-', 'linewidth', .25);
        set(gca, pref);
        if i == 1
            set(gca, 'ytick', yt, 'yticklabel', ytl,...
                'fontname', 'arial',...
                'fontsize', 14);
            hy = ylabel('Cumulative Percentage of Sequences',...
                'fontname', 'arial',...
                'fontsize', 18);
            set(hy, 'position', [-.8 50 1]);
        end
    end

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');

return


    % draw set of bezier ribbons
    function ribbons(SA, SB, clin, cm)
       
        for w = 1:length(clin)
            ia = find(SA.lin == clin(w));
            ib = find(SB.lin == clin(w));
            [xa, xb] = hbez(...
                0, SA.y(ia), SA.y(ia + 1),...
                1, SB.y(ib), SB.y(ib + 1), .5);
            ic = mod(w, size(cm, 1)) + 1;
            h = fill(xa, xb, cm(ic, :),...
                'edgecolor', cm(ic, :), 'linewidth', .25);
            alpha(h, .5);
        end
        
    return



