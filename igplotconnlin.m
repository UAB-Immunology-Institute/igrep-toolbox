function igplotconnlin(D, pops, varargin)

% IGPLOTCONNLIN Create a 2-population shared lineage size plot
%
%    IGPLOTCONNLIN(datasource, pop) creates a scatterplot comparing lineage
%    sizes between two populations in a database (both common and specific
%    to a single population).  Position of a circle corresponds to the
%    number of sequences in each of the two populations.  If there is a
%    non-zero number of sequences for both populations, the circle is red;
%    otherwise, if the lineage involves only one of the two populations,
%    then the circle is gray and is positioned along x = 0 or y = 0.  Size
%    of the circle (log) reflects the number of occurrences (lineages) with
%    the specific numbers of sequences in the two populations.  Axes are
%    log10(1 + number of sequences).  The aspect ratio of the plot is
%    forced to be square and x- and y- ranges are the same.
%
%    IGPLOTCONNLIN(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        alternate names to display on plot (instead of
%                       population names) - a two element cell vector of
%                       char arrays
%        'maxxy'        specify the maximum extent of the x- and y-axes -
%                       this parameter only has effect if the specified
%                       value is greater than the autoscale value.
%  
%    Usage:
%
%        IGPLOTCONNLIN(D, {'p7', 'p11'})
%
%        IGPLOTCONNLIN(D, {'p7', 'p11'}, 'label', {'mem', 'PC'},...
%            'maxxy', 5)
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
    addRequired(pa, 'D',...
        @(x) validateattributes(x, {'igrep'}, {'nonempty'}));
    addRequired(pa, 'pops',...
        @(x) validatetext(x, 2));
    addParameter(pa, 'label', pops,...
        @(x) validatetext(x, 2));
    addParameter(pa, 'maxxy', 0,...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0}));
    parse(pa, D, pops, varargin{:});    
    pops = cellstr(pops);
    label = cellstr(pa.Results.label);
    maxxy = pa.Results.maxxy;

    % get data
    m0 = table2array(D.lintable(pops));
    m = m0(:, 2:end);
    n1 = length(find(m(:, 1) > 0));
    n2 = length(find(m(:, 2) > 0));
    nc = length(find(prod(m, 2) > 0));
    
    % count unique pairs of lineage sizes (add 1 b/c these must be indices)
% 	[p1, p2, n] = find(accumarray(m + 1, 1)); % fails for very large clones
    m1 = mat2str(m(:, 1));
    m2 = mat2str(m(:, 2));
    m1 = strsplit(m1(2:(end - 1)), ';')';
    m2 = strsplit(m2(2:(end - 1)), ';')';
    sep = repmat({'~'}, length(m1), 1);
    z = tabulate(strcat(m1, sep, m2));
    comb = cellfun(@(x) strsplit(x, '~'), z(:, 1), 'UniformOutput', false);
    comb = 1 + cellfun(@str2num, reshape([comb{:}], 2, length(comb))');
    n = cell2mat(z(:, 2));
    p1 = comb(:, 1);
    p2 = comb(:, 2);

    % figure layout
    pw = 500;
    ph = pw;
    left = 50;
    bot = 50;
    top = 15;
    right = 15;
    fw = left + pw + right;
    fh = bot + ph + top;
    pos = [left / fw, bot / fh, pw / fw, ph / fh];

    % make window, axes
    figure('position', [20 20 fw fh]);
    ax = subplot('position', pos, 'nextplot', 'add');

    % plot data
    x = log10(p1);
    y = log10(p2);
    offset = 25;
    scale = 200;
    s = offset + (log10(n) * scale);
    ic = find(prod([x y], 2) ~= 0);
    c = ones(length(x), 3) * .4;
    c(ic, :) = repmat([1 0 0], length(ic), 1);
    s = scatter(x, y, s, c);
    
    % configure datatips
    k = get(ax, 'children');
    xdata = (10 .^ k.XData) - 1;
    ydata = (10 .^ k.YData) - 1;
    sdata = 10 .^ ((k.SizeData - offset) / scale);
    s.DataTipTemplate.DataTipRows(1) = dataTipTextRow('x seqs', xdata);
    s.DataTipTemplate.DataTipRows(2) = dataTipTextRow('y seqs', ydata);
    s.DataTipTemplate.DataTipRows(3) = dataTipTextRow('n', sdata);
    
    % force same axis limits for x and y; autoscale or specified
    b = max([xlim ylim]);    
    if maxxy ~= 0 && maxxy >= b
        b = maxxy;
    end
    a = -(b * .05);
    
    % format plot
    set(gca,...
        'xlim', [a b],...
        'ylim', [a b],...
        'ticklength', [0 0],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'linewidth', 1,...
        'box', 'off');
    tx = [label{1} ' lineages (' num2str(nc) ' / ' num2str(n1) ' common)'];
    ty = [label{2} ' lineages (' num2str(nc) ' / ' num2str(n2) ' common)'];
    xlabel(tx, 'fontname', 'arial', 'fontsize', 16, 'interpreter', 'none');        
    ylabel(ty, 'fontname', 'arial', 'fontsize', 16, 'interpreter', 'none');

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return




