function igplotmh(AU, varargin)

% IGPLOTMH Create a pairwise Morisita-Horn overlap heat map
%
%    IGPLOTMH(datasource) creates an nxn heat map of Morisita-Horn overlap
%    indices for every pair of populations in a data source.  The first
%    argument is a IGREP object.  Note: hard "0" is black; this function
%    uses alpha so resulting figure should be saved as PDF, not EPS.
%
%    IGPLOTMH(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%    parameter name/value pairs.
%        'poporder'     alternate ordering (or subset) of populations -
%                       supply list of indicies of populations in the
%                       desired order (relative to the original order
%                       without using this option); the largest index 
%                       cannot exceed the number of populations
%
%    Usage:
%
%        IGPLOTMH(D)
%
%        IGPLOTMH(D, 'poporder', [5 8 7 2])
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
    addParameter(pa, 'poporder', AU.upop,...
        @(x) validatetext(x));
        
    parse(pa, AU, varargin{:});    
    pop = pa.Results.poporder;
    
    % if specifying popultions, make sure unique and in database
    if ~isempty(setdiff(pop, AU.upop))
        error('one or more populations not in database');
    end
    if length(pop) ~= length(unique(pop))
        error('specified population list not unique');
    end

    % all-by-all morisita-horn overlap matrix
    o = AU.mhoverlap;

    % population names
    p = o.Properties.RowNames;

    % matrix
    m = table2array(o);
    
    % number of populations
    n = length(p);

    % get indices of retained populations if subset specified
    if nargin > 1
        k = nan(length(pop), 1);
        for i = 1:length(pop), k(i) = find(strcmp(pop{i}, AU.upop)); end        
        p = p(k);
        m = m(k, k);
        n = length(k);
    end
    
    % colormap
    ncol = 2048;    
    cmap = parula(ncol);
    darkness = mean(cmap, 2);
    
    % column major order of cell darknesses (to determine black/white font)
    idx = round((m * (ncol - 1)) + 1);
    dmap = darkness(idx);
    dmap(dmap > .5) = 1;
    dmap(dmap <= .5) = 0;
    dmap = dmap(:);   
    
    % column major order of cell values converted to string
    h = strsplit(sprintf('%.2f|', m), '|');
    h = h(1:(end - 1))';    

    % figure constants
    pix = 25;
    pw = n * pix;
    sp = 6;
    left = 100;
    top = 100;
    right = 50;
    bot = 40;
    fw = left + sp + pw + right;
    fh = bot + pw + sp + top;
    pref = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    fpref = struct(...
        'fontname', 'arial',...
        'fontsize', 12,...
        'interpreter', 'none');
    tpref = struct(...
        'fontname', 'arial',...
        'fontsize', 8,...
        'horizontalalignment', 'center');
    
    % create figure
    figure('position', [20 20 fw fh]);

    % row labels
    subplot('position', [0, bot / fh, left / fw, pw / fh]);
    ht = text(ones(n, 1), 1:n, p, fpref, 'horizontalalignment', 'right');
    set(gca, pref, 'xlim', [0 1], 'ylim', [.5 n + .5], 'ydir', 'reverse');

    % determine max text extent and adjust figure dimensions accordingly
    te = cell2mat(get(ht, 'extent'));
    maxte = max(te(:, 3));
    left = left * maxte * 1.02;
    if left < 100, left = 100; end
    top = left * .9;
    right = top * cos(pi / 3);
    fw = left + pw + right;
    fh = bot + pw + top;
    set(gcf, 'position', [20 20 fw fh]);
    set(gca, 'position', [0, bot / fh, left / fw, pw / fh]);    

    % column labels
    subplot('position', [(left + sp) / fw, (bot + pw + sp) / fh, pw / fw, top / fh]);
    text(1:n, zeros(n, 1), p, fpref, 'rotation', 60);
    set(gca, pref, 'xlim', [.5 n + .5], 'ylim', [0 1]);

    % database name
    subplot('position', [(left + sp) / fw, 0, pw / fw, bot / fh]);
    text(.5, .5, AU.dbname,...
            'fontname', 'courier',...
            'fontsize', 14,...
            'horizontalalignment', 'center',...
            'interpreter', 'none');
    set(gca, pref, 'xlim', [0 1], 'ylim', [0 1]);

    % manually create heatmap
    subplot('position', [(left + sp) / fw, bot / fh, pw / fw, pw / fh],...
        'nextplot', 'add');
    
    % find hard zeros and make those cells completely transparent (alpha)
    iz = ones(size(m));
    iz(m == 0) = 0;
    imagesc(m, 'alphadata', iz, [0 1])    
    colormap(cmap);
    plot([1 1]' * (.5:1:(n + .5)), [.5 n + .5], '-',...
        'color', [1 1 1] * .5, 'linewidth', .5);
    plot([.5 n + .5], [1 1]' * (.5:1:(n + .5)), '-',...
        'color', [1 1 1] * .5, 'linewidth', .5);
    
    % print m-h values in each cell
    [i, j] = find(ones(n));
    text(j(dmap == 0), i(dmap == 0), h(dmap == 0), tpref, 'color', 'w');
    text(j(dmap == 1), i(dmap == 1), h(dmap == 1), tpref, 'color', 'k');
     
    % clean up axes
    set(gca, pref,...
        'xlim', [.5 n + .5],...
        'ylim', [.5 n + .5],...
        'color', 'k',...
        'ydir', 'reverse');
       
    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return



        