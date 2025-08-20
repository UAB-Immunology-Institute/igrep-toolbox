function igplotclonalityiso(D, poplist, iso, s, varargin)

% IGPLOTCLONALITYISO Create many clonality plots highlighting isotypes
%
%    IGPLOTCLONALITYISO(D, pops, isotypes, col) generates multiple
%    clonality plots on a single figure for populations in the specified
%    data source and colorizes lineages by percent composition of the 
%    specified isotype.  The first argument is a IGREP object.  The second
%    argument is a list of population names to plot - they are plotted in
%    the order that they appear in this list, from top to bottom.  The
%    third argument is the isotype to highlight. The fourth argument is the
%    color to use for 100% composition of the specified isotype - it is
%    used to generate a gradient of the specified color (100%) to white
%    (0%) - this should be a three-element vector rather than the
%    single-letter color designation available to other functions.
%
%    IGPLOTCLONALITYISO(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        an n-element (where n is the number of populations)
%                       cell vector of character arrays or list of strings
%                       of alternate population names - if this is not
%                       specified, the database population names are used
%        'minx'         the minimum value along the x-axis t display
%                       (default is 0) - must be between 0 and 100
%
%    Usage:
%
%        IGPLOTCLONALITYISO(D, {'pop1', 'pop2', 'pop3'}, 'G', [.8 0 0])
%
%        IGPLOTCLONALITYISO(D, {'pop1', 'pop2', 'pop3'}, 'G', [.8 0 0],...
%            'label', {'Naive', 'Memory', 'ASC'}, 'minx', 20)
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
    addRequired(pa, 'poplist',...
        @(x) validatetext(x));
    addRequired(pa, 'iso',...
        @(x) validateattributes(x, {'char'}, {'row'}));
    addRequired(pa, 's',...
        @(x) validateattributes(x, {'numeric'},...
        {'row', 'numel', 3, '>=', 0, '<=', 1}));
    n = length(poplist);
    addParameter(pa, 'label', poplist,...
        @(x) validatetext(x, n));
    addParameter(pa, 'minx', 0,...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0, '<', 100}));
    parse(pa, D, poplist, iso, s, varargin{:});    
    poplist = cellstr(poplist);
    t = cellstr(pa.Results.label);
    xmin = pa.Results.minx;
        
    % flip poplist so plots start from top, not bottom
    poplist = flipud(poplist(:))';
    t = flipud(t(:))';
            
    % colormap - 0 to 100    
    if s(1) < 1, r = (1:-(1 - s(1))/99:s(1))'; else, r = ones(100, 1); end
    if s(2) < 1, g = (1:-(1 - s(2))/99:s(2))'; else, g = ones(100, 1); end
    if s(3) < 1, b = (1:-(1 - s(3))/99:s(3))'; else, b = ones(100, 1); end
    cmap = [r g b];
        
    % populate histograms
    H = struct;
    for i = 1:n
        a = D.rankliniso(poplist{i}); 
        
        k = find(strcmp(iso, a.Properties.VariableNames(3:end)), 1);
        
%         k = find(strcmp(iso, a(1, :)));
        if isempty(k)
            error(['isotype ' iso ' not in database']);
        end
%         H(i).lins  = cell2mat(a(2:end, 2));
%         H(i).linsi = cell2mat(a(2:end, k));
        H(i).lins  = a.total;
        H(i).linsi = a.(iso);
        H(i).color = cmap(1 + round((H(i).linsi ./ H(i).lins) * 99), :);        
    end
    
    % figure constants
    left = 12;
    pw = 1200;
    ph = 30;
    bot = 100;
    sp = 15; 
    hsp = 7;
    top = 2;
    right = 10;
    cbarw = 250;
    cbarh = 30;
    cbot = 5;
    fw = left + pw + hsp + right + hsp;
    fh = bot + (n * ph) + ((n - 1) * sp) + top;
    yoff = bot + ((0:(n - 1)) * (ph + sp))';
    pos1 = [...
        repmat((left + pw + hsp) / fw, n, 1),...
        yoff / fh,...
        repmat(right / fw, n, 1),...
        repmat(ph / fh, n, 1)];
    pos2 = [...
        repmat(left / fw, n, 1),...
        yoff / fh,...
        repmat(pw / fw, n, 1),...
        repmat(ph / fh, n, 1)];
    pos3 = [...
        (left + pw - cbarw) / fw,...
        cbot / fh,...
        cbarw / fw,...
        cbarh / fh];

    % threshold for drawing isotype colors (norm lin size - %) - 1/2 pixel
    thr = (100 - xmin) * .5 / pw;
        
    % subplot preferences
    pref = struct(...
        'ticklength', [0 0],...
        'yticklabel', {{}},...
        'xlim', [xmin 100],...
        'ylim', [0 1],...
        'box', 'on',...
        'linewidth', 1.5,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'layer', 'top');
    preft = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off',...
        'xlim', [0 1],...
        'ylim', [0 1]);
    prefl = struct(...
        'fontname', 'arial',...
        'fontsize', 16,...
        'interpreter', 'none');
    prefc = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'box', 'on',...
        'linewidth', .5,...
        'xlim', [.5 100.5],...
        'ylim', [0 1],...
        'layer', 'top');
    
    % figure window
    figure('position', [20 20 fw fh]);
    
    % handles for population labels and axes
    ht = nan(n, 1);
    ha1 = nan(n, 1);
    ha2 = nan(n, 1);
    
    % make each clonality bar
    for i = 1:n
        
        % population label
        ha1(i) = subplot('position', pos1(i, :));
        ht(i) = text(0, .5, t{i}, prefl);
        set(gca, preft);
            
        % clonality bar    
        ha2(i) = subplot('position', pos2(i, :), 'nextplot', 'add');
        normlin = 100 * H(i).lins ./ sum(H(i).lins);
        c = cumsum(normlin);
        ibig = find(normlin > thr);
        if ~isempty(ibig)
            for q = 2:length(ibig)
                lup = c(ibig(q));
                llo = c(ibig(q) - 1);
                xf = [llo llo lup lup llo];
                yf = [0 1 1 0 0];
                fill(xf, yf, H(i).color(ibig(q), :), 'edgecolor', 'none');
            end    
        end
        nc = length(c);
        y = reshape((ones(nc, 1) * [0 1 0])', (nc * 3), 1);
        x = reshape([c c c]', (nc * 3), 1);
        plot(x, y, '-', 'color', [1 1 1] * .2, 'linewidth', .25);
        set(gca, pref);
        if i > 1, set(gca, 'xticklabel', {}); end
        if i == 1
            xlabel({' '; 'Cumulative Percentage of Sequences'},...
                'fontname', 'arial',...
                'fontsize', 16);
        end
        
    end
            
    % make colorbar
    ha3 = subplot('position', pos3, 'nextplot', 'add');
    for i = 1:100
        rectangle('position', [i - .5, 0, 1, 1],...
            'facecolor', cmap(i, :),...
            'edgecolor', cmap(i, :));
    end
    set(gca, prefc);
    title(['Percent ' iso], 'fontname', 'arial', 'fontsize', 16);
    
    % adjust figure position and axis positions based on label extents
    te = cell2mat(get(ht, 'extent'));
    maxte = max(te(:, 3));
    right = right * maxte;
    fw = left + pw + hsp + right + hsp;
    fh = bot + (n * ph) + ((n - 1) * sp) + top;
    set(gcf, 'position', [20 20 fw fh]);
    pos1 = [...
        repmat((left + pw + hsp) / fw, n, 1),...
        yoff / fh,...
        repmat(right / fw, n, 1),...
        repmat(ph / fh, n, 1)];
    pos2 = [...
        repmat(left / fw, n, 1),...
        yoff / fh,...
        repmat(pw / fw, n, 1),...
        repmat(ph / fh, n, 1)];
    pos3 = [...
        (left + pw - cbarw) / fw,...
        cbot / fh,...
        cbarw / fw,...
        cbarh / fh];
    for i = 1:n
        set(ha1(i), 'position', pos1(i, :));
        set(ha2(i), 'position', pos2(i, :));
    end
    set(ha3, 'position', pos3);   
    
    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w');     
    
return    

