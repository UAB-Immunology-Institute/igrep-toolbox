function igplotvj(AU, pop, varargin)

% IGPLOTVJ Create a VJ-usage heatmap for a selected population
%
%    IGPLOTVJ(datasource, pop) creates a VJ-usage heat map for a specified
%    population in the data source.  Columns are v-genes and rows are
%    j-genes. The first argument is a IGREP object.  The second argument is
%    the desired population in the data source to plot.  The frequency of a
%    particular v-gene and j-gene combination is represented by color using
%    the parula color scale (black = 0, blue = low frequency, yellow = high
%    frequency).  Note: hard "0" is black; this function uses alpha so
%    resulting figure should be saved as PDF, not EPS.
%
%    IGPLOTVJ(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies optional
%    parameter name/value pairs.
%        'label'        alternate name to display on plot (instead of
%                       population names)
%        'maxp'         maximum percentage of occurrence of V-J usage (by
%                       sequence or lineage) - Autoscale if not specified
%        'method'       specifies how to count gene usage: must be either
%                       'sequence' (default) or 'lineage'
%
%    Usage:
%
%        IGPLOTVJ(D, 'p2')
%
%        IGPLOTVJ(D, 'p2', 'method', 'lineage', 'maxp', 20, 'label', 'PB')
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
        @(x) validateattributes(x, {'igrep'},...
        {'nonempty'}));
    addRequired(pa, 'pop',...
        @(x) validateattributes(x, {'char'},...
        {'nonempty'}));
    addParameter(pa, 'label', pop,...
        @(x) validatetext(x, 1));
    addParameter(pa, 'maxp', 0,...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>', 0}));
    addParameter(pa, 'method', 'sequence',...
        @(x) validateenum(x, {'sequence', 'lineage'}));
    parse(pa, AU, pop, varargin{:});
    maxp = pa.Results.maxp;
    label = char(pa.Results.label);
    method = char(pa.Results.method);

    % get v- and j-gene usage counts
    a = AU.nvj(pop, method);
    v = a.Properties.RowNames;
    j = a.Properties.VariableNames;
    c = table2array(a);
    nv = length(v);
    nj = length(j);
    frq = (100 * c / sum(c(:)))';
            
    % set max percentage or autoscale to max
    if maxp == 0, maxp = max(frq(:)); end
        
    % if alternate population name empty, use database population name
    if isempty(label), label = pop; end
    
    % figure constants
    pix = 20;
    pw = nv * pix;
    ph = nj * pix;
    bot = 80;
    left = 80;
    sp = 35;
    cw = 20;
    right = 15;
    top = 30;
    fw = left + pw + sp + cw + right;
    fh = bot + ph + top;
    posvlab = [left / fw, 0, pw / fw, bot / fh];
    posjlab = [0, bot / fh, left / fw, ph / fh];
    posmain = [left / fw, bot / fh, pw / fw, ph / fh];
    postop = [left / fw, (bot + ph) / fh, pw / fw, top / fh];
    poscb = [(left + pw + sp) / fw, bot / fh, cw / fw, ph / fh];
    fpr = struct(...
        'fontname', 'arial',...
        'fontsize', 12,...
        'interpreter', 'none');
    pr = struct(...
        'ticklength', [0 0],...
        'xticklabel', {{}},...
        'yticklabel', {{}},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'on',...
        'ydir', 'normal');
                
    % make figure
    figure('position', [20, 20, fw, fh]);

    % figure label
    subplot('position', postop)
    text(.5, .5, [label ' (by ' method ')'], 'fontname', 'arial',...
        'fontsize', 18, 'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca, pr, 'xlim', [0 1], 'ylim', [0 1]);

    % v-gene labels
    subplot('position', posvlab);   
    text(1:nv, .95 * ones(nv, 1), v, fpr, 'rotation', -60);
    set(gca, pr, 'xlim', [.5 nv + .5], 'ylim', [0 1]);

    % j-gene labels
    subplot('position', posjlab);
    text(.95 * ones(nj, 1), 1:nj, j, fpr, 'horizontalalignment', 'right');
    set(gca, pr, 'xlim', [0 1], 'ylim', [.5 nj + .5]);

    % heat map - make "0" black using alpha
    subplot('position', posmain, 'nextplot', 'add');
    iz = ones(size(frq));
    iz(frq == 0) = 0;
    imagesc(frq, 'alphadata', iz, [0 maxp])    
    colormap(parula(2048));
    plot([1 1]' * (.5:1:(nv + .5)), [.5 nj + .5], '-',...
        'color', [1 1 1] * .5, 'linewidth', .5);
    plot([.5 nv + .5], [1 1]' * (.5:1:(nj + .5)), '-',...
        'color', [1 1 1] * .5, 'linewidth', .5);
    set(gca, pr,...
        'xlim', [.5 nv + .5],...
        'ylim', [.5 nj + .5],...
        'color', 'k');

    % heat map colorbar
    subplot('position', poscb);
    imagesc(flipud(linspace(0, 1)'), [0 1]);
    colormap(parula(2048));
    set(gca, 'ticklength', [0 0], 'xticklabel', {}, 'yticklabel', {},...
        'box', 'on');
    title([num2str(sprintf('%.1f', maxp)) ' %'], 'fontname', 'arial',...
        'fontsize', 12, 'fontweight', 'normal');
       
    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return
    
        
        

