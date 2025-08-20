function igplotmuthist(D, binsz, maxf, poplist, cols, varargin)

% IGPLOTMUTHIST mutation frequency histograms for an analysis unit.
%
%    IGPLOTMUTHIST(datasourse, binsize, maxf, pops, cols) plots the
%    distributions of V-gene mutation frequencies for the specified 
%    population in the analysis unit using the specified colors.  The first
%    argument is a IGREP object containing the data.  The second argument
%    is the histogram bin size (in percentages) = "1" is a typical value to
%    use.  The third argument is the maximum value for the histogram (e.g. 
%    40%).  The fourth argument is a cell vector of population names.  The
%    fifth argument is a vector of colors with the same number of elements
%    as the number of populations.
%
%    IGPLOTMUTHIST(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'normalize'    must be either 'yes' or 'no' (default) - specifies
%                       whether to normalize each histogram by the
%                       corresponding number of sequences
%        'label'        an n-element (where n is the number of populations)
%                       cell vector of character arrays or list of strings
%                       of alternate population names - if this is not
%                       specified, the database population names are used
%        'isotype'      a vector of isotypes - one per population - if no
%                       isotype filtering is desired for a particular
%                       population, use '' for that element
%        'region'       a vector of regions for mutation frequency
%                       computation - one per population - valid values are
%                       "v", "fr" and "cdr"; default is "v"
%        'mtype'        a vector of mutation types - one per population -
%                       valid values are "total", "nonsilent" and "silent";
%                       default is "total"
%
%    Usage:
%
%        IGPLOTMUTHIST(D, 1, 25, {'p1', 'p2', 'p3', 'p4'},...
%            {'r', 'g', 'b', 'm'})
%
%        IGPLOTMUTHIST(D, 1, 25, {'p1', 'p1', 'p2', 'p2'},...
%            {'r', [1 .5 .5], 'b', [.5 .5 1]},...
%            'region', {'fr', 'cdr', 'fr', 'cdr'}, 'normalize', 'yes')
%
%        IGPLOTMUTHIST(D, 1, 25, {'p1', 'p1', 'p1'},...
%            {'r', 'b', 'k'},...
%            'region', {'fr', 'fr', 'fr'},...
%            'isotype', {'A', 'G', ''},...
%            'mtype', {'nonsilent', 'nonsilent', 'nonsilent'},...
%            'normalize', 'yes')
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
    addRequired(pa, 'binsz',...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>', 0}));
    addRequired(pa, 'maxf',...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>', 0}));
    addRequired(pa, 'poplist', @(x) validatetext(x));
    n = length(poplist);
    addRequired(pa, 'color', @(x) validatecolor(x, n));    
    addParameter(pa, 'label', poplist, @(x) validatetext(x, n));            
    addParameter(pa, 'normalize', 'no',...
        @(x) validateenum(x, {'yes', 'no'}));    
    addParameter(pa, 'isotype', repmat({''}, 1, n),...
        @(x) validatetext(x, n));
    addParameter(pa, 'mtype', repmat({'total'}, 1, n),...
        @(x) validatetext(x, n));
    addParameter(pa, 'region', repmat({'v'}, 1, n),...
        @(x) validatetext(x, n));    
    parse(pa, D, binsz, maxf, poplist, cols, varargin{:});      
    poplist = cellstr(pa.Results.poplist);
    isolist = cellstr(pa.Results.isotype);
    mtype = cellstr(pa.Results.mtype);
    region = cellstr(pa.Results.region);
    cols = pa.Results.color;
    label = cellstr(pa.Results.label);    
    if strcmp(pa.Results.normalize, 'no')
        nrm = 0;
    else
        nrm = 1;
    end
    
    % retrieve data
    nbin = floor(maxf / binsz);
    y = nan(nbin, n); 
    lab = cell(n, 1);
    for i = 1:n             
        temp = [label{i} ': ' region{i} ', ' mtype{i}];     
        if ~isempty(isolist{i}), temp = [temp ' (' isolist{i} ')']; end
        lab(i) = cellstr(temp);
        if isempty(isolist{i})
            m = D.muthist(poplist{i}, region{i}, mtype{i}, binsz, maxf);
        else
            m = D.muthist(poplist{i}, region{i}, mtype{i}, binsz, maxf,...
                {isolist{i}, '', ''});
        end
        if i == 1, x = m.center; end
        y(:, i) = m{:, 4};
    end
   
    % normalize each column by total number of sequences?
    ylab = 'Number of Sequences';
    if nrm == 1
        y = y ./ sum(y);
        ylab = 'Fraction of Sequences';
    end
    
    % figure constants and set up figure
    left = 70;
    right = 20;
    bottom = 50;
    top = 10;
    pw = 600;
    ph = 400;
    fw = left + pw + right;
    fh = bottom + ph + top;
    figure('position', [30 30 fw fh]);
    subplot('position', [left / fw, bottom / fh, pw / fw, ph / fh],...
        'nextplot', 'add');
    
    % convert into column indices to find which histograms to plot
    h = nan(n, 1);
    for i = 1:n
        h(i) = plot(x, y(:, i), '-', 'color', cols{i}, 'linewidth', 1.5);      
    end
    set(gca,...
        'ticklength', [0 0],...
        'linewidth', 1,...
        'fontname', 'arial',...
        'fontsize', 12,...
        'yticklabel', num2str(get(gca, 'ytick')'),...
        'box', 'off');
    xlabel('Mutation Frequency (%)',...
        'fontname', 'arial',...
        'fontsize', 16);
    ylabel(ylab,...
        'fontname', 'arial',...
        'fontsize', 16);
    
    % make legend
    hleg = legend(h, lab, 'location', 'northeast');
    set(hleg,...
        'box', 'off',...
        'fontname', 'arial',...
        'fontsize', 14,...
        'interpreter', 'none');

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return
    
