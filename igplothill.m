function dtab = igplothill(D, q, poplist, isolist, cols, varargin)

% IGPLOTHILL Create Hill diversity profile plots over a range of orders
%
%    IGPLOTHILL(datasource, orders, pop, iso, colors) creates
%    a diversity profile plot over the specified range of orders.  The
%    first argument is a IGREP object.  The second argument is a vector of
%    n orders across which to evaluate the diversity index.  The third
%    argument is a flag: a value of 1 plots evenness (diversity divided by
%    zero-order-diversity); any other value plots diversity.  The next
%    three arguments are all m-element cell vectors of strings
%    corresponding to m desired plots.  The fourth argument is a list of
%    population names, the fifth is a list of isotypes (use '' if to
%    include all), and the sixth is a list of line colors.  The return
%    value is an n-order by m-population table of data that corresponds to
%    the plots.
%
%    IGPLOTHILL(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'evenness'     if this value is 'yes' evenness will be plotted; if
%                       the value is 'no' diversity will be plotted
%        'sequences'    only consider random subsample of n sequences per
%                       population/isotype
%        'label'        an n-element (where n is the number of populations)
%                       cell vector of character arrays or list of strings
%                       of alternate population names - if this is not
%                       specified, the database population names are used
%
%    When plotting diversity, the y-scale is log.  When plotting evenness,
%    the y-scale is linear.
%
%    Usage:
%
%        IGPLOTHILL(D, 0:.1:10, {'pop1', 'pop2'}, {'', 'A'}, {'r', 'b'})
%
%        IGPLOTHILL(D, 0:.1:10, {'pop1', 'pop2'}, {'', 'A'},...
%            {'r', [.5 .5 .5]}, 'evenness', 'yes')
%
%        IGPLOTHILL(D, 0:.1:10, {'pop1', 'pop2'}, {'', 'A'}, {'r', 'b'},...
%            'sequences', 10000, 'label', {'Plasmablast', 'Memory'})
%
%        (note, colors can also be specified as [r g b])
%
%    See:
%        - Hill (1973) Diversity and evenness: a unifying notation and its
%          consequences. Ecology 54:427-432.
%        - Stern et al (2014) B Cells populating the multiple sclerosis
%          brain mature in the draining cervical lymph nodes. Sci Trans Med
%          6(248):248ra107, DOI: 10.1126/scitranslmed.3008879.
%        - Grieff et al (2015) A bioinformatic framework for immune 
%          repertoire diversity profiling enables detection of 
%          immunological status. Genome Med 7:49,
%          DOI: 10.1186/s13073-015-0169-8.
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
    addRequired(pa, 'q',...
        @(x) isnumeric(x));
    addRequired(pa, 'poplist',...
        @(x) validatetext(x));
    n = length(poplist);
    addRequired(pa, 'isolist',...
        @(x) validatetext(x, n));
    addRequired(pa, 'cols',...
        @(x) validatecolor(x, n));
    addParameter(pa, 'evenness', 'no',...
        @(x) validateenum(x, {'yes', 'no'}));
    addParameter(pa, 'label', poplist, @(x) validatetext(x, n));        
    addParameter(pa, 'sequences', [],...
        @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>', 0}));
    parse(pa, D, q, poplist, isolist, cols, varargin{:});    
    poplist = cellstr(poplist);
    isolist = cellstr(isolist);
    label = pa.Results.label(:);    
    nsamp = pa.Results.sequences;
    if strcmp(pa.Results.evenness, 'no')
        ev = 0;
    else
        ev = 1;
    end
    
    % initialize diversity data
    d = nan(length(q), n);
    if isempty(nsamp)
        for i = 1:n
            if ~isempty(isolist{i})
                t = D.ranklin(poplist{i}, isolist{i});
            else
                t = D.ranklin(poplist{i});
            end
            d(:, i) = hilldiv(t{:, 'size'}, q);
        end
    else
        for i = 1:n
            if ~isempty(isolist{i})
                t = D.ranklin(poplist{i}, isolist{i}, nsamp);
            else
                t = D.ranklin(poplist{i}, nsamp);
            end
            if sum(t{:, 2}) < nsamp
                w = poplist{i};
                if ~isempty(isolist{i})
                    w = [w ' (' isolist{i} ')'];
                end
                warning([w ': < ' num2str(nsamp) ' sequences']);
            end
            d(:, i) = hilldiv(t{:, 2}, q);
        end
    end

    % evenness or diversity
    ylab = 'Diversity';
    sc = 'log';
    if ev == 1
        d = d ./ repmat(d(1, :), length(q), 1);
        ylab = 'Evenness';
        sc = 'linear';        
    end
    
    % return value as table
    dtab = array2table([q' d], 'VariableNames', [{'order'}, poplist(:)']);
    
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

    % add plots
    h = nan(n, 1);
    for i = 1:n
        h(i) = plot(q, d(:, i), '-', 'color', cols{i}, 'linewidth', 1.5);
    end
    set(gca,...
        'yscale', sc,...
        'ticklength', [0 0],...
        'linewidth', 1,...
        'fontname', 'arial',...
        'fontsize', 12,...
        'box', 'off');
    xlabel('Order', 'fontname', 'arial', 'fontsize', 16);
    ylabel(ylab, 'fontname', 'arial', 'fontsize', 16);
    
    % make legend
    leg = cell(n, 1);
    for i = 1:n
        if isempty(isolist{i})
            leg(i) = cellstr(['  ' label{i}]);
        else
            leg(i) = cellstr(['  ' label{i} ' (' isolist{i} ')']);
        end
    end
    hleg = legend(h, leg, 'location', 'northeast');
    set(hleg,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'interpreter', 'none',...
        'box', 'off');
    
    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return


