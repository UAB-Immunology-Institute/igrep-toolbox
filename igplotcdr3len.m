function igplotcdr3len(AU, pop, seq, method, varargin)

% IGPLOTCDR3LEN CDR3 length distribution histogram for specific population
%
%    IGPLOTCDR3LEN(datasource, pop, f) generates a histogram of CDR3
%    lengths for a particular population. The first argument is a IGREP
%    object containing the data. The second argument the population to
%    consider.  The third argument is either 'nt' or 'aa' for nucleotide or
%    amino acid, respectively.  The fourth argument is either 'sequence' or
%    'lineage' to specify what to use to generate the histogram.
%
%    IGPLOTCDR3LEN(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        an alternate name for the population - if this is 
%                       is not specified, the database population names are
%                       used
%        'facecolor'    face color for histogram bars specified as either a
%                       single-character or a rgb triplet.
%        'filter'       filter sequences by isotype, V-Gene or J-Gene -
%                       this is a 3-element cell array of strings where the
%                       first element is the isotype to filter, the second
%                       is the V-Gene and the third is the J-Gene - an
%                       empty character array prevents filtering for that
%                       item, and if the V-Gene or J-Gene is preceded by
%                       "!" then allow ecverything but
%
%    Usage:
%
%        IGPLOTCDR3LEN(D, 'mypop', 'nt', 'lineage')
%
%        IGPLOTCDR3LEN(D, 'mypop', 'aa', 'lineage',...
%            'label', 'Plasmablast', 'facecolor', [.2 .2 .8])
%
%        IGPLOTCDR3LEN(D, 'mypop', 'aa', 'lieage',...
%            'label', 'Plasmablast', 'facecolor', 'b',...
%            'filter', {'A', '', ''})
%
%        IGPLOTCDR3LEN(D, 'mypop', 'aa', 'sequence',...
%            'label', 'Plasmablast', 'facecolor', 'b',...
%            'filter', {'G', 'IGHV4-34', 'IGHJ4'})
%
%        IGPLOTCDR3LEN(D, 'mypop', 'aa', 'sequence',...
%            'label', 'Plasmablast', 'facecolor', 'b',...
%            'filter', {'M', '!IGHV1-2', ''})
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
        @(x) validateattributes(x, {'char', 'string'},...
        {'nonempty'}));
    addRequired(pa, 'seq',...
        @(x) validateattributes(x, {'char'},...
        {'nonempty'}));
    addRequired(pa, 'method',...
        @(x) validateenum(x, {'lineage', 'sequence'}));
    addParameter(pa, 'label', pop,...
        @(x) validatetext(x, 1));
    addParameter(pa, 'filter', {},...
        @(x) validateattributes(x, {'cell'},...
        {'nonempty', 'vector', 'numel', 3}));
    addParameter(pa, 'facecolor', [.7 .7 .7],...
        @(x) validatecolor(x));     
    parse(pa, AU, pop, seq, method, varargin{:});
    label     = pa.Results.label;
    filt      = pa.Results.filter;
    facecolor = pa.Results.facecolor;

    % if filtering sequences
    if isempty(filt)
        c = AU.cdr3len(pop, seq);
        t = char(label);
    else
        c = AU.cdr3len(pop, seq, filt);    
        t = {char(label)};    
        if ~isempty(filt{1}), t = [t ['Isotype = ' filt{1}]]; end
        if ~isempty(filt{2}), t = [t ['V-Gene = '  filt{2}]]; end
        if ~isempty(filt{3}), t = [t ['J-Gene = '  filt{3}]]; end    
        t = strjoin(t, ', ');
    end
    m = c{:, ['cdr3len_' seq]};
    n = c{:, 'n'};

    % if tabulating by lineage, use first occurrence of lin to get cdr len
    if isequal(method, 'lineage')    
        m = c{:, ['cdr3len_' seq]};
        val = m;
        num = [num2str(length(val)) ' Lineages'];
    else
        val = repelem(m, n);
        num = [num2str(length(val)) ' Sequences'];        
    end

    % nucleotides or amino acids
    if isequal(seq, 'nt')
        ed = 1.5:3:97.5;
        bc = ed(2:end) - 1.5;
        xlab = '# Nucleotides in CDR3';
    else
        ed = .5:1:32.5;
        bc = ed(2:end) - .5;
        xlab = '# Amino Acids in CDR3';
    end

    % figure constants
    pw = 650;
    ph = 250;
    left = 70;
    right = 10;
    bot = 60;
    top = 60;
    fw = left + pw + right;
    fh = bot + ph + top;
    pos  = [left / fw, bot / fh, pw / fw, ph / fh];
    pos1 = [left / fw, (bot + ph) / fh, pw / fw, top / fh];

    % draw figure window
    figure('position', [20 20 fw fh]);

    % axes for title text
    subplot('position', pos1);
    text(.5, .8, t,...
        'fontname', 'arial',...
        'fontsize', 16,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    text(.5, .4, num,...
        'fontname', 'arial',...
        'fontsize', 14,...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca,...
        'ticklength', [0 0],...
        'xtick', [],...
        'ytick', [],...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');

    % axes for histogram
    ax = subplot('position', pos);
    % histogram(m, ed, 'facecolor', facecolor, 'facealpha', 1)
    histogram(val, ed, 'facecolor', facecolor, 'facealpha', 1)
    set(ax,...
        'xlim', [min(ed) max(ed)],...
        'xtick', bc,...
        'ticklength', [0 0],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'ylim', [0 max(ylim)],...
        'linewidth', 1,...
        'box', 'off');
    ax.YAxis.Exponent = 0;
    ht = xlabel(xlab, 'fontname', 'arial', 'fontsize', 16);
    ht.Position(2) = -.12 * max(ylim);

    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 

return


