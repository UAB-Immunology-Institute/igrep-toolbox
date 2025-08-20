function igplotisotype(AU, poplist, varargin)

% IGPLOTISOTYPE isotype distribution barplots for populations
%
%    IGPLOTISOTYPE(datasource, pops) generates a figure with stacked bar
%    plots showing the proportional distribution of isotypes for each of
%    the specified populations. The first argument is a IGREP object
%    containing the data. The second argument is a list of populations to
%    consider.  This argument should be empty if you wish to show all
%    populations in the database.
%
%    IGPLOTISOTYPE(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        an n-element (where n is the number of populations)
%                       cell vector of character arrays or list of strings
%                       of alternate population names - if this is not
%                       specified, the database population names are used
%        'keepu'        include unassigned isotypes ("U") - valid choices
%                       are 'yes' (default) and 'no'
%        'color'        use specifed color order as opposed to
%                       Matlab-defined color order - this must be a cell
%                       vector where each cell can be either a
%                       single-character color designation or an rgb
%                       3-tuple - the number of colors must equal the
%                       number of isotypes in the database (and one less if
%                       you are omitting "U"s, and colors in this list will
%                       be ordered on the bars from bottom to top
%
%    Usage:
%
%        IGPLOTISOTYPE(D, {})
%
%        IGPLOTISOTYPE(D, {'pop1', 'pop2'})
%
%        IGPLOTISOTYPE(D, {'pop1', 'pop2'}, 'label', {'MEM', 'ASC'},...
%            'keepu', 'no', 'color', {'r', [0 .6 0], 'b'})
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

    % manually test first argument to see if igrep object
    if ~isequal(class(AU), 'igrep')
        error('first argument must be an IGREP object');
    end
        
    % get list of isotypes and populations in database
    ui = string(AU.sequences.Properties.VariableNames(3:end));
    ni = length(ui);
    pops = AU.upop;
    
    % manually test poplist to get number of populations
    if isempty(poplist)
        poplist = AU.upop;        
    else
        if ~isstring(poplist) || min(size(poplist)) > 1
            error('population list must be a cell vector of strings');
        end
    end
    [~, idx] = ismember(poplist, pops);
    if ~isempty(find(idx == 0, 1))
        error('one or more populations not in database');
    end
    np = length(poplist);
            
    % attribute-value pair option handling
    pa = inputParser;
    addParameter(pa, 'label', poplist, @(x) validatetext(x, np));    
    addParameter(pa, 'color', '', @(x) validatecolor(x));
    addParameter(pa, 'keepu', 'yes', @(x) validateenum(x, {'yes', 'no'}));
    parse(pa, varargin{:});
    label = pa.Results.label(:);
    color = pa.Results.color;
               
    % set up data - test whether including "U"
    if strcmp(pa.Results.keepu, 'no')
        seqs0 = removevars(AU.sequences, 'U');   
        seqs0 = seqs0{:, 3:end};
        ni = ni - 1;   
        % ui = ui(~strcmp(ui, 'U'));
        ui = ui(ui ~= "U");
    else
        seqs0 = AU.sequences{:, 3:end};
    end
    if ~isempty(color) && length(color) ~= ni
        error(['expecting ' num2str(ni) ' colors']);
    end
    
    % data matrix for plotting
    seqs = seqs0(idx, :);
    [np, ni] = size(seqs);
    if np == 1
        total = sum(seqs); 
        comp = 100 * seqs / total;
    else
        total = sum(seqs, 2);
        comp = 100 * seqs ./ repmat(total, 1, ni);
    end
        
    % figure constants 
    ticpix = 7;
    barw = 60;
    pw = barw * np;
    ph = 300;
    left = 65;
    bot = 120;
    right = 120;
    top = 25;
    fw = left + pw + right;
    fh = bot + ph + top;
    pos = [left / fw, bot / fh, pw / fw, ph / fh];
    poslab = [left / fw, 0, pw / fw, bot / fh];
    postop = [left / fw, (bot + ph) / fh, pw / fw, top / fh];
    pref = struct(...
        'box', 'off',...
        'ticklength', [0 0],...
        'xtick', [],...
        'ytick', [],...
        'xcolor', 'w',...
        'ycolor', 'w');
    tic = ticpix / (((pw > ph) * pw) + ((pw <= ph) * ph));
    
    % draw figure window
    figure('position', [20 20 fw fh]);

    % axes for population names
    subplot('position', poslab);
    text(1:np, ones(np, 1) * .9, label,...
        'fontname', 'arial',...
        'fontsize', 16,...
        'interpreter', 'none',...
        'rotation', -45);
    set(gca, pref, 'xlim', [.25 np + .25], 'ylim', [0 1]);
    
    % axes for sequence counts
    subplot('position', postop);
    text(1:np, ones(np, 1) * .5, strip(cellstr(num2str(total(:)))),...
        'fontname', 'arial',...
        'fontsize', 12,...
        'horizontalalignment', 'center');
    set(gca, pref, 'xlim', [.25 np + .25], 'ylim', [0 1]);
    
    % axes for bar plot
    subplot('position', pos);        
    b = bar(1:np, comp, .5, 'stacked');
    set(gca,...
        'box', 'off',...
        'linewidth', 1,...
        'xcolor', 'w',...
        'xtick', [],...
        'xlim', [.25 np + .25],...
        'tickdir', 'out',...
        'ticklength', [tic 0],...
        'fontname', 'arial',...
        'fontsize', 12,...
        'ylim', [0 100]);
    if ~isempty(color), for i = 1:ni, b(i).FaceColor = color{i}; end, end
    ylabel('Isotype Composition', 'fontname', 'arial', 'fontsize', 16);    
    hleg = legend(fliplr(b), fliplr(ui));
    set(hleg,...
        'box', 'off',...
        'fontname', 'arial',...
        'fontsize', 16);
    hleg.Position(1) = (left + pw + (barw / 2)) / fw;
    hleg.Position(2) = bot / fh;
    
    % clean up
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
    
return
