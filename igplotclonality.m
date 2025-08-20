function r = igplotclonality(AU, pop, varargin)

% IGPLOTCLONALITY Create a clonality plot for a sample/population
%
%    IGPLOTCLONALITY(datasource, pop) creates a clonality plot for the
%    specified population in a specified data source.  Lineages are ordered
%    by size (from bottom to top) and the normalized lineage size (100 *
%    number of sequences in lineage / total number of sequences) is shown
%    on the x-axis.  The first argument is a IGREP object.  The second
%    argument is a string specifying the population to plot.  The return
%    value is a three-element vector containing the number of sequences,
%    the number of lineages, and the number of "expanded" lineages (i.e.
%    the number of lineages that are above the "knee" in the clonality
%    curve.
%
%    IGPLOTCLONALITY(..., 'PARAM1', val1, 'PARAM2', val2, ...) specifies
%    optional parameter name/value pairs.
%        'label'        alternate name to display on plot (instead of
%                       population name)
%        'maxx'         maximum value to display on x-axis (autoscales if
%                       unspecified or if value is "0")
%        'miny'         minimum y-value to plot (range is from 0 to 100 if
%                       unspecified)
%        'linecolor'    1x3 numeric RGB color specification for lines - 
%                       each element must be between 0 and 1; can also use
%                       Matlab's standard single character color
%                       designations
%        'knee'         show knee on clonality plot, above which lineages
%                       could be considered "expanded" - this is computed
%                       as follows: a line is drawn from the upper right to
%                       lower left on the clonality curve, and then
%                       distances from each point on the curve to the line
%                       (and orthogonal to that line) are found - the knee
%                       is the value corresponding to the largest distance.
%        'v'            highlight lineages that have specified v-gene
%        'j'            highlight lineages that have specified j-gene
%        'id'           highlight specified lineage (numeric lineage ID)
%        'facecolor'    1x3 numeric RGB color specification for highlighted
%                       lineages - each element must be between 0 and 1;
%                       can also use Matlab's standard single character 
%                       color designations
%
%    Usage:
%
%        IGPLOTCLONALITY(D, '42PB')
%
%        IGPLOTCLONALITY(D, '42PB, 'label', 'Plasmablast', 'miny', 20,...
%            'maxx', 25, 'linecolor', [.8 0 0], 'knee', 'on')
%
%        IGPLOTCLONALITY(D, '42PB', 'v', 'IGHV4-34', 'j', 'IGHJ3')
%
%        IGPLOTCLONALITY(D, '42PB', 'v', 'IGHV1-2', 'facecolor', [1 .6 0])
%
%        IGPLOTCLONALITY(D, '42PB', 'id', 1234, 'facecolor', 'r')
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
        @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
    addParameter(pa, 'v', '', @(x) validatetext(x, 1));
    addParameter(pa, 'j', '', @(x) validatetext(x, 1));
    addParameter(pa, 'id', [], @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>', 0}));        
    addParameter(pa, 'label', pop, @(x) validatetext(x, 1));
    addParameter(pa, 'linecolor', 'b', @(x) validatecolor(x));
    addParameter(pa, 'facecolor', [1 .5 .5], @(x) validatecolor(x));
    addParameter(pa, 'miny', 0, @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0, '<', 100}));
    addParameter(pa, 'maxx', 0, @(x) validateattributes(x, {'numeric'},...
        {'nonempty', 'scalar', '>=', 0}));
    addParameter(pa, 'knee', 'off',...
        @(x) validateenum(x, {'on', 'off'}));
    parse(pa, AU, pop, varargin{:});
    miny = pa.Results.miny;
    maxx = pa.Results.maxx;
    linecolor = pa.Results.linecolor;
    label = pa.Results.label;
    knee = char(pa.Results.knee);
    vfilt = char(pa.Results.v);
    jfilt = char(pa.Results.j);
    lfilt = pa.Results.id;
    hlcolor = pa.Results.facecolor;
    
    % if specifying filter for v or j gene, check if valid gene
    if ~isempty(vfilt) && ~ismember(vfilt, AU.uv)
        error([vfilt ' vgene not in database']);
    end
    if ~isempty(jfilt) && ~ismember(jfilt, AU.uj)
        error([jfilt ' jgene not in database']);
    end
    
    % get lineage sizes for specified pop
    if ~isempty(vfilt) || ~isempty(jfilt)        
        t = flipud(AU.lintable({pop}, 2));
        v = t.v;
        j = t.j;
    else
        t = flipud(AU.lintable({pop}));
    end    
    sz = t.(pop);
    lid = t.lineage;
    sz = sz';
    
    % compute cumulative lineage sizes
    cumseq = cumsum(sz);
    nseq = max(cumseq);
    nlin = length(sz);
    cpercent = 100 * cumseq / nseq;
    normlinsize = 100 * sz / nseq;
        
    % coords for cumulative plot
    x = normlinsize';          
    y = cpercent';  
    
    % normalized lineage size threshold to plot
    thr = .01;
    
    % compute trapezoids for highlighting
    ftext = {};
    if ~isempty(vfilt) || ~isempty(jfilt) || ~isempty(lfilt)
        k = find(x > thr);
        if ~isempty(vfilt)
            k = intersect(find(strcmp(v, vfilt)), k);
            ftext = [ftext, {['V = ' vfilt]}];
        end
        if ~isempty(jfilt)
            k = intersect(find(strcmp(j, jfilt)), k);
            ftext = [ftext, {['J = ' jfilt]}];
        end
        if ~isempty(lfilt)
            k = intersect(find(lfilt == lid), k);
            ftext = [ftext, {['Lineage = ' num2str(lfilt)]}];
        end
        xt = [zeros(length(k), 1), x(k - 1), x(k), zeros(length(k), 2)]';
        yt = [y(k - 1), y(k - 1), y(k), y(k), y(k - 1)]';        
    end
    ftext = strjoin(ftext, ', ');
    
    % draw horizontal clonality bars as single line - way faster
    z = zeros(length(x), 1);
    xline = reshape([z x z]', length(x) * 3, 1);
    yline = reshape([y y y]', length(y) * 3, 1);   
    
    % figure out where "knee" in curve is
    m = 100 / x(end);
    b = y + (x / m);
    x0 = b / (m + (1 / m));
    y0 = m * x0;
    dd = sqrt((x - x0).^2 + (y - y0).^2);
    sgn = x0 - x;
    ineg = find(sgn < 0);
    dd(ineg) = -1 * dd(ineg);            
    idmax = find(dd == max(dd));
    
    % return values: # seq, # lineages, # expanded lineages (above knee)
    r = [nseq, nlin, length(x) - idmax];    

    % if necessary, autoscale x-axis
    if maxx == 0, maxx = max(normlinsize) * 1.002; end

    % figure constants (in pixels)
    sw = 200;
    sh = 550;
    spleft = 65;
    spright = 15;
    spbot = 60;
    sptop = 10;    
    th = 100;
    
    % font sizes
    fonttic = 16;
    fonttitle = 20;
    fontsubt = 16;
    fontax = 18;
 
    % calculate figure dimensions
    fw = spleft + sw + spright;
    fh = spbot + sh + sptop + th;
    
    % plot positions
    posMain  = [spleft / fw, spbot / fh, sw / fw, sh / fh];
    posTitle = [spleft / fw, (spbot + sh + sptop) / fh, sw / fw, th / fh];            
       
    % create figure
    figure('position', [60 60 fw fh]);

    % plot to hold figure title
    subplot('position', posTitle);
    text(0, .55, 'Sequences:',  'fontname', 'arial', 'fontsize', fontsubt);
    text(0, .35, 'Lineages:', 'fontname', 'arial', 'fontsize', fontsubt);  
    text(1, .55, num2str(nseq), 'fontname', 'arial',...
        'fontsize', fontsubt, 'horizontalalignment', 'right');
    text(1, .35, num2str(nlin), 'fontname', 'arial',...
        'fontsize', fontsubt, 'horizontalalignment', 'right');
    text(.5, .10, ftext, 'fontname', 'arial',...
        'fontsize', fontsubt, 'horizontalalignment', 'center');
    text(.5, .8, label,...
        'fontname', 'arial',...
        'fontsize', fonttitle,...
        'fontweight', 'bold',...
        'horizontalalignment', 'center',...
        'interpreter', 'none');
    set(gca,...
        'xlim', [0 1],...
        'ylim', [0 1],...
        'ticklength', [0 0],...
        'xticklabel', {},...
        'yticklabel', {},...
        'xcolor', 'w',...
        'ycolor', 'w',...
        'box', 'off');
    
    % main plot
    hax = subplot('position', posMain, 'nextplot', 'add');    
    if ~isempty(vfilt) || ~isempty(jfilt) || ~isempty(lfilt)    
        fill(xt, yt, hlcolor);    
    end
    plot(xline, yline, '-', 'linewidth', .75, 'color', linecolor);
    plot(x, y, '-', 'linewidth', 3.5, 'color', linecolor);        
    plot([1 1] * x(end), [0 100], '--', 'linewidth', 1.5, 'color',...
        linecolor);
    if isequal(knee, 'on')
        plot(x(idmax), y(idmax), 'o', 'color', linecolor,...
            'linewidth', 1, 'markersize', 15);
    end
    set(gca,...
        'linewidth', 1,...
        'ticklength', [0 0],...
        'box', 'off',...
        'fontname', 'arial',...
        'fontsize', fonttic,...
        'xlim', [0 maxx],...
        'ylim', [miny 100 + ((100 - miny) * .001)]);  
    ylabel('Cumulative Percentage of Sequences', 'fontname', 'arial',...
        'fontsize', fontax)
    xlabel('Norm. Lineage Size (%)', 'fontname', 'arial',...
        'fontsize', fontax);
    
    % clean up
    set(hax, 'layer', 'top');
    set(gcf, 'inverthardcopy', 'off', 'color', 'w'); 
        
return
                      
    
    