function [muth, stat] = muthist(obj, pop, r, t, bin, maxf, varargin)
    % IGREP/MUTHIST Mutation frequency histogram
    %
    % [mut, s] = MUTHIST(obj, pop, r, t, bin, maxf) returns a histogram of
    % mutation frequencies for a particular population (first argument).
    % The second argument is one of "v", "fr", "cdr" or "all" and specifies
    % whether to compute the mutation frequencies over the full v-region or
    % to break it down by FR (FR1 + FR2 + FR3) and CDR (CDR1 + CDR2), or
    % show all three methods.  The third argument is one of "total",
    % "nonsilent", "silent" or "all" and specifies which types of mutations
    % to consider in the mutation frequency calculation ("all" shows it all
    % three ways).  The fourth argument is the bin size (in percent, "1" is
    % a typically reasonable value) and the fifth argument is the maximum
    % frequency to consider (e.g. "40%").  The first output is a Matlab  
    % table where the first row contains the column header.  The first 
    % three columns are the bin (start, end and center), and the subsequent 
    % columns are the binned mutation frequencies of the desired
    % population, regions and types.  The second output argument is a
    % table of means and medians for the desired population, regions and
    % types.
    %
    % [mut, s] = MUTHIST(obj, pop, bin, maxf, filt) only considers
    % sequences meeting the filter constraints as described in the
    % igrep/mut method.
    %
    % Dot notation usage:
    %     [m, s] = D.MUTHIST('pop', 'all', 'all', 1, 40)
    %     [m, s] = D.MUTHIST('pop', 'all', 'nonsilent', 1, 40)
    %     [m, s] = D.MUTHIST('pop', 'fr', 'all', 1, 40)
    %     [m, s] = D.MUTHIST('pop', 'cdr', 'silent', 1, 40)
    %     [m, s] = D.MUTHIST('pop', 'v', 'all', 1, 40,...
    %                  {'G', '', ''})
    %     [m, s] = D.MUTHIST('pop', 'v', 'all', 1, 40,...
    %                  {'', '!IGHV1-2', ''}) 
    
    % Chris Fucile
    % Alex Rosenberg
    % 5 August 2025
    % University of Alabama at Birmingham
    % Department of Biomedical Informatics and Data Science
    % and UAB Immunology Institute
    % Copyright (c) 2025. All rights reserved.
    % This software is offered with no guarantees of any kind.
    
    validateattributes(bin, {'numeric'},...
        {'nonempty', 'scalar', '>', 0});
    validateattributes(maxf, {'numeric'},...
        {'nonempty', 'scalar', '>', 0});
    validatepop(obj, pop);  
    if nargin > 6
        m = mut(obj, pop, r, t, varargin{1});
    else
        m = mut(obj, pop, r, t);
    end    
    e = 0:bin:maxf;
    bctr = e(1:(end - 1)) + (bin / 2);     
    h0 = m.Properties.VariableNames;    
    h = h0(cellfun('isempty', regexp(h0, 'nucleotides')));
    h = h(5:end);
    head = [{'start', 'end', 'center'}, h];    
    muth = nan(length(e) - 1, length(head));    
    muth(:, 1:3) = [e(1:(end - 1))', e(2:end)', bctr'];    
    for i = 1:length(h), muth(:, i + 3) = histcounts(m.(h{i}), e)'; end    
    muth = array2table(muth, 'VariableNames', head);
    s = [{'average'; 'median'},...
        num2cell([mean(m{:, h}, 'omitnan'); median(m{:, h}, 'omitnan')])];
    stat = cell2table(s, 'VariableNames', [{'stat'}, h]);
end  
