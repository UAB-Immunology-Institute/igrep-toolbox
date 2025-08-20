function fig2pdf(fh, name)

% FIG2PDF save figure to PDF
%
%    FIG2PDF(gcf, filename) saves the current figure (gcf) to a PDF file.
%    If you specify a particular file handle, that figure will be saved.
%
%    Usage:
%
%        FIG2PDF(gcf, '~/Desktop/myfile.pdf')

% Chris Fucile
% Alex Rosenberg
% 5 August 2025
% University of Alabama at Birmingham
% Department of Biomedical Informatics and Data Science
% and UAB Immunology Institute
% Copyright (c) 2025. All rights reserved.
% This software is offered with no guarantees of any kind.

    fh.PaperPositionMode = 'auto';
    pos = fh.PaperPosition;
    fh.PaperSize = [pos(3) pos(4)];
    print(fh, name, '-dpdf');

return