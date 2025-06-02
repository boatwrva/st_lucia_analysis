function set_figure_LUCIA(fig)
% Apply custom style settings a MATLAB fig
% Attempt to match "packages\python\student_led_cruise\student_led_cruise\gom_subduction.mplstyle"
% The equivalent fields from the mpl stylesheet are indicted below as
% py: fieldname
% Some matplotlib fields have no direct equivalent. They are:
%   figure.constrained_layout   by default matplotlib is not using
%       constrained_layout
%   figure.constrained_layout.h_pad
%   figure.constrained_layout.w_pad
%   figure.constrained_layout.hspace
%   figure.constrained_layout.wspace
%   figure.dpi                  this needs to be set when saving the figure
%   legend.handlelength         this sets the length of the token in the
%       legend. If you need to do this in matlab see
%       https://www.mathworks.com/matlabcentral/answers/95161-how-can-i-modify-the-lengh-of-the-lines-in-a-legend
%   mathtext.default            matplotlib is using serif and sitx
%   mathtext.fontset
%   mathtext.rm
%   mathtext.it
%   savefig.dpi                 again must be set when saving the figure
%
% Jamie Hilditch July 2024

% py: axes.formatter.use_mathtext
set(fig, 'DefaultAxesTickLabelInterpreter', 'latex');
set(fig, 'DefaultLegendInterpreter', 'latex');
set(fig, 'DefaultTextInterpreter', 'latex');

% py: axes.labelsize, axes.titlesize
set(fig, 'DefaultAxesFontSize', 16);
set(fig, 'DefaultAxesTitleFontSizeMultiplier', 9/7); % titlesize / labelsize -> titlesize = 9

% py: figure.facecolor
set(fig, 'DefaultFigureColor', 'white');

% fonts: ideally we want computer modern roman to match the latex
% py: font.family, font.serif, font.size
set(fig, 'DefaultAxesFontName', 'Times New Roman');
set(fig, 'DefaultTextFontName', 'Times New Roman');
set(fig, 'DefaultTextFontSize', 16);

% py: legend.fontsize, legend.frameon, legend.edgecolor
set(fig, 'DefaultLegendFontSize', 16);
set(fig, 'DefaultLegendBox', 'on');
set(fig, 'DefaultLegendEdgeColor', [1 1 1]);

% py: lines.linewidth, lines.markersize
set(fig, 'DefaultLineLineWidth', 0.5);
set(fig, 'DefaultLineMarkerSize', 5);

% py: xtick.direction, ytick.direction
set(fig, 'DefaultAxesTickDir', 'in');

% Save figure properties
set(fig, 'DefaultFigurePaperPositionMode', 'auto');
set(fig, 'DefaultFigurePosition', [100, 100, 800, 600]); % adjust as needed
set(fig, 'DefaultFigureInvertHardcopy', 'off');
set(fig, 'DefaultFigurePaperType', 'A4');
set(fig, 'DefaultFigurePaperUnits', 'centimeters');
set(fig, 'DefaultFigurePaperSize', [21, 29.7]); % A4 size

end
