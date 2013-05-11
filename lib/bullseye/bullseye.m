% NOTE: Modified slightly by Derek Macklin (macklin@stanford.edu)
% BULLSEYE plots patch data in polar coordinates.  The DATA matrix is
% mapped to a bullseye plot as follows.  The matrix rows are mapped to
% the circumference and matrix columns are mapped to radial direction.  For
% example, the top-left (first row, first column) data point in the matrix is
% mapped to the inner-most radial patch and starts from 0 degrees (far
% right of bullseye, not top) and extends counter-clockwise.
%
% If you use the hard coded cardiac anatomic lables then this first
% row is within the INFEROTLATERAL segment.  The last row is within the
% INFERIOR segment.
%
% SYNTAX: h = bullseye(DATA, varargin)
%
% color       - matrix of data (ROWSxCOLUMNS==CIRCUMFERENCExRADIUS)
% param.N     - number of points used to define each patch
% param.rho   - 2x1 or 1x2 vector of inner and outer radius
% param.tht   - 2x1 or 1x2 vector of beginning and ending circumference in DEGREES
% param.tht0  - 1x1 value in degrees indicating the circumferential location of the first segment (>0 is CCW)
% param.labels- 0 or 1, flag that turns on anatomic labeling
% param.lines - 1x2 or 2x1, determines number of RADIAL x CIRCUMFERENTIAL lines to use for anatomic sector emphasis
%
% h     - handle for the patch object
%
% Example: 
% bullseye(rand(3,5), 'N', 100, 'tht', [45 180]);
%
% Feedback:  Rate this file @ <a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/16458', '-browser')">MATLAB Central File Exchange</a>
%
% DBE 2004/07/21
% DBE 2005/10/24 Minor bug fixes

% http://www.mathworks.com/matlabcentral/fileexchange/16458

function [h, t] = bullseye(data, varargin)

data = flipud(data);  % The matrices are defined from left to right...graph is CCW

param.N          = [];
param.rho        = [];
param.tht        = [];
param.tht0       = [];
param.label      = [];
param.lines      = [];
param.axesHandle = [];

i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'n'
            param.N     = varargin{i+1};   i = i + 2;
        case 'rho'
            param.rho   = varargin{i+1};   i = i + 2;
        case 'tht'
            param.tht   = varargin{i+1};   i = i + 2;
        case 'tht0'
            param.tht0  = varargin{i+1};   i = i + 2;
        case {'label', 'labels'}
            param.label = varargin{i+1};   i = i + 2;
        case 'lines'
            param.lines = varargin{i+1};   i = i + 2;
        case 'axeshandle'
            param.axesHandle = varargin{i+1}; i = i + 2;
        otherwise
            error(['Input argument ', varargin{i}, ' not valid.']);
    end
end

if isempty(param.N   ),       param.N     = 5;        end
if isempty(param.rho ),       param.rho   = [1 5];    end
if isempty(param.tht ),       param.tht   = [0 360];  end
if isempty(param.tht0),       param.tht0  = 0;        end
if isempty(param.label),      param.label = 0;        end
if isempty(param.lines),      param.lines = [];       end
if isempty(param.axesHandle), param.axesHandle = gca; end

[X, Y] = gen_xy(data, param);
h = patch(X, Y, lin(data')', 'EdgeColor', 'none', 'Parent', param.axesHandle);

if ~isempty(param.lines)
    [X, Y] = gen_xy(zeros(param.lines), param);
    h = patch(X, Y, lin(zeros(param.lines)')', 'Parent', param.axesHandle);
    set(h, 'LineWidth', 3, 'FaceColor', 'none');
end

axis equal off;

if param.label, t = bull_labels(param); end

% c=colorbar;

function [X, Y] = gen_xy(data, param)
ind = 1;
X = zeros(2 * param.N, size(data, 2));
Y = zeros(2 * param.N, size(data, 2));
for j = 1:size(data, 1)
    for k = 1:size(data, 2)
        dtht = -diff(param.tht) / size(data, 1);
        tht1 = dtht * (j - 1) + param.tht(1) - param.tht0 - 90;
        tht2 = dtht *  j      + param.tht(1) - param.tht0 - 90;
        
        dr = diff(param.rho) / size(data, 2);
        rho1 = dr * (k - 1) + param.rho(1);
        rho2 = dr *  k      + param.rho(1);
        
        ang = linspace(tht1 / 180 * pi, tht2 / 180 * pi, param.N);
        [arc1.x, arc1.y] = pol2cart(ang, rho1);
        [arc2.x, arc2.y] = pol2cart(ang, rho2);
        X(:, ind) = [arc1.x arc2.x(end:-1:1)];
        Y(:, ind) = [arc1.y arc2.y(end:-1:1)];
        ind = ind + 1;
    end
end

function t = bull_labels(param)

% Labels are added CCW...this alone determines that the first ROW is INFEROLATERAL etc...
% label{1} = 'inferoseptal';
% label{2} = 'inferior';
% label{3} = 'inferolateral';
% label{4} = 'anterolateral';
% label{5} = 'anterior';
% label{6} = 'anteroseptal';
label{1} = 'anterior';
label{2} = 'anteroseptal';
label{3} = 'inferoseptal';
label{4} = 'inferior';
label{5} = 'inferolateral';
label{6} = 'anterolateral';

t = zeros(length(label), 1);
for k = 1:length(label)
    ang = (60 * (k - 1) - param.tht0 - 60) * (pi / 180);
    [x, y] = pol2cart(ang, 0.875 * max(param.rho));
    t(k) = text(x, y, label{k}, 'HorizontalAlignment', 'Center', 'FontSize', 14, 'Parent', param.axesHandle);
    if rem((ang * 180 / pi), 360) >= 180
        set(t(k), 'rotation', ang * (180 / pi) + 90);
    elseif (ang * 180 / pi) == 0 || (ang * 180 / pi) == 360
        set(t(k), 'rotation', ang * (180 / pi) + 90);
    else
        set(t(k), 'rotation', ang * (180 / pi) - 90);
    end
end

%  This function linearizes a matrix (m) of any dimension (eg M=m(:)).
%  If an index vector (ind) is given then the the ind entries of m(:) are
%  returned.
%
% SYNTAX: m = lin(m);
%         m = lin(m, ind);
%
% DBE 12/22/03

function m = lin(m, ind)

m = m(:);

if nargin == 2
    m = m(ind);
end