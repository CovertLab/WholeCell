function tick2text(varargin)
%TICK2TEXT Changes tick labels to text objects
%
% tick2text
% tick2text(ax)
% tick2text(ax, param1, val1, ...)
%
% This function creates text objects that mimic x, y, and z tick labels,
% while the original tick labels are removed.  The text object format
% allows more modification (such as color, rotation, etc) than the original
% tick labels, and also allows TeX strings to be utilized in tick labels
% (see xformat and yformat inputs for more details).    
%
% The handles of these new ticks are saved in the application data of the
% axis as 'XTickText' and/or 'YTickText'.  Using these handles, the text
% objects can be modified as a group or individually (this method is used,
% rather than simply returning the handles of the new text objects, since
% zooming and panning will cause new objects to be created and old ones to
% be deleted). 
%
% The modifications to the tick labels are maintained during zooming and
% panning as well.  Note that there will be a delay in the relabeling when
% panning, i.e. the old labels will drag with the mouse until you release
% the button, at which point the new labels will appear.  Also note that
% the new tick label format cannot be updated if you reset tick
% locations or axis limits manually, so apply tick2text to an axis after
% these adjustments have already been made.
%
% Input variables:
%
%   ax:         handle of axis to be modified (default to current axis)
%
%   param/val:  optional parameter/value pairs
%
%               axis:           string containing characters 'x', 'y',
%                               and/or 'z', indicating whether to replace
%                               xticks, yticks, and/or zticks,
%                               respectively. Default is 'xy'.
%
%               xformat:        formatting for x tick label.  This can be
%                               either a formatting string (see sprintf) or
%                               the handle to a function that takes as
%                               input a scalar number and returns a
%                               character array.                        
%
%               yformat:        formatting for y tick label.  This can be
%                               either a formatting string (see sprintf) or
%                               the handle to a function that takes as
%                               input a scalar number and returns a
%                               character array.
%
%               zformat:        formatting for z tick label.  This can be
%                               either a formatting string (see sprintf) or
%                               the handle to a function that takes as
%                               input a scalar number and returns a
%                               character array.
%
%               ytickoffset:    distance that y labels are set back off the
%                               y-axis, relative to width of axis.  Default
%                               is 0.04 and works well for a single 2D axis
%                               in the default-sized figure.  Depending on
%                               the size of the axis, length of new tick
%                               labels, and presence of axis labels, this
%                               may need to be adjusted.
%
%               xtickoffset:    distance that x labels are set back off the
%                               y-axis, relative to width of axis.  Default
%                               is 0.04 and works well for a single 2D axis
%                               in the default-sized figure.  Depending on
%                               the size of the axis, length of new tick
%                               labels, and presence of axis labels, this
%                               may need to be adjusted. 
%
%               ztickoffset:    distance that z labels are set back off the
%                               back plane (x = ymax plane), relative to
%                               y-span of axis.  Default is 0.08 and works
%                               well for a single 3D axis in the
%                               default-sized figure.  Depending on the
%                               size of the axis, length of new tick
%                               labels, and presence of axis labels, this
%                               may need to be adjusted.
%
%               panon:          logical scalar.  If true (default), set the
%                               customized pan function to update tick
%                               labels when panning.  Because custom
%                               panning can only be set for a figure, not
%                               an individual axis, you may want to disable
%                               this if you are applying customized ticks
%                               to multiple axes in a figure.
%
% Example:
%   
%   This example creates two axes.  The top axis shows Matlab's default
%   tick labels.  The tick2text function is applied to the bottom axis,
%   changing the x axis to show ticks as multiples of pi, and the y axis to
%   label with full fixed-point values.  To demonstrate how to retrieve the
%   new tick handles, the x labels are then rotated 340 degrees.
%
%     x = linspace(0,2*pi);
%     y = sin(x) + 100000;
% 
%     ax(1) = subplot(2,1,1);
%     plot(x,y);
%     title('Original tick marks');
% 
%     ax(2) = subplot(2,1,2);
%     plot(x,y);
%     title('Modified tick marks');
% 
%     set(ax, 'xlim', [0 2*pi]);
% 
%     tick2text(ax(2), 'yformat', '%.2f', ...
%                      'xformat', @(x) sprintf('%.2g\\pi', x/pi), ...
%                      'ytickoffset', 0.07, ...
%                      'xtickoffset', 0.06);
% 
%     hx = getappdata(ax(2), 'XTickText');
%     set(hx, 'Rotation', 340, 'fontsize', 12);

% Copyright 2009 Kelly Kearney

%----------------------------
% Parse and check input
%----------------------------

% Get axis handle

if ishandle(varargin{1})
    ax = varargin{1};
    pv = varargin(2:end);
else
    ax = gca;
    pv = varargin;
end

if ~isscalar(ax) || ~ishandle(ax) || ~strcmp(get(ax, 'type'), 'axes')
    error('Handle must be to a single axis');
end

% Get parent figure (for panning)

fig = get(ax, 'parent');
while ~strcmp(get(fig, 'type'), 'figure')
    fig = get(fig, 'parent');
end

% Parse and check optional variables

Options = struct('xformat', [], 'yformat', [], 'zformat', [], ...
                 'ytickoffset', .04, 'xtickoffset', .04, ...
                 'ztickoffset', .08, 'axis', 'xy', 'panon', true);
             
Options = parse_pv_pairs(Options, pv);

if ~isscalar(Options.ytickoffset) || ~isscalar(Options.xtickoffset) || ~isscalar(Options.ztickoffset)
    error('Offset values must be numerical scalars');
end

if ~ischar(Options.axis) || ~isempty(regexp(Options.axis, '[^xyz]'))
    error('Axis must be string including x, y, and/or z');
end

% Set formatting functions

if isempty(Options.xformat)
    xformatfun = @num2str;
elseif ischar(Options.xformat)
    xformatfun = @(x) num2str(x, Options.xformat);
elseif strcmp(class(Options.xformat), 'function_handle')
    xformatfun = Options.xformat;
else 
    error('xformat must be a formatting string or function handle');
end

if isempty(Options.yformat)
    yformatfun = @num2str;
elseif ischar(Options.yformat)
    yformatfun = @(x) num2str(x, Options.yformat);
elseif strcmp(class(Options.yformat), 'function_handle')
    yformatfun = Options.yformat;
else 
    error('yformat must be a formatting string or function handle');
end

if isempty(Options.zformat)
    zformatfun = @num2str;
elseif ischar(Options.zformat)
    zformatfun = @(x) num2str(x, Options.zformat);
elseif strcmp(class(Options.zformat), 'function_handle')
    zformatfun = Options.zformat;
else 
    error('zformat must be a formatting string or function handle');
end

%----------------------------
% Change exisiting tick 
% labels
%----------------------------

texttick(ax, xformatfun, yformatfun, zformatfun, Options.ytickoffset, Options.xtickoffset, Options.ztickoffset, Options.axis);

% Add customized zoom
    
try
    hzoom = zoom(ax);
    set(hzoom, 'ActionPostCallback', {@newzoompan, ax, xformatfun, yformatfun, zformatfun, Options.ytickoffset, Options.xtickoffset, Options.ztickoffset, Options.axis});
catch
    warning('Unable to apply customized zoom');
end
    
% Add customized pan

if Options.panon
    hpan = pan(fig);
    set(hpan, 'ActionPostCallback', {@newzoompan, ax, xformatfun, yformatfun, zformatfun, Options.ytickoffset, Options.xtickoffset, Options.ztickoffset, Options.axis});
end

%----------------------------
% Update tick labels when
% zooming and panning
%----------------------------

function newzoompan(h, ed, varargin)
texttick(varargin{:})

%----------------------------
% Change tick labels to text 
% objects
%----------------------------

function texttick(haxis, xformat, yformat, zformat, ytickoffset, xtickoffset, ztickoffset, tickax)

% Get current values of ticks and axis limits

temp = get(haxis, {'XTick', 'YTick', 'ZTick', 'XLim', 'YLim', 'ZLim', 'XScale', 'YScale', 'ZScale'});
[xtick, ytick, ztick, xlim, ylim, zlim, xscale, yscale, zscale] = deal(temp{:});

% Create new labels

xticklab = arrayfun(xformat, xtick, 'uni', 0);
yticklab = arrayfun(yformat, ytick, 'uni', 0);
zticklab = arrayfun(zformat, ztick, 'uni', 0);

% Calculate positions of new tick labels

if strcmp(xscale, 'linear')
    ytickoffset = ytickoffset * diff(xlim);
else
    ytickoffset = 10^(ytickoffset .* diff(log10(xlim)));
end
  
if strcmp(yscale, 'linear')
    xtickoffset = xtickoffset * diff(ylim);
else
    xtickoffset = 10^(xtickoffset .* diff(log10(ylim)));
end

if strcmp(yscale, 'linear')
    ztickoffset = ztickoffset * diff(ylim);
else
    ztickoffset = 10^(ztickoffset .* diff(log10(ylim)));
end

xlabely = repmat(ylim(1) - xtickoffset, size(xtick));
xlabelz = repmat(zlim(1), size(xtick));
ylabelx = repmat(xlim(1) - ytickoffset, size(ytick));
ylabelz = repmat(zlim(1), size(ytick));
zlabelx = repmat(xlim(1), size(ztick));
zlabely = repmat(ylim(2) + ztickoffset, size(ztick));

% Create or update x ticks

if ~isempty(strfind(tickax, 'x'))
    set(haxis, 'XTickLabel', []);
    if isappdata(haxis, 'XTickText')

        hxold = getappdata(haxis, 'XTickText');

        % Get properties of current ticks

        props = get(hxold);
        params = fieldnames(props(1));
        vals = struct2cell(props(1));

        % Ignore properties that always change between ticks or are read only

        [tf,loc] = ismember({'Extent', 'String', 'BeingDeleted', 'Children', 'Position', 'Type'}, params);
        params(loc) = [];
        vals(loc) = [];
        paramval = [params'; vals'];

        % Create new ticks with same properties

        delete(hxold);
        hx = text(xtick, xlabely, xlabelz, xticklab, paramval{:}, 'parent', haxis);      
        setappdata(haxis, 'XTickText', hx);

    else
        hx = text(xtick, xlabely, xlabelz, xticklab, 'horizontalalignment', 'center', 'parent', haxis);      
        setappdata(haxis, 'XTickText', hx);
    end
end
   
% Create or update y ticks

if ~isempty(strfind(tickax, 'y'))
    set(haxis, 'YTickLabel', []);
    if isappdata(haxis, 'YTickText')
        hyold = getappdata(haxis, 'YTickText');

        % Get properties of current ticks

        props = get(hyold);
        params = fieldnames(props(1));
        vals = struct2cell(props(1));

        % Ignore properties that always change between ticks or are read only

        [tf,loc] = ismember({'Extent', 'String', 'BeingDeleted', 'Children', 'Position', 'Type'}, params);
        params(loc) = [];
        vals(loc) = [];
        paramval = [params'; vals'];

        % Create new ticks with same properties

        delete(hyold);
        hy = text(ylabelx, ytick, ylabelz, yticklab, paramval{:}, 'parent', haxis);      
        setappdata(haxis, 'YTickText', hy);

    else
        hy = text(ylabelx, ytick, ylabelz, yticklab, 'horizontalalignment', 'center', 'parent', haxis);
        setappdata(haxis, 'YTickText', hy);
    end
end

% Create or update z ticks

if ~isempty(strfind(tickax, 'z'))
    set(haxis, 'ZTickLabel', []);
    if isappdata(haxis, 'ZTickText')
        hzold = getappdata(haxis, 'ZTickText');

        % Get properties of current ticks

        props = get(hzold);
        params = fieldnames(props(1));
        vals = struct2cell(props(1));

        % Ignore properties that always change between ticks or are read only

        [tf,loc] = ismember({'Extent', 'String', 'BeingDeleted', 'Children', 'Position', 'Type'}, params);
        params(loc) = [];
        vals(loc) = [];
        paramval = [params'; vals'];

        % Create new ticks with same properties

        delete(hzold);
        hz = text(zlabelx, zlabely, ztick, zticklab, paramval{:}, 'parent', haxis);      
        setappdata(haxis, 'ZTickText', hz);

    else
        hz = text(zlabelx, zlabely, ztick, zticklab, 'horizontalalignment', 'center', 'parent', haxis);
        setappdata(haxis, 'ZTickText', hz);
    end
end

%----------------------------
% Parse optional input
%----------------------------

% This subfunction is copied directly from John D'Errico's
% parse_pv_pairs.m (FEX file ID #9082).

function params=parse_pv_pairs(params,pv_pairs)

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
  error 'Property/value pairs must come in PAIRS.'
end
if n<=0
  % just return the defaults
  return
end

if ~isstruct(params)
  error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
  p_i = lower(pv_pairs{2*i-1});
  v_i = pv_pairs{2*i};
  
  ind = strmatch(p_i,lpropnames,'exact');
  if isempty(ind)
    ind = find(strncmp(p_i,lpropnames,length(p_i)));
    if isempty(ind)
      error(['No matching property found for: ',pv_pairs{2*i-1}])
    elseif length(ind)>1
      error(['Ambiguous property name: ',pv_pairs{2*i-1}])
    end
  end
  p_i = propnames{ind};
  
  % override the corresponding default in params
  params = setfield(params,p_i,v_i);
  
end