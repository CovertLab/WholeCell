classdef PlotUtil
    methods (Static = true)
        function compartmentData = selectCompartmentsForPlot(sim, compartments, data, sortOrder)
            sz = size(data);
            sz(2) = length(compartments);
            sz(end+1) = 1;
            compartmentData = zeros(sz);
            
            for j=1:length(compartments)
                switch compartments(j)
                    case -3,   compartmentData(:,j,:) = sum(data, 2);
                    case -2,   compartmentData(:,j,:) = mean(data, 2);
                    case -1,   compartmentData(:,j,:) = median(data, 2);
                    case 0,    compartmentData(:,j,:) = std(data, 2);
                    otherwise, compartmentData(:,j,:) = data(:,compartments(j),:);
                end
            end
            
            if(~exist('sortOrder','var')); sortOrder='Compartment, Property'; end;
            switch sortOrder
                case 'None',
                case 'Property, Compartment', compartmentData = reshape(permute(compartmentData, [2 1 3:length(sz)]), [sz(1)*sz(2) sz(3:end)]);
                case 'Compartment, Property', compartmentData = reshape(compartmentData, [sz(1)*sz(2) sz(3:end)])';
            end
        end
        
        function labelCompartmentPlot(sim, plotHandles, compartments, varargin)
            if length(plotHandles) == 1 || (length(plotHandles) > length(compartments) && nargin == 3); return, end;
            if nargin == 4 && length(plotHandles) == length(varargin{1})
                compartmentLabels=varargin{1};
            else
                compartmentLabels = edu.stanford.covert.cell.sim.util.PlotUtil.compartmentPlotLabels(sim, compartments, varargin{:});
                compartmentLabels = reshape(compartmentLabels, 1, []);
            end
            legend(plotHandles, compartmentLabels);
        end
        
        function compartmentLabels = compartmentPlotLabels(sim, compartments, labels)
            if exist('labels', 'var')
                numProperties = length(labels);
            else
                numProperties = 1;
            end
            compartmentLabels = cell(numProperties, length(compartments));
            
            for j = 1:length(compartments)
                switch compartments(j)
                    case -3,   compartmentLabel = 'Sum';
                    case -2,   compartmentLabel = 'Mean';
                    case -1,   compartmentLabel = 'Median';
                    case 0,    compartmentLabel = 'Std Dev';
                    otherwise, compartmentLabel = sim.compartment.names{compartments(j)};
                end
                
                for i = 1:numProperties
                    if exist('labels', 'var')
                        compartmentLabels{i,j} = [compartmentLabel ', ' labels{i}];
                    else
                        compartmentLabels{i,j} = compartmentLabel;
                    end
                end
            end
        end
        
        function [axesHandle, figHandle] = newAxesHandle(paperSize)
            figHandle = figure();

            if nargin >= 1
                set(figHandle, 'PaperUnits', 'centimeters');
                set(figHandle, 'PaperSize', paperSize);
                set(figHandle, 'PaperPositionMode', 'manual');
                set(figHandle, 'PaperPosition', [0 0 paperSize]);
                set(figHandle, 'Units', 'normalized');
            end
            
            axesHandle = subplot(1, 1, 1, 'Parent', figHandle);
        end
        
        function [axesHandles, xAxesHandles, titleHandles, offSetAxesHandles] = multiElementPlot(figHandle, axesSizes, xSpan, options)            
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if nargin < 4
                options = struct;
            end
            
            if ~iscell(axesSizes)
                [axesHandles, xAxesHandles, titleHandles, offSetAxesHandles] = PlotUtil.multiElementPlotHelper(figHandle, axesSizes, xSpan, options);
            else
                axesHandles = cell(size(axesSizes));
                xAxesHandles = cell(size(axesSizes));
                titleHandles = cell(size(axesSizes));
                offSetAxesHandles = cell(size(axesSizes));
                if ~iscell(xSpan)
                    xSpan = repmat({xSpan}, size(axesSizes));
                end
            
                if ~isfield(options, 'position');
                    options.position = [0.02 0.12 0.95 0.8];
                end
                
                H = 0;
                for i = 1:numel(axesSizes)
                    H = max(H, sum(axesSizes{i}) + numel(axesSizes{i}) - 1);
                end
                
                if isfield(options, 'colWidths')
                    colWidths = options.colWidths;
                else
                    colWidths = ones(size(axesSizes));
                end
                if isfield(options, 'yAxesLabelWidths')
                    yAxesLabelWidths = options.yAxesLabelWidths;
                else
                    yAxesLabelWidths = 0.075 * ones(size(axesSizes));
                end
                
                for i = 1:numel(axesSizes)
                    axesSizes{i} = axesSizes{i} * (H - numel(axesSizes{i}) + 1) / sum(axesSizes{i});
                    
                    tmpOptions = struct;
                    
                    if isfield(options, 'titleStr')
                        if iscell(options.titleStr)
                            tmpOptions.titleStr = options.titleStr{i};
                        elseif i == 1;
                            tmpOptions.titleStr = options.titleStr;
                        end
                    end
                    
                    if isfield(options, 'position') && iscell(options.position)
                        tmpOptions.position = options.position{i};
                    else
                        x0 = options.position(1) + sum(colWidths(1:i-1)) * options.position(3) / sum(colWidths) + yAxesLabelWidths(i);
                        w = colWidths(i) * options.position(3) / sum(colWidths) - yAxesLabelWidths(i);
                        y0 = options.position(2);
                        h = options.position(4);
                        tmpOptions.position = [x0 y0 w h];
                    end
                    
                    if isfield(options, 'colorOrder')
                        if iscell(options.colorOrder)
                            tmpOptions.colorOrder = options.colorOrder{i};
                        else
                            tmpOptions.colorOrder = options.colorOrder;
                        end
                    end
                    
                    if isfield(options, 'xlabelStr')
                        if iscell(options.xlabelStr)
                            tmpOptions.xlabelStr = options.xlabelStr{i};
                        else
                            tmpOptions.xlabelStr = options.xlabelStr;
                        end
                    end
                    
                    if isfield(options, 'xdata')
                        if iscell(options.xdata)
                            tmpOptions.xdata = options.xdata{i};
                        else
                            tmpOptions.xdata = options.xdata;
                        end
                    end
                    
                    if isfield(options, 'ydata')
                        tmpOptions.ydata = options.ydata{i};
                    end
                    
                    if isfield(options, 'ylabelStr')
                        if iscell(options.ylabelStr)
                            tmpOptions.ylabelStr = options.ylabelStr{i};
                        else
                            tmpOptions.ylabelStr = options.ylabelStr;
                        end
                    end
                    
                    if isfield(options, 'plotFunc')
                        if iscell(options.plotFunc)
                            tmpOptions.plotFunc = options.plotFunc{i};
                        else
                            tmpOptions.plotFunc = options.plotFunc;
                        end
                    end
                    
                    [axesHandles{i}, xAxesHandles{i}, titleHandles{i}, offSetAxesHandles{i}] = PlotUtil.multiElementPlotHelper(figHandle, axesSizes{i}, xSpan{i}, tmpOptions);
                end
            end
        end
        
        function [axesHandles, xAxesHandle, titleHandle, offSetAxesHandles] = multiElementPlotHelper(figHandle, axesSizes, xSpan, options)            
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if exist('options', 'var')
                if isfield(options, 'titleStr')
                    titleStr = options.titleStr;
                end
                if isfield(options, 'position')
                    position = options.position;
                end
                if isfield(options, 'colorOrder')
                    colorOrder = options.colorOrder;
                end
                if isfield(options, 'xlabelStr')
                    xlabelStr = options.xlabelStr;
                end
            end
            
            if ~exist('position', 'var')
                position = [0.13 0.12 0.8 0.8];
            end
            if ~exist('xlabelStr', 'var')
                xlabelStr = 'Time (h)';
            end
            
            %value plots
            nPlots = numel(axesSizes);
            axesHandles = zeros(nPlots, 1);
            axesX = position(1);
            axesW = position(3);
            axesY0 = position(2);
            axesTotH = position(4);
            a = 0.2;
            b = 0.2;
            c = 0.2;
            for i = 1:nPlots
                axesHandle = axes('Parent', figHandle);
                axesHandles(i) = axesHandle;
                
                if i == 1 && exist('titleStr', 'var')
                    titleHandle = title(axesHandle, titleStr, 'FontSize', 10);
                else
                    titleHandle = NaN;
                end
                
                set(axesHandle, ...
                    'FontSize', 7, ...
                    'TickDir', 'out', ...
                    'YMinorTick', 'off', ...
                    'XTick', [], ...
                    'XColor', get(figHandle, 'Color'));
                box(axesHandle, 'off');
                ylabel(axesHandle, '', ...
                    'FontSize', 8, ...
                    'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'center');
                xlim(axesHandle, xSpan);
                hold(axesHandle, 'on');
                
                if exist('colorOrder', 'var')
                    set(axesHandle, 'ColorOrder', colorOrder);
                end
                
                %layout
                axesH = axesTotH * a * axesSizes(i) / (a * sum(axesSizes) + b * (nPlots - 1) + c);
                axesY = axesY0 + axesTotH - axesTotH * (a * sum(axesSizes(1:i)) + b * (i - 1)) / ...
                    (a * sum(axesSizes) + b * (nPlots - 1) + c);
                set(axesHandle, 'Position', [axesX axesY axesW axesH]);
            end
            
            %x axis
            xAxesHandle = axes('Parent', figHandle);
            plot(xAxesHandle, xSpan(:)', [0 0], 'LineStyle', 'none', 'Marker', 'none');
            if range(xSpan) >= 10
                n = 10^floor(log10(range(xSpan)));
                set(xAxesHandle, 'XTick', unique(n * (ceil(xSpan(1) / n):floor(xSpan(2) / n))));
            elseif range(xSpan) >= 5
                n = 5;
                set(xAxesHandle, 'XTick', unique(n * (ceil(xSpan(1) / n):floor(xSpan(2) / n))));
            end
            
            set(xAxesHandle, ...
                'FontSize', 7, ...
                'TickDir', 'out', ...
                'YTick', [], ...
                'XMinorTick', 'off', ...
                'YColor', get(figHandle, 'Color'));
            box(xAxesHandle, 'off');
            xlabel(xAxesHandle, xlabelStr, ...
                'FontSize', 8, ...
                'VerticalAlignment', 'top', ...
                'HorizontalAlignment', 'center');
            xlim(xAxesHandle, xSpan);
            axesH = 1e-6;
            axesY = axesY0;
            set(xAxesHandle, 'Position', [axesX axesY axesW axesH]);
            
            %y axis labels
            if isfield(options, 'ylabelStr')
                if ~iscell(options.ylabelStr)
                    options.ylabelStr = {options.ylabelStr};
                end
                for i = 1:numel(options.ylabelStr)
                    ylabel(axesHandles(i), options.ylabelStr{i});
                end
            end
            
            %plot data
            if isfield(options, 'ydata')
                if ~iscell(options.ydata)
                    options.ydata = {options.ydata};
                end
                if ~isfield(options, 'plotFunc')
                    options.plotFunc = repmat({@PlotUtil.plotLine}, size(options.ydata));
                elseif ~iscell(options.plotFunc)
                    options.plotFunc = repmat({options.plotFunc}, size(options.ydata));
                end
                
                for i = 1:numel(options.ydata)
                    if ~isempty(options.ydata{i})
                        options.plotFunc{i}(axesHandles(i), options.xdata, options.ydata{i}, false, true, false);
                    end
                end
                
                PlotUtil.alignYAxesLabels(axesHandles);
                if ~isfield(options, 'offsetYAxes')
                    options.offsetYAxes = 0.03;
                end
                if options.offsetYAxes
                    offSetAxesHandles = PlotUtil.offsetYAxes(axesHandles, options.offsetYAxes, xAxesHandle);
                else
                    offSetAxesHandles = zeros(size(axesHandles));
                end
            else
                offSetAxesHandles = zeros(size(axesHandles));
            end
        end
        
        function handle = plotLine(axesHandle, x, y, setXLim, setYLim, buttonDownCallback)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            handle = plot(axesHandle, x, y);
            
            if nargin < 6 || buttonDownCallback
                for i = 1:length(handle)
                    s = sprintf('fprintf(''%d\\n'')', i);
                    set(handle(i), 'ButtonDownFcn', s);
                end
            end
            
            if nargin < 4 || setXLim
                xlim([min(x) max(x)]);
            end
            
            if nargin < 5 || setYLim
                minY = min(y(:));
                maxY = max(y(:));
                PlotUtil.formatYAxes(axesHandle, minY, maxY);
            end
        end
        
        function handle = plotScatter(axesHandle, x, y, msize, color, mstyle)
            handle = scatter(axesHandle, x, y, msize, color, mstyle);
%             handle = get(handle, 'Children');
%             for i = 1:length(handle)
%                 s = sprintf('fprintf(''%d\\n'')', i);
%                 set(handle(i), 'ButtonDownFcn', s);
%             end
  
            if range(x)
                xlim([min(x) max(x)]);
            end
            
            minY = min(y(:));
            maxY = max(y(:));
            
            if maxY - minY > eps
                ylim(axesHandle, [minY maxY]);
                set(axesHandle, 'YTick', [minY mean([minY maxY]) maxY]);
            else
                set(axesHandle, 'YTick', y(1));
            end
        end
        
        function handle = plotBoundedLineFromRawData(axesHandle, x, y)
            quantiles = zeros(2, length(x));
            for i = 1:length(x)
                tfs = ~isnan(y(:, i));
                quantiles(:, i) = prctile(y(tfs, i), [25 75], 1);
            end
            
            qmean = mean(quantiles,1);
            b = abs(quantiles-repmat(qmean, 2, 1))';
            handle = boundedline(x, qmean, ...
                b, 'cmap', winter(4), 'transparency', 0.5);
            set(handle(1),'Visible','Off');
            hold on;
            handle = plot(x,[nanmean(y, 1);quantiles]);
            set(handle(1),'Color',[0 0 0],'LineWidth',2);
            for i = 2:3
                set(handle(i), 'Color', [0.0 0.0 0.0], 'LineStyle', '--');
            end
            hold off;

            xlim([min(x) max(x)]);

            minY = min(quantiles(1,:));
            maxY = max(quantiles(2,:));

            if(maxY-minY>eps)
                ylim(axesHandle, [minY maxY]);
                set(axesHandle, 'YTick', [minY mean([minY maxY]) maxY]);
            else
                set(axesHandle, 'YTick', y(1));
            end
        end
        
        function handles = plotColoredTicks(axesHandle, ticLocs, colorOrder)
            % Loop backwards so that if two simulations end at nearly the
            % same time, you'll see the color of the first one
            hold(axesHandle, 'on');
            handles = zeros(length(ticLocs), 1);
            for i = length(ticLocs):-1:1
                % Only plot if not a NaN
                if ~isnan(ticLocs(i))
                    handles(i) = line([ticLocs(i) ticLocs(i)], [0 1], ...
                        'Parent', axesHandle, ...
                        'Color', colorOrder(i, :));
                end
            end
        end
        
        function handle = plotStackedArea(axesHandle, xdata, ydata, setXLim, setYLim, buttonDownCallback)
            if size(xdata, 2) == size(ydata, 2)
                xdata = xdata';
                ydata = ydata';
            end
            handle = area(axesHandle, xdata, ydata, 'EdgeColor', 'none');
        end
        
        function handle = plotNormalizedStackedArea(axesHandle, xdata, ydata, varargin)
            if size(xdata, 2) == size(ydata, 2)
                xdata = xdata';
                ydata = ydata';
            end
            ydata = ydata ./ repmat(sum(ydata, 2) / 100, 1, size(ydata, 2));
            ydata(isnan(ydata)) = 0;
            handle = edu.stanford.covert.cell.sim.util.PlotUtil.plotStackedArea(axesHandle, xdata, ydata, varargin{:});
        end
        
        function handle = plotBubbles(axesHandle, xdata, ydata, varargin)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if size(xdata, 2) == 1
                xdata = [xdata ones(size(xdata))];
                ydata = [ydata ones(size(ydata))];
            end
            
            n = 360;
            theta = (1:n)/n * 2 * pi;
            x = xdata(:, 1 * ones(n, 1)) + xdata(:, 2 * ones(n, 1)) .* cos(theta(ones(size(xdata, 1), 1), :));
            y = ydata(:, 1 * ones(n, 1)) + ydata(:, 2 * ones(n, 1)) .* sin(theta(ones(size(xdata, 1), 1), :));
            colors = (1:size(xdata, 1))';
            handle = fill(x', y', colors', ...
                'EdgeColor', 'k', ...
                'Parent', axesHandle);           
            
            colormap(PlotUtil.getRedGreenColorOrder(colors));
        end
        
        function handle = plotXAxis(axesHandle, x, xLabel)
            handle = edu.stanford.covert.cell.sim.util.PlotUtil.plotLine(axesHandle, x, ones(1,length(x)));
            
            set(handle, 'Color', 'w');
            ycolor = get(axesHandle, 'YColor');
            set(axesHandle, 'XColor', ycolor);
            set(axesHandle, 'YColor', get(get(axesHandle,'Parent'), 'Color'));
            set(axesHandle, 'YTick', []);
            box(axesHandle, 'off');
            set(axesHandle, 'TickDir', 'out');
            xlabel(axesHandle, xLabel, 'fontSize', 12);
        end
        
        function labelSubplots(axesHandles, labels, x, y)
            if nargin >= 2 && ~iscell(labels)
                if numel(axesHandles) > 1
                    labels = num2cell(labels);
                else
                    labels = {labels};
                end
            elseif nargin < 2
                labels = num2cell(char(65 + (1:numel(axesHandles)) - 1));
            end
            if nargin < 3
                x = -0.025;
            end
            if nargin < 4
                y = 1.10;
            end
            
            for i = 1:numel(axesHandles)
                text(x, y, labels{i}, ...
                    'FontSize', 10, ...
                    'Parent', axesHandles(i), ...
                    'Units', 'normalized', ...
                    'HorizontalAlignment', 'right', ...
                    'VerticalAlignment', 'bottom');
            end
        end
        
        function formatYAxes(axesHandle, minY, maxY, fontSize)
            if nargin < 4
                fontSize = 7;
            end
            
            if maxY - minY > eps
                ylim(axesHandle, [minY maxY]);
                ytick = [minY maxY];
            else
                ytick = minY;
            end
            set(axesHandle, 'YTick', ytick);
            if all(ceil(ytick) == ytick) && all(ytick < 1e7)
                set(axesHandle, 'YTickLabel', cellfun(@num2str, num2cell(ytick), 'UniformOutput', false));
            else
                tick2text(axesHandle, 'axis', 'y', 'ytickoffset', 0.06)
                h = getappdata(axesHandle, 'YTickText');
                set(h, 'Interpreter', 'tex', 'FontSize', fontSize, 'HorizontalAlignment', 'right')
                for i = 1:numel(ytick)
                    if ytick(i) == 0 || all(ytick == ceil(ytick) & ytick <= 1e3)
                        str = sprintf('%d', ytick(i));
                    elseif abs(ytick(i)) < 1e-2 || abs(ytick(i)) > 1e3
                        tmp = sprintf('%.0e', ytick(i));
                        str = sprintf('%s{\\times}10^{%d}', tmp(1:find(tmp == 'e')-1), str2double(tmp(find(tmp=='e')+1:end)));
                    elseif abs(ytick(i)) < 1e3 && abs(ytick(i)) >= 1e-2
                        str = sprintf(['%1.' num2str(ceil(-log10(diff(ytick)))+1) 'f'], ytick(i));
                    else
                        tmp = sprintf('%.1e', ytick(i));
                        str = sprintf('%s{\\times}10^{%d}', tmp(1:find(tmp == 'e')-1), str2double(tmp(find(tmp=='e')+1:end)));
                    end
                    set(h(i), 'String', str);
                end
            end
        end
        
        function axesHandles2 = offsetYAxes(axesHandles, offset, varargin)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~iscell(axesHandles)
                axesHandles2 = PlotUtil.offsetYAxesHelper(axesHandles, offset, varargin{:});
            else
                axesHandles2 = cell(size(axesHandles));
                for i = 1:numel(axesHandles)
                    axesHandles2{i} = PlotUtil.offsetYAxesHelper(axesHandles{i}, offset, varargin{:});
                end
            end
        end
        
        function axesHandles2 = offsetYAxesHelper(axesHandles, offset, xAxesHandle)
            axesHandles2 = zeros(size(axesHandles));
            for i = 1:numel(axesHandles)
                axesHandle = axesHandles(i);
                axesHandle2 = axes;
                axesHandles2(i) = axesHandle2;
                props = get(axesHandle);
                fields = fieldnames(props);
                for j = 1:numel(fields)
                    try %#ok<TRYNC>
                        set(axesHandle2, fields{j}, props.(fields{j}));
                    end
                end
                set(axesHandle2, 'XLim', get(axesHandle, 'XLim'))
                set(axesHandle2, 'YLim', get(axesHandle, 'YLim'))
                     
                %copy properties
                pos = get(axesHandle, 'Position');
                tmp = offset * pos(3);
                if offset > 0
                    pos(1) = pos(1) - tmp;
                else
                    set(axesHandle, 'Position', [pos(1)-tmp pos(2) pos(3)+tmp pos(4)]);                    
                end
                set(axesHandle2, 'Position', pos)
                uistack(axesHandle2, 'bottom');
                
                %copy title properties
                props = get(get(axesHandle, 'Title'));
                fields = fieldnames(props);
                for j = 1:numel(fields)
                    try %#ok<TRYNC>
                        set(get(axesHandle2, 'Title'), fields{j}, props.(fields{j}));
                    end
                end
                pos = get(get(axesHandle, 'Title'), 'Position');
                pos(1) = pos(1) + tmp;
                set(get(axesHandle2, 'Title'), 'Position', pos)
                
                %copy xlabel properties
                props = get(get(axesHandle, 'XLabel'));
                fields = fieldnames(props);
                for j = 1:numel(fields)
                    try %#ok<TRYNC>
                        set(get(axesHandle2, 'XLabel'), fields{j}, props.(fields{j}));
                    end
                end
                pos = get(get(axesHandle, 'XLabel'), 'Position');
                pos(1) = pos(1) + tmp;
                set(get(axesHandle2, 'XLabel'), 'Position', pos)
                
                %copy ylabel properties
                props = get(get(axesHandle, 'YLabel'));
                fields = fieldnames(props);
                for j = 1:numel(fields)
                    try %#ok<TRYNC>
                        set(get(axesHandle2, 'YLabel'), fields{j}, props.(fields{j}));
                    end
                end
                
                set(axesHandle, 'Visible', 'off');
            end
            
            if exist('xAxesHandle', 'var')                
                xpos = get(xAxesHandle, 'Position');
                tmp = offset * xpos(3);
                set(xAxesHandle, 'Position', [xpos(1)-tmp xpos(2) xpos(3)+tmp xpos(4)]);
            end
        end
        
        function alignYAxesLabels(axesHandles, varargin)
            import edu.stanford.covert.cell.sim.util.PlotUtil;
            
            if ~iscell(axesHandles)
                PlotUtil.alignYAxesLabelsHelper(axesHandles, varargin{:});
            else
                for i = 1:numel(axesHandles)
                    PlotUtil.alignYAxesLabelsHelper(axesHandles{i}, varargin{:});
                end
            end
        end
        
        function alignYAxesLabelsHelper(axesHandles, xpos)
            if nargin < 2
                xpos = Inf;
                for i = 1:numel(axesHandles)
                    pos = get(get(axesHandles(i), 'YLabel'), 'Position');
                    xpos = min(xpos, pos(1));
                end
            end
            
            for i = 1:numel(axesHandles)
                pos = get(get(axesHandles(i), 'YLabel'), 'Position');
                set(get(axesHandles(i), 'YLabel'), 'Position', [xpos pos(2:end)]);
            end
        end
        
        function colorOrder = getRedGreenColorOrder(lineOrder)
            nLines = numel(lineOrder);
            nLines1 = floor(nLines/2);
            nLines2 = ceil(nLines/2);
            
            colorOrder = [[(nLines1:-1:1)'/nLines1; zeros(nLines2, 1)] [zeros(nLines1, 1); (1:nLines2)'/nLines2] zeros(nLines, 1)];
            colorOrder(lineOrder, :) = colorOrder;
        end
        
        function handle = drawCircle(centerX, centerY, radius, axesHandle)
            if nargin < 3
                throw(MException('PlotUtil:drawCircle:error', 'Insufficient arguments'));
            end
            
            handle = rectangle('Position', [centerX-radius centerY-radius 2*radius 2*radius], 'Curvature', [1 1]);
            
            if nargin == 4 && ~isempty(axesHandle)
                set(handle, 'Parent', axesHandle)
            end
        end
        
        function h = labelSubFigure(string, position, figHandle, fontSize, units)
            if nargin < 4
                fontSize = 8;
            end
            if nargin < 5
                units = 'normalized';
            end
            
            h = annotation(figHandle, 'textbox', [0 0 0 0], ...
                'String', string, ...
                'FontWeight', 'bold', ...
                'FontSize', fontSize, ...
                'HorizontalAlignment', 'Left', ...
                'VerticalAlignment', 'middle', ...
                'EdgeColor', 'none', ...
                'FitBoxToText', 'on', ...
                'Margin', -25);
            set(h, 'Position', position, 'units', units);
        end
    end
end
