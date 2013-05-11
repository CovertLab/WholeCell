classdef PrintUtil
    %print a cell array to various places
    methods (Static)                    
        function printToStdIO(content, colLabels, options)
            if exist('options', 'var') && isstruct(options) && isfield(options, 'indentation')
                content = edu.stanford.covert.cell.sim.util.PrintUtil.applyIndentation(content, options.indentation);
            end
            
            %write file
            colWidths = edu.stanford.covert.cell.sim.util.PrintUtil.calcColWidths(content, colLabels);
            colTerminators = [repmat({'\t'}, 1, numel(colLabels)-1) '\n'];
            
            for j = 1:numel(colLabels)
                fprintf(['%-' num2str(colWidths(j)) 's' colTerminators{j}], colLabels{j});
            end
            for j = 1:numel(colLabels)
                fprintf(['%-' num2str(colWidths(j)) 's' colTerminators{j}], repmat('=', 1, colWidths(j)));
            end
            
            for i = 1:size(content, 1)
                for j = 1:numel(colLabels)
                    if isempty(content{i,j})
                        fprintf(['%-' num2str(colWidths(j)) 's' colTerminators{j}], content{i,j});
                    elseif isnumeric(content{i,j})
                        if ceil(content{i,j}) == content{i,j}
                            fprintf(['%' num2str(colWidths(j)) 'd' colTerminators{j}], content{i,j});
                        else
                            if ispc
                                fprintf(['%.' num2str(colWidths(j)-2-5) 'e' colTerminators{j}], content{i,j});
                            else
                                fprintf(['%.' num2str(colWidths(j)-2-4) 'e' colTerminators{j}], content{i,j});
                            end
                        end
                    elseif islogical(content{i,j})
                        fprintf(['%' num2str(colWidths(j)) 'd' colTerminators{j}], content{i,j});
                    elseif ischar(content{i,j})
                        fprintf(['%-' num2str(colWidths(j)) 's' colTerminators{j}], content{i,j});
                    else
                        throw(MException('PrintUtil:error', 'unsupported type %s', class(content{i,j})));
                    end
                end
            end
            
            fprintf('\n');
        end
        
        function printToFile(content, colLabels, fileName, varargin)
            %open file
            if ismember(fileName(find(fileName == '.', 1, 'last')+1 : end), {'xls','xlsx'})
                if ispc
                    edu.stanford.covert.cell.sim.util.PrintUtil.printToExcel(content, colLabels, fileName, varargin{:});
                    return;
                else
                    warning('WholeCell:warning', 'Excel not available, saving in csv format');
                end
            end
            fid = fopen(fileName, 'w');
            
            %write file
            colWidths = edu.stanford.covert.cell.sim.util.PrintUtil.calcColWidths(content, colLabels);
            colTerminators = [repmat({'\t'}, 1, numel(colLabels)-1) '\n'];
            
            for j = 1:numel(colLabels)
                fprintf(fid, ['%-' num2str(colWidths(j)) 's' colTerminators{j}], colLabels{j});
            end
            
            for i = 1:size(content, 1)
                for j = 1:numel(colLabels)
                    if isempty(content{i,j})
                        fprintf(fid, ['%-' num2str(colWidths(j)) 's' colTerminators{j}], content{i,j});
                    elseif isnumeric(content{i,j})
                        if ispc
                            fprintf(fid, ['%.' num2str(colWidths(j)-2-5) 'e' colTerminators{j}], content{i,j});
                        else
                            fprintf(fid, ['%.' num2str(colWidths(j)-2-4) 'e' colTerminators{j}], content{i,j});
                        end                        
                    elseif islogical(content{i,j})
                        fprintf(fid, ['%' num2str(colWidths(j)) 'd' colTerminators{j}], content{i,j});
                    elseif ischar(content{i,j})
                        fprintf(fid, ['%-' num2str(colWidths(j)) 's' colTerminators{j}], content{i,j});
                    else
                        throw(MException('PrintUtil:error', 'unsupported type %s', class(content{i,j})));
                    end
                end
            end
            
            %close file
            fclose(fid);
        end
        
        function printToExcel(content, colLabels, fileName, sheet, options)
            if exist('options', 'var') && isstruct(options) && isfield(options, 'indentation')
                content = edu.stanford.covert.cell.sim.util.PrintUtil.applyIndentation(content, options.indentation);
            end
            if exist('sheet', 'var')
                xlswrite(fileName, [colLabels; content], sheet);
            else
                xlswrite(fileName, [colLabels; content]);
            end
        end
        
        function colWidths = calcColWidths(content, colLabels)
            elWidths = zeros(size(content));
            for i = 1:size(content, 1)
                for j = 1:size(content, 2)
                    if isnumeric(content{i,j})
                        elWidths(i,j) = max(9, length(num2str(content{i,j})));
                    elseif islogical(content{i,j})
                        elWidths(i,j) = 1;
                    elseif ischar(content{i,j})
                        elWidths(i,j) = length(content{i,j});
                    else
                        throw(MException('PrintUtil:error', 'unsupported type %s', class(content{i,j})));
                    end
                end
            end            
            colWidths = max([elWidths; cellfun(@length, colLabels)], [], 1);
        end
        
        function content = applyIndentation(content, indentation)
            content(:, 1) = cellfun(@(content, indentation) [repmat('    ', 1, indentation) content], content(:,1), num2cell(indentation), 'UniformOutput', false);
        end
    end
end
