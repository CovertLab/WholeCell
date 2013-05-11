%Analysis of process runtimes
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/12/2011
classdef RunTime
    methods (Static)
        function run(fileName)
            profSummary = edu.stanford.covert.cell.sim.analysis.RunTime.calcProcessRunTimes();
            subtotals = ...
                [profSummary.evolveState]' + ...
                [profSummary.calcResourceRequirements_Current]' + ...
                [profSummary.copyToState]' + ...
                [profSummary.copyFromState]';
            total = sum(subtotals);
            
            fid = fopen(fileName, 'w');
            fprintf(fid, '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">\n');
            fprintf(fid, '<html>\n');
            fprintf(fid, '\t<head>\n');
            fprintf(fid, '\t\t<title>%s</title>\n', 'Profile Summary');
            fprintf(fid, '\t\t<style>\n');
            fprintf(fid, 'body{ font-family:arial; font-size: 10pt;}\n');
            fprintf(fid, 'th{text-align:left; font-weight:bold; background-color:#CCCCCC;}\n');
            fprintf(fid, 'td,th{padding-right:10px;}\n');
            fprintf(fid, 'td, tfoot th{text-align:right}\n');
            fprintf(fid, 'td:first-child, tfoot th:first-child{text-align:left;}\n');
            fprintf(fid, 'td:last-child,th:last-child{padding-right:0px;}\n');
            fprintf(fid, '\t\t</style>\n');
            fprintf(fid, '\t</head>\n');
            fprintf(fid, '\t<body>\n');
            fprintf(fid, '\t\t<table>\n');
            fprintf(fid, '\t\t\t<thead>\n');
            fprintf(fid, '\t\t\t\t<tr>\n');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Process');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Evolve State');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Substrate Requirements');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Copy from State');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Copy to State');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Percent');
            fprintf(fid, '\t\t\t\t</tr>\n');
            fprintf(fid, '\t\t\t</thead>\n');
            fprintf(fid, '\t\t\t<tbody>\n');
            for i = 1:numel(profSummary)
                fprintf(fid, '\t\t\t\t<tr>\n');
                fprintf(fid, '\t\t\t\t\t<td>%s</td>\n', profSummary(i).name);
                fprintf(fid, '\t\t\t\t\t<td>%.3f</td>\n', profSummary(i).evolveState);
                fprintf(fid, '\t\t\t\t\t<td>%.3f</td>\n', profSummary(i).calcResourceRequirements_Current);
                fprintf(fid, '\t\t\t\t\t<td>%.3f</td>\n', profSummary(i).copyFromState);
                fprintf(fid, '\t\t\t\t\t<td>%.3f</td>\n', profSummary(i).copyToState);
                fprintf(fid, '\t\t\t\t\t<td>%.1f</td>\n', subtotals(i)/total*100);
                fprintf(fid, '\t\t\t\t</tr>\n');
            end
            fprintf(fid, '\t\t\t</tbody>\n');
            fprintf(fid, '\t\t\t<tfoot>\n');
            fprintf(fid, '\t\t\t\t\t<th>%s</th>\n', 'Percent');
            fprintf(fid, '\t\t\t\t\t<th>%.1f</th>\n', sum([profSummary.evolveState])/total*100);
            fprintf(fid, '\t\t\t\t\t<th>%.1f</th>\n', sum([profSummary.calcResourceRequirements_Current])/total*100);
            fprintf(fid, '\t\t\t\t\t<th>%.1f</th>\n', sum([profSummary.copyFromState])/total*100);
            fprintf(fid, '\t\t\t\t\t<th>%.1f</th>\n', sum([profSummary.copyToState])/total*100);
            fprintf(fid, '\t\t\t\t\t<th>%.1f</th>\n', 100);
            fprintf(fid, '\t\t\t</tfoot>\n');
            fprintf(fid, '\t\t</table>\n');
            fprintf(fid, '\t</body>\n');
            fprintf(fid, '</html>\n');
            fclose(fid);
        end
        
        function times = calcProcessRunTimes(profData)
            if ~exist('profData', 'var')
                profData = profile('info');
            end
            
            functionNames = {profData.FunctionTable.FunctionName};           
            
            names = {};
            for i = 1:numel(profData.FunctionTable)
                processClass = regexp(profData.FunctionTable(i).FunctionName, '(.*?)>.*?\.evolveState$', 'tokens');
                if isempty(processClass); continue;
                end                
                
                names = [names; processClass{1}(1)]; %#ok<AGROW>
            end
            names = [{'Process'}; sort(names)];
            
            times = struct;            
            for i = 1:numel(names)
                times(i).name = names{i};
                
                %evolveState
                idx = find(strcmp(functionNames, [times(i).name '>' times(i).name '.evolveState']), 1);
                if ~isempty(idx)
                    times(i).evolveState = profData.FunctionTable(idx).TotalTime;
                else
                    times(i).evolveState = 0;
                end
                
                %copyFromState
                idx = find(strcmp(functionNames, [times(i).name '>' times(i).name '.copyFromState']), 1);
                if ~isempty(idx)
                    times(i).copyFromState = profData.FunctionTable(idx).TotalTime;
                else
                    times(i).copyFromState = 0;
                end
                
                %copyToState
                idx = find(strcmp(functionNames, [times(i).name '>' times(i).name '.copyToState']), 1);
                if ~isempty(idx)
                    times(i).copyToState = profData.FunctionTable(idx).TotalTime;
                else
                    times(i).copyToState = 0;
                end
                
                %calcResourceRequirements_Current
                idx = find(strcmp(functionNames, [times(i).name '>' times(i).name '.calcResourceRequirements_Current']), 1);
                if ~isempty(idx)
                    times(i).calcResourceRequirements_Current = profData.FunctionTable(idx).TotalTime;
                else
                    times(i).calcResourceRequirements_Current = 0;
                end
            end
        end
    end
end