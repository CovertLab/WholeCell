%Documentation
% Calls m2html to generate documentation of m code
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/8/2011
classdef Documentation
    methods (Static)
        function generate()
			import edu.stanford.covert.cell.sim.util.Documentation;
		
            [~, folder] = fileparts(pwd);
			try
				Documentation.cleanup(['doc' filesep 'm2html']);				
			end
			if ~exist(['doc' filesep 'm2html'], 'dir')
				mkdir(['doc' filesep 'm2html']);
			end
            directories = [{folder}; cellfun(@(x) [folder filesep x], [
                edu.stanford.covert.util.IOUtil.getDirectoryNamesRecursively('src')], ...
                'UniformOutput', false)];
            cd('..');
            try
                m2html(...
                    'mFiles', directories, ...
                    'htmlDir', [folder filesep 'doc' filesep 'm2html'], ...
                    'recursive', 'off',...
                    'graph', 'on', ...
                    'todo', 'on', ...
                    'verbose', 'off');
				cd(folder);
            catch exception
				cd(folder);
				exception.rethrow();
			end            
        end
		
		function cleanup(folder)
			import edu.stanford.covert.cell.sim.util.Documentation;
			
			files = dir(folder);
			for i = 1:numel(files)
				if isequal(files(i).name, '.') || isequal(files(i).name, '..') || isequal(files(i).name, '.svn')
					continue;
				end
				if isdir([folder filesep files(i).name])
					Documentation.cleanup([folder filesep files(i).name]);
				else
					delete([folder filesep files(i).name]);
				end
			end
		end
    end
end