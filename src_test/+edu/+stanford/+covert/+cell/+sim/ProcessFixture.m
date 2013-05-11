%ProcessFixture
% Class for generating fixtures for Process classes.
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef ProcessFixture

    properties (Constant)
        pathPrefix = edu.stanford.covert.cell.sim.ProcessFixture.determinePathPrefix();
    end

    methods (Static)
        function store(process, filename)
            import edu.stanford.covert.cell.sim.ProcessFixture;
            if ~exist('filename', 'var')
                filename = ProcessFixture.defaultFilename(process);
            end
            process.verbosity = 0;
            process.seed = 1;
            process.seedRandStream();

            fixture = process; %#ok<NASGU>
            filename = [ProcessFixture.pathPrefix filename];
            save('-v7', filename, 'fixture');
            edu.stanford.covert.io.editMatFileHeader(filename);
        end

        function process = load(process, filename, loadStates)
            import edu.stanford.covert.cell.sim.ProcessFixture;
            if ~exist('filename', 'var') || isempty(filename)
                filename = ProcessFixture.defaultFilename(process);
            end
            if ~exist('loadStates', 'var') || isempty(loadStates)
                loadStates = true;
            end
            
            tmp = load([ProcessFixture.pathPrefix filename]);
            process = tmp.fixture;
            
            if ~loadStates
                metaClass = metaclass(process);
                process.states = [];
                for i = 1:numel(metaClass.Properties)
                    name = metaClass.Properties{i}.Name;
                    if ~metaClass.Properties{i}.Dependent && isa(process.(name), 'edu.stanford.covert.cell.sim.CellState')
                        process.(name) = [];
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function value = determinePathPrefix()
            loc = which('edu.stanford.covert.cell.sim.ProcessFixture');
            i = strfind(loc, filesep);
            value = [loc(1:i(end)) '+process' filesep 'fixtures' filesep];
        end

        function value = defaultFilename(process)
            className = class(process);
            i = strfind(className, '.');
            value = [className(i(end)+1:end) '.mat'];
        end

        function value = determinePropertyNames(process)
            mc = metaclass(process);
            value = cell(numel(mc.Properties), 1);
            for i = 1:length(value)
                p = mc.Properties{i};
                if strcmp(p.GetAccess,'public') && ...
                   strcmp(p.SetAccess,'public') && ...
                   ~p.Dependent
                    value{i} = p.Name;
                end
            end
            value = value(~cellfun('isempty', value));
        end
    end
end
