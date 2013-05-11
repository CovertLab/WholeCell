%CellStateFixture
% Class for generating fixtures for CellState classes.
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef CellStateFixture

    properties (Constant)
        pathPrefix = edu.stanford.covert.cell.sim.CellStateFixture.determinePathPrefix();
    end

    methods (Static)
        function store(state, filename)
            import edu.stanford.covert.cell.sim.CellStateFixture;
            if ~exist('filename', 'var')
                filename = CellStateFixture.defaultFilename(state);
            end
            state.verbosity = 0;
            state.seed = 1;
            state.seedRandStream();            

            fixture = state; %#ok<NASGU>
            filename = [CellStateFixture.pathPrefix filename];
            save('-v7', filename, 'fixture');
            edu.stanford.covert.io.editMatFileHeader(filename);
        end

        function state = load(state, filename, loadStates)
            import edu.stanford.covert.cell.sim.CellStateFixture;
            if ~exist('filename', 'var') || isempty(filename)
                filename = CellStateFixture.defaultFilename(state);
            end
            if ~exist('loadStates', 'var') || isempty(loadStates)
                loadStates = true;
            end

            tmp = load([CellStateFixture.pathPrefix filename]);
            state = tmp.fixture;
            
            if ~loadStates
                names = fieldnames(state);
                for i = 1:numel(names)
                    if isa(state.(names{i}), 'edu.stanford.covert.cell.sim.CellState')
                        state.(names{i}) = [];
                    end
                end
            end
        end
    end

    methods (Static, Access = private)
        function value = determinePathPrefix()
            loc = which('edu.stanford.covert.cell.sim.CellStateFixture');
            i = strfind(loc, filesep);
            value = [loc(1:i(end)) '+state' filesep 'fixtures' filesep];
        end

        function value = defaultFilename(state)
            className = class(state);
            i = strfind(className, '.');
            value = [className(i(end)+1:end) '.mat'];
        end

        function value = determinePropertyNames(state)
            mc = metaclass(state);
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
