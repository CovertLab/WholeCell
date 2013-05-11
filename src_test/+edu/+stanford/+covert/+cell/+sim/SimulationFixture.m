%SimulationFixture
% Class for generating fixtures for Simuation classes.
%
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef SimulationFixture
    properties (Constant)
        pathPrefix = edu.stanford.covert.cell.sim.SimulationFixture.determinePathPrefix();
        defaultProps = edu.stanford.covert.cell.sim.SimulationFixture.determineDefaultPropertyNames();
        defaultFilename = 'Simulation.mat';
    end

    methods (Static)
        function store(sim, filename)
            import edu.stanford.covert.cell.sim.SimulationFixture;
            if ~exist('filename', 'var')
                filename = SimulationFixture.defaultFilename;
            end
            
            fixture = sim; %#ok<NASGU>
            filename = [SimulationFixture.pathPrefix filename];
            save('-v7', filename, 'fixture');
            edu.stanford.covert.io.editMatFileHeader(filename);
        end

        function sim = load(filename, processes)
            import edu.stanford.covert.cell.sim.SimulationFixture;
            if ~exist('filename', 'var') || isempty(filename)
                filename = SimulationFixture.defaultFilename;
            end

            tmp = load([SimulationFixture.pathPrefix filename]);
            sim = tmp.fixture;

            if exist('processes', 'var')                    
                if iscell(processes)
                    processes = cellfun(@(x) ['Process_' x], processes, 'UniformOutput', false);

                    tfs = false(size(sim.processes));
                    for i = 1:numel(sim.processes)
                        tfs(i) = ismember(sim.processes{i}.wholeCellModelID, processes);
                    end
                    sim.setForTest('processes', sim.processes(tfs));

                    tfs = false(size(sim.processesInInitOrder));
                    for i = 1:numel(sim.processesInInitOrder)
                        tfs(i) = ismember(sim.processesInInitOrder{i}.wholeCellModelID, processes);
                    end
                    sim.setForTest('processesInInitOrder', sim.processesInInitOrder(tfs));
                    [~, idxs] = ismember(...
                        cellfun(@(x) x.wholeCellModelID, sim.processesInInitOrder, 'UniformOutput', false), ...
                        cellfun(@(x) x.wholeCellModelID, sim.processes, 'UniformOutput', false));
                    sim.setForTest('processInitOrderIndexs', idxs);

                    tfs = false(size(sim.processesInEvalOrder));
                    for i = 1:numel(sim.processesInEvalOrder)
                        tfs(i) = ismember(sim.processesInEvalOrder{i}.wholeCellModelID, processes);
                    end
                    sim.setForTest('processesInEvalOrder', sim.processesInEvalOrder(tfs));
                    [~, idxs] = ismember(...
                        cellfun(@(x) x.wholeCellModelID, sim.processesInEvalOrder, 'UniformOutput', false), ...
                        cellfun(@(x) x.wholeCellModelID, sim.processes, 'UniformOutput', false));
                    sim.setForTest('processEvalOrderIndexs', idxs);
                end
            else
                metaClass = metaclass(sim);
                for i = 1:numel(metaClass.Properties)
                    if ...
                            strcmp(metaClass.Properties{i}.GetAccess, 'public') && ...
                            ~metaClass.Properties{i}.Dependent && ...
                            isa(sim.(metaClass.Properties{i}.Name), 'edu.stanford.covert.cell.sim.Process')
                        sim.setForTest(metaClass.Properties{i}.Name, []);
                    end
                end
                sim.setForTest('processes', {});                    
                sim.setForTest('processesInInitOrder', {});
                sim.setForTest('processesInEvalOrder', {});
                sim.setForTest('processEvalOrderIndexs', []);
                sim.setForTest('processInitOrderIndexs', []);
            end
            
            %set options
            sim.applyOptions('verbosity', 0, 'seed', 1);
            sim.seedRandStream();
            
            %seed random generators
            for i = 1:length(sim.states)
                state = sim.states{i};
                state.seed = 1;
                state.seedRandStream();
            end
            
            for i = 1:length(sim.processes)
                process = sim.processes{i};
                process.seed = 1;
                process.seedRandStream();
            end
            
            %rebuild function handles in enzyme properties (not preserved well
            %through saving in MATLAB <= 2009b)
            m = sim.process('DNASupercoiling');
            if ~isempty(m)
                m.enzymeProperties = m.buildEnzymeProperties();
            end
            
            %% set references in processes to state objects
            for i = 1:length(sim.states)
                sim.states{i}.storeObjectReferences(sim);
            end
            for i = 1:length(sim.processes)
                sim.processes{i}.storeObjectReferences(sim);
            end
            
            %% set dynamic properties of process that come from communication with simulation
            for i = 1:length(sim.processes)
                sim.processes{i}.copyFromState();
            end
            
            %% state caches
            s = sim.state('Chromosome');
            if ~isempty(s)
                s.invalidate();
            end
            
            s = sim.state('Mass');
            if ~isempty(s)
                s.calcMass();
            end
            
            s = sim.state('Geometry');
            if ~isempty(s)
                s.calculateVolume();
            end
        end
    end

    methods (Static, Access = private)
        function value = determinePathPrefix()
            loc = which('edu.stanford.covert.cell.sim.SimulationFixture');
            i = strfind(loc, filesep);
            value = [loc(1:i(end)) 'fixtures' filesep];
        end

        function value = determineDefaultPropertyNames()
            metaData = metaclass(edu.stanford.covert.cell.sim.Simulation);
            value = cell(numel(metaData.Properties), 1);
            for i = 1:length(value)
                p = metaData.Properties{i};
                if strcmp(p.GetAccess,'public') && ~p.Dependent && ~p.Constant
                    value{i} = p.Name;
                end
            end
            value = value(~cellfun('isempty', value));
            value = setdiff(value, {
                'processes'
                'processEvalOrderIndexs'
                'processInitOrderIndexs'
                'processesInInitOrder'
                'processesInEvalOrder'
                });
        end
    end
end
