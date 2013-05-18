% Instantiate simulation from knowledge base.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/20/2012
classdef SimulationRunnerDecay1 < handle
    %default options
    properties (SetAccess = protected)
        kbCache          = 'data/knowledgeBase.mat'
        simCache         = 'data/Simulation_fitted%s.mat'
        simJsonCache     = 'data/%s.json'
        
        outDir           = ''
        logToDisk        = false
        logToDb          = false
        childDirs        = {}
        
        useCachedKb      = true
        useCachedSim     = true
        cacheKb          = false
        cacheSimulation  = false
    end
    
    properties (Constant = true)
        DEFAULT_DISPLACED_BY = {
            'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE'
            'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
            'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
            'MG_001_DIMER'
            'MG_001_MONOMER'
            'MG_073_206_421_TETRAMER'
            'MG_094_HEXAMER'
            'MG_097_MONOMER'
            'MG_105_OCTAMER'
            'MG_184_DIMER'
            'MG_235_MONOMER'
            'MG_244_DIMER'
            'MG_254_MONOMER'
            'MG_262_MONOMER'
            'MG_339_MONOMER'
            'MG_352_DIMER'
            'MG_358_359_10MER'
            'MG_438_MONOMER'
            'MG_498_MONOMER'
            'RNA_POLYMERASE_HOLOENZYME'
            'RNA_POLYMERASE'
            };
        DEFAULT_DISPLACES = {
            'DNA_GYRASE'
            'MG_122_MONOMER'
            'MG_203_204_TETRAMER'
            'MG_213_214_298_6MER_ADP'
            'MG_469_1MER_ADP'
            'MG_469_1MER_ATP'
            };
    end
    
    methods
        function this = SimulationRunnerDecay1(varargin)
            if nargin >=1 && isstruct(varargin{1})
                options = varargin{1};
            else
                options = struct();
                for i = 1:2:numel(varargin)
                    options.(varargin{i}) = varargin{i+1};
                end
            end
            
            if isfield(options, 'childDir1') && isfield(options, 'childDir2')
                options.childDirs = {options.childDir1; options.childDir2};
            end
            
            metaClass = metaclass(this);
            for i = 1:numel(metaClass.Properties)
                if isfield(options, metaClass.Properties{i}.Name)
                    this.(metaClass.Properties{i}.Name) = options.(metaClass.Properties{i}.Name);
                end
            end
        end
        
        function [sim, kb] = run(this)
            metaClass = metaclass(this);
            simName = strrep(metaClass.Name, 'edu.stanford.covert.cell.sim.runners.', '');
            
            %construct KB and simulation
            fprintf('Constructing %s simulation ...\n', simName);
            [sim, kb] = this.constructKbAndSimulation();
            
            %run simulation
            fprintf('Running %s simulation ...\n', simName);
            this.runSimulation(sim);
        end
        
        function [sim, kb] = constructKbAndSimulation(this)
            import edu.stanford.covert.cell.sim.Simulation;
            import edu.stanford.covert.cell.sim.util.FitConstants;
            
            %% load default KB and simulation
            %load default network structure and parameters
            [kb, kbWid] = this.loadDefaultNetwork();
            
            %initialize simulation from default KB
            defSim = this.loadDefaultSimulation(kb, kbWid);
            
            %% modify network
            if isequal('edu.stanford.covert.cell.sim.runners.SimulationRunnerDecay1', class(this))
                sim = defSim;
            else
                %% modify structure
                this.modifyNetworkStructure(kb);
                
                %% create simulation with modified structure
                modSim = Simulation(kb.states, kb.processes);
                modSim.initializeConstants(kb);
                
                %% merge fitted parameters into new simulation
                defFitter = FitConstants(defSim, struct('verbosity', 1));
                modFitter = FitConstants(modSim, struct('verbosity', 1));
                
                %get parameter values
                defParams = defFitter.constructParameterVectorFromSimulation();
                modParams = modFitter.initializeFittedConstants();
                
                [defRnaExp, defNmpComp, defAaComp, defRnaWtFracs, defRnaDecayRates, ...
                    defBioComp, defBioProd, defByprod, defUnaccECons] = ...
                    defFitter.extractParameterVector(defParams);
                [modRnaExp, ~, ~, ~, modRnaDecayRates, ...
                    modBioComp, modBioProd, modByprod, ~] = ...
                    modFitter.extractParameterVector(modParams);
                
                %merge default parameter values into modified simulation
                defMet = defSim.state('Metabolite');
                modMet = modSim.state('Metabolite');
                [tfs, idxs] = ismember(modMet.wholeCellModelIDs, defMet.wholeCellModelIDs);
                modBioComp(tfs) = defBioComp(idxs(tfs));
                modBioProd(tfs) = defBioProd(idxs(tfs));
                modByprod(tfs) = defByprod(idxs(tfs));
                modBioComp(~tfs) = 0;
                modBioProd(~tfs) = 0;
                modByprod(~tfs) = 0;
                
                defRna = defSim.state('Rna');
                modRna = modSim.state('Rna');
                [tfs, idxs] = ismember(modRna.wholeCellModelIDs(modRna.matureIndexs), defRna.wholeCellModelIDs(defRna.matureIndexs));
                modRnaExp(tfs) = defRnaExp(idxs(tfs));
                modRnaDecayRates(tfs) = defRnaDecayRates(idxs(tfs));
                modRnaExp(~tfs) = 0;
                
                modNmpComp = defNmpComp;
                modAaComp = defAaComp;
                modRnaWtFracs = defRnaWtFracs;
                modUnaccECons = defUnaccECons;
                
                %set parameters
                modParams = modFitter.constructParameterVector(...
                    modRnaExp, modNmpComp, modAaComp, modRnaWtFracs, modRnaDecayRates, ...
                    modBioComp, modBioProd, modByprod, modUnaccECons);
                modFitter.applyParameterVectorToSimulation(modParams);
                
                %% modify parameter values
                this.modifyNetworkParameters(modSim);
                
                %% return new simulation
                sim = modSim;
            end
        end
        
        function [val, growths] = isSimulationFitted(~, sim, nSample)
            import xunit.utils.compareFloats;
            if nargin < 3
                nSample = 100;
            end
            
            %handle to MetabolicReaction state
            mr = sim.state('MetabolicReaction');
            
            %save filter
            initialGrowthFilterWidth = mr.initialGrowthFilterWidth;
            mr.initialGrowthFilterWidth = Inf;
            
            %sample initial conditions
            growths = zeros(nSample, 1);
            for i = 1:nSample
                sim.applyOptions('seed', i);
                sim.initializeState();
                growths(i) = mr.growth;
            end
            
            %restore filter
            mr.initialGrowthFilterWidth = initialGrowthFilterWidth;
            
            %test
            val = ...
                compareFloats(mr.meanInitialGrowthRate, mean(growths(growths > 1e-6)),   'elementwise', 'relative', 1e-1, sqrt(eps)) && ...
                compareFloats(mr.meanInitialGrowthRate, median(growths(growths > 1e-6)), 'elementwise', 'relative', 1e-1, sqrt(eps));
        end
        
        % write simulation data
        function cacheSimulationParametersJson(this, sim)
            data = {'options', 'parameters', 'fixedConstants', 'fittedConstants'};
            for  i = 1:numel(data)
                fid = fopen(sprintf(this.simJsonCache, data{i}), 'w');
                tmp = sim.(['get' upper(data{i}(1)) data{i}(2:end)])();
                if ismember(data{i}, {'fixedConstants', 'fittedConstants'})
                    for j = 1:numel(sim.states)
                        stateID = sim.states{j}.wholeCellModelID(7:end);
                        fields = fieldnames(tmp.states.(stateID));
                        for k = 1:numel(fields)
                            tmp.states.(stateID).(fields{k}) = [];
                        end
                    end
                    
                    for j = 1:numel(sim.processes)
                        processID = sim.processes{j}.wholeCellModelID(9:end);
                        fields = fieldnames(tmp.processes.(processID));
                        for k = 1:numel(fields)
                            tmp.processes.(processID).(fields{k}) = [];
                        end
                    end
                end
                fwrite(fid, edu.stanford.covert.io.jsonFormat(tmp));
                fclose(fid);
            end
        end
        
        %% write simulation test fixture
        function cacheSimulationTestFixtures(~, sim)
            %import
            import edu.stanford.covert.cell.sim.CellStateFixture;
            import edu.stanford.covert.cell.sim.ProcessFixture;
            import edu.stanford.covert.cell.sim.SimulationFixture;
            
            %store simulation fixture
            SimulationFixture.store(sim, 'Simulation.mat');
            
            %% unset mass reference in geometry state to reduce fixture sizes
            mass = sim.state('Mass');
            sim.state('Geometry').mass = [];
            
            %% write state test fixtures
            for i = 1:length(sim.states)
                CellStateFixture.store(sim.states{i});
            end
            
            %% write process test fixtures
            for i = 1:length(sim.processes)
                ProcessFixture.store(sim.processes{i});
            end
            
            %% reset mass reference in geometry state, and regenerate fixtures
            sim.state('Geometry').mass = mass;
            CellStateFixture.store(sim.state('Geometry'));
            ProcessFixture.store(sim.process('Metabolism'));
        end
        
        function runSimulation(this, sim)
            %% import
            import edu.stanford.covert.cell.sim.constant.Condition;
            import edu.stanford.covert.cell.sim.util.ConditionSet;
            import edu.stanford.covert.cell.sim.util.DatabaseLogger;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            import edu.stanford.covert.cell.sim.util.SummaryLogger;
            import edu.stanford.covert.db.MySQLDatabase;
            
            %% process options, setup loggers
            summaryLogger = SummaryLogger(1, 2);
            loggers = {summaryLogger};
            ic = [];
            
            %load perturbations
            if ~isempty(this.outDir) && exist([this.outDir filesep 'conditions.xml'], 'file')
                data = ConditionSet.parseConditionSet(sim, [this.outDir filesep 'conditions.xml']);
                data.metadata.knowledgeBaseWID = knowledgeBaseWID;
                sim.applyOptions(data.options);
                sim.applyOptions(data.perturbations);
                sim.applyParameters(data.parameters);
                summaryLogger.setOptions(struct('verbosity', sim.verbosity, 'outputDirectory', this.outDir));
            end
            
            %load initial conditions
            if ~isempty(this.outDir) && exist([this.outDir filesep 'initialConditions.mat'], 'file')
                ic = load([this.outDir filesep 'initialConditions.mat']);
            end
            
            %setup disk logger
            if this.logToDisk
                diskLogger = DiskLogger(this.outDir, 100);
                if exist('data', 'var')
                    diskLogger.addMetadata(data.metadata);
                else
                    diskLogger.addMetadata(struct(...
                        'shortDescription', '', ...
                        'longDescription', '', ...
                        'email', '', ...
                        'firstName', '', ...
                        'lastName', '', ...
                        'affiliation', '', ...
                        'knowledgeBaseWID', '', ...
                        'revision', '', ...
                        'differencesFromRevision', '', ...
                        'userName', '', ...
                        'hostName', '', ...
                        'ipAddress', '' ...
                        ));
                end
                loggers = [loggers; {diskLogger}];
            end
            
            %setup database logger
            if this.logToDb
                databaseLogger = DatabaseLogger(MySQLDatabase(config), 1000);
                databaseLogger.addMetadata(data.metadata);
                loggers = [loggers; {databaseLogger}];
            end
            
            %% simulate
            try
                decayMultiplier = 50;
                pm = sim.state('ProteinMonomer');
                matchA = strcmp(pm.wholeCellModelIDs, 'MonomerA');
                for i = 1:length(matchA)
                    if matchA(i) == 1
                        pm.decayRates(i) = decayMultiplier * pm.decayRates(i);
                    end
                end
                pc = sim.state('ProteinComplex');
                matchAA = strcmp(pc.wholeCellModelIDs, 'ComplexAA');
                for i = 1:length(matchAA)
                    if matchAA(i) == 1
                        pc.decayRates(i) = decayMultiplier * pc.decayRates(i);
                    end
                end

                % cell cycle up to, but not including division
                currentTime = clock;
                newSeed = 10000 * (60 * currentTime(5) + currentTime(6));
                sim.applyOptions('seed', newSeed);
                sim.run(ic,  loggers);
                
                % division
                if ~isempty(this.childDirs) && ...
                        sim.state('Chromosome').ploidy == 2 && ...
                        sim.state('Chromosome').segregated && ...
                        sim.state('Geometry').pinchedDiameter == 0
                    [~, daughters] = sim.divideState();
                    
                    daughter1 = daughters(1); %#ok<NASGU>
                    daughter2 = daughters(2); %#ok<NASGU>
                    
                    save([this.childDirs{1} filesep 'initialConditions.mat'], '-struct', 'daughter1');
                    save([this.childDirs{2} filesep 'initialConditions.mat'], '-struct', 'daughter2');
                end
            catch exception
                if ~isempty(this.outDir)
                    errFile = ['tmp' filesep 'simulation_error_' datestr(now, 'yyyy_mm_dd_HH_MM_SS_FFF') '.mat'];
                else
                    errFile = [this.outDir filesep 'err.mat'];
                end
                
                try
                    save(errFile, 'sim');
                catch saveException
                    exception.addCause(MException('Simulation:runTimeError', ...
                        'Simulation run-time error at %d. Unable to save simulation and exiting.', ...
                        sim.state('Time').values)).addCause(saveException).rethrow();
                end
                
                exception.addCause(MException('Simulation:runTimeError', ...
                    'Simulation run-time error at %d. Saving simulation at ''%s'' and exiting.', ...
                    sim.state('Time').values, errFile)).rethrow();
            end
        end
    end
    
    methods (Access = protected)
        function [knowledgeBase, knowledgeBaseWID] = loadDefaultNetwork(this)
            import edu.stanford.covert.cell.kb.KnowledgeBase;
            import edu.stanford.covert.cell.kb.KnowledgeBaseUtil;
            import edu.stanford.covert.db.MySQLDatabase;
            
            if ~this.useCachedKb || ~exist(this.kbCache, 'file')
                % connect to database
                dbConnectionParameters = config();
                database = MySQLDatabase(dbConnectionParameters);
                
                % construct latest knowledge base from database
                knowledgeBaseWID = KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
                knowledgeBase = KnowledgeBase(database, knowledgeBaseWID);
                
                %serialize
                knowledgeBase.serializeLinks();
                
                % save
                if this.cacheKb
                    save(this.kbCache, 'knowledgeBase', 'knowledgeBaseWID');
                end
                
                %deserialize
                knowledgeBase.deserializeLinks();
                
                % cleanup
                database.close();
            else
                %load cached KB
                tmp = load(this.kbCache);
                knowledgeBase = tmp.knowledgeBase;
                knowledgeBaseWID = tmp.knowledgeBaseWID;
                knowledgeBase.deserializeLinks();
            end
        end
        
        function sim = loadDefaultSimulation(this, kb, kbWid)
            import edu.stanford.covert.cell.sim.Simulation;
            import edu.stanford.covert.cell.sim.util.CachedSimulationObjectUtil;
            import edu.stanford.covert.cell.sim.util.FitConstants;
            
            if ~this.useCachedKb || ~this.useCachedSim || ~exist(sprintf(this.simCache, ''), 'file')
                sim = Simulation(kb.states, kb.processes);
                sim.initializeConstants(kb);
                sim.applyOptions('seed', 1);
                
                %fit simulation
                fitter = FitConstants(sim, struct('verbosity', 1));
                fitter.run();
                
                %choose an initial simulation state with desired growth rate
                sim.applyOptions('seed', 1);
                mr = sim.state('MetabolicReaction');
                initGrowthFilterW = mr.initialGrowthFilterWidth;
                mr.initialGrowthFilterWidth = 1e-2;
                sim.initializeState();
                mr.initialGrowthFilterWidth = initGrowthFilterW;
                
                %cache simulation object
                if this.cacheSimulation
                    CachedSimulationObjectUtil.store(sim, kbWid, this.simCache);
                end
            else
                sim = CachedSimulationObjectUtil.load();
                sim.constructRandStream();
            end
        end
        
        function modifyNetworkStructure(~, kb)
            kb.invalidate();
            kb.calcIndices();
            kb.computeDependentProperties();
        end
        
        function newRxns = createDefaultDisplacementReactions(this, kb, newProts)
            import edu.stanford.covert.cell.kb.Reaction;
            
            cytosol = findobj(kb.compartments, 'wholeCellModelID', 'c');
            nucleoid = findobj(kb.compartments, 'wholeCellModelID', 'd');
            
            chrState = findobj(kb.states, 'wholeCellModelID', 'State_Chromosome');
            
            newRxns = Reaction.empty(0, 1);
            for i = 1:numel(newProts)
                for j = 1:numel(this.DEFAULT_DISPLACED_BY)
                    tmp = Reaction(kb, NaN, ...
                        {sprintf('DNABoundProteinRelease_%s_%s', this.DEFAULT_DISPLACED_BY{j}, newProts(i).wholeCellModelID)}, ...
                        {sprintf('DNA bound protein release (%s releases %s)', this.DEFAULT_DISPLACED_BY{j}, newProts(i).name)}, ...
                        {''}, {''}, {'N'}, {'F'}, NaN, NaN, ...
                        {''}, NaN, NaN, {'1/min'}, ...
                        {''}, NaN, NaN, {'1/min'},  ...
                        NaN, NaN, ...
                        {''}, {''}, ...
                        0, 1000, {'mmol/gDCW/h'});
                    newRxns = [newRxns; tmp]; %#ok<AGROW>
                    
                    tmp.state = chrState;
                    
                    tmpMon = findobj(kb.proteinMonomers, 'wholeCellModelID', this.DEFAULT_DISPLACED_BY{j});
                    tmpCpx = findobj(kb.proteinComplexs, 'wholeCellModelID', this.DEFAULT_DISPLACED_BY{j});
                    if ~isempty(tmpMon)
                        tmp.proteinMonomers = [tmpMon; tmpMon];
                        tmp.proteinMonomerCompartments = [nucleoid; cytosol];
                        tmp.proteinMonomerCoefficients = [1; -1];
                        
                        tmp.proteinComplexs = [newProts(i); newProts(i)];
                        tmp.proteinComplexCompartments = [cytosol; nucleoid];
                        tmp.proteinComplexCoefficients = [1; -1];
                    else
                        tmp.proteinComplexs = [newProts(i); newProts(i); tmpCpx; tmpCpx];
                        tmp.proteinComplexCompartments = [cytosol; nucleoid; nucleoid; cytosol];
                        tmp.proteinComplexCoefficients = [1; -1; 1; -1];
                    end
                end
                
                for j = 1:numel(this.DEFAULT_DISPLACES)
                    tmpMon = findobj(kb.proteinMonomers, 'wholeCellModelID', this.DEFAULT_DISPLACES{j});
                    tmpCpx = findobj(kb.proteinComplexs, 'wholeCellModelID', this.DEFAULT_DISPLACES{j});
                    if ~isempty(tmpMon)
                        tmp = Reaction(kb, NaN, ...
                            {sprintf('DNABoundProteinRelease_%s_%s', newProts(i).wholeCellModelID), tmpMon.wholeCellModelID}, ...
                            {sprintf('DNA bound protein release (%s releases %s)', newProts(i).name, tmpMon.name)}, ...
                            {''}, {''}, {'N'}, {'F'}, NaN, NaN, ...
                            {''}, NaN, NaN, {'1/min'}, ...
                            {''}, NaN, NaN, {'1/min'},  ...
                            NaN, NaN, ...
                            {''}, {''}, ...
                            0, 1000, {'mmol/gDCW/h'});
                        
                        if isa(newProts, 'edu.stanford.covert.cell.kb.ProteinMonomer')
                            tmp.proteinMonomers = [newProts(i); newProts(i); tmpMon; tmpMon];
                            tmp.proteinMonomerCompartments = [cytosol; nucleoid; nucleoid; cytosol];
                            tmp.proteinMonomerCoefficients = [-1; 1; -1; 1];
                        else
                            tmp.proteinMonomers = [tmpMon; tmpMon];
                            tmp.proteinMonomerCompartments = [nucleoid; cytosol];
                            tmp.proteinMonomerCoefficients = [-1; 1];
                            
                            tmp.proteinComplexs = [newProts(i); newProts(i)];
                            tmp.proteinComplexCompartments = [cytosol; nucleoid];
                            tmp.proteinComplexCoefficients = [-1; 1];
                        end
                    else
                        tmp = Reaction(kb, NaN, ...
                            {sprintf('DNABoundProteinRelease_%s_%s', newProts(i).wholeCellModelID), tmpCpx.wholeCellModelID}, ...
                            {sprintf('DNA bound protein release (%s releases %s)', newProts(i).name, tmpCpx.name)}, ...
                            {''}, {''}, {'N'}, {'F'}, NaN, NaN, ...
                            {''}, NaN, NaN, {'1/min'}, ...
                            {''}, NaN, NaN, {'1/min'},  ...
                            NaN, NaN, ...
                            {''}, {''}, ...
                            0, 1000, {'mmol/gDCW/h'});
                        
                        if isa(newProts, 'edu.stanford.covert.cell.kb.ProteinMonomer')
                            tmp.proteinMonomers = [newProts(i); newProts(i)];
                            tmp.proteinMonomerCompartments = [cytosol; nucleoid];
                            tmp.proteinMonomerCoefficients = [-1; 1];
                            
                            tmp.proteinComplexs = [tmpCpx; tmpCpx];
                            tmp.proteinComplexCompartments = [nucleoid; cytosol];
                            tmp.proteinComplexCoefficients = [-1; 1];
                        else
                            tmp.proteinComplexs = [newProts(i); newProts(i); tmpCpx; tmpCpx];
                            tmp.proteinComplexCompartments = [cytosol; nucleoid; nucleoid; cytosol];
                            tmp.proteinComplexCoefficients = [-1; 1; -1; 1];
                        end
                    end
                    
                    tmp.state = chrState;
                    
                    newRxns = [newRxns; tmp]; %#ok<AGROW>
                end
            end
            
            chrState.reactions = [chrState.reactions; newRxns];
            kb.reactions = [kb.reactions; newRxns];
        end
        
        function [rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                bioComp, bioProd, byprod, unaccECons] = ...
                modifyNetworkParameters(this, sim, ...
                rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                bioComp, bioProd, byprod, unaccECons) %#ok<INUSL,MANU>
        end
    end
end


