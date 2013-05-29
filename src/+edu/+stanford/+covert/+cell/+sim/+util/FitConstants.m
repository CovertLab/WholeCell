%FitConstants
%
% Introduction
% ================================================
% This class is the last stage in fitting the simulation. The purpose of
% this class is to resolve conflicts among the experimental data separately
% used to parameterize each of the processes, as well as to calculate the
% biomass composition / production used to parameterize the flux-balance
% analysis (FBA) metabolic model. Furthermore, the goal of this class is to
% resolve these conflicts in a way which keeps the values of the parameters
% as close to their experimental observed values as possible.
%
% Fitting workflow
% ================================================
% This class is the last, of several stages, of simulation fitting:
% 1. Collect parameter values from literature and databases
% 2. Fit free parameters over individual processes
%    - Fit chromosome condensation parameters to recapitulate observed
%      average SMC spacing
%    - Fit repliciation initiation cooperativity constants to
%      recapitulate desired repliciation initiation duration (observed
%      cell cycle length - calculated replication duration - simulation
%      cytokinesis duration)
% 3. Fit free parameters over groups of processes
% 4. Fit dry weight fractions to accomodate a full chromosome
%    replicated by the end of the replication phase (which is before the
%    end of the cytokinesis phase)
%    - Increase DNA dry weight fraction
%    - Decrease all other dry weight fractions
% 5. Fit parameters over entire model. Find set of parameter values
%    closest to their experimentally measured values which satisfy
%    several joint parameter constraints. Until convergence,
%    - Calculate biomass composition, production, and unaccounted enegry
%      consumption ("dark enegry")
%    - Calculate constraints using chosen geneExpressionRobustness
%    - Satisfy constriants heuristically
%    - Solve non-linear constrained optimization problem
%
% Remaining unconstrained free parameters
% ================================================
% The fitting work flow identifies and/or adjusts the value of every
% parameter in the simulation with 5 exceptions:
% - proteinMisfoldingRate
% - tmRNABindingProbability
% - geneExpressionRobustness
% - initialFractionNTPsInRNAs
% - initialFractionAAsInMonomers
%
% These are the only truly free parameters across the entire simulation.
%
% Parameter conflict resolution
% ================================================
% In particular, the goal of this class is to resolve conflicts among
% several pieces of experimental data, and to do so in a least squares
% sense:
% - gene expression
% - protein expression
% - RNA weight fractions
% - NMP composition
% - AA composition
% - gene sequences
% - protein sequences
% - transcription unit structure
% - RNA half lives
% - tRNA synthetase, transferase rates
% - RNA polymerase elongation rate
% - Ribosome elongation rate
% - RNA polymerase state expectations
% - cell cycle length
% - cell cycle phase lengths
% - genetic code
% - protein complex composition
% - enzyme kinetics (enzymeBounds)
% - transport rates (reactionBounds)
%
% Stated more formally, the goal of this class is to the identify the set
% of gene expression, NMP composition, AA composition, RNA weight
% fractions, and RNA decay rates which satisfy several linear and
% non-linear constraints and minimze their sum of squares deviation from
% their experimentally observed values. Stated mathmetatically,
%
% Minimize
%   ||W*(z - experimentalZ)||_2 =
%   z'*W*W*z - 2*z'*W*W*experimentalZ + experimentalZ'*W*W*experimentalZ
%
%           (rnaExpression)
% Where z = (nmpComposition)
%           (aaComposition)
%           (rnaWeightFractions)
%           (rnaDecayRates)
%
% Subject to:
% - normalization                                            ones * rnaExp  = 1
%                                                           ones * nmpComp  = 1
%                                                            ones * aaComp  = 1
%                                                        ones * rnaWtFracs  = 1
% - RNA type distribution                                          mRNAExp  = I_mRNA * rnaExp
%                                                                  rRNAExp  = I_rRNA * rnaExp
%                                                                  sRNAExp  = I_sRNA * rnaExp
%                                                                  tRNAExp  = I_tRNA * rnaExp
%                                      mRNAMWs * mRNAExp / rnaMWs * rnaExp  = rnaWtFracs(mRNAWtFracIdxs)
%                                      rrnaMWs * rRNAExp / rnaMWs * rnaExp  = rnaWtFracs(rRNAWtFracIdxs)
%                                      srnaMWs * sRNAExp / rnaMWs * rnaExp  = rnaWtFracs(sRNAWtFracIdxs)
%                                      trnaMWs * tRNAExp / rnaMWs * rnaExp  = rnaWtFracs(tRNAWtFracIdxs)
% - Monomer expression                                              monExp  = matureRNAGeneComp(mRNAIdxs, :) * rnaExp /
%                                                                             (ones * matureRNAGeneComp(mRNAGeneIdxs, :) * rnaExp)
% - NMP Composition                rnaBaseCnts * rnaExp / rnaLens * rnaExp  = nmpComp
% - AA Composition                   monAACnts * monExp / monLens * monExp  = aaComp
% - Doubling lower bounds                                     ribosome exp >= min exp for protein doubling
%                                                       RNA polymerase exp >= min exp for RNA doubling
%                                                 transcription factor exp >= min exp for RNA doubling
%                                                   translation factor exp >= min exp for protein doubling
%                                                                 tRNA exp >= min exp for protein doubling
%                                                      tRNA synthetase exp >= min exp for protein doubling
%                                                                         ...
% - FtsZ expression: held to value used to calculate cytokinesis duration
% - DnaA expression: held to value used to fit replication initiation
%   duration
% - Topoisomerase I / gyrase expression: constrained to produce observed
%   steady-state superhelical density
%
% That is we pose the problem as one of non-linear constrianed
% optimization, and use the MATLAB fmincon routine to identify the optimal
% parameter set, z. Linear constraints are implemented as pairs of
% matrices, A, and right-hand sides representing equality and inequality
% constriants. Dependent linear constraints are automatically removed.
% Non-linear constraints are implemented as class methods, are calculated
% in part using the calcResourceRequirements_LifeCycle methods of the
% processes through the calcResourceRequirements method of this class.
%
% However, because fmincon has troubling identifying solutions
% which satisfy all of the non-linear constraints, before executing fmincon
% we first hueristically identify a consistent set of parameter values:
% - transcription unit expression <- average expression of genes in
%   transcription unit; RNA weight fractions ./ RNA fraction
%   molecular weights
% - transcription unit decay rate <- average decay rate of genes in
%   transcription unit
% - NMP composition <- transcription unit sequences *
%   transcription unit expression
% - AA composition <- protein monomer sequences *
%   transcription unit composition(mRNAs, :) * transcription unit
%   expression
%
% The initial heuristic procedure modifies RNA expression, NMP and AA
% composition, and rRNA half lives. The initial heuristic procedure doesn't
% modify RNA weight fractions or m/s/tRNA half lives.
%
% Biomass composition, production calculation
% ================================================
% During, and following fitting we calculated the biomass composition,
% production, byproduct secretion, and unaccounted energy consumption. See
% documentation above calcResourceRequirements method for units and
% normalization conditions. Biomass production - byproducts forms the FBA
% objective. Biomass composition is used by to initialize the cell prior to
% simulation.
%
% Assumptions
% ================================================
% - uniform protein half lives
% - RNA half lives are for free species
% - RNA expression includes free and bound RNA species
% - bound RNA, protein (eg. by RNA polymerase, DNA) has same half life
%   as free RNA, protein
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 3/22/2011
classdef FitConstants < handle
    %options
    properties
        method = 'heuristic';
        verbosity = 0;
        maxIter = 100;
        tolerance = 2e-2;
        fminconOptions = optimset(...
            'TolFun', 1e-3, ...
            'TolCon', 1e-6, ...
            'TolX', 1e-3, ...
            'MaxTime', 10 * 60, ...
            'MaxFunEvals', 1e6, ...
            'ScaleProblem', 'none', ...
            'Algorithm', 'interior-point', ...
            'GradObj', 'on', ...
            'GradConstr', 'on', ...
            'Hessian', 'user-supplied', ...
            'Display', 'off');
    end
    
    %references
    properties
        simulation
    end
    
    %counts
    properties
        nGenes
        nRNAs
        nNMPs
        nAAs
        nRNAWtFracs
        nVars
    end
    
    %indices
    properties
        rnaExpIdxs
        nmpIdxs
        aaIdxs
        rnaDecayRateIdxs
        rnaWtFracIdxs
        
        rnaIndex_dnaA
        rnaIndex_ftsZ
        
        monomerIndex_gyrase
        monomerIndex_topoisomeraseI
        monomerIndex_DnaA
        monomerIndex_FtsZ
    end
    
    %constants
    properties
        gryStoichiometry
        dnaACnt
        ftsZCnt
    end
    
    methods
        function this = FitConstants(sim, options)
            this.simulation = sim;
            if nargin > 1
                if isfield(options, 'method'), this.method = options.method; end
                if isfield(options, 'verbosity'), this.verbosity = options.verbosity; end
                if isfield(options, 'maxIter'), this.maxIter = options.maxIter; end
                if isfield(options, 'tolerance'), this.tolerance = options.tolerance; end
                if isfield(options, 'fminconOptions'), this.fminconOptions = options.fminconOptions; end
            end
            
            %references
            g = sim.gene;
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            geom = sim.state('Geometry');
            ring = sim.state('FtsZRing');
            time = sim.state('Time');
            
            %counts
            this.nGenes = length(g.wholeCellModelIDs);
            this.nRNAs = length(r.matureIndexs);
            this.nNMPs = length(m.nmpIndexs);
            this.nAAs = length(m.aminoAcidIndexs);
            this.nRNAWtFracs = 6;
            this.nVars = this.nRNAs + this.nNMPs + this.nAAs + this.nRNAs + this.nRNAWtFracs;
            
            %indices
            this.rnaExpIdxs       = 1:this.nRNAs;
            this.nmpIdxs          = this.nRNAs + 1:this.nRNAs + this.nNMPs;
            this.aaIdxs           = this.nRNAs + this.nNMPs + 1:this.nRNAs + this.nNMPs + this.nAAs;
            this.rnaDecayRateIdxs = this.nRNAs + this.nNMPs + this.nAAs + 1:this.nRNAs + this.nNMPs + this.nAAs + this.nRNAs;
            this.rnaWtFracIdxs    = this.nRNAs + this.nNMPs + this.nAAs + this.nRNAs + 1:this.nRNAs + this.nNMPs + this.nAAs + this.nRNAs + this.nRNAWtFracs;
            
            [~, this.monomerIndex_gyrase] = ismember({'MG_003_MONOMER', 'MG_004_MONOMER'}, pm.wholeCellModelIDs(pm.matureIndexs));
            [~, this.monomerIndex_topoisomeraseI] = ismember({'MG_122_MONOMER'}, pm.wholeCellModelIDs(pm.matureIndexs));
            [~, this.monomerIndex_DnaA] = ismember({'MG_469_MONOMER'}, pm.wholeCellModelIDs(pm.matureIndexs));
            [~, this.monomerIndex_FtsZ] = ismember({'MG_224_MONOMER'}, pm.wholeCellModelIDs(pm.matureIndexs));
            
            this.rnaIndex_dnaA = find(r.matureRNAGeneComposition(g.mRNAIndexs(this.monomerIndex_DnaA), :));
            this.rnaIndex_ftsZ = find(r.matureRNAGeneComposition(g.mRNAIndexs(this.monomerIndex_FtsZ), :));
            
            %stoichiometries
            this.gryStoichiometry = max(max(max(pc.proteinComplexComposition(g.mRNAIndexs(this.monomerIndex_gyrase), :, :))));
            
            %fix DnaA copy number at value used in fitDuration method of
            %the Replication initation unit test class to fit replication
            %initiation duration
            this.dnaACnt = 54;
            
            %fix FtsZ copy number at value used in Cytokinesis medium test
            %to calculate cytokinesis duration
            numEdges = ring.calcNumEdges(...
                geom.calculateWidth(mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) / geom.density), ...
                ring.filamentLengthInNm);
            this.ftsZCnt = ceil(4.00 / exp(log(2) * (1 - time.cytokinesisDuration / time.cellCycleLength)) * ...
                numEdges * ring.numFtsZSubunitsPerFilament);
        end
        
        function this = run(this)
            %import classes
            import edu.stanford.covert.cell.sim.constant.Condition;
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            
            %toggle off warnings
            warningStatus = warning('query', 'WholeCell:warning');
            warning('off', 'WholeCell:warning');
            
            %seed random number generator
            seed = sim.seed;
            stateSeeds = zeros(size(sim.states));
            processSeeds = zeros(size(sim.processes));
            sim.applyOptions('seed', 0);
            sim.seedRandStream();
            for i = 1:numel(sim.states)
                o = sim.states{i};
                stateSeeds(i) = o.seed;
                o.seed = 0;
                o.seedRandStream();
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                processSeeds(i) = o.seed;
                o.seed = 0;
                o.seedRandStream();
            end
            
            %initialize fitted constants
            paramVec = this.initializeFittedConstants();
            this.applyParameterVectorToSimulation(paramVec);
            
            %fit constants (note: constants are applied to simulation
            %in the execution of these methods since some of the process
            %methods which are invoked depend on these parameters)
            switch this.method
                case 'analytic',  paramVec = this.fitAnalytically(paramVec);
                case 'heuristic', paramVec = this.fitHeuristically(paramVec);
                otherwise,        throw(MException('FitConstants:error', 'unknown method %s', this.method));
            end
            
            %apply fitted constants to simulation
            this.applyParameterVectorToSimulation(paramVec);
            
            %reset seeds
            sim.applyOptions('seed', seed);
            sim.seedRandStream();
            for i = 1:numel(sim.states)
                o = sim.states{i};
                o.seed = stateSeeds(i);
                o.seedRandStream();
            end
            for i = 1:numel(sim.processes)
                o = sim.processes{i};
                o.seed = processSeeds(i);
                o.seedRandStream();
            end
            
            %reload warning state
            warning(warningStatus.state, 'WholeCell:warning');
        end
        
        % initial guess of parameter values
        %- transcription unit expression: imputed, averaged gene expression
        %- transcription unit decay rates: imputed, averaged gene decay rates
        function paramVec = initializeFittedConstants(this)
            %% import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = this.simulation;
            t = sim.state('Time');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            mass = sim.state('Mass');
            mr = sim.state('MetabolicReaction');
            met = sim.process('Metabolism');
            
            mr.growth0 = 1 / mass.timeAveragedCellWeight / t.cellCycleLength;
            
            %% constants
            cIdx = sim.compartment.cytosolIndexs;
            
            %% biomass composition, production, byproducts
            biomassComposition = m.experimentalBiomassComposition;
            biomassProduction = max(0,  m.experimentalBiomassComposition);
            byproducts        = max(0, -m.experimentalBiomassComposition);
            
            biomassProduction(m.ntpIndexs, :) = biomassComposition(m.ntpIndexs, :) + biomassComposition(m.nmpIndexs, :);
            biomassProduction(m.nmpIndexs, :) = 0;
            biomassProduction(m.dntpIndexs, :) = biomassComposition(m.dntpIndexs, :) + biomassComposition(m.dnmpIndexs, :);
            biomassProduction(m.dnmpIndexs, :) = 0;
            
            unaccE = met.growthAssociatedMaintenance - sum(biomassProduction(m.ntpIndexs, cIdx));
            biomassProduction(m.ntpIndexs(1), cIdx) = biomassProduction(m.ntpIndexs(1), cIdx) + unaccE;
            biomassProduction(m.waterIndexs,  cIdx) = biomassProduction(m.waterIndexs,  cIdx) + unaccE;
            byproducts(m.ndpIndexs(1),        cIdx) = byproducts(m.ndpIndexs(1),        cIdx) + unaccE;
            byproducts(m.phosphateIndexs,     cIdx) = byproducts(m.phosphateIndexs,     cIdx) + unaccE;
            byproducts(m.hydrogenIndexs,      cIdx) = byproducts(m.hydrogenIndexs,      cIdx) + unaccE;
            
            normalization = mass.cellInitialDryWeight * ConstantUtil.nAvogadro / ...
                (sum(biomassComposition, 2)' * m.molecularWeights);
            biomassComposition = normalization * biomassComposition;
            biomassProduction  = normalization * biomassProduction;
            byproducts         = normalization * byproducts;
            unaccE             = normalization * unaccE;
            
            waterComp = mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * mass.fractionWetWeight  * ...
                ConstantUtil.nAvogadro / m.molecularWeights(m.waterIndexs);
            biomassComposition(m.waterIndexs, cIdx) = biomassComposition(m.waterIndexs, cIdx) + waterComp;
            biomassProduction(m.waterIndexs, cIdx) = biomassProduction(m.waterIndexs, cIdx) + waterComp;
            
            %% expression
            geneExp = ComputationUtil.imputeMissingData(r.expectedGeneExpression(:, 1));
            geneExp = geneExp / sum(geneExp);
            
            rnaExp = ComputationUtil.invertCompositionMatrix(r.matureRNAGeneComposition) * geneExp;
            rnaExp = rnaExp / sum(rnaExp);
            
            %% decay rates
            matureRNAGeneComposition    = r.matureRNAGeneComposition;
            invMatureRNAGeneComposition = ComputationUtil.invertCompositionMatrix(matureRNAGeneComposition);
            expectedGeneDecayRates      = ComputationUtil.imputeMissingData(r.expectedGeneDecayRates);
            rnaDecayRates               = invMatureRNAGeneComposition * expectedGeneDecayRates;
            
            %% match expression, decay rates
            rnaExp = (r.nascentRNAMatureRNAComposition * ...
                ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition) * ...
                (rnaExp .* (log(2) / t.cellCycleLength + rnaDecayRates))) ./ (log(2) / t.cellCycleLength + rnaDecayRates);
            
            %% construct parameter vector
            paramVec = this.constructParameterVector(...
                rnaExp, ...
                m.experimentalNMPComposition, ...
                m.experimentalAAComposition, ...
                r.expectedWeightFractions, ...
                rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, ...
                0.5 * unaccE);
        end
        
        function unscaledParamVec = fitAnalytically(this, unscaledParamVec0)
            import edu.stanford.covert.util.ComputationUtil;
            
            %% 1. Use approximate solution as initial guess for optimization problem
            paramVec0 = unscaledParamVec0;
            [~, ~, ~, ~, ~, paramVec0{1}(this.rnaExpIdxs)] = this.calcMacromolecularCounts(paramVec0);
            
            %heuristic fit
            unscaledParamVec = this.fitHeuristically(unscaledParamVec0);
            paramVec = unscaledParamVec;
            [~, ~, ~, ~, ~, rnaCnts] = this.calcMacromolecularCounts(paramVec);
            paramVec{1}(this.rnaExpIdxs) = rnaCnts;
            
            %% 2. Setup optimization problem
            
            %objective
            %1. Minimize change in parameters
            %2. Weight transcription unit expression, NMP production, AA production equally
            W = ones(size(paramVec0{1}));
            H = diag(W);
            objFcn = @(paramVec) this.objFunc(paramVec, paramVec0{1}, H);
            
            %equality constraints
            [Aineq, bineq, Aeq, beq] = this.linConstraintFunc(paramVec);
            
            %simple bounds
            [~, ~, ~, ~, rnaCntLBs, ~, rnaCntUBs, ~] = this.calcResourceRequirements(unscaledParamVec);
            
            lb = 1e-12 * ones(size(paramVec{1})); %don't use zero for lower bounds because zeros makes the problem ill conditioned
            ub = ones(size(paramVec{1}));
            lb(this.rnaExpIdxs) = min(rnaCntLBs, rnaCnts);
            ub(this.rnaExpIdxs) = max(rnaCntUBs, rnaCnts);
            
            %% 3. Check that initial guess satisfies constraints and bounds
            if ~all(Aineq * paramVec{1} < bineq)
                throw(MException('Simulation:error', 'Initial guess does not satisfy linear equality constraints'));
            end
            
            if max(abs(Aeq * paramVec{1} - beq)) > 1e-10
                throw(MException('Simulation:error', 'Initial guess does not satisfy linear equality constraints'));
            end
            
            [cineq, ceq] = this.nlinConstraintFunc(paramVec{1});
            if any(cineq > 0)
                throw(MException('Simulation:error', 'Initial guess does not satisfy non-linear inequality constraints'));
            end
            if max(abs(ceq)) > 2e-8
                throw(MException('Simulation:error', 'Initial guess does not satisfy non-linear equality constraints'));
            end
            
            if any(rnaCnts < lb(this.rnaExpIdxs)) || any(rnaCnts > ub(this.rnaExpIdxs))
                throw(MException('Simulation:error', 'Initial guess does not satisfy simple bound constraints'));
            end
            
            %% 4. Solve optimization problem
            %options
            fminconOpts = this.fminconOptions;
            fminconOpts.TypicalX = paramVec0{1};
            fminconOpts.HessFcn = @(paramVec, lambda) this.hessFunc(paramVec, lambda, H);
            
            %solve
            [paramVec{1}, ~, exitFlag, output] = fmincon(...
                objFcn, paramVec{1}, Aineq, bineq, Aeq, beq, lb, ub, @this.nlinConstraintFunc, fminconOpts);
            if ~strcmpi(fminconOpts.Display, 'off')
                fprintf('\n');
            end
            if exitFlag <= 0
                throw(MException('Simulation:error', 'Nonlinear programming error:\n%s', output.message));
            end
            if exitFlag > 0
                warning('WholeCell:warning', output.message);
            end
            
            %% 5. Heuristic fit
            unscaledParamVec = paramVec;
            unscaledParamVec{1}(this.rnaExpIdxs) = paramVec{1}(this.rnaExpIdxs) / sum(paramVec{1}(this.rnaExpIdxs));
            unscaledParamVec = this.fitHeuristically(unscaledParamVec);
        end
        
        function [val, gradient, hessian] = objFunc(~, paramVec, paramVec0, H)
            val = paramVec' * H * paramVec   -   2 * paramVec0' * H * paramVec;
            
            if nargout >= 2
                gradient = 2 * H * paramVec   -   2 * H * paramVec0;
            end
            
            if nargout >= 3
                hessian = 2 * H;
            end
        end
        
        function [Aineq, bineq, Aeq, beq] = linConstraintFunc(this, paramVec)
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            sim = this.simulation;
            g = sim.gene;
            m = sim.state('Mass');
            r = sim.state('Rna');
            sc = sim.process('DNASupercoiling');
            
            topTf = any(r.matureRNAGeneComposition(g.mRNAIndexs(this.monomerIndex_topoisomeraseI), :), 1);
            gyrTf = any(r.matureRNAGeneComposition(g.mRNAIndexs(this.monomerIndex_gyrase), :), 1);
            [rnaExp, ~, ~, rnaWtFracs, rnaDecayRates] = this.extractParameterVector(paramVec);
            
            %% inequality constraints
            %1. topoisomerase IV / gyrase
            %   fit topoisomerase I and gyrase to achieve the observed
            %   steady-state superhelical density of
            %   chromosome.equilibriumSuperhelicalDensity. topoisomerase IV
            %   doesn't contribute to the steady-state because it is
            %   inactive at the observed steady-state superhelical density.
            Aineq = zeros(0, this.nVars);
            bineq = zeros(0, 1);
            
            %% equality constraints
            %1. RNA expression normalized
            %2. NMP composition normalized
            %3. AA composition normalized
            %4. RNA weight fractions normalized
            %5. FtsZ expression equals value already used to fit
            %   cytokinesis duration. Here we implement this approximately
            %   as a linear constraint.
            %6. DnaA expression equals value already used to fit
            %   replication initiation duration. Here we implement this
            %   approximately as a linear constraint.
            
            %initialize
            Aeq = zeros(7+numel(this.rnaWtFracIdxs), this.nVars);
            beq = zeros(7+numel(this.rnaWtFracIdxs), 1);
            
            %RNA expression = RNA weight
            Aeq(1, this.rnaExpIdxs) = r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro;
            beq(1) = m.cellInitialDryWeight * m.dryWeightFractionRNA;
            
            %NMP composition is normalized (relative NMP composition sums to 1)
            Aeq(2, this.nmpIdxs) = 1;
            beq(2) = 1;
            
            %AA composition is normalized (relative AA composition sums to 1)
            Aeq(3, this.aaIdxs) = 1;
            beq(3) = 1;
            
            %RNA weight fractions (mRNA, rRNA, sRNA, tRNA) are constrained
            Aeq(3+(1:numel(this.rnaWtFracIdxs)), this.rnaWtFracIdxs) = eye(numel(this.rnaWtFracIdxs));
            beq(3+(1:numel(this.rnaWtFracIdxs))) = rnaWtFracs;
            
            %individual RNA decay rates may be adjusted, but total decay
            %rate is constrained, roughly equivalent to constraining the
            %average RNA decay rate
            Aeq(end-3, this.rnaDecayRateIdxs) = 1;
            beq(end-3) = sum(rnaDecayRates);
            
            %DnaA expression is constrained to a particular value so that
            %replication initiation duration doesn't need to be refit
            Aeq(end-2, this.rnaExpIdxs(this.rnaIndex_dnaA)) = 1;
            beq(end-2) = rnaExp(this.rnaIndex_dnaA);
            
            %FtsZ expression is constrained to a particular value so that
            %cytokinesis duration doesn't need to be recalculated
            Aeq(end-1, this.rnaExpIdxs(this.rnaIndex_ftsZ)) = 1;
            beq(end-1) = rnaExp(this.rnaIndex_ftsZ);
            
            %net supercoiling activity is zero (topoisomerase I positive
            %supercoiling activity is balanced by gyrase negative supercoiling activity)
            Aeq(end, this.rnaExpIdxs(topTf)) = sc.topoIDeltaLK  * sc.topoIActivityRate;
            Aeq(end, this.rnaExpIdxs(gyrTf)) = sc.gyraseDeltaLK * sc.gyraseActivityRate / this.gryStoichiometry * 1.3; %Factor of 1.3 added from observing simulation, trying to match initial and final conditions
            beq(end) = 0;
            
            %remove dependent linear constraints
            [Aeq, beq] = ComputationUtil.calcIndependentLinearConstraints(Aeq, beq);
        end
        
        function [c, ceq, gradientC, gradientCeq] = nlinConstraintFunc(this, paramVec)
            sim = this.simulation;
            g = sim.gene;
            t = sim.state('Time');
            r = sim.state('Rna');
            p = sim.state('ProteinMonomer');
            
            rnaLens       = r.lengths(r.matureIndexs);
            rnaBaseCounts = this.getRnaNMPCounts();
            rnaMWs        = r.molecularWeights(r.matureIndexs);
            monLens       = p.lengths(p.matureIndexs);
            monAACounts   = this.getMonomerAACounts();
            monDecayRates = p.decayRates(p.matureIndexs);
            
            matureRNAGeneComp = r.matureRNAGeneComposition;
            
            %extract
            [rnaExp, nmpComp, aaComp, rnaWtFracs] = this.extractParameterVector(paramVec);
            
            rnaWts = this.formulateRnaWtFractionConstraints();
            
            %inequality constraints and gradient
            c = [];
            if nargout >= 3
                gradientC = [];
            end
            
            %equality constraints
            %1. RNA weight fractions / RNA expression
            %2. RNA base counts, expresion / NMP composition
            %3. Protein AA counts, expression / AA composition
            ceq = [
                (rnaWts - rnaWtFracs * rnaMWs') * rnaExp
                (rnaBaseCounts' - nmpComp * rnaLens') * rnaExp
                (monAACounts' - aaComp * monLens') * ((matureRNAGeneComp(g.mRNAIndexs, :) * rnaExp) ./ (log(2) / t.cellCycleLength + monDecayRates))
                ];
            
            %equality constraints gradient
            if nargout >= 4
                gradientCeqRNAExp = [
                    (rnaWts - rnaWtFracs * rnaMWs')
                    (rnaBaseCounts' - nmpComp * rnaLens')
                    (monAACounts' - aaComp * monLens') * (matureRNAGeneComp(g.mRNAIndexs, :) ./ repmat(log(2) / t.cellCycleLength + monDecayRates, 1, this.nRNAs))
                    ]';
                
                gradientCeqNMPComp = [
                    zeros(this.nRNAWtFracs, this.nNMPs)
                    -eye(this.nNMPs) * (rnaLens' * rnaExp)
                    zeros(this.nAAs, this.nNMPs)
                    ]';
                
                gradientCeqAAComp = [
                    zeros(this.nRNAWtFracs, this.nAAs)
                    zeros(this.nNMPs, this.nAAs)
                    -eye(this.nAAs) * (monLens' * ((matureRNAGeneComp(g.mRNAIndexs, :) * rnaExp) ./ (log(2) / t.cellCycleLength + monDecayRates)))
                    ]';
                
                gradientCeqRNADecayRates = [
                    zeros(this.nRNAWtFracs, this.nRNAs)
                    zeros(this.nNMPs, this.nRNAs)
                    zeros(this.nAAs, this.nRNAs)
                    ]';
                
                gradientCeqRnaWtFracs = [
                    -eye(this.nRNAWtFracs) * (rnaMWs' * rnaExp)
                    zeros(this.nNMPs, this.nRNAWtFracs)
                    zeros(this.nAAs, this.nRNAWtFracs)
                    ]';
                
                gradientCeq = [
                    gradientCeqRNAExp
                    gradientCeqNMPComp
                    gradientCeqAAComp
                    gradientCeqRNADecayRates
                    gradientCeqRnaWtFracs
                    ];
            end
        end
        
        %Hessian of non-linear objective function and constraints.
        function hessian = hessFunc(this, ~, lambda, H)
            sim = this.simulation;
            
            g = sim.gene;
            t = sim.state('Time');
            r = sim.state('Rna');
            p = sim.state('ProteinMonomer');
            
            rnaLens = r.lengths(r.matureIndexs);
            rnaMWs = r.molecularWeights(r.matureIndexs);
            monLens  = p.lengths(p.matureIndexs);
            monDecayRates = p.decayRates(p.matureIndexs);
            
            %objective function
            hessian = 2 * H;
            
            %RNA weight fractions
            for i = 1:this.nRNAWtFracs
                tempHessian = zeros(size(hessian));
                tempHessian(this.rnaExpIdxs, this.rnaWtFracIdxs(i)) = -rnaMWs;
                tempHessian(this.rnaWtFracIdxs(i), this.rnaExpIdxs) = -rnaMWs;
                
                hessian = hessian + lambda.eqnonlin(i) * tempHessian;
            end
            
            %NMPs
            for i = 1:this.nNMPs
                tempHessian = zeros(size(hessian));
                tempHessian(this.rnaExpIdxs, this.nmpIdxs(i)) = -rnaLens;
                tempHessian(this.nmpIdxs(i), this.rnaExpIdxs) = -rnaLens;
                
                hessian = hessian + lambda.eqnonlin(i + this.nRNAWtFracs) * tempHessian;
            end
            
            %AAs
            for i = 1:this.nAAs
                tempHessian = zeros(size(hessian));
                tempHessian(this.rnaExpIdxs, this.aaIdxs(i)) = -(r.matureRNAGeneComposition(g.mRNAIndexs, :) ...
                    ./ repmat(log(2) / t.cellCycleLength + monDecayRates, 1, this.nRNAs))' * monLens;
                tempHessian(this.aaIdxs(i), this.rnaExpIdxs) = -(r.matureRNAGeneComposition(g.mRNAIndexs, :) ...
                    ./ repmat(log(2) / t.cellCycleLength + monDecayRates, 1, this.nRNAs))' * monLens;
                
                hessian = hessian + lambda.eqnonlin(i + this.nRNAWtFracs + this.nNMPs) * tempHessian;
            end
        end
        
        function W = hessMultFunc(~, Hinfo, Y)
            W = Hinfo * Y;
        end
        
        %Modifies RNA expression, NMP and AA composition, and rRNA half
        %lives. Doesn't modify RNA weight fractions or m/s/tRNA half lives.
        %Assumes uniform protein half lives.
        function paramVec = fitHeuristically(this, paramVec)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %references
            sim = this.simulation;
            g = sim.gene;
            t = sim.state('Time');
            mass = sim.state('Mass');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            sc = sim.process('DNASupercoiling');
            
            %constants
            rnaMWs = r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro;
            monMWs = pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro;
            matureRNAGeneComposition = r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs);
            matureRNAGeneComposition(matureRNAGeneComposition == 0) = NaN;
            
            %0. extract
            [rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, ...
                unaccountedEnergyConsumption] = ...
                this.extractParameterVector(paramVec);
            
            %1. match half lives, expression
            rnaDecayRates(r.matureRRNAIndexs) = mean(rnaDecayRates(r.matureRRNAIndexs));
            rnaExp = (r.nascentRNAMatureRNAComposition * ...
                ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition) * ...
                (rnaExp .* (log(2) / t.cellCycleLength + rnaDecayRates))) ./ (log(2) / t.cellCycleLength + rnaDecayRates);
            
            %2. match gene expression, weight fractions
            rnaCnt = zeros(size(r.matureIndexs));
            rnaCnt(r.matureMRNAIndexs) = ...
                rnaExp(r.matureMRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.mRNAWeightFractionIndexs) / (rnaMWs(r.matureMRNAIndexs)' * rnaExp(r.matureMRNAIndexs));
            rnaCnt(r.matureRRNAIndexs) = ...
                rnaExp(r.matureRRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                sum(rnaWtFracs(r.rRNAWeightFractionIndexs)) / (rnaMWs(r.matureRRNAIndexs)' * rnaExp(r.matureRRNAIndexs));
            rnaCnt(r.matureSRNAIndexs) = ...
                rnaExp(r.matureSRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.sRNAWeightFractionIndexs) / (rnaMWs(r.matureSRNAIndexs)' * rnaExp(r.matureSRNAIndexs));
            rnaCnt(r.matureTRNAIndexs) = ...
                rnaExp(r.matureTRNAIndexs) * ...
                mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                rnaWtFracs(r.tRNAWeightFractionIndexs) / (rnaMWs(r.matureTRNAIndexs)' * rnaExp(r.matureTRNAIndexs));
            
            monCnt = (r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * rnaCnt(r.matureMRNAIndexs)) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            monCnt = mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * monCnt / (monMWs' * monCnt);
            
            %3. Find gene expression, RNA decay rates, RNA weight fractions,
            %   nmp composition, aa composition which  satisfy gene expression lower bounds
            for i = 1:this.maxIter
                %print status
                if this.verbosity > 0
                    fprintf('Heuristic Iter: %d ...', i);
                end
                
                %transcription unit expression bounds
                rnaExp = rnaCnt / sum(rnaCnt);
                [biomassComposition, biomassProduction, byproducts, unaccountedEnergyConsumption, ...
                    minRNACnt, minMonCnt, maxRNACnt, maxMonCnt] = ...
                    this.calcResourceRequirements(this.constructParameterVector(...
                    rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, biomassComposition, biomassProduction, byproducts, ...
                    unaccountedEnergyConsumption));
                
                %gryase/topoisomerase I equilibrium
                maxSupercoilingActivity = max(...
                    max(monCnt(this.monomerIndex_topoisomeraseI), minMonCnt(this.monomerIndex_topoisomeraseI)) * (sc.topoIDeltaLK * sc.topoIActivityRate), ...
                    min(max(monCnt(this.monomerIndex_gyrase), minMonCnt(this.monomerIndex_gyrase))) * (-sc.gyraseDeltaLK * sc.gyraseActivityRate) / this.gryStoichiometry);
                maxMonCnt(this.monomerIndex_topoisomeraseI) = maxSupercoilingActivity / (sc.topoIDeltaLK * sc.topoIActivityRate);
                maxMonCnt(this.monomerIndex_gyrase) = maxSupercoilingActivity / (-sc.gyraseDeltaLK * sc.gyraseActivityRate) * this.gryStoichiometry;
                minMonCnt(this.monomerIndex_topoisomeraseI) = maxMonCnt(this.monomerIndex_topoisomeraseI);
                minMonCnt(this.monomerIndex_gyrase) = maxMonCnt(this.monomerIndex_gyrase);
                
                %stricter bounds
                maxRNACnt = min(maxRNACnt, mass.cellInitialDryWeight * mass.dryWeightFractionRNA     ./ rnaMWs);
                maxMonCnt = min(maxMonCnt, mass.cellInitialDryWeight * mass.dryWeightFractionProtein ./ monMWs);
                
                %check that upper and lower bounds don't conflict
                if any(minMonCnt ./ maxMonCnt > 1 + this.tolerance)
                    throw(MException('FitConstants:error', 'Minimum expression exceeds maximum expression'));
                end
                
                %apply robustness
                minRNACnt(r.matureRRNAIndexs) = this.calcResourceRequirementsRobustness(minRNACnt(r.matureRRNAIndexs), r.geneExpressionRobustness);
                minRNACnt(r.matureSRNAIndexs) = this.calcResourceRequirementsRobustness(minRNACnt(r.matureSRNAIndexs), r.geneExpressionRobustness);
                minRNACnt(r.matureTRNAIndexs) = max(r.minTRnaCnt, minRNACnt(r.matureTRNAIndexs));
                minRNACnt(r.matureTRNAIndexs) = this.calcResourceRequirementsRobustness(minRNACnt(r.matureTRNAIndexs), 2 * r.geneExpressionRobustness);
                minMonCnt(minMonCnt > 0) = this.calcResourceRequirementsRobustness(minMonCnt(minMonCnt > 0), r.geneExpressionRobustness, pm.minimumAverageExpression);
                minMonCnt(minMonCnt > 0) = max(minMonCnt(minMonCnt > 0), max(sum(pc.proteinComplexComposition(g.mRNAIndexs(minMonCnt > 0), :, :), 3), [], 2) * pc.minimumAverageExpression);
                
                minRNACnt = min(minRNACnt, maxRNACnt);
                minMonCnt = min(minMonCnt, maxMonCnt);
                
                %ribosomal proteins
                minMonCnt = max(minMonCnt, sum(pc.proteinComplexComposition(g.mRNAIndexs, pc.ribosome70SIndexs, :), 3) * 200);
                
                %DnaA, FtsZ
                minMonCnt(this.monomerIndex_DnaA) = this.dnaACnt;
                minMonCnt(this.monomerIndex_FtsZ) = this.ftsZCnt;
                maxMonCnt(this.monomerIndex_DnaA) = this.dnaACnt;
                maxMonCnt(this.monomerIndex_FtsZ) = this.ftsZCnt;
                
                if i > 1 && ...
                        ~any((rnaCnt ./ minRNACnt) < (1 - this.tolerance)) && ...
                        ~any((monCnt ./ minMonCnt) < (1 - this.tolerance))
                    if this.verbosity > 0
                        fprintf(' converged.\n');
                    end
                    break;
                end
                
                if this.verbosity > 0
                    fprintf(' RNA constraint violations: %d, Protein constraint violations: %d.\n', ...
                        sum((rnaCnt ./ minRNACnt) < (1 - this.tolerance)), ...
                        sum((monCnt ./ minMonCnt) < (1 - this.tolerance)));
                    
                    if this.verbosity > 1
                        ratio = rnaCnt ./ minRNACnt;
                        idxs = find(ratio < 1 - this.tolerance);
                        for j = 1:numel(idxs)
                            fprintf('\t%s\t%.3e\n', r.wholeCellModelIDs{r.matureIndexs(idxs(j))}, ratio(idxs(j)));
                        end
                        
                        ratio = monCnt ./ minMonCnt;
                        idxs = find(ratio < 1 - this.tolerance);
                        for j = 1:numel(idxs)
                            fprintf('\t%s\t%.3e\n', pm.wholeCellModelIDs{pm.matureIndexs(idxs(j))}(1:6), ratio(idxs(j)));
                        end
                    end
                end
                
                %monomers
                adjMinMonCnt = (r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * ...
                    max(matureRNAGeneComposition .* ...
                    repmat(minMonCnt .* (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs)), 1, numel(r.matureMRNAIndexs)) , [], 1)') ./ ...
                    (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
                adjMaxMonCnt = (r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * ...
                    min(matureRNAGeneComposition .* ...
                    repmat(maxMonCnt .* (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs)), 1, numel(r.matureMRNAIndexs)), [], 1)') ./ ...
                    (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
                adjMaxMonCnt = max(adjMinMonCnt, adjMaxMonCnt);
                
                if any(adjMinMonCnt > adjMaxMonCnt) && any(adjMinMonCnt ./ minMonCnt < 1 - this.tolerance) && any(adjMaxMonCnt ./ maxMonCnt > 1 + this.tolerance)
                    throw(MException('FitConstants:error', 'Cannot satisfy constriants'));
                end
                
                if mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein < adjMinMonCnt' * monMWs || ...
                        mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein > adjMaxMonCnt' * monMWs
                    throw(MException('FitConstants:error', 'Minimum monomer expression greater than protein weight fraction'))
                end
                if any(monCnt < adjMinMonCnt)
                    monCnt = adjMinMonCnt + ...
                        max(0, monCnt - adjMinMonCnt) * ...
                        (mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein - adjMinMonCnt' * monMWs) / (monMWs' * max(0, monCnt - adjMinMonCnt));
                end
                if any(monCnt > adjMaxMonCnt)
                    monCnt = adjMaxMonCnt - ...
                        min(0, monCnt - adjMaxMonCnt) * ...
                        (adjMaxMonCnt' * monMWs - mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein) / (monMWs' * min(0, monCnt - adjMaxMonCnt));
                end
                
                %mRNA
                mrnaExp = ComputationUtil.invertCompositionMatrix(r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs)) * ...
                    (monCnt .* (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs)));
                rnaCnt(r.matureMRNAIndexs) = mrnaExp * ...
                    mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * ...
                    rnaWtFracs(r.mRNAWeightFractionIndexs) / (rnaMWs(r.matureMRNAIndexs)' * mrnaExp);
                
                %rRNA
                if any(minRNACnt(r.matureRRNAIndexs) > rnaCnt(r.matureRRNAIndexs))
                    throw(MException('FitConstants:error', 'Minimum rRNA expression greater than rRNA weight fraction'));
                end
                
                %s/tRNA
                %TODO: handle maxRNACnt
                minMatureSRNAProd = minRNACnt(r.matureSRNAIndexs) .* ...
                    (log(2) / t.cellCycleLength + rnaDecayRates(r.matureSRNAIndexs));
                minMatureSRNACnt = ...
                    (r.nascentRNAMatureRNAComposition(r.matureSRNAIndexs, :) * ...
                    (max(r.nascentRNAMatureRNAComposition(r.matureSRNAIndexs, :) .* ...
                    minMatureSRNAProd(:, ones(size(r.nascentRNAMatureRNAComposition, 2),1)), [], 1))') ./ ...
                    (log(2) / t.cellCycleLength + rnaDecayRates(r.matureSRNAIndexs));
                matureSRNACnt = rnaCnt(r.matureSRNAIndexs);
                if mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * rnaWtFracs(r.sRNAWeightFractionIndexs) < minMatureSRNACnt' * rnaMWs(r.matureSRNAIndexs)
                    throw(MException('FitConstants:error', 'Minimum sRNA expression greater than sRNA weight fraction'))
                end
                matureSRNACnt = minMatureSRNACnt + ...
                    max(0, matureSRNACnt - minMatureSRNACnt) * ...
                    (mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * rnaWtFracs(r.sRNAWeightFractionIndexs) - minMatureSRNACnt' * rnaMWs(r.matureSRNAIndexs)) / ...
                    (rnaMWs(r.matureSRNAIndexs)' * max(0, matureSRNACnt - minMatureSRNACnt));
                rnaCnt(r.matureSRNAIndexs) = matureSRNACnt;
                
                minMatureTRNAProd = minRNACnt(r.matureTRNAIndexs) .* ...
                    (log(2) / t.cellCycleLength + rnaDecayRates(r.matureTRNAIndexs));
                minMatureTRNACnt = ...
                    (r.nascentRNAMatureRNAComposition(r.matureTRNAIndexs, :) * ...
                    (max(r.nascentRNAMatureRNAComposition(r.matureTRNAIndexs, :) .* ...
                    minMatureTRNAProd(:, ones(size(r.nascentRNAMatureRNAComposition, 2),1)), [], 1))') ./ ...
                    (log(2) / t.cellCycleLength + rnaDecayRates(r.matureTRNAIndexs));
                matureTRNACnt = rnaCnt(r.matureTRNAIndexs);
                if mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * rnaWtFracs(r.tRNAWeightFractionIndexs) < minMatureTRNACnt' * rnaMWs(r.matureTRNAIndexs)
                    throw(MException('FitConstants:error', 'Minimum tRNA expression greater than tRNA weight fraction'))
                end
                matureTRNACnt = minMatureTRNACnt + ...
                    max(0, matureTRNACnt - minMatureTRNACnt) * ...
                    (mass.cellInitialDryWeight * mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * rnaWtFracs(r.tRNAWeightFractionIndexs) - minMatureTRNACnt' * rnaMWs(r.matureTRNAIndexs)) / ...
                    (rnaMWs(r.matureTRNAIndexs)' * max(0, matureTRNACnt - minMatureTRNACnt));
                rnaCnt(r.matureTRNAIndexs) = matureTRNACnt;
                
                %NMP composition
                nmpComp = this.getRnaNMPCounts()' * rnaCnt;
                nmpComp = nmpComp / sum(nmpComp);
                
                %AA composition
                aaComp = this.getMonomerAACounts()' * monCnt;
                aaComp = aaComp / sum(aaComp);
                
                %calculate total monomers, and check within range
                [~, ~, ~, ~, ~, ~, totMons] = this.calcMacromolecularCounts(...
                    this.constructParameterVector(rnaCnt / sum(rnaCnt), nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                    biomassComposition, biomassProduction, byproducts, unaccountedEnergyConsumption));
                
                monExp = (r.matureRNAGeneComposition(g.mRNAIndexs, r.matureMRNAIndexs) * mrnaExp) ./ ...
                    (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
                monCnt = monExp * mass.cellInitialDryWeight * mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein / (monExp' * monMWs);
                
                if any(totMons ./ adjMaxMonCnt > 1 + this.tolerance) || any(totMons ./ adjMinMonCnt < 1 - this.tolerance)
                    throw(MException('FitConstants:error', 'Expression doesn''t obey constraints'))
                end
            end
            
            if i == this.maxIter
                throw(MException('Simulation:error', 'Gene expression bounds cannot be satisfied'));
            end
            
            %commit
            paramVec = this.constructParameterVector(...
                rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, unaccountedEnergyConsumption);
        end
        
        %Calculates:
        %Quantity              Units                                  Dimensions                     Normalization
        %===================   ====================================   ============================   ===========================================================
        %biomass composition   molecules / cell at life cycle start   nMetabolites x nCompartments           bmCmmp' * effMetMWs = cell mass at life cycle start
        %biomass production    molecules / cell life cycle            nMetabolites x nCompartments   (bmProd - byProd)' * metMWs = cell mass at life cycle start
        %by products           molecules / cell life cycle            nMetabolites x nCompartments   (bmProd - byProd)' * metMWs = cell mass at life cycle start
        %
        %Also calculates lower bounds on expression of key RNAs and
        %proteins For each functional gene product unit (eg. complex,
        %monomer, tRNA)
        %  1. Calculate minimum number of functional units
        %     needed per initial cell
        %  2. Lookup stoichiometry of gene products (sRNA, tRNA, rRNA,
        %     protein monomers) in functional unit
        %  3. Calculate RNA and protein monomer expression necessary to
        %     achieve number calculated in step (1).
        function [lclBmComp, lclBmProd, lclByProd, unaccECons, minRNAExp, minMonExp, maxRNAExp, maxMonExp] = ...
                calcResourceRequirements(this, paramVec)
            %% import classes
            import edu.stanford.covert.cell.sim.util.MetaboliteUtil
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = this.simulation;
            g = sim.gene;
            mass = sim.state('Mass');
            t = sim.state('Time');
            geom = sim.state('Geometry');
            m = sim.state('Metabolite');
            c = sim.state('Chromosome');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            mr = sim.state('MetabolicReaction');
            dnaDamage = sim.process('DNADamage');
            met = sim.process('Metabolism');
            
            %% store simulation constants
            paramVec0 = this.constructParameterVectorFromSimulation();
            
            %% constants
            pcComp = sum(pc.proteinComplexComposition, 3);
            [~, ~, ~, ~, rnaDecayRates] = this.extractParameterVector(paramVec);
            
            %% initialize state
            this.applyParameterVectorToSimulation(paramVec);
            initialGrowthFilterWidth = mr.initialGrowthFilterWidth;
            mr.initialGrowthFilterWidth = Inf;
            sim.initializeState();
            mr.initialGrowthFilterWidth = initialGrowthFilterWidth;
            
            %% Macromolecule counts
            [freeRnas, freeMons, freeCpxs, rnaInCpxs, monInCpxs, ...
                ~, ~, ~, rnaDecays, monDecays, cpxDecays, ...
                freeNTPs, freeAAs] = ...
                this.calcMacromolecularCounts(paramVec);
            
            %% biomass composition
            bmComp = zeros(numel(m.wholeCellModelIDs), 1);
            
            %water
            waterComp = mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) * mass.fractionWetWeight  * ...
                ConstantUtil.nAvogadro / m.molecularWeights(m.waterIndexs);
            bmComp(m.waterIndexs) = waterComp;
            
            %DNA
            if mass.dryWeightFractionDNA
                dnmpComp = [getBaseCounts(c.sequence); nnz(c.damagedBases == m.m6ADIndexs)];
                dnmpComp(1) = dnmpComp(1) - nnz(c.damagedBases == m.m6ADIndexs);
                chrWt = dnmpComp' * (m.molecularWeights([m.dnmpIndexs; m.getIndexs({'m6dAMP'})]) - ...
                    ConstantUtil.elements.O - ConstantUtil.elements.H) / ConstantUtil.nAvogadro;
                bmComp([m.dnmpIndexs; m.getIndexs({'m6dAMP'})]) = ...
                    bmComp([m.dnmpIndexs; m.getIndexs({'m6dAMP'})]) + ...
                    dnmpComp * min(1, mass.dryWeightFractionDNA * mass.cellInitialDryWeight / chrWt);
                
                dnmpComp = getBaseCounts(c.sequence);
                bmComp(m.dntpIndexs) = bmComp(m.dntpIndexs) + ...
                    dnmpComp * max(0, mass.dryWeightFractionDNA * mass.cellInitialDryWeight - chrWt) / ...
                    (dnmpComp' * m.molecularWeights(m.dntpIndexs)) * ConstantUtil.nAvogadro;
            end
            
            %RNAs + Protein
            macroComp = ...
                + r.baseCounts(r.matureIndexs, :)' * freeRnas ...
                + pm.baseCounts(pm.matureIndexs, :)' * freeMons ...
                + pc.baseCounts(pc.matureIndexs, :)' * freeCpxs;
            bmComp = bmComp + macroComp;
            
            %Small molecules -- lipids, polyamines, carbohydrates,
            %vitamins, ions, nucleotides
            expBioCompMolFracs = m.experimentalBiomassCompositionMolFractions;
            expBioCompMolFracs([m.diphosphateIndexs; m.phosphateIndexs; m.hydrogenIndexs], :) = 0;
            
            V = mass.cellInitialDryWeight / (1 - mass.fractionWetWeight) / geom.density; %volume (L)
            
            smBmComp = zeros(numel(m.wholeCellModelIDs), 1);
            smBmComp(m.ntpIndexs) = freeNTPs;
            smBmComp(m.aminoAcidIndexs) = freeAAs;
            smBmComp(m.lipidIndexs) = ...
                expBioCompMolFracs(m.lipidIndexs) * ...
                mass.dryWeightFractionLipid / (expBioCompMolFracs(m.lipidIndexs)' * m.molecularWeights(m.lipidIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            smBmComp(m.polyamineIndexs) = ...
                expBioCompMolFracs(m.polyamineIndexs) * ...
                mass.dryWeightFractionPolyamine / (expBioCompMolFracs(m.polyamineIndexs)' * m.molecularWeights(m.polyamineIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            smBmComp(m.carbohydrateIndexs) = ...
                expBioCompMolFracs(m.carbohydrateIndexs) * ...
                mass.dryWeightFractionCarbohydrate / (expBioCompMolFracs(m.carbohydrateIndexs)' * m.molecularWeights(m.carbohydrateIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            smBmComp(m.vitaminIndexs) = ...
                expBioCompMolFracs(m.vitaminIndexs) * ...
                mass.dryWeightFractionVitamin / (expBioCompMolFracs(m.vitaminIndexs)' * m.molecularWeights(m.vitaminIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            smBmComp(m.ionIndexs) = ...
                expBioCompMolFracs(m.ionIndexs) * ...
                mass.dryWeightFractionIon / (expBioCompMolFracs(m.ionIndexs)' * m.molecularWeights(m.ionIndexs)) * ...
                mass.cellInitialDryWeight * ConstantUtil.nAvogadro;
            
            if mass.dryWeightFractionNucleotide > 0
                nucComp = zeros(size(smBmComp));
                nucComp(m.ntpIndexs([1 3]))  = 1e-3 * ConstantUtil.nAvogadro * V * m.meanNTPConcentration;
                nucComp(m.ndpIndexs([1 3]))  = 1e-3 * ConstantUtil.nAvogadro * V * m.meanNDPConcentration;
                nucComp(m.nmpIndexs)         = 1e-3 * ConstantUtil.nAvogadro * V * m.meanNMPConcentration;
                nucComp(m.phosphateIndexs)   = sum(nucComp(m.ndpIndexs));
                nucComp(m.diphosphateIndexs) = sum(nucComp(m.nmpIndexs));
                nucComp(m.hydrogenIndexs)    = sum(nucComp(m.ndpIndexs)) + sum(nucComp(m.nmpIndexs));
                smBmComp = smBmComp + nucComp;
                
                assert(abs(mass.dryWeightFractionNucleotide * mass.cellInitialDryWeight - ...
                    nucComp' * m.molecularWeights / ConstantUtil.nAvogadro) / ...
                    (mass.dryWeightFractionNucleotide * mass.cellInitialDryWeight) < 1e-2);
            end
            
            bmComp = bmComp + smBmComp;
            
            %% Production, byproducts, minimum/maximum gene expression
            %initialize
            bmProd = zeros(numel(m.wholeCellModelIDs), 1);
            byProd = zeros(numel(m.wholeCellModelIDs), 1);
            minRNAExp = zeros(size( r.matureIndexs));
            minMonExp = zeros(size(pm.matureIndexs));
            minCpxExp = zeros(size(pc.matureIndexs));
            maxRNAExp = Inf(size( r.matureIndexs));
            maxMonExp = Inf(size(pm.matureIndexs));
            maxCpxExp = Inf(size(pc.matureIndexs));
            m.processBiomassProduction = zeros(numel(m.wholeCellModelIDs), numel(sim.processes));
            m.processByproduct = zeros(numel(m.wholeCellModelIDs), numel(sim.processes));
            
            %constants
            constants = edu.stanford.covert.util.StructUtil.catstruct(...
                sim.getFixedConstants(), sim.getFittedConstants());
            constants.processes.DNADamage.substrateGlobalIndexs = dnaDamage.substrateGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalIndexs = dnaDamage.substrateStimulusGlobalIndexs;
            constants.processes.DNADamage.substrateStimulusCompartmentIndexs = dnaDamage.substrateStimulusCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusGlobalCompartmentIndexs = dnaDamage.substrateStimulusGlobalCompartmentIndexs;
            constants.processes.DNADamage.substrateStimulusLocalIndexs = dnaDamage.substrateStimulusLocalIndexs;
            
            %integrated state
            states = struct;
            states.rnas0               = freeRnas;
            states.rnas                = freeRnas * t.cellCycleLength / log(2);
            states.rnaDecays           = rnaDecays;
            states.rnaDecays0          = rnaDecays * log(2) / t.cellCycleLength;
            states.rnaProductions      = freeRnas + rnaInCpxs + rnaDecays;
            states.rnaProductions0     = (freeRnas + rnaInCpxs) * log(2) / t.cellCycleLength + (freeRnas .* rnaDecayRates);
            states.rnaProductions0(setdiff(1:end, r.matureMRNAIndexs)) = ...
                states.rnaProductions0(setdiff(1:end, r.matureMRNAIndexs)) + ...
                pcComp(setdiff(1:end, g.mRNAIndexs), :) * (freeCpxs .* pc.decayRates(pc.matureIndexs));
            
            states.monomers0           = freeMons;
            states.monomers            = freeMons * t.cellCycleLength / log(2);
            states.monomerDecays       = monDecays;
            states.monomerDecays0      = monDecays * log(2) / t.cellCycleLength;
            states.monomerProductions  = freeMons + monInCpxs + monDecays;
            states.monomerProductions0 = (freeMons + monInCpxs) * log(2) / t.cellCycleLength + (freeMons .* pm.decayRates(pm.matureIndexs));
            states.monomerProductions0 = ...
                states.monomerProductions0 + ...
                pcComp(g.mRNAIndexs, :) * (freeCpxs .* pc.decayRates(pc.matureIndexs));
            
            states.complexs0           = freeCpxs;
            states.complexs            = freeCpxs * t.cellCycleLength / log(2);
            states.complexDecays       = cpxDecays;
            states.complexDecays0      = cpxDecays * log(2) / t.cellCycleLength;
            states.complexProductions  = freeCpxs + cpxDecays;
            states.complexProductions0 = freeCpxs * log(2) / t.cellCycleLength + (freeCpxs .* pc.decayRates(pc.matureIndexs));
            
            %processes
            for i = 1:length(sim.processes)
                p = sim.processes{i};
                
                [tmpBiomassProd, tmpByproducts, tmpMinEnzExp, tmpMaxEnzExp] = p.calcResourceRequirements_LifeCycle(constants, states);
                
                gblIdxs = p.substrateMetaboliteGlobalIndexs(:, 1);
                lclIdxs = p.substrateMetaboliteLocalIndexs(:, 1);
                
                bmProd(gblIdxs) = bmProd(gblIdxs) + sum(tmpBiomassProd(lclIdxs, :), 2);
                byProd(gblIdxs) = byProd(gblIdxs) + sum(tmpByproducts(lclIdxs, :), 2);
                
                if ~isempty(tmpBiomassProd) && ~isempty(tmpByproducts)
                    m.processBiomassProduction(gblIdxs, i) = tmpBiomassProd(lclIdxs, 1);
                    m.processByproduct(gblIdxs, i) = tmpByproducts(lclIdxs, 1);
                end
                
                if ~isempty(p.enzymeRNALocalIndexs)
                    minRNAExp(p.enzymeRNAGlobalIndexs(:, 1)) = max(minRNAExp(p.enzymeRNAGlobalIndexs(:, 1)), tmpMinEnzExp(p.enzymeRNALocalIndexs));
                    maxRNAExp(p.enzymeRNAGlobalIndexs(:, 1)) = min(maxRNAExp(p.enzymeRNAGlobalIndexs(:, 1)), tmpMaxEnzExp(p.enzymeRNALocalIndexs));
                end
                if ~isempty(p.enzymeMonomerLocalIndexs)
                    minMonExp(p.enzymeMonomerGlobalIndexs(:, 1)) = max(minMonExp(p.enzymeMonomerGlobalIndexs(:, 1)), tmpMinEnzExp(p.enzymeMonomerLocalIndexs));
                    maxMonExp(p.enzymeMonomerGlobalIndexs(:, 1)) = min(maxMonExp(p.enzymeMonomerGlobalIndexs(:, 1)), tmpMaxEnzExp(p.enzymeMonomerLocalIndexs));
                end
                if ~isempty(p.enzymeComplexLocalIndexs)
                    minCpxExp(p.enzymeComplexGlobalIndexs(:, 1)) = max(minCpxExp(p.enzymeComplexGlobalIndexs(:, 1)), tmpMinEnzExp(p.enzymeComplexLocalIndexs));
                    maxCpxExp(p.enzymeComplexGlobalIndexs(:, 1)) = min(maxCpxExp(p.enzymeComplexGlobalIndexs(:, 1)), tmpMaxEnzExp(p.enzymeComplexLocalIndexs));
                end
            end
            
            %water
            bmProd(m.waterIndexs) = bmProd(m.waterIndexs) + waterComp;
            
            %small molecules
            bmProd = bmProd + smBmComp;
            
            %dark energy
            %- energy known to be necessary for growth, but
            %- unaccounted in model
            metProd = bmProd - byProd;
            accountedEnergy = sum([
                + 2 * metProd([m.ntpIndexs; m.dntpIndexs])
                + 1 * metProd([m.ndpIndexs; m.dndpIndexs])
                + 0 * metProd([m.nmpIndexs; m.dnmpIndexs])
                + 1 * metProd([m.diphosphateIndexs])
                + 0 * metProd([m.phosphateIndexs])
                ]);
            
            unaccECons = met.growthAssociatedMaintenance * ...
                ConstantUtil.nAvogadro / 1000 * mass.cellInitialDryWeight ...
                - accountedEnergy;
            
            bmProd(m.atpIndexs)       = unaccECons + bmProd(m.atpIndexs);
            bmProd(m.waterIndexs)     = unaccECons + bmProd(m.waterIndexs);
            byProd(m.adpIndexs)       = unaccECons + byProd(m.adpIndexs);
            byProd(m.phosphateIndexs) = unaccECons + byProd(m.phosphateIndexs);
            byProd(m.hydrogenIndexs)  = unaccECons + byProd(m.hydrogenIndexs);
            
            %% sort biomass composition, production, byproducts to correct compartments
            lclBmComp = MetaboliteUtil.localizeMetabolites(sim, bmComp);
            lclBmProd = MetaboliteUtil.localizeMetabolites(sim, bmProd);
            lclByProd = MetaboliteUtil.localizeMetabolites(sim, byProd);
            
            %% monomers, complexs in metabolism as substrates
            monIdxs = met.substrateMonomerGlobalIndexs(:, 1);
            cpxIdxs = met.substrateComplexGlobalIndexs(:, 1);
            minMonExp(monIdxs) = max(minMonExp(monIdxs), 1);
            minCpxExp(cpxIdxs) = max(minCpxExp(cpxIdxs), 1);
            
            %% calculate minimum, maximum RNA, monomer expression
            monInCpxs1 = find(any(pcComp(g.mRNAIndexs, isfinite(minCpxExp)),2));
            minRNAExp(setdiff(1:end, r.matureMRNAIndexs)) = max(...
                minRNAExp(setdiff(1:end, r.matureMRNAIndexs)), ...
                pcComp(setdiff(1:end, g.mRNAIndexs), :) * minCpxExp);
            minMonExp(monInCpxs1) = max(...
                minMonExp(monInCpxs1), ...
                pcComp(g.mRNAIndexs(monInCpxs1), isfinite(minCpxExp)) * minCpxExp(isfinite(minCpxExp)));
            
            monInCpxs2 = find(any(pcComp(g.mRNAIndexs, isfinite(maxCpxExp)),2));
            maxRNAExp(setdiff(1:end, r.matureMRNAIndexs)) = min(...
                maxRNAExp(setdiff(1:end, r.matureMRNAIndexs)), ...
                pcComp(setdiff(1:end, g.mRNAIndexs), :) * maxCpxExp);
            maxMonExp(monInCpxs2) = min(...
                maxMonExp(monInCpxs2), ...
                pcComp(g.mRNAIndexs(monInCpxs2), isfinite(maxCpxExp)) * maxCpxExp(isfinite(maxCpxExp)));
            
            maxMonExp(union(monInCpxs1, monInCpxs2)) = max(...
                minMonExp(union(monInCpxs1, monInCpxs2)),...
                maxMonExp(union(monInCpxs1, monInCpxs2)));
            
            %% Hold DnaA, FtsZ expression constant to avoid refitting
            minMonExp(this.monomerIndex_DnaA) = this.dnaACnt;
            minMonExp(this.monomerIndex_FtsZ) = this.ftsZCnt;
            
            maxMonExp(this.monomerIndex_DnaA) = this.dnaACnt;
            maxMonExp(this.monomerIndex_FtsZ) = this.ftsZCnt;
            
            %% reset simulation constants
            this.applyParameterVectorToSimulation(paramVec0);
        end
        
        function [freeRnas, freeMons, freeCpxs, ...
                rnaInCpxs, monInCpxs, ...
                totRnas, totMons, totCpxs, ...
                rnaDecays, monDecays, cpxDecays, freeNTPs, freeAAs] = ...
                calcMacromolecularCounts(this, paramVec)
            import edu.stanford.covert.util.ComputationUtil;
            import edu.stanford.covert.util.ConstantUtil;
            
            %% references
            sim = this.simulation;
            g = sim.gene;
            t = sim.state('Time');
            mass = sim.state('Mass');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            pc = sim.state('ProteinComplex');
            
            %% constants
            pcComp = sum(pc.proteinComplexComposition, 3);
            ignoreCpxIdxs = [
                find(pc.formationProcesses(pc.matureIndexs) == sim.processIndex('Metabolism'));
                pc.getIndexs({'RIBOSOME_30S'; 'RIBOSOME_30S_IF3'; 'RIBOSOME_50S'; 'RNA_POLYMERASE_HOLOENZYME'})];
            
            %% fitted constants
            [rnaExp, ~, ~, ~, rnaDecayRates] = this.extractParameterVector(paramVec);
            
            %% expected numbers of macromolecules
            totRnas = mass.initialFractionNTPsInRNAs * mass.dryWeightFractionRNA * mass.cellInitialDryWeight * rnaExp / ...
                (rnaExp' * r.molecularWeights(r.matureIndexs) / ConstantUtil.nAvogadro);
            freeNTPs = (1 - mass.initialFractionNTPsInRNAs) * mass.dryWeightFractionRNA * mass.cellInitialDryWeight * ...
                r.baseCounts(r.matureIndexs, m.nmpIndexs)' * rnaExp / ...
                ((r.baseCounts(r.matureIndexs, m.nmpIndexs)' * rnaExp)' * m.molecularWeights(m.ntpIndexs) / ConstantUtil.nAvogadro);
            
            monExp = (r.matureRNAGeneComposition(g.mRNAIndexs, :) * rnaExp) ./ ...
                (log(2) / t.cellCycleLength + pm.decayRates(pm.matureIndexs));
            totMons = mass.initialFractionAAsInMonomers * mass.dryWeightFractionProtein * mass.cellInitialDryWeight * ...
                monExp / (monExp' * pm.molecularWeights(pm.matureIndexs) / ConstantUtil.nAvogadro);
            freeAAs = (1 - mass.initialFractionAAsInMonomers) * mass.dryWeightFractionProtein * mass.cellInitialDryWeight * ...
                pm.baseCounts(pm.matureIndexs, m.aminoAcidIndexs)' * monExp / ...
                ((pm.baseCounts(pm.matureIndexs, m.aminoAcidIndexs)' * monExp)' * m.molecularWeights(m.aminoAcidIndexs) / ConstantUtil.nAvogadro);
            
            subunits = zeros(size(g.wholeCellModelIDs));
            subunits(g.mRNAIndexs) = totMons;
            subunits(g.rRNAIndexs) = totRnas(r.matureRRNAIndexs);
            subunits(g.sRNAIndexs) = totRnas(r.matureSRNAIndexs);
            subunits(g.tRNAIndexs) = totRnas(r.matureTRNAIndexs);
            totCpxs = min(repmat(subunits, 1, numel(pc.matureIndexs)) ./ pcComp, [], 1)';
            totCpxs(ignoreCpxIdxs) = 0;
            
            [subs2Nets, cpxs2Nets, nets] = edu.stanford.covert.util.findNonInteractingRowsAndColumns(...
                sum(pc.proteinComplexComposition, 3));
            for i = 1:numel(nets)
                if size(nets{i}, 2) <= 1
                    continue;
                end
                
                tmpSubIdxs = find(subs2Nets == i);
                tmpCpxIdxs = find(cpxs2Nets == i);
                tmpPcComp = pcComp(tmpSubIdxs, tmpCpxIdxs);
                tmpSubs = subunits(tmpSubIdxs);
                tmpCpxs = zeros(size(tmpCpxIdxs));
                while true
                    if all(sum(tmpCpxIdxs, 1) == sum(tmpCpxIdxs(:, 1)))
                        tmpRates = prod(repmat(tmpSubs / mean(tmpSubs), 1, numel(tmpCpxIdxs)) .^ tmpPcComp, 1)';
                    else
                        tmpRates = prod(repmat(tmpSubs / sum(totMons), 1, numel(tmpCpxIdxs)) .^ tmpPcComp, 1)';
                    end
                    tmpRates(ismember(tmpCpxIdxs, ignoreCpxIdxs)) = 0;
                    if ~any(tmpRates) || max(min(repmat(tmpSubs, 1, numel(tmpCpxIdxs)) ./ tmpPcComp, [], 1)) < 1e-3
                        break;
                    end
                    tmpCpxs = tmpCpxs + tmpRates * min(tmpSubs ./ (tmpPcComp * tmpRates));
                    tmpSubs = subunits(tmpSubIdxs) - tmpPcComp * tmpCpxs;
                end
                totCpxs(tmpCpxIdxs) = tmpCpxs;
            end
            
            totCpxs(pc.getIndexs({'MG_213_214_298_6MER_ADP'; 'MG_213_214_298_6MER'})) = ...
                + totCpxs(pc.getIndexs({'MG_213_214_298_6MER_ADP'; 'MG_213_214_298_6MER'})) ...
                + [1; -1] * totCpxs(pc.getIndexs('MG_213_214_298_6MER'));
            
            idxs = [sim.process('FtsZPolymerization').enzymeComplexGlobalIndexs; pc.getIndexs('MG_224_9MER_GDP')];
            idx = sim.process('FtsZPolymerization').enzymeGlobalIndexs(sim.process('FtsZPolymerization').enzymeIndexs_FtsZ_GTP);
            totCpxs(idx) = sum(totCpxs(idxs));
            totCpxs(setdiff(idxs, idx)) = 0;
            
            idxs = sim.process('ReplicationInitiation').enzymeComplexGlobalIndexs;
            idx = sim.process('ReplicationInitiation').enzymeGlobalIndexs(sim.process('ReplicationInitiation').enzymeIndexs_DnaA_1mer_ATP);
            totCpxs(idx) = sum(totCpxs(idxs));
            totCpxs(setdiff(idxs, idx)) = 0;
            
            %RNA, monomers free / in complexes
            rnaInCpxs = zeros(size(r.matureIndexs));
            rnaInCpxs(setdiff(1:end, r.matureMRNAIndexs)) = ...
                pcComp(setdiff(1:end, g.mRNAIndexs), :) * totCpxs;
            assert(all(totRnas - rnaInCpxs > -1e-8));
            freeRnas = max(0, totRnas - rnaInCpxs);
            
            monInCpxs = pcComp(g.mRNAIndexs, :) * totCpxs;
            assert(all(totMons - monInCpxs > -1e-8));
            freeMons = max(0, totMons - monInCpxs);
            
            freeCpxs = totCpxs;
            
            f = max(0, ...
                (totRnas' * r.molecularWeights(r.matureIndexs) + totMons' * pm.molecularWeights(pm.matureIndexs)) / ...
                (freeRnas' * r.molecularWeights(r.matureIndexs) + freeMons' * pm.molecularWeights(pm.matureIndexs) + freeCpxs' * pc.molecularWeights(pc.matureIndexs)));
            freeRnas = f * freeRnas;
            freeMons = f * freeMons;
            freeCpxs = f * freeCpxs;
            rnaInCpxs = f * rnaInCpxs;
            monInCpxs = f * monInCpxs;
            totRnas = f * totRnas;
            totMons = f * totMons;
            totCpxs = f * totCpxs;
            
            %decays
            cpxDecays = freeCpxs .* pc.decayRates(pc.matureIndexs) * t.cellCycleLength / log(2);
            rnaDecays = freeRnas .* rnaDecayRates * t.cellCycleLength / log(2);
            rnaDecays(setdiff(1:end, r.matureMRNAIndexs)) = ...
                rnaDecays(setdiff(1:end, r.matureMRNAIndexs)) + ...
                pcComp(setdiff(1:end, g.mRNAIndexs), :) * cpxDecays;
            monDecays = freeMons .* pm.decayRates(pm.matureIndexs) * t.cellCycleLength / log(2) + ...
                pcComp(g.mRNAIndexs, :) * cpxDecays;
        end
        
        %Robustness
        %  Set lower bound on average expression equal to calculate
        %  lowwer bound plus nSigma standard deviations so that
        %  probability expression falls below calculated lower bound
        %  (lb) is small (~5%).
        %
        %  Since distribution is approximately poisson, assume variance
        %  approximately equal to mean. Solve for mean such that
        %  probability that quantity falls less than lb is small (~5%).
        %  That is solve for avg:
        %     lb = avg - nSigma * sqrt(avg)
        %    avg = lb + nSigma^2/2 + sqrt(1/4*(2*lb+nSigma^2)^2 - lb^2)
        function exp = calcResourceRequirementsRobustness(~, exp, robustness, minExp)
            exp = exp + robustness^2/2 + sqrt(1/4*(2*exp + robustness^2).^2 - exp.^2);
            
            if nargin >= 4
                exp = max(exp, minExp);
            end
        end
        
        function [rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, ...
                unaccountedEnergyConsumption] = ...
                extractParameterVector(this, vec)
            if nargout <= 5 && isnumeric(vec)
                vec = {vec, [], [], [], []};
            end
            rnaExp        = vec{1}(this.rnaExpIdxs, 1);
            nmpComp       = vec{1}(this.nmpIdxs, 1);
            aaComp        = vec{1}(this.aaIdxs, 1);
            rnaWtFracs    = vec{1}(this.rnaWtFracIdxs, 1);
            rnaDecayRates = vec{1}(this.rnaDecayRateIdxs, 1);
            biomassComposition = vec{2};
            biomassProduction = vec{3};
            byproducts = vec{4};
            unaccountedEnergyConsumption = vec{5};
        end
        
        function vec = constructParameterVector(this, ...
                rnaExp, nmpComp, aaComp, rnaWtFracs, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, unaccountedEnergyConsumption)
            vec = cell(5, 1);
            vec{1} = zeros(this.nVars, 1);
            vec{1}(this.rnaExpIdxs, 1)       = rnaExp;
            vec{1}(this.nmpIdxs, 1)          = nmpComp;
            vec{1}(this.aaIdxs, 1)           = aaComp;
            vec{1}(this.rnaWtFracIdxs, 1)    = rnaWtFracs;
            vec{1}(this.rnaDecayRateIdxs, 1) = rnaDecayRates;
            vec{2} = biomassComposition;
            vec{3} = biomassProduction;
            vec{4} = byproducts;
            vec{5} = unaccountedEnergyConsumption;
        end
        
        function paramVec = constructParameterVectorFromSimulation(this)
            %references
            sim = this.simulation;
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            met = sim.process('Metabolism');
            
            %data
            paramVec = this.constructParameterVector(...
                r.expression(r.matureIndexs), ...
                sum(m.biomassComposition(m.nmpIndexs, :), 2) / sum(sum(m.biomassComposition(m.nmpIndexs, :), 2)), ...
                sum(m.biomassComposition(m.aminoAcidIndexs, :), 2) / sum(sum(m.biomassComposition(m.aminoAcidIndexs, :), 2)), ...
                this.formulateRnaWtFractionConstraints * r.expression(r.matureIndexs) / sum(this.formulateRnaWtFractionConstraints * r.expression(r.matureIndexs)), ...
                r.decayRates(r.matureIndexs), ...
                m.biomassComposition, ...
                m.biomassProduction, ...
                m.byproducts, ...
                met.unaccountedEnergyConsumption);
        end
        
        function applyParameterVectorToSimulation(this, paramVec)
            %import classes
            import edu.stanford.covert.util.ComputationUtil;
            
            %references
            sim = this.simulation;
            g = sim.gene;
            time = sim.state('Time');
            c = sim.state('Chromosome');
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pc = sim.state('ProteinComplex');
            t = sim.process('Transcription');
            met = sim.process('Metabolism');
            tr = sim.process('TranscriptionalRegulation');
            sc = sim.process('DNASupercoiling');
            
            %data
            invNascentRNAMatureRNAComp = ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition);
            
            %extract
            [rnaExp, ~, ~, ~, rnaDecayRates, ...
                biomassComposition, biomassProduction, byproducts, ...
                unaccountedEnergyConsumption] = ...
                this.extractParameterVector(paramVec);
            
            % allocate object space
            r.decayRates = zeros(numel(r.wholeCellModelIDs), 1);
            r.expression = zeros(numel(r.wholeCellModelIDs), 1);
            
            %RNA decay rates
            r.decayRates(r.nascentIndexs)       = invNascentRNAMatureRNAComp * rnaDecayRates;
            r.decayRates(r.processedIndexs)     = rnaDecayRates;
            r.decayRates(r.intergenicIndexs)    = realmax;
            r.decayRates(r.matureIndexs)        = rnaDecayRates;
            r.decayRates(r.boundIndexs)         = 0;
            r.decayRates(r.misfoldedIndexs)     = rnaDecayRates;
            r.decayRates(r.damagedIndexs)       = realmax;
            r.decayRates(r.aminoacylatedIndexs) = rnaDecayRates;
            
            %RNA expression
            r.expression(r.nascentIndexs)       = 0;
            r.expression(r.processedIndexs)     = 0;
            r.expression(r.intergenicIndexs)    = 0;
            r.expression(r.matureIndexs)        = rnaExp;
            r.expression(r.boundIndexs)         = 0;
            r.expression(r.misfoldedIndexs)     = 0;
            r.expression(r.damagedIndexs)       = 0;
            r.expression(r.aminoacylatedIndexs) = 0;
            
            %compute RNA polymerase binding probability from gene-wise expression and
            %decay rates, and average across genes in each transcription unit
            t.transcriptionUnitBindingProbabilities = ...
                invNascentRNAMatureRNAComp * (rnaExp .* (log(2) / time.cellCycleLength + rnaDecayRates));
            
            %Deconvolve from RNA polymerase binding probability effects of
            %- active transcriptional regulators
            %- supercoiling
            %- transcription unit copy number variation over cell cycle
            if ~isempty(c.polymerizedRegions) && ~isempty(tr)
                tfFoldChange = tr.calcBindingProbabilityFoldChange(tr.bindTranscriptionFactors());
            else
                tfFoldChange = ones(size(c.transcriptionUnitStartCoordinates, 1), 2);
            end
            if ~isempty(sc)
                scProbFoldChange = sc.calcRNAPolymeraseBindingProbFoldChange([1 1], ...
                    c.sequenceLen, ...
                    c.sequenceLen / c.relaxedBasesPerTurn * (1 + c.equilibriumSuperhelicalDensity));
            else
                scProbFoldChange = ones(size(c.transcriptionUnitStartCoordinates, 1), 2);
            end
            
            t.transcriptionUnitBindingProbabilities = ...
                t.transcriptionUnitBindingProbabilities ...
                ./ tfFoldChange(:, 1) ...
                ./ scProbFoldChange(:, 1) ...
                ./ this.calcEffectiveMeanTranscriptionUnitCopyNumbers(this.simulation);
            
            %normalize RNA polymerase binding probability
            t.transcriptionUnitBindingProbabilities = ...
                t.transcriptionUnitBindingProbabilities / ...
                sum(t.transcriptionUnitBindingProbabilities);
            
            %update Metabolism process metabolismProduction and
            %fbaReactionStoichiometryMatrix properties
            m.biomassComposition = biomassComposition;
            m.biomassProduction  = biomassProduction;
            m.byproducts         = byproducts;
            met.formulateFBA(biomassProduction - byproducts, unaccountedEnergyConsumption);
            
            %calculate RNA polymerase state transition probabilities
            [freeRnas, ~, freeCpxs, rnaInCpxs, ~, ~, ~, totCpxs] = this.calcMacromolecularCounts(paramVec);
            rnaProds = (freeRnas + rnaInCpxs) * log(2) / time.cellCycleLength + (freeRnas .* rnaDecayRates);
            rnaProds(setdiff(1:end, r.matureMRNAIndexs)) = ...
                rnaProds(setdiff(1:end, r.matureMRNAIndexs)) + ...
                sum(pc.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3) * (freeCpxs .* pc.decayRates(pc.matureIndexs));
            ntpProd = (ComputationUtil.invertCompositionMatrix(r.nascentRNAMatureRNAComposition) * rnaProds)' * r.lengths(r.nascentIndexs);
            
            nPols = totCpxs(t.enzymeGlobalIndexs(t.enzymeIndexs_rnaPolymerase));
            t.stateTransitionProbabilities = t.calcStateTransitionProbabilities(...
                nPols, ntpProd, t.transcriptionUnitBindingProbabilities);
        end
        
        function baseCounts = getRnaNMPCounts(this)
            sim = this.simulation;
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            baseCounts = r.baseCounts(r.processedIndexs, m.nmpIndexs);
        end
        
        function aaCounts = getMonomerAACounts(this)
            sim = this.simulation;
            m = sim.state('Metabolite');
            p = sim.state('ProteinMonomer');
            pmod = sim.process('ProteinModification');
            
            aaCounts = p.baseCounts(p.processedIIIndexs, m.aminoAcidIndexs);
            aaCounts(:, m.getIndexs('CYS') == m.aminoAcidIndexs) = ...
                aaCounts(:, m.getIndexs('CYS') == m.aminoAcidIndexs) + ...
                p.baseCounts(p.processedIIIndexs, m.getIndexs('diacylglycerolCys'));
            [tfs, idxs] = ismember(pmod.substrateWholeCellModelIDs, m.wholeCellModelIDs(m.aminoAcidIndexs));
            aaCounts(:, idxs(tfs)) = ...
                aaCounts(:, idxs(tfs)) + ...
                pmod.reactionModificationMatrix' * max(0, -pmod.reactionStoichiometryMatrix(tfs, :))';
        end
        
        function tuConstraints = formulateTranscriptionUnitConstraints(this)
            rnaGeneComp = this.simulation.state('Rna').nascentRNAGeneComposition;
            tuConstraints  = zeros(this.nGenes - this.nRNAs, this.nGenes);
            k = 0;
            for j = 1:this.nRNAs
                iGene = find(rnaGeneComp(:, j));
                if numel(iGene) <= 1
                    continue;
                end
                
                tuConstraints(k+1:k+length(iGene)-1, iGene(1)) = 1;
                tuConstraints((k+1:k+length(iGene)-1) + (iGene(2:end)-1) * nTUConstraints) = -1;
                k = k + length(iGene) - 1;
            end
        end
        
        function rnaWts = formulateRnaWtFractionConstraints(this)
            sim = this.simulation;
            r = sim.state('Rna');
            
            rnaMWs = r.molecularWeights(r.matureIndexs);
            
            numRNAs   = length(r.matureIndexs);
            numMRNAs  = length(r.matureMRNAIndexs);
            numRRNAs  = length(r.matureRRNAIndexs);
            numSRNAs  = length(r.matureSRNAIndexs);
            numTRNAs  = length(r.matureTRNAIndexs);
            
            indicatorMRNA = zeros(numMRNAs, numRNAs);
            indicatorSRNA = zeros(numSRNAs, numRNAs);
            indicatorTRNA = zeros(numTRNAs, numRNAs);
            indicatorRibosomalRRNA = zeros(numRRNAs, numRNAs);
            indicatorMRNA((1:numMRNAs)' + (r.matureMRNAIndexs - 1) * numMRNAs) = 1;
            indicatorSRNA((1:numSRNAs)' + (r.matureSRNAIndexs - 1) * numSRNAs) = 1;
            indicatorTRNA((1:numTRNAs)' + (r.matureTRNAIndexs - 1) * numTRNAs) = 1;
            indicatorRibosomalRRNA((1:numRRNAs)' + (r.matureRibosomalRRNAIndexs - 1) * numRRNAs) = 1;
            
            %extract
            rnaWts = [
                rnaMWs(r.matureMRNAIndexs)'                * indicatorMRNA;
                diag(rnaMWs(r.matureRibosomalRRNAIndexs))  * indicatorRibosomalRRNA;
                rnaMWs(r.matureSRNAIndexs)'                * indicatorSRNA
                rnaMWs(r.matureTRNAIndexs)'                * indicatorTRNA];
        end
    end
    
    methods (Static = true)
        function value = calcEffectiveMeanTranscriptionUnitCopyNumbers(sim)
            c = sim.state('Chromosome');
            time = sim.state('Time');
            
            value = (...
                + 1 * time.replicationInitiationDuration ...
                + (1 + abs(c.transcriptionUnitStartCoordinates + c.transcriptionUnitLengths - 1 - c.sequenceLen/2) / (c.sequenceLen/2)) * time.replicationDuration ...
                + 2 * time.cytokinesisDuration ...
                ) / time.cellCycleLength;
            
            iTU = sim.process('Transcription').transcriptionUnitIndexs_DnaAR12345Boxes;
            value(iTU) = 1/3 * (value(iTU) - time.replicationInitiationDuration / time.cellCycleLength); %1/3 fit based on simulations of multiple generations
        end
    end
end
