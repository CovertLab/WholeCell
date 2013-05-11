%FtsZPolymerization
%
% @wholeCellModelID Process_FtsZPolymerization
% @name             FtsZPolymerization
% @description
%
% In the cytosol, FtsZ can exist in one of multiple states: inactivated monomer,
% deactivated monomer (GDP bound), activated monomer (GTP bound), nucleated
% (dimer of two activated monomers), or elongated (polymer of three or more
% activated monomers). FtsZ molecules move between these states at rates
% obtained from Chen et al. 2005 and Surovtsev et al. 2008.
%
% The maximum polymer length is a fittable parameter. We currently use 9 (40nm),
% the lower end of the range given for E. coli in Anderson 2004.
%
% FtsZ polymerization is modeled using a set of differential equations
% modified from that described in Surovtsev et al. 2008, involving the
% activation, nucleation, and elongation of FtsZ polymers.
% The main modifications were that the equations were simplified to not
% include annealing and cyclization of FtsZ polymers.
%
% Solving the equations results in a real-valued distribution of monomers and
% filament lengths at each time step. This process discretizes the distribution
% at each time step for compatibility with the rest of the simulation.
%
% References
% ===============
% 1. Surovtsev, I.V., Morgan, J.J., Lindahl, P.A. (2008). Kinetic Modelling of the
%    Assembly, Dynamic Steady State, and Contraction of the FtsZ Ring in
%    Prokaryotic Cytokinesis. Plos CB 4: 1-19. [PUB_0164]
% 2. Chen, Y., Bjornson, K., Redick, S.D., Erickson, H.P. (2005). A rapid
%    fluorescence assay for ftsZ assembly indicates cooperative assembly with a
%    dimer nucleus. Biophysical journal 88: 505-514. [PUB_0200]
% 3. Li, Z., Trimble, M.J., Brun, Y.V., Jensen, G.J. (2007). The structure of
%    FtsZ filaments in vivo suggests a force-generating role in cell division.
%    EMBO 26: 4694-4708. [PUB_0611]
% 4. Anderson, D.E., Gueiros-Filho, F.J., Erickson, H.P. (2004). Assembly
%    Dynamics of FtsZ Rings in Bacillus subtilis and Escherichia coli and
%    Effects of FtsZ-Regulating Proteins. Journal of Bacteriology 186:
%    5775-5781. [PUB_0217]
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/9/2010

classdef FtsZPolymerization < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'maxPolymerLength';
            'activationFwd';
            'activationRev';
            'exchangeFwd';
            'exchangeRev';
            'nucleationFwd';
            'nucleationRev';
            'elongationFwd';
            'elongationRev';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};     %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = {     %whole cell model IDs of substrates
            'GDP'; 'GTP'; 'PI'; 'H2O'; 'H'};
        substrateIndexs_gdp       = 1;
        substrateIndexs_gtp       = 2;
        substrateIndexs_phosphate = 3;
        substrateIndexs_water     = 4;
        substrateIndexs_hydrogen  = 5;

        enzymeWholeCellModelIDs = {     %whole cell model IDs of enzymes
            'MG_224_MONOMER';           %cell division protein FtsZ
            'MG_224_MONOMER_GDP';       %cell division protein FtsZ deactivated
            'MG_224_MONOMER_GTP';       %cell division protein FtsZ activated
            'MG_224_2MER_GTP';          %cell division protein FtsZ activated 2mer
            'MG_224_3MER_GTP';          %cell division protein FtsZ activated 3mer
            'MG_224_4MER_GTP';          %cell division protein FtsZ activated 4mer
            'MG_224_5MER_GTP';          %cell division protein FtsZ activated 5mer
            'MG_224_6MER_GTP';          %cell division protein FtsZ activated 6mer
            'MG_224_7MER_GTP';          %cell division protein FtsZ activated 7mer
            'MG_224_8MER_GTP';          %cell division protein FtsZ activated 8mer
            'MG_224_9MER_GTP'};         %cell division protein FtsZ activated 9mer
        enzymeIndexs_FtsZ           = 1;
        enzymeIndexs_FtsZ_GDP       = 2;
        enzymeIndexs_FtsZ_GTP       = 3;
        enzymeIndexs_FtsZ_dimer     = 4;
        enzymeIndexs_FtsZ_9mer      = 11;
        enzymeIndexs_FtsZ_activated = (3:11)';
    end

    %fixed biological constants
    properties
        maxPolymerLength      %number of activated FtsZ subunits in longest polymer
        activationFwd         % 1.1    1/s     kact1  [PUB_0200]
        activationRev         % 0.01   1/s     kact2  [PUB_0200]
        exchangeFwd           % 1e4    1/(Ms)  kexf   [PUB_0164]
        exchangeRev           % 5e3    1/(Ms)  kexr   [PUB_0164]
        nucleationFwd         % 4.2e6  1/(Ms)  knuc1  [PUB_0200]
        %nucleationRev        % 7700   1/s     knuc2  [PUB_0200]
        nucleationRev         % 40     1/s     knuc2  [PUB_0164]
        elongationFwd         % 5.1e6  1/(Ms)  kel1   [PUB_0200]
        elongationRev         % 2.9    1/s     kel2   [PUB_0200]
    end

    %constructor
    methods
        function this = FtsZPolymerization(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(varargin{:});
            this.maxPolymerLength = length(this.enzymeIndexs_FtsZ_activated);
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, constants, states)
            import import edu.stanford.covert.cell.sim.process.Cytokinesis;
            
            %initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %substrate and byproducts: FtsZ polymerization uses GTP to make
            %FtsZ-GTP filaments for cytokinesis. Hydrolysis of the bound GTP is
            %coupled to contraction in cytokinesis. Here we account for the
            %bound GTP and GDP. The free metabolites involved in the hydrolysis
            %of bound GTP to GDP are accounted for in Cytokinesis.
            [~, gtpCost] = Cytokinesis.calcRequiredPinchingCycles(this.geometry.width, ...
                constants.states.FtsZRing.numFtsZSubunitsPerFilament, constants.states.FtsZRing.filamentLengthInNm);
            gdpCost = gtpCost;
            
            gtpCost = ...
                + gtpCost ...
                + (1:9) * states.complexProductions(this.enzymeGlobalIndexs(this.enzymeIndexs_FtsZ_activated));
            gdpCost = ...
                + gdpCost ...
                - states.complexProductions(this.enzymeGlobalIndexs(this.enzymeIndexs_FtsZ_GDP)) ...
                - 9 * states.complexProductions(this.complex.getIndexs('MG_224_9MER_GDP'));
            
            bmProd(this.substrateIndexs_gtp) = gtpCost;
            byProd(this.substrateIndexs_gdp) = gdpCost;
            
            %enzymes: accounted for in Cytokinesis
        end    

        %initialization
        function initializeState(this)
            %Approach steady state, stop when converged
            maxSteps = 50; %max number of steps to take to try to reach a steady state
            tol = max(1, 0.002 * sum(this.enzymes));  %tolerance
            substrates = this.substrates;
            this.substrates(:) = 1e6;
            
            enzymes = this.enzymes;
            for iStep = 1:maxSteps
                this.evolveState();
                
                %check if enzymes converged
                if norm(enzymes - this.enzymes, Inf) < tol
                    break;
                end
                
                enzymes = this.enzymes;
            end
            
            this.substrates = substrates;
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_gtp) = ...
                + this.enzymes(this.enzymeIndexs_FtsZ) ...
                + this.enzymes(this.enzymeIndexs_FtsZ_GDP);
        end

        %simulation
        function evolveState(this)
            if ~any(this.enzymes)
                return;
            end
            
            %solve ODE
            [~, odeSolutions] = this.integrateODEs(this.moleculesToConcentration(this.enzymes));
            
            %preserve mass balance by counting fractional polymers as monomers
            enzymes = this.discretizeEnzymes(odeSolutions(:, find(all(odeSolutions >= 0, 1), 1, 'last')));
            
            %update substrate counts
            [this.enzymes, this.substrates] = this.applySubstrateLimits(enzymes, this.substrates);
        end
    end
    
    %model helper functions
    methods
        %integrates ODE equations using modified ODE23S method
        %
        %requirements
        %- tspan is increasing real 2 element vector with first element 0
        %- d^2Y/dt^2 = df/ft = 0
        function [tout, yout] = integrateODEs(this, y0)
            params = [
                this.activationFwd
                this.activationRev
                this.exchangeFwd
                this.exchangeRev
                this.nucleationFwd
                this.nucleationRev
                this.elongationFwd
                this.elongationRev
                this.moleculesToConcentration(this.substrates);
                ];

            % Handle solver arguments
            neq = numel(y0);
            tfinal = this.stepSizeSec;
            hmax = 0.1 * this.stepSizeSec;
            threshold = 0.1 / edu.stanford.covert.util.ConstantUtil.nAvogadro / this.geometry.volume;
            rtol = 1e-2;
            
            % Initialize method parameters.
            pow = 1/3;
            d = 1 / (2 + sqrt(2));
            e32 = 6 + sqrt(2);
            
            % Compute the initial slope yp.
            t = 0;
            y = y0;
            f0 = diff(y0, params);
            dfdy = jacobian(y, params);
            
            % Compute an initial step size h using yp = y'(t).
            wt = max(abs(y), threshold);
            
            % Compute y''(t) and a better initial step size.
            DfDt = dfdy*f0;
            rh = 1.25 * sqrt(0.5 * max(abs(DfDt ./ wt))) / rtol^pow;
            absh = 1/rh;
            
            % Allocate memory if we're generating output.
            chunk = 100;
            tout = zeros(1, chunk);
            yout = zeros(neq, chunk);
            nout = 1;
            tout(nout) = t;
            yout(:,nout) = y;
            
            % THE MAIN LOOP
            done = false;
            notFirstStep = false;
            while ~done
                hmin = 16*eps*t;
                absh = min(hmax, max(hmin, absh));
                h = absh;
                
                % Stretch the step if within 10% of tfinal-t.
                if 1.1*absh >= abs(tfinal - t)
                    h = tfinal - t;
                    absh = abs(h);
                    done = true;
                end
                
                if notFirstStep % J is already computed on first step
                    dfdy = jacobian(y, params);
                end
                
                % LOOP FOR ADVANCING ONE STEP.
                nofailed = true; % no failed attempts
                while true % Evaluate the formula.
                    Miter = eye(neq) - (h*d)*dfdy;  % sparse if dfdy is sparse
                    k1aux = f0;
                    
                    [L,U,p] = lu(Miter,'vector');
                    k1 = U \ (L \ k1aux(p));
                    f1 = diff(y + 0.5*h*k1, params);
                    k2 = (U \ (L \ (f1(p) - k1(p)))) + k1;
                    
                    tnew = t + h;
                    if done
                        tnew = tfinal; % Hit end point exactly.
                    end
                    h = tnew - t; % Purify h.
                    
                    ynew = y + h*k2;
                    f2 = diff(ynew, params);
                    k3aux = (f2 - e32*(k2 - f1) - 2*(k1 - f0));
                    k3 = U \ (L \ k3aux(p));
                    
                    % Estimate the error.
                    err = (absh/6) * ...
                        max(abs((k1-2*k2+k3) ./ max(max(abs(y),abs(ynew)),threshold)));
                    
                    % Accept the solution only if the weighted error is no more than the
                    % tolerance rtol.  Estimate an h that will yield an error of rtol on
                    % the next step or the next try at taking this step, as the case may be,
                    % and use 0.8 of this value to avoid failures.
                    if err > rtol % Failed step
                        if absh <= hmin
                            warning('MATLAB:ode23s:IntegrationTolNotMet',['Failure at t=%e.  ' ...
                                'Unable to meet integration tolerances without reducing ' ...
                                'the step size below the smallest value allowed (%e) ' ...
                                'at time t.'],t,hmin);
                            
                            tout = tout(1:nout);
                            yout = yout(:, 1:nout);
                            return;
                        end
                        
                        nofailed = false;
                        absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                        h = absh;
                        done = false;
                    else % Successful step
                        break;
                    end
                end % while true
                notFirstStep = true;
                
                oldnout = nout;
                nout = nout + 1;
                if nout > length(tout)
                    tout = [tout, zeros(1,chunk)];  %#ok<AGROW> % requires chunk >= refine
                    yout = [yout, zeros(neq,chunk)]; %#ok<AGROW>
                end
                idx = oldnout+1:nout;
                tout(idx) = tnew;
                yout(:,idx) = ynew;
                
                if done
                    break;
                end
                
                % If there were no failures compute a new h.
                if nofailed
                    % Note that absh may shrink by 0.8, and that err may be 0.
                    temp = 1.25*(err/rtol)^pow;
                    if temp > 0.2
                        absh = absh / temp;
                    else
                        absh = 5.0*absh;
                    end
                end
                
                % Advance the integration one step.
                t = tnew;
                y = ynew;
                f0 = f2; % because formula is FSAL
                
            end % while ~done
            
            tout = tout(1:nout);
            yout = yout(:, 1:nout);
        end
        
        function enzymes = discretizeEnzymes(this, enzymeConcentrations)
            enzymes = this.randStream.stochasticRound(this.concentrationToMolecules(enzymeConcentrations));
            nMonomers = [1 1 1:this.maxPolymerLength];
            while true
                dMonomer = nMonomers * (enzymes - this.enzymes);
                if dMonomer == 0
                    break;
                elseif dMonomer > 0
                    idx = ceil(numel(nMonomers) * this.randStream.rand());
                    if enzymes(idx) > 0
                        enzymes(idx) = enzymes(idx) - 1;
                        if idx > 2
                            enzymes(idx - 1) = enzymes(idx - 1) + 1;
                        end
                    end
                elseif dMonomer < 0
                    idx = ceil(numel(nMonomers) * this.randStream.rand());
                    if idx > 2
                        if enzymes(idx - 1) > 0
                            enzymes(idx - 1) = enzymes(idx - 1) - 1;
                            enzymes(idx) = enzymes(idx) + 1;
                        end
                    else
                        enzymes(idx) = enzymes(idx) + 1;
                    end
                end
            end
        end
        
        function [enzymes, substrates] = applySubstrateLimits(this, enzymes, substrates)
            nGTP = [0 0 1:this.maxPolymerLength];
            nGDP = [0 1 zeros(1, this.maxPolymerLength)];
            
            while nGTP * (enzymes - this.enzymes) > substrates(this.substrateIndexs_gtp)
                idx = find(enzymes(this.enzymeIndexs_FtsZ_activated), 1, 'first');
                enzymes(this.enzymeIndexs_FtsZ_activated(idx)) = ...
                    enzymes(this.enzymeIndexs_FtsZ_activated(idx)) - 1;
                enzymes(this.enzymeIndexs_FtsZ) = ...
                    enzymes(this.enzymeIndexs_FtsZ) + idx;
            end
            
            while nGDP * (enzymes - this.enzymes) > ...
                    substrates(this.substrateIndexs_gdp) + ...
                    min(substrates(this.substrateIndexs_water), ...
                    substrates(this.substrateIndexs_gtp) - nGTP * (enzymes - this.enzymes))
                enzymes(this.enzymeIndexs_FtsZ_GDP) = ...
                    enzymes(this.enzymeIndexs_FtsZ_GDP) - 1;
                enzymes(this.enzymeIndexs_FtsZ) = ...
                    enzymes(this.enzymeIndexs_FtsZ) + 1;
            end
            
            substrates(this.substrateIndexs_gtp) = ...
                + substrates(this.substrateIndexs_gtp) ...
                - nGTP * (enzymes - this.enzymes);
            substrates(this.substrateIndexs_gdp) = ...
                + substrates(this.substrateIndexs_gdp) ...
                - nGDP * (enzymes - this.enzymes);
            
            substrates = substrates ...
                + [1; -1; 1; -1; 1] * max(0, -substrates(this.substrateIndexs_gdp));
        end
        
        function result = moleculesToConcentration(this, count)
            N = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            result = count / (N * this.geometry.volume);
        end
        
        function result = concentrationToMolecules(this, conc)
            N = edu.stanford.covert.util.ConstantUtil.nAvogadro;
            result = conc * (N * this.geometry.volume);
        end
    end
end

%Computes dy/dt given y, the vector of FtsZ molecule concentrations.
function dydt = diff(y, p)
% activationFwd = p(1);
% activationRev = p(2);
% exchangeFwd = p(3);
% exchangeRev = p(4);
% nucleationFwd = p(5);
% nucleationRev = p(6);
% elongationFwd = p(7);
% elongationRev = p(8);
% gdpConcentration = p(9);
% gtpConcentration = p(10);

dydt = zeros(size(y));

%inactive monomer
dydt(1, :) = ...
    - p(1) * y(1, :) ...
    + p(2) * y(3, :);
%deactivated monomer
dydt(2, :) = ...
    - p(3) * y(2, :) * p(10) ...
    + p(4) * y(3, :) * p(9);
%activated monomer
dydt(3, :) = ...
    + p(1) * y(1, :) ...
    - p(2) * y(3, :) ...
    + p(3) * y(2, :) * p(10) ...
    - p(4) * y(3, :) * p(9) ...
    - 2 * p(5) * y(3, :) .^2 ...
    + 2 * p(6) * y(4, :) ...
    - p(7) * y(3, :) .* sum(y(4:end-1, :), 1) ...
    + p(8) * sum(y(5:end, :), 1);
%activated dimer
dydt(4, :) = ...
    + p(5) * y(3, :) .^2 ...
    - p(6) * y(4, :) ...
    - p(7) * y(3, :) .* y(4, :) ...
    + p(8) * y(5, :);
%activated polymer
for i = 5:size(y, 1)-1
    dydt(i, :) = ...
        + p(7) * y(3, :) .* y(i-1, :) ...
        - p(8) * y(i, :) ...
        - p(7) * y(3, :) .* y(i, :) ...
        + p(8) * y(i+1, :);
end
dydt(end, :) = ...
    + p(7) * y(3, :) .* y(end-1, :)...
    - p(8) * y(end, :);
end

function J = jacobian(y, p)
% activationFwd = p(1);
% activationRev = p(2);
% exchangeFwd = p(3);
% exchangeRev = p(4);
% nucleationFwd = p(5);
% nucleationRev = p(6);
% elongationFwd = p(7);
% elongationRev = p(8);
% gdpConcentration = p(9);
% gtpConcentration = p(10);

J = zeros(numel(y));

%inactive monomer
J(1, 1) = - p(1);
J(1, 3) = + p(2);
%deactivated monomer
J(2, 2) = - p(3) * p(10);
J(2, 3) = + p(4) * p(9);
%activated monomer
J(3, 1) = + p(1);
J(3, 2) = + p(3) * p(10);
J(3, 3) = ...
    - p(2) ...
    - p(4) * p(9) ...
    - 4 * p(5) * y(3) ...
    - p(7) * sum(y(4:end-1));
J(3, 4) = ...
    + 2 * p(6);
J(3, 5:end) = ...
    + p(8);
J(3, 4:end-1) = J(3, 4:end-1) ...
    - p(7) * y(3);
%activated dimer
J(4, 3) = ...
    + 2 * p(5) * y(3) ...
    - p(7) * y(4);
J(4, 4) = ...
    - p(6) ...
    - p(7) * y(3);
J(4, 5) = ...
    + p(8);
%activated polymer
for i = 5:size(y, 1)-1
    J(i, 3) = ...
        + p(7) * y(i-1) ...
        - p(7) * y(i);
    J(i, i-1) = ...
        + p(7) * y(3);
    J(i, i) = ...
        - p(8) ...
        - p(7) * y(3);
    J(i, i+1) = ...
        + p(8);
end
J(end, 3)     = + p(7) * y(end-1);
J(end, end-1) = + p(7) * y(3);
J(end, end)   = - p(8);
end