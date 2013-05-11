%DNASupercoiling
%
% @wholeCellModelID Process_DNASupercoiling
% @name             DNA Supercoiling
% @description
%   Biology
%   ===========
%   DNA gyrase and topoisomerase IV each use 2 ATP to induce 2
%   negative
%   supercoils each time they act. Occasionally, topoisomerase IV can also
%   induce positive supercoils but that is not considered in our model.
%   Topoisomerase I acts to induce positive supercoils. Gyrase, topoisomerase
%   IV, and topoisomerase I act at rates of 1.2, 2.5, and 1 strand passing
%   events per second respectively.
% 
%   We use a calculation of the DNA's linking number (LK) in order to track
%   the supercoiling of the DNA.  The &delta;LK is the difference between the
%   current level of DNA supercoiling and the relaxed level of DNA
%   supercoiling. The LKrelaxed is defined as number of base pairs/10.5, in which
%   10.5 is the number of bases per turn in a relaxed double helix. As the
%   replication loops move, the LKcurrent deviates from this relaxed state,
%   and gyrases and topoisomerases help bring the DNA back to the relaxed
%   state. The superhelical density, or specific linking number density,
%   &sigma;<sub>sp</sub>, is defined as the:
%      (LKcurrent - LKrelaxed)/LKrelaxed. 
%   The activity of gyrases and topoisomerases depends on the &sigma;<sub>sp</sub>
%   of the DNA. Although there is likely a more complex relationship between the
%   DNA supercoiling and enzyme activity, we model the enzyme activity as a
%   combination of step functions and logistic functions. 
%   Topoisomerase IV can only act if &sigma;<sub>sp</sub> is
%   higher than 0. Topoisomerase I can only act if &sigma;<sub>sp</sub> is lower
%   than zero. Gyrase can only act if &sigma;<sub>sp</sub> is higher than -0.1.
%   Within, the regions where gyrase and topoI are allowed to act, we apply
%   a logistic function describing the probability of activity, centered
%   around equilibrium sigma, -0.06. For TopoI, the probability of activity
%   approaches 1 as the sigma gets more positive. For gyrase, the
%   probability of activity approaches 1 as
%   
%   We model up to three regions on the chromosomes and track their LKs
%   separately. Before replication, "unreplicated DNA" is the only region
%   present. The enzymes may act on the chromosome affecting the LKcurrent,
%   and experimentally it has been shown that the enzymes would obtain an
%   LKcurrent such that the steady state &sigma;<sub>sp</sub> is =-0.06. As the
%   replication loop progresses, "unreplicated DNA" is the region downstream of
%   the two replication loops.  The number of bases in this region decreases
%   during replication, meaning that the LKrelaxed decreases. The LKcurrent then,
%   is too high, and must be brought down towards the LKrelaxed by inducing
%   negative supercoils.
%   
%   The second and third regions are the replicated DNA upstream of the
%   replication loops on each of the two chromosomes. As new uncoiled DNA is
%   formed, it is already coiled, but gyrases and topoisomerases can continue
%   to act on this DNA, ideally maintaining the steady state &sigma;<sub>sp</sub>.
%   After replication is complete, these two regions are the only two that exist.
%   
%   It is essential that the &detlta;LK in the region downstream of the replication
%   loops be brought down to 0 by the end of replication.
%   
%   Another consideration is the processivity of the enzymes when bound to the
%   DNA. Topoisomerase IV is highly processive and will stay bound to the DNA
%   as long as the &sigma;<sub>sp</sub> is greater than zero. Topoisomerase I is
%   not highly processive, and essentially acts on the DNA and falls back off
%   right away gyrase will stay bound to the DNA for about 30-60 seconds. We model
%   gyrase processivity as a poisson distribution with &lambda; = 45 seconds.
%   
%   In this process, we go through all free gyrases and topoisomerases in
%   random order,  determine what regions they can bind in, and randomly bind
%   them to a large enough open position on the DNA. We adjust the linking
%   number based on all enzyme actions, and account for the usage of ATP. We
%   also track the processivity of gyrase and topoisomerase IV, and unbind
%   them from the DNA when appropriate.
%   
%   If replication is in progress, we first knock off any enzymes that the
%   replication loop collides into.
%
%   There is an effect of supercoiling on the probabilities of gene 
%   transcription (Peter 2004). While fold change at differnt sigmas have
%   been calculated for many E. coli genes, here, for simplicity, we are
%   only including the effects on the 5 supercoiling genes: gyrB, gyrA,
%   parC, parE, and topA. These genes exist in 3 transcription units. Peter
%   2004 has data for the fold change of expression of each of these genes
%   are various values of sigma (ranging from sigma -0.06-0.02) at various 
%   experimental conditions. Due to the limited data and large variation within
%   the data, we have decided to use a linear fit of the fold change data and 
%   extrapolate the linear fit within the sigma range of -0.08 to 0.07 where 
%   the linear fit seems reasonable. Outside of this range, we estimate a
%   constant fold change of expression. While there is a separate set of
%   data for each of the 5 genes, our model requires a single probability
%   of transcription for each transcription unit. The data for the first
%   gene in each transcription unit is used (gyrB, and parE). TopA is
%   transcribed with genes that are not supercoiling related. The
%   expression of those genes (MG_119, 120, 121) will also get affected by 
%   this process. The general trend is that gyrase and topoIV will have an
%   increased expression when sigma is higher than the equilibrium, and
%   that topoI will have a higher expression when sigma is lower than the
%   equilibrium. Fold changes for gyrase will vary between 0.14 and 6.57.
%   Fold changes for topoIV will vary between 0.98 and 1.14. Fold changes
%   for topoI will vary between 0.042 and 1.15. 
%   
%   References
%   ===============
%   1. Ullsperger, C., Cozzarelli, N.R. (1996). Contrasting enzymatic activities
%      of topoisomerase IV and DnA gyrase from Escherichia coli. Journal of Bio
%      Chem 271: 31549-31555. [PUB_0236]
%   2. Dekker, N.H., Viard, T., Bouthier de la Tour, C., Duguet, M., Bensimon,
%      D., Croquette, V. (2003). Thermophilic Topoisomerase I on a single DNA
%      molecule. Journal of molecular biology 329: 271-282. [PUB_0502]
%   3. Gore, J., Bryant, Z., Stone, M.D., Nollmann, M., Cozzarelli, N.R.,
%      Bustamante, C. (2006). Mechanochemical analysis of DNA gyrase using rotor
%      bead tracking. Nature 439: 100-104. [PUB_0751]
%   4. Bates, A. (2006). DNA Topoisomerases: Single Gyrase Caught in the Act.
%      Current Biology 16: 204-206. [PUB_0752]
%   5. Peng, H., Marians, K.J. (1995). The Interaction of Escherichia coli
%      Topoisomerase IV with DNA. Journal of Biological Chemistry 42:
%      25286–25290. [PUB_0694]
%   6. Wang, J. (1996) DNA Topoisomerases. Annual Reviews 65: 635-692. [PUB_0693]
%   7. Peter, B.J., Arsuaga, J., Breier, A.M., Khodursky, A.B., Brown,
%      P.O., Cozzarelli, N.R. (2004) Genomic transcriptional response to loss
%      of chromosomal supercoiling in Escherichia coli. Genome Biology 5:
%      1-13. [PUB_0920]
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/18/2010
classdef DNASupercoiling < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'gyraseMeanDwellTime';
            'gyraseSigmaLimit';
            'topoISigmaLimit';
            'topoIVSigmaLimit';
            'gyraseDeltaLK';
            'topoIDeltaLK';
            'topoIVDeltaLK';
            'gyraseActivityRate';
            'topoIActivityRate';
            'topoIVActivityRate';
            'gyraseATPCost';
            'topoIATPCost';
            'topoIVATPCost';
            'foldChangeSlopes'; 
            'foldChangeIntercepts';
            'foldChangeLowerSigmaLimit';
            'foldChangeUpperSigmaLimit'; 
            'numTranscriptionUnits';
            'tuIndexs';
            'tuCoordinates';
            };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};   %whole cell model IDs of stimuli
        substrateWholeCellModelIDs = {   %whole cell model IDs of substrates
            'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'};
        substrateIndexs_atp       = 1;
        substrateIndexs_adp       = 2;
        substrateIndexs_phosphate = 3;
        substrateIndexs_water     = 4;
        substrateIndexs_hydrogen  = 5;

        enzymeWholeCellModelIDs = {      %whole cell model IDs of enzymes
            'DNA_GYRASE';                %DNA gyrase
            'MG_203_204_TETRAMER';       %DNA topoisomerase IV
            'MG_122_MONOMER'};           %DNA topoisomerase I
        enzymeIndexs_gyrase = 1;
        enzymeIndexs_topoIV = 2;
        enzymeIndexs_topoI  = 3;
    end
    
    %fixed biological constants
    properties
        gyraseMeanDwellTime       %Mean dwell time of gyrase to bind DNA [45s; PUB_0753]
        
        gyraseSigmaLimit          %superhelical density (\sigma_{sp}) above which gyrase acts [-0.1; PUB_0753]
        topoISigmaLimit           %superhelical density (\sigma_{sp}) above which topoisomerase I acts [0; PUB_0502]
        topoIVSigmaLimit          %superhelical density (\sigma_{sp}) below which topoisomerase IV acts [0; PUB_0753]
        
        gyraseDeltaLK             %number of supercoils introduced by each gyrase strand passing event [-2; PUB_0695]
        topoIDeltaLK              %number of supercoils introduced by each topoisomerase I strand passing event [1; PUB_0502]
        topoIVDeltaLK             %number of supercoils introduced by each topoisomerase IV strand passing event [-2]
        
        gyraseActivityRate        %strand passing events per second per enzyme [1.2; PUB_0753]
        topoIActivityRate         %strand passing events per second per enzyme [1.0; PUB_0502]
        topoIVActivityRate        %strand passing events per second per enzyme [2.5; PUB_0750]
        
        gyraseATPCost             %ATP hydrolyzed per gyrase-catalyzed strand passing event [2]
        topoIATPCost              %ATP hydrolyzed per topoisomerase I-catalyzed strand passing event [0]
        topoIVATPCost             %ATP hydrolyzed per topoisomerase IV-catalyzed strand passing event [2]
        
        topoILogisiticConst       %fittable parameter controlling the steapness of the the logistic function of topoI activity around equilibrium sigma [100]
        gyrLogisiticConst         %fittable parameter controlling the steapness of the the logistic function of gyrase activity around equilibrium sigma [-100]
        
        enzymeProperties          %struct built by buildEnzymeProperties
        
        foldChangeSlopes          %slopes of linear fit of sigma-gene expression data [gyrase, topo IV, topo I] [PUB_0920] 
        foldChangeIntercepts      %y-intercepts of linear fit of sigma-gene expression data [gyrase, topo IV, topo I] [PUB_0920]
        foldChangeLowerSigmaLimit %lower bound of sigma for applying linear fit (fittable parameter)
        foldChangeUpperSigmaLimit %upper bound of sigma for applying linear fit (fittable parameter)
        
        numTranscriptionUnits     %number of transcription units
        tuIndexs                  %transcription units whose transcription probs will be affected by supercoiling [gyrase, topo IV, topo I]
        tuCoordinates             %start coords of transcription units whose transcription probs will be affected by supercoiling [gyrase, topo IV, topo I]
    end
    
    %references to cell state
    properties
        rnaPolymerase             %reference to RNA polymerase state
    end

    %constructor
    methods
        function this = DNASupercoiling(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.states = [this.states; {this.rnaPolymerase}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, simulation, varargin{:});            
            
            c = this.chromosome;
            
            this.enzymeProperties = this.buildEnzymeProperties();
            
            %get transcription unit indices
            [~, geneIndexs] = ismember({'MG_003', 'MG_203', 'MG_122'}, this.gene.wholeCellModelIDs);
            [this.tuIndexs, ~] = find(this.rna.nascentRNAGeneComposition(geneIndexs, :)');            
            this.numTranscriptionUnits = numel(c.transcriptionUnitStartCoordinates);
            this.tuCoordinates = c.transcriptionUnitStartCoordinates(this.tuIndexs);
        end
    end

    %model
    methods        
        %Calculate sufficient metabolites and enzymes to achieve the observed
        %mean superhelical density after replication of the daughter chromosome:
        %
        %  (LK after replication) = 2 * (LK before replication)
        %                         = (LK before replication)
        %                            + (links introduced by DNA pol)
        %                            - (links pushed ahead of replisome and resolved by gyrase)
        %                            - (links introduced behind replisome by gryase to return to steady-state superhelical density)
        %               2*LK_{ss} = LK_{ss} + 2*LK_{0} - LK_{ss} - 2*(LK_{0} - LK_{ss})
        %
        %That is calculate:
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        %
        %FitConstants accounts for topoisomerase I / gyrase expression balance
        %required to achieve the observed steady state superhelical density.
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, constants, states)
            %initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %references
            c = this.chromosome;
            t = constants.states.Time;
            
            %number of links to be resolved by gyrase, and number of gyrase
            %link resolution reactions
            %- links that need to be resolved because of replication
            %- links that need to be resolved to counter-balance the effect
            %  of topoisomerase I
            nRxn = 0;
            
            LK_0  =  c.sequenceLen / c.relaxedBasesPerTurn;
            LK_ss = (1 + c.equilibriumSuperhelicalDensity) * LK_0;
            dLK = - LK_ss - 2*(LK_0 - LK_ss);
            nRxn = nRxn + dLK / this.gyraseDeltaLK;
            
            topoIProps = this.enzymeProperties(this.enzymeIndexs_topoI);
            gyrProps = this.enzymeProperties(this.enzymeIndexs_gyrase);
            nRxn = nRxn + 0.25 * ... %factor of 0.25 added based on observing simulations
                abs(states.monomers(this.enzymeGlobalIndexs(this.enzymeIndexs_topoI)) * ...
                topoIProps.activityRate * topoIProps.deltaLK / gyrProps.deltaLK);
            
            %substrates and byproducts
            nATP = nRxn * this.enzymeProperties(this.enzymeIndexs_gyrase).atpCost;
            bmProd(this.substrateIndexs_atp)       = nATP;
            bmProd(this.substrateIndexs_water)     = nATP;
            byProd(this.substrateIndexs_adp)       = nATP;
            byProd(this.substrateIndexs_phosphate) = nATP;
            byProd(this.substrateIndexs_hydrogen)  = nATP;
            
            %Gyrase
            minEnzExp(this.enzymeIndexs_gyrase) = ...
                nRxn / t.replicationDuration / this.gyraseActivityRate / ...
                exp(log(2) * t.replicationInitiationDuration / t.cellCycleLength);
        end

        %initialization:
        %- Simulation initializeState method initializes 1
        %  chromosome with linkingNumber such that the superhelical density is
        %  equal to the equilibriumSuperhelicalDensity
        %- Here we bind all gyrase to the  
        function initializeState(this)
            %bind gyrase
            nGyraseBinding = this.randStream.stochasticRound(...
                this.enzymes(this.enzymeIndexs_gyrase) * ...
                (1 - 1/this.gyraseMeanDwellTime/this.stepSizeSec));
            this.bindProteinToChromosomeStochastically(this.enzymeIndexs_gyrase, ...
                nGyraseBinding);
            
            %Calculate effects of supercoiling on probability of transcription
            this.rnaPolymerase.supercoilingBindingProbFoldChange = this.calcRNAPolymeraseBindingProbFoldChange();
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            %chromosome reference
            c = this.chromosome;
            
            %initialize
            result = zeros(size(this.substrates));
            
            %ATP needed for gyrase activity
            [~, sigmas] = find(c.superhelicalDensity);
            if isempty(sigmas)
                sigmas = 0;
            end
            
            enzProps = this.enzymeProperties(this.enzymeIndexs_gyrase);
            nATP = enzProps.atpCost * ceil(...
                (this.enzymes(this.enzymeIndexs_gyrase) + this.boundEnzymes(this.enzymeIndexs_gyrase)) * ...
                enzProps.activityRate *  max(enzProps.probOfActivityFunc(sigmas, c)));
            
            result(this.substrateIndexs_atp)   = nATP;
            result(this.substrateIndexs_water) = nATP;
        end

        %simulation
        function evolveState(this)
            import edu.stanford.covert.util.CircularSparseMat;
            import edu.stanford.covert.util.ComputationUtil;
            
            c = this.chromosome;
            
            %get enzyme properties
            enzProps = this.enzymeProperties;

            % Calculate sigmas
            [tmpPosStrands, lengths] = find(c.doubleStrandedRegions);
            tmpIdxs = find(mod(tmpPosStrands(:, 2), 2));
            dsPosStrands = tmpPosStrands(tmpIdxs, :);
            lengths = lengths(tmpIdxs, 1);
            
            [~, linkingNumbers] = find(c.linkingNumbers);
            linkingNumbers = linkingNumbers(tmpIdxs, 1);
            sigmas = (linkingNumbers - lengths/c.relaxedBasesPerTurn) ./ (lengths/c.relaxedBasesPerTurn);

            % Determine which regions accommodate which proteins.
            legal = [...
                sigmas > this.gyraseSigmaLimit ...
                sigmas > this.topoIVSigmaLimit ...
                sigmas < this.topoISigmaLimit];

            % Unbind topoIV where it can't act.
            this.releaseProteinFromChromosome(...
                this.enzymeIndexs_topoIV, Inf,...
                dsPosStrands(legal(:, this.enzymeIndexs_topoIV), :),... %protected regions
                lengths(legal(:, this.enzymeIndexs_topoIV)));
            
            % Probabilistically unbind gyrase.
            this.releaseProteinFromChromosome(this.enzymeIndexs_gyrase, 1 / this.gyraseMeanDwellTime, [], []);
            
            % Bind enzymes in regions with accommodating sigmas.
            order = this.randStream.randperm(length(this.enzymes));
            transientlyBinding = zeros(size(legal));
            for i = 1:length(order)
                j = order(i);
                
                if enzProps(j).meanDwellTime == 0
                    if numel(lengths) == 1 && lengths == c.sequenceLen
                        [~, mons] = find(c.monomerBoundSites);
                        [~, cpxs] = find(c.complexBoundSites);
                        space = lengths ...
                            - sum(c.monomerDNAFootprints(mons, 1)) ...
                            - sum(c.complexDNAFootprints(cpxs, 1)) ...
                            - nnz(c.damagedSites);
                        transientlyBinding(1, j) = min(this.enzymes(j, 1), space / this.enzymeDNAFootprints(j, 1));
                    elseif numel(lengths) == 2 && all(lengths == c.sequenceLen)
                        [monPosStrnds, mons] = find(c.monomerBoundSites);
                        [cpxPosStrnds, cpxs] = find(c.complexBoundSites);
                        space = zeros(2, 1);
                        space(1) = lengths(1) ...
                            - sum(c.monomerDNAFootprints(mons(monPosStrnds(:, 2) <=2 ), 1)) ...
                            - sum(c.complexDNAFootprints(cpxs(cpxPosStrnds(:, 2) <=2 ), 1)) ...
                            - nnz(c.damagedSites(:, 1:2));
                        space(2) = lengths(2) ...
                            - sum(c.monomerDNAFootprints(mons(monPosStrnds(:, 2) >2 ), 1)) ...
                            - sum(c.complexDNAFootprints(cpxs(cpxPosStrnds(:, 2) >2 ), 1)) ...
                            - nnz(c.damagedSites(:, 3:4));
                        transientlyBinding(:, j) = min(this.enzymes(j, 1), sum(space) / this.enzymeDNAFootprints(j, 1)) * ...
                            space / sum(space);
                        if this.randStream.rand < 0.5
                            transientlyBinding(1, j) = ComputationUtil.roundHalfUp(transientlyBinding(1, j));
                            transientlyBinding(2, j) = ComputationUtil.roundHalfDown(transientlyBinding(2, j));
                        else
                            transientlyBinding(2, j) = ComputationUtil.roundHalfUp(transientlyBinding(2, j));
                            transientlyBinding(1, j) = ComputationUtil.roundHalfDown(transientlyBinding(1, j));
                        end
                    else
                        [monPosStrnds, mons] = find(c.monomerBoundSites);
                        [cmpPosStrnds, cmps] = find(c.complexBoundSites);
                        dmgPosStrnds = find(c.damagedSites);
                        monPos = monPosStrnds(:, 1);
                        cmpPos = cmpPosStrnds(:, 1);
                        dmgPos = dmgPosStrnds(:, 1);
                        monChrs = ceil(monPosStrnds(:, 2) / 2);
                        cmpChrs = ceil(cmpPosStrnds(:, 2) / 2);
                        dmgChrs = ceil(dmgPosStrnds(:, 2) / 2);
                        
                        space = zeros(size(lengths));
                        for k = 1:numel(lengths)
                            if ~legal(k, j)
                                continue;
                            end
                            
                            monTfs = monPos >= dsPosStrands(k, 1) & monPos <= dsPosStrands(k, 1) + lengths(k, 1) - 1;
                            cmpTfs = cmpPos >= dsPosStrands(k, 1) & cmpPos <= dsPosStrands(k, 1) + lengths(k, 1) - 1;
                            dmgTfs = dmgPos >= dsPosStrands(k, 1) & dmgPos <= dsPosStrands(k, 1) + lengths(k, 1) - 1;
                            monTfs(monTfs) = monChrs(monTfs, 1) == dsPosStrands(k, 2);
                            cmpTfs(cmpTfs) = cmpChrs(cmpTfs, 1) == dsPosStrands(k, 2);
                            dmgTfs(dmgTfs) = dmgChrs(dmgTfs, 1) == dsPosStrands(k, 2);
                            
                            space(k, 1) = ...
                                + lengths(k, 1) ...
                                - sum(c.monomerDNAFootprints(mons(monTfs, 1), 1)) ...
                                - sum(c.monomerDNAFootprints(cmps(cmpTfs, 1), 1)) ...
                                - sum(dmgTfs);
                        end
                        
                        if numel(lengths) == 1
                            transientlyBinding(:, j) = min(this.enzymes(j, 1), space / this.enzymeDNAFootprints(j, 1));
                        else
                            transientlyBinding(:, j) = space / sum(space) * min(this.enzymes(j, 1), sum(space) / this.enzymeDNAFootprints(j, 1));
                        end
                    end
                elseif any(legal(:, j))
                    this.bindProteinToChromosomeStochastically(...
                        j, [], dsPosStrands(legal(:,j),:), lengths(legal(:,j)));
                end
            end

            % Calculate new linking numbers and do substrate accounting.
            enzProps = enzProps(this.randStream.randperm(length(enzProps)));
            for i = 1:size(dsPosStrands, 1)
                len = lengths(i);
                for j = 1:length(enzProps)
                    if legal(i, enzProps(j).idx)
                        if enzProps(j).meanDwellTime == 0
                            nBound = transientlyBinding(i, enzProps(j).idx);
                        else
                            nBound = size(this.findProteinInRegion(...
                                dsPosStrands(i, 1), dsPosStrands(i, 2), len, enzProps(j).idx), 1);
                        end
                        % use logistic function to determine the
                        % probability of enzyme activity
                        % Use sigmas from beginning of timestep, before any
                        % enzymes have acted
                        probOfActivity = enzProps(j).probOfActivityFunc(sigmas(i),c);
                        nStrandPassingEvents = min([
                            this.randStream.stochasticRound(nBound * enzProps(j).activityRate *  probOfActivity);
                            fix(this.substrates([this.substrateIndexs_atp; this.substrateIndexs_water]) / enzProps(j).atpCost)
                            ]);
                        %adjust linking number
                        linkingNumbers(i) = linkingNumbers(i) + enzProps(j).deltaLK * nStrandPassingEvents;
                        
                        %substrate accounting
                        nATP = nStrandPassingEvents * enzProps(j).atpCost;
                        this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - nATP;
                        this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - nATP;
                        this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + nATP;
                        this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + nATP;
                        this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + nATP;
                    end
                end
            end
            c.linkingNumbers = CircularSparseMat(...
                [dsPosStrands; dsPosStrands(:, 1) dsPosStrands(:, 2)+1], [linkingNumbers; linkingNumbers], ...
                [c.sequenceLen c.nCompartments], 1);    
            %Calculate effects of supercoiling on probability of transcription
            this.rnaPolymerase.supercoilingBindingProbFoldChange = this.calcRNAPolymeraseBindingProbFoldChange(dsPosStrands, lengths, linkingNumbers);
        end
        
        function foldChange = calcRNAPolymeraseBindingProbFoldChange(this, dsPosStrands, lengths, linkingNumbers)
            c = this.chromosome;
            
            if nargin == 1
                [tmpPosStrands, lengths] = find(c.doubleStrandedRegions);
                tmpIdxs = find(mod(tmpPosStrands(:, 2),2));
                dsPosStrands = tmpPosStrands(tmpIdxs, :);
                lengths = lengths(tmpIdxs, 1);
                
                [~, linkingNumbers] = find(c.linkingNumbers);
                linkingNumbers = linkingNumbers(tmpIdxs, 1);
            end
            
            foldChange = ones(size(c.transcriptionUnitStartCoordinates, 1), 2); %initialize fold change
            sigmas = (linkingNumbers - lengths/c.relaxedBasesPerTurn) ./ (lengths/c.relaxedBasesPerTurn);
            for i = 1:length(this.tuCoordinates)
                %find double stranded region where transcription unit lies
                rgnIdxs = find(...
                    dsPosStrands(:,1)               <= this.tuCoordinates(i) & ...
                    dsPosStrands(:,1) + lengths - 1 >= this.tuCoordinates(i));
                
                for j = 1:numel(rgnIdxs)
                    %piecewise linear function of sigma of double stranded region
                    %where transcription unit lies
                    thresholdedSigma = min(...
                        max(sigmas(rgnIdxs(j)), this.foldChangeLowerSigmaLimit),...
                        this.foldChangeUpperSigmaLimit);
                    foldChange(...
                        this.tuIndexs(i), ceil(dsPosStrands(rgnIdxs(j), 2)/2)) = ...
                        this.foldChangeIntercepts(i) + this.foldChangeSlopes(i) * thresholdedSigma;
                end
            end
        end
        
        function value = buildEnzymeProperties(this)
            value = [...
                struct(...
                    'idx', this.enzymeIndexs_gyrase,...
                    'sigmaLimit', this.gyraseSigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas > this.gyraseSigmaLimit, ...
                    'atpCost', this.gyraseATPCost,...
                    'deltaLK', this.gyraseDeltaLK,...
                    'activityRate', this.gyraseActivityRate,...
                    'meanDwellTime', this.gyraseMeanDwellTime,...
                    'probOfActivityFunc', @(sigma, c) 1 ./ (1 + exp(this.gyrLogisiticConst * (sigma - c.equilibriumSuperhelicalDensity))));
                struct(...
                    'idx', this.enzymeIndexs_topoIV,...
                    'sigmaLimit', this.topoIVSigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas > this.topoIVSigmaLimit, ...
                    'atpCost', this.topoIVATPCost,...
                    'deltaLK', this.topoIVDeltaLK,...
                    'activityRate', this.topoIVActivityRate,...
                    'meanDwellTime', NaN,...
                    'probOfActivityFunc', @(sigma, c) ones(size(sigma)));
                struct(...
                    'idx', this.enzymeIndexs_topoI,...
                    'sigmaLimit', this.topoISigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas < this.topoISigmaLimit, ...
                    'atpCost', this.topoIATPCost,...
                    'deltaLK', this.topoIDeltaLK,...
                    'activityRate', this.topoIActivityRate,...
                    'meanDwellTime', 0,...
                    'probOfActivityFunc', @(sigma, c) 1 ./ (1 + exp(this.topoILogisiticConst * (sigma - c.equilibriumSuperhelicalDensity))))];
        end
    end
end
