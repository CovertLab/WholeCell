%Translocate Proteins
%
% @wholeCellModelID Process_ProteinTranslocation
% @name             Protein translocation
% @description
%   Biology
%   ===============
%   Proteins are produced in the cytoplasm, and integral membrane proteins,
%   lipoproteins, and extracellular proteins must be translocated into and
%   through the cell membrane to reach their intended destination. These three
%   types of proteins are recognized by the translocation machinery through
%   N-terminal signal sequences. Integral membrane proteins are first
%   recognized by the signal recognition particle (SRP) which in turn delivers
%   them to the the preprotein translocase and pore. Lipoproteins and
%   extracellular proteins are recognized directly by the preprotein translocase
%   and pore. After association with the preprotein translocase, proteins are
%   pushed through the preprotein translocation pore by an ATP-dependent,
%   step-wise mechanism. After insertion into/through the membrane lipoprotein
%   and extracellular protein signal peptides are cleaved (see Protein
%   Processing II process).
%
%   This process simulates protein translocation into and through the cell
%   membrane involving four enzymes:
%   - signal recognition particle
%   - signal recognition particle receptor
%   - translocase ATPase
%   - translocase pore
%
%   Knowledge Base
%   ===============
%   We assigned each protein monomer and complex to one of the localizations:
%   - Integral membrane
%   - Lipoprotein
%   - Cytoplasmic
%   - Extracellular
%   - Terminal oganelle, cytoplasmic
%   - Terminal organelle, integral membrane
%
%   The localization of each protein, and the signal peptide length of
%   lipoproteins and secreted proteins was compiled from several sources:
%   - Computational prediction of membrane spanning domains and signal peptides
%     - Phobius [PUB_0262]
%     - PrediSi [PUB_0255]
%     - SignalP-HMM [PUB_0263]
%     - SignalP-NN [PUB_0263]
%     - SOSUI [PUB_0261, PUB_0264]
%     - SPdb database of observed signal peptides [PUB_0253]
%   - Mass-Spec determination of the N-terminal residue of each protein [PUB_0280]
%   - Databases of protein localization: BRENDA [PUB_0570], DBSubLoc [PUB_0573],
%     EchoBase [PUB_0574], GenoBase [PUB_0386], PSortDB [PUB_0572], and UniProt
%     [PUB_0096]
%   - Primary literature of the composition of the terminal organelle [PUB_0088,
%     PUB_0089, PUB_0091, PUB_0092, PUB_0093, PUB_0406, PUB_0407, PUB_0408, PUB_0409]
%   - Primary literature [PUB_0284, PUB_0303]
%
%   Additionally, we tried unsuccessfully to include computational predictions
%   from these sources:
%   - SecretomeP – didn't predicted any secreted peptides [PUB_0252]
%   - LipPred    – had CGI and bad request errors [PUB_0254]
%   - SIG-Pred   – found no signal sequences [PUB_0255]
%   - sigcleave  – unclear what it returns [PUB_0257]
%   - TatFind    – identified no Tat signal peptides [PUB_0258]
%   - PilFind    – identified no type IV pilin-like signal peptides [PUB_0259]
%   - SPEPLip    – provides no easy way to query on genome-scale [PUB_0260]
%
%   Roughly 25% of M. genitalium proteins require translocation.
%     Localization         No. Monomers
%     =================    ============
%     Cytosol              363
%     Integral membrane     82
%     Lipoprotein           14
%     Secreted              20
%     -----------------    ------------
%     Total                479
%
%   Representation
%   ===============
%   The substrates, enzymes, and monomers properties represent the counts of
%   metabolites, translocation enzymes, and protein. The substrates and enzymes
%   properties have compartment dimension length 1, meaning that the process only
%   accesses the counts of each metabolite and enzyme in the relevant
%   compartment. The monomers property has compartment dimension length 5,
%   meaning it accesses protein monomers in all compartments of the simulation.
%   Although it is known that protein translocation proceeds in discrete steps,
%   this process doesn't not represent any intermediate states of protein
%   translocation. The process treats protein translocation as an all-or-nothing
%   event.
%
%   monomerSRPPathways is a boolean which represents whether or not each protein
%   monomer requires the signal recognition particle (SRP) and its receptor to
%   translocate. This variable is true for integral membrane proteins, and false
%   otherwise. monomerLengths represents the number of amino acids in each
%   protein monomer species.
%   ceil(monomerLengths/preproteinTranslocase_aaTranslocatedPerATP) represents
%   the ATP cost of the translocase to translocate each monomer.
%   SRP_GTPUsedPerMonomer represents the GTP cost to translocate each integral
%   membrane protein.
%
%   Initialization
%   ===============
%   All protein monomers are initialized to the mature state, and in their
%   correct localization. This is achieved by the simulation class
%   initializeState method.
%
%   Simulation
%   ===============
%   In a randomized order, for each protein monomer across all species that
%   requires translocation
%   1. Calculate amount of ATP, GTP, SRP, and translocase required to
%      translocate the single monomer.
%   2. Terminate if insufficient resources exist to translocate the monomer.
%   3. Update the counts of ATP, GTP, ADP, GDP, Pi, H2O, and H+.
%   4. Decrement the counts of available SRP and translocase
%   5. Increment the count of the monomer in its localized compartment.
%      Decrement the count of the monoemr in the cytosol compartment.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/30/2010
classdef  ProteinTranslocation < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'translocaseSpecificRate';
            'preproteinTranslocase_aaTranslocatedPerATP';
            'SRP_GTPUsedPerMonomer';
            'monomerLengths';
            'monomerCompartments';
            'monomerSRPPathways';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'monomers'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ATP';'GTP';'ADP';'GDP';'H';'H2O';'PI'};
        substrateIndexs_atp       = 1; %index of ATP within substrateWholeCellModelIDs
        substrateIndexs_gtp       = 2; %index of ATP within substrateWholeCellModelIDs
        substrateIndexs_adp       = 3; %index of ADP within substrateWholeCellModelIDs
        substrateIndexs_gdp       = 4; %index of ADP within substrateWholeCellModelIDs
        substrateIndexs_hydrogen  = 5; %index of hydrogen within substrateWholeCellModelIDs
        substrateIndexs_water     = 6; %index of water within substrateWholeCellModelIDs
        substrateIndexs_phosphate = 7; %index of inorganic phosphate within substrateWholeCellModelIDs

        enzymeWholeCellModelIDs = {    %enzyme whole cell model ids
            'MG_0001_048';                   %signal recognition particle
            'MG_297_MONOMER';                %signal recognition particle receptor
            'MG_072_DIMER';                  %preprotein translocase, SecA subunit
            'MG_055_170_277_464_476_20MER'}; %preprotein translocase pore
        enzymeIndexs_signalRecognitionParticle         = 1; %index of signal recognition particle within enzymeWholeCellModelIDs
        enzymeIndexs_signalRecognitionParticleReceptor = 2; %index of signal recognition particle receptor within enzymeWholeCellModelIDs
        enzymeIndexs_translocaseATPase                 = 3; %index of translocase ATPase within enzymeWholeCellModelIDs
        enzymeIndexs_translocasePore                   = 4; %index of translocase pore within enzymeWholeCellModelIDs
                
        monomerWholeCellModelIDs         %whole cell model IDs of nascent protein monomers
        monomerIndexs_translocating      %indices of monomers that need to be translocated within monmers
        monomerIndexs_nontranslocating   %indices of monomers that don't need to be translocated within monmers
    end
    
    %fixed biological constants
    properties
        translocaseSpecificRate                     %number amino acids translocated per translocase per second (2.710e12) [PUB_0001]
        preproteinTranslocase_aaTranslocatedPerATP  %number amino acids translocated by translocase per ATP hydrolysis (35) [PUB_0001, PUB_0002, PUB_0003]
        SRP_GTPUsedPerMonomer                       %GTPs required for SRP to shuttle each monomer (2) [PUB_0003, PUB_0010]

        monomerLengths                              %lengths of nascent monomers
        monomerCompartments                         %indices of monomer compartments within compartment.wholeCellModelIDs
        monomerSRPPathways                          %true ==> monomer requires SRP to translocate; false ==> monomer doesn't require SRP
    end

    %global state (stored locally for convenience)
    properties
        monomers                    %numbers of nascent protein monomers
    end

    %constructor
    methods
        function this = ProteinTranslocation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
        
            this.monomerWholeCellModelIDs        = this.monomer.wholeCellModelIDs(this.monomer.processedIIndexs);
            this.monomerLengths                  = this.monomer.lengths(this.monomer.processedIIndexs);
            this.monomerCompartments             = double(knowledgeBase.proteinMonomerCompartments);
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;
            this.monomerIndexs_translocating     = find(this.monomerCompartments ~= this.compartment.cytosolIndexs);
            this.monomerIndexs_nontranslocating  = find(this.monomerCompartments == this.compartment.cytosolIndexs);
            this.monomerSRPPathways              = ismember({knowledgeBase.proteinMonomers.signalSequenceType}',{'lipoprotein','secretory'}');
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            this.monomers = this.monomer.counts(this.monomer.processedIIndexs, :, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.monomer.counts(this.monomer.processedIIndexs, :, :) = this.monomers;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.monomers = zeros(length(this.monomerWholeCellModelIDs), this.compartment.count, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            %preprotein translocase ATP cost
            translocaseEnergy = states.monomerProductions(this.monomerIndexs_translocating)' * ...
                ceil(this.monomerLengths(this.monomerIndexs_translocating) / this.preproteinTranslocase_aaTranslocatedPerATP);
            
            bmProd(this.substrateIndexs_atp)       = bmProd(this.substrateIndexs_atp)       + translocaseEnergy;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + translocaseEnergy;
            byProd(this.substrateIndexs_adp)       = byProd(this.substrateIndexs_adp)       + translocaseEnergy;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + translocaseEnergy;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + translocaseEnergy;
            
            %GTP cost of SRP for monomers with N-terminal signal sequences
            srpEnergy = this.SRP_GTPUsedPerMonomer * this.monomerSRPPathways' * states.monomerProductions;
            
            bmProd(this.substrateIndexs_gtp)       = bmProd(this.substrateIndexs_gtp)       + srpEnergy;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + srpEnergy;
            byProd(this.substrateIndexs_gdp)       = byProd(this.substrateIndexs_gdp)       + srpEnergy;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + srpEnergy;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + srpEnergy;
            
            %% enzymes
            %preprotein translocase and pore
            minEnzExp([this.enzymeIndexs_translocaseATPase; this.enzymeIndexs_translocasePore]) = ...
                2 * states.monomerProductions0(this.monomerIndexs_translocating)' * ...
                ceil(this.monomerLengths(this.monomerIndexs_translocating) / this.preproteinTranslocase_aaTranslocatedPerATP) / ...
                this.translocaseSpecificRate;
            minEnzExp(this.enzymeIndexs_translocasePore) = minEnzExp(this.enzymeIndexs_translocaseATPase);

            %SRP for monomers with N-terminal signal sequences
            minEnzExp([this.enzymeIndexs_signalRecognitionParticle; this.enzymeIndexs_signalRecognitionParticleReceptor]) = ...
                2 * this.monomerSRPPathways' * states.monomerProductions0;
        end

        %initialization: monomers initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            
            translocases = ceil(min(this.enzymes([this.enzymeIndexs_translocaseATPase this.enzymeIndexs_translocasePore])) * this.translocaseSpecificRate * this.stepSizeSec);
            SRPs         = ceil(min(this.enzymes([this.enzymeIndexs_signalRecognitionParticle this.enzymeIndexs_signalRecognitionParticleReceptor])) * this.stepSizeSec);
            
            %ATP cost of preprotein translocase
            result(this.substrateIndexs_atp) = min(...
                translocases, ...
                this.monomers(this.monomerIndexs_translocating, this.compartment.cytosolIndexs)' * ...
                ceil(this.monomerLengths(this.monomerIndexs_translocating) / this.preproteinTranslocase_aaTranslocatedPerATP));
            
            %GTP cost of SRP for monomers with N-terminal signal sequences
            result(this.substrateIndexs_gtp) = this.SRP_GTPUsedPerMonomer * min(...
                SRPs, ...
                this.monomerSRPPathways' * this.monomers(:, this.compartment.cytosolIndexs));
            
            %water
            result(this.substrateIndexs_water) = ...
                + result(this.substrateIndexs_atp) ...
                + result(this.substrateIndexs_gtp);
        end

        %simulation
        function evolveState(this)
            %randomly select protein monomers to translocate
            monomersNeedingTranslocation = cumsum(this.monomers(this.monomerIndexs_translocating, this.compartment.cytosolIndexs)); %cumulative sum of monomers that need to be translocated

            %stop if no monomers need to be translocated
            if isempty(monomersNeedingTranslocation) || monomersNeedingTranslocation(end) == 0
                return;
            end
            
            %resources required for translocation
            %-substrates: ATP, water
            %-translocase pore, ATPase
            %-signal recognition particle, receptor
            atp          = this.substrates(this.substrateIndexs_atp);
            gtp          = this.substrates(this.substrateIndexs_gtp);
            water        = this.substrates(this.substrateIndexs_water);
            translocases = min(this.enzymes([this.enzymeIndexs_translocaseATPase this.enzymeIndexs_translocasePore])) * this.translocaseSpecificRate * this.stepSizeSec / this.preproteinTranslocase_aaTranslocatedPerATP;
            SRPs         = min(this.enzymes([this.enzymeIndexs_signalRecognitionParticle this.enzymeIndexs_signalRecognitionParticleReceptor])) * this.stepSizeSec;

            %randomly select protein monomers to translocate
            randomMonomerOrder = this.randStream.randperm(monomersNeedingTranslocation(end)); %random ordering over monomers that need to be translocated

            for i = 1:monomersNeedingTranslocation(end)
                %indices
                idx_translocatingMonomer = find(randomMonomerOrder(i) <= monomersNeedingTranslocation, 1, 'first'); %index of randomly selected monomer within monomerIndexs_translocating
                idx_monomer = this.monomerIndexs_translocating(idx_translocatingMonomer); %index of randomly selected monomer within monomers

                %cost of translocating this monomer
                monomerTranslocaseCost = ceil(this.monomerLengths(idx_monomer)/this.preproteinTranslocase_aaTranslocatedPerATP); %translocase cost of translocating the randomly selected monomer (1/translocaseSpecificRate translocase enzymes per preproteinTranslocase_aaTranslocatedPerATP amino acids)
                monomerSRPCost         = this.monomerSRPPathways(idx_monomer);              %SRP cost of translocating the randomly selected monomer (1 signal recognition particle and receptor per monomer)
                monomerATPCost         = monomerTranslocaseCost;                            %ATP cost of translocating the randomly selected monomer (translocase: 1 ATP per preproteinTranslocase_aaTranslocatedPerATP amino acids)
                monomerGTPCost         = monomerSRPCost * this.SRP_GTPUsedPerMonomer;       %GTP cost of translocating the randomly selected monomer (SRP: 2 GTP per monomer)
                monomerWaterCost       = monomerATPCost + monomerGTPCost;                   %water cost of translocating the randomly selected monomer (SRP: 2 water per monomer; translocase: 1 water per preproteinTranslocase_aaTranslocatedPerATP amino acids)

                %stop if insufficient resources available for translocation
                if ...
                        translocases < monomerTranslocaseCost || ...
                        SRPs         < monomerSRPCost         || ...
                        atp          < monomerATPCost         || ...
                        gtp          < monomerGTPCost         || ...
                        water        < monomerWaterCost
                    break;
                end

                %update protein monomers
                this.monomers(idx_monomer, this.monomerCompartments(idx_monomer)) = ...
                    this.monomers(idx_monomer, this.monomerCompartments(idx_monomer)) + 1;
                this.monomers(idx_monomer, this.compartment.cytosolIndexs) = ...
                    this.monomers(idx_monomer, this.compartment.cytosolIndexs) - 1;

                %update substrates
                atp          = atp          - monomerATPCost;
                gtp          = gtp          - monomerGTPCost;
                water        = water        - monomerWaterCost;
                translocases = translocases - monomerTranslocaseCost;
                SRPs         = SRPs         - monomerSRPCost;

                this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - monomerATPCost;
                this.substrates(this.substrateIndexs_gtp)       = this.substrates(this.substrateIndexs_gtp)       - monomerGTPCost;
                this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - monomerWaterCost;
                this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + monomerATPCost;
                this.substrates(this.substrateIndexs_gdp)       = this.substrates(this.substrateIndexs_gdp)       + monomerGTPCost;
                this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + monomerWaterCost;
                this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + monomerWaterCost;
            end
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.monomers, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    this.monomer.molecularWeights(this.monomer.processedIIndexs)' * sum(this.monomers, 2) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + ...
                    permute(this.monomer.molecularWeights(this.monomer.processedIIndexs)' * permute(sum(this.monomers, 2),[1 3 2]),[1 3 2]) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
