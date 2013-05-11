%Assemble Terminal Organelle
%
% @wholeCellModelID Process_TerminalOrganelleAssembly
% @name             Terminal Organelle Assembly
% @description
%   Biology
%   ===============
%   The terminal organelle is a 300 x 80 nm membrane-bound bleb believed to
%   be involved in several processes including motility, adhesion,
%   replication, and cytokinesis. The terminal organelle is electron lucent 
%   and non-fibrillar and enriched for adhesins and adhesin accessory proteins 
%   which help the adhesins localize to the terminal organelle. The terminal 
%   organelle has been investigated in several studies by Mitchell Balish at
%   Miami University of Ohio and Duncan Krause of the University of Georgia. In
%   particular, they have examined the composition and assembly of the electron
%   dense (protein) portion of the terminal organelle, and have examined the
%   duplication and migration to the opposite pole of the daughter terminal
%   organelle during cell division. Pich et al have studied the regulation of
%   terminal organelle duplication.
%
%   The terminal organelle is composed of 8 proteins.
%   - HMW1-3
%   - MgPa
%   - P32, P65, P110, P200
%
%   These 8 proteins assemble in a stereotyped order; that is, there is a
%   heirarchical pattern to the assembly of the terminal organelle: several
%   proteins require other proteins to be already localized in the terminal
%   organelle to be incorporated into the terminal organelle, and two proteins
%   HMW1 and HMW2 require the other either to already been incorporated in the
%   terminal organelle, or free in the cytoplasm to be incorporated into the
%   terminal organelle. No kinetic information has been reported on terminal
%   organelle assembly.
%
%   During replication a second organelle is believed to duplicate from the
%   first organelle and migrate toward the opposite pole.
%
%   Knowledge Base
%   ===============
%   The hierarchical pattern of terminal organelle assembly among the 8 proteins
%   is encoded in 10 reactions in the knowledge base.
%
%   Representation
%   ===============
%   Proteins are incorporated according to one of two patterns:
%   - cytoplasm -> terminal organelle cytoplasm
%   - membrane  -> terminal organelle membrane
%
%   The counts of terminal organelle proteins are represented by the substrates
%   property. Because proteins are only incorporated in two patterns, the
%   substrates property has two compartments: one to represent the
%   unincorporated compartment of each protein (cytoplasm/membrane) and a
%   second for the incorporated compartment of each protein (terminal organelle
%   cytoplasm/membrane).
%
%   The hierarchical pattern of terminal organelle assembly is represented in
%   this process by localizationReactions and localizationSubstrates.
%   localizationSubstrates represents the protein translocated by each of the
%   10 reactions. localizationReactions represents the proteins (and the
%   compartment in which the protein must be located) required for each
%   translocation reaction to proceed.
%
%   Initialization
%   ===============
%   Terminal organelle proteins are all initialized into the terminal organelle
%   compartments.
%
%   Simulation
%   ===============
%     1. Compute which reactions can proceed based on the amounts of
%        unincorporated and incorporated proteins.
%     2. Compute which proteins can localize to the terminal organelle based on
%        step (1).
%     3. Incorporate the proteins that can localize into the terminal organelle.
%     4. Repeat steps 1-3 until no additional proteins can localize.
%
%   References
%   ===============
%   1. Structure, function, and assembly of the terminal organelle of
%      Mycoplasma pneumoniae (2001). FEMS Microbiology Letters. 198: 1-7.
%   2. Cellular engineering in a minimal microbe: structure and assembly of
%      the terminal organelle of Mycoplasma pneumoniae (2004). Mol
%      Microbiol. 51(4): 917-24.
%   3. Razin S, Jacobs E (1992). Mycoplasma adhesion. J Gen Microbiol.
%      138(3): 407-22. [PUB_0088].
%   4. Balish (2006). Subcellular structures of mycoplasmas. Front Biosci.
%      11: 2017-27. [PUB_0407]
%   5. Chaudhry R, Varshney AK, Malhotra P (2007). Adhesion proteins of
%      Mycoplasma pneumoniae. Front Biosci. 12: 690-9. [PUB_ 0406]
%   6. Balish MF, Krause DC (2006). Mycoplasmas: A Distinct Cytoskeleton
%      for Wall-Less Bacteria. J Mol Microbiol Biotechnol. 11(3-5): 244-55.
%      [PUB_0091]
%   7. Pich OQ, Burgos R, Querol E, Pinol J (2009). P110 and P140 
%      cytadherence-related proteins are negative effectors of terminal 
%      organelle duplication in Mycoplasma genitalium. PLoS One. 4 (10):
%      e7452. [PUB_0794].
%   8. Pich OQ, Burgos R, Ferrer-Navarro M, Querol E, Pinol J (2008). 
%      Role of Mycoplasma genitalium MG218 and MG317 cytoskeletal proteins in 
%      terminal organelle organization, gliding motility and cytadherence. 
%      Microbiology. 154 (Pt 10): 3188-98. [PUB_0803]
%   9. Boonmee A, Ruppert T, Herrmann R (2009). The gene mpn310 (hmw2) from 
%      Mycoplasma pneumoniae encodes two proteins, HMW2 and HMW2-s, which
%      differ in size but use the same reading frame. FEMS Microbiol Lett.
%      290(2): 174-81. [PUB_0804]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef  TerminalOrganelleAssembly < edu.stanford.covert.cell.sim.ReactionProcess
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'localizationReactions';
            'localizationSubstrates';
            'localizationThreshold';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        compartmentIndexs_incorporated   = 1; %index of terminal organelle cytosol and membrane compartments of substrates
        compartmentIndexs_unincorporated = 2; %index of cytosol and membrane compartments of substrates

        substrateWholeCellModelIDs__={ %whole cell model IDs of substrates
            'MG_191_MONOMER';  %MgPa adhesin
            'MG_192_MONOMER';  %P110 protein
            'MG_217_MONOMER';  %P65 adhesin
            'MG_218_MONOMER';  %high molecular weight cytadherence accessory protein 2
            'MG_312_MONOMER';  %high molecular weight cytadherence accessory protein 1
            'MG_317_MONOMER';  %high molecular weight cytadherence accessory protein 3
            'MG_318_MONOMER';  %P32 adhesin
            'MG_386_MONOMER'}; %P200 protein
        substrateIndexs_MgPa = 1; %index within substrates of MgPa adhesin
        substrateIndexs_P110 = 2; %index within substrates of P110 protein
        substrateIndexs_P65  = 3; %index within substrates of P65 adhesin
        substrateIndexs_HMW2 = 4; %index within substrates of high molecular weight cytadherence accessory protein 2
        substrateIndexs_HMW1 = 5; %index within substrates of high molecular weight cytadherence accessory protein 1
        substrateIndexs_HMW3 = 6; %index within substrates of high molecular weight cytadherence accessory protein 3
        substrateIndexs_P32  = 7; %index within substrates of P32 adhesin
        substrateIndexs_P200 = 8; %index within substrates of P200 protein

        enzymeIndexs_HMW1 %index within enzymes of high molecular weight cytadherence accessory protein 1
        enzymeIndexs_HMW2 %index within enzymes of high molecular weight cytadherence accessory protein 2
        enzymeIndexs_HMW3 %index within enzymes of high molecular weight cytadherence accessory protein 3
        enzymeIndexs_P32  %index within enzymes of P32 adhesin
    end

    %fixed biological constants
    properties
        localizationReactions       %adjacency matrix between proteins and the reactions that incorporate them into the terminal organelle [reactions X substrates]
        localizationSubstrates      %adjacency matrix between reactions and the proteins required to carry out that reaction [reactions X substrates X compartment (incorporation/unincorporated)]
        localizationThreshold       %number of proteins required to carry out each reaction [reactions X 1], equals sum(sum(localizationSubstrates,3),2)
    end

    %constructor
    methods
        function this = TerminalOrganelleAssembly(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.ReactionProcess(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            options.retainEnzymeCompartments       = true;
            options.retainModificationCompartments = true;

            this.initializeConstants@edu.stanford.covert.cell.sim.ReactionProcess(...
                knowledgeBase, simulation, options);

            %enzymes
            this.enzymeIndexs_HMW1 = this.enzymeIndexs({'MG_312_MONOMER'}); %high molecular weight cytadherence accessory protein 1
            this.enzymeIndexs_HMW2 = this.enzymeIndexs({'MG_218_MONOMER'}); %high molecular weight cytadherence accessory protein 2
            this.enzymeIndexs_HMW3 = this.enzymeIndexs({'MG_317_MONOMER'}); %high molecular weight cytadherence accessory protein 3
            this.enzymeIndexs_P32  = this.enzymeIndexs({'MG_318_MONOMER'}); %P32 adhesin

            %proteins
            this.substrateWholeCellModelIDs = sort(this.substrateWholeCellModelIDs__);
            [this.substrateNames...
             this.substrateStimulusLocalIndexs...
             this.substrateStimulusGlobalIndexs...
             this.substrateStimulusCompartmentIndexs...
             this.substrateStimulusGlobalCompartmentIndexs...
             this.substrateMetaboliteLocalIndexs...
             this.substrateMetaboliteGlobalIndexs...
             this.substrateMetaboliteCompartmentIndexs...
             this.substrateMetaboliteGlobalCompartmentIndexs...
             this.substrateRNALocalIndexs...
             this.substrateRNAGlobalIndexs...
             this.substrateRNACompartmentIndexs...
             this.substrateRNAGlobalCompartmentIndexs...
             this.substrateMonomerLocalIndexs...
             this.substrateMonomerGlobalIndexs...
             this.substrateMonomerCompartmentIndexs...
             this.substrateMonomerGlobalCompartmentIndexs...
             this.substrateComplexLocalIndexs...
             this.substrateComplexGlobalIndexs...
             this.substrateComplexCompartmentIndexs...
             this.substrateComplexGlobalCompartmentIndexs...
             ~, ~, ~, ...
             this.substrateGlobalIndexs...
             this.substrateMolecularWeights] = this.initializeConstantsHelper(...
                this.substrateWholeCellModelIDs, false);

            this.substrateStimulusGlobalIndexs        = repmat(this.substrateStimulusGlobalIndexs,1,2);
            this.substrateMetaboliteGlobalIndexs      = repmat(this.substrateMetaboliteGlobalIndexs,1,2);
            this.substrateRNAGlobalIndexs             = repmat(this.substrateRNAGlobalIndexs,1,2);
            this.substrateMonomerGlobalIndexs         = repmat(this.substrateMonomerGlobalIndexs,1,2);
            this.substrateComplexGlobalIndexs         = repmat(this.substrateComplexGlobalIndexs,1,2);

            this.substrateStimulusCompartmentIndexs   = repmat(this.substrateStimulusCompartmentIndexs,1,2);
            this.substrateMetaboliteCompartmentIndexs = repmat(this.substrateMetaboliteCompartmentIndexs,1,2);
            this.substrateRNACompartmentIndexs        = repmat(this.substrateRNACompartmentIndexs,1,2);
            this.substrateMonomerCompartmentIndexs    = repmat(this.substrateMonomerCompartmentIndexs,1,2);
            this.substrateComplexCompartmentIndexs    = repmat(this.substrateComplexCompartmentIndexs,1,2);

            this.substrateMonomerCompartmentIndexs(this.substrateMonomerCompartmentIndexs(:,this.compartmentIndexs_unincorporated) == this.compartment.terminalOrganelleCytosolIndexs, this.compartmentIndexs_unincorporated)  = this.compartment.cytosolIndexs;
            this.substrateMonomerCompartmentIndexs(this.substrateMonomerCompartmentIndexs(:,this.compartmentIndexs_unincorporated) == this.compartment.terminalOrganelleMembraneIndexs, this.compartmentIndexs_unincorporated) = this.compartment.membraneIndexs;
            this.substrateComplexCompartmentIndexs(this.substrateComplexCompartmentIndexs(:,this.compartmentIndexs_unincorporated) == this.compartment.terminalOrganelleCytosolIndexs, this.compartmentIndexs_unincorporated)  = this.compartment.cytosolIndexs;
            this.substrateComplexCompartmentIndexs(this.substrateComplexCompartmentIndexs(:,this.compartmentIndexs_unincorporated) == this.compartment.terminalOrganelleMembraneIndexs, this.compartmentIndexs_unincorporated) = this.compartment.membraneIndexs;
            
            this.substrateStimulusGlobalCompartmentIndexs = sub2ind(...
                [numel(this.stimulus.wholeCellModelIDs) this.compartment.count],...
                this.substrateStimulusGlobalIndexs,...
                this.substrateStimulusCompartmentIndexs);
            this.substrateMetaboliteGlobalCompartmentIndexs = sub2ind(...
                [numel(this.metabolite.wholeCellModelIDs) this.compartment.count],...
                this.substrateMetaboliteGlobalIndexs,...
                this.substrateMetaboliteCompartmentIndexs);
            this.substrateRNAGlobalCompartmentIndexs = sub2ind(...
                [numel(this.rna.wholeCellModelIDs) this.compartment.count],...
                this.rna.matureIndexs(this.substrateRNAGlobalIndexs),...
                this.substrateRNACompartmentIndexs);
            this.substrateMonomerGlobalCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count],...
                this.monomer.matureIndexs(this.substrateMonomerGlobalIndexs),...
                this.substrateMonomerCompartmentIndexs);
            this.substrateComplexGlobalCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count],...
                this.complex.matureIndexs(this.substrateComplexGlobalIndexs),...
                this.substrateComplexCompartmentIndexs);

            %Localization requirements (see references)
            %localizationSubstrates(i,j): unincorporated terminalOrganelleProtein j localized in reaction i
            %localizationReactions(i,j,k): protein j in compartment k required for localization reaction i
            %localizationThreshold: number of protein species that must be present for localization of protein; sum of dims 2, 3 of localizationReactions
            substrateIndexs = zeros(size(this.substrateWholeCellModelIDs));
            substrateIndexs(this.substrateMonomerLocalIndexs) = this.gene.mRNAIndexs(this.substrateMonomerGlobalIndexs(:,1));
            substrateIndexs(this.substrateComplexLocalIndexs) = this.substrateComplexGlobalIndexs(:,1) + knowledgeBase.numGenes;
            this.localizationSubstrates = sum(isnan(this.reactionModificationMatrix(:, substrateIndexs, :)),3);

            numProteins  = length(this.substrateWholeCellModelIDs__);
            numReactions = size(this.reactionCatalysisMatrix,1);
            this.localizationReactions = zeros(numReactions, numProteins, 2);

            for i = 1:numReactions
                [enzymeLocalIndex, enzymeGlobalCompartmentIndex] = find(permute(...
                    this.reactionCatalysisMatrix(i,:,:), [2 3 1]));

                if isempty(enzymeLocalIndex); continue; end;

                if any(this.enzymeMonomerLocalIndexs == enzymeLocalIndex)
                    enzymeGlobalIndex = this.enzymeGlobalIndexs(enzymeLocalIndex);
                    proteinLocalIndex = this.substrateMonomerLocalIndexs(this.substrateMonomerGlobalIndexs(:,1) == enzymeGlobalIndex);
                    proteinLocalCompartmentIndex = find(this.substrateMonomerCompartmentIndexs(proteinLocalIndex,:) == enzymeGlobalCompartmentIndex);
                else
                    enzymeGlobalIndex = this.enzymeGlobalIndexs(enzymeLocalIndex);
                    proteinLocalIndex = this.substrateComplexLocalIndexs(this.substrateComplexGlobalIndexs(:,1) == enzymeGlobalIndex);
                    proteinLocalCompartmentIndex = find(this.substrateComplexCompartmentIndexs(proteinLocalIndex,:) == enzymeGlobalCompartmentIndex);
                end

                this.localizationReactions(i, proteinLocalIndex, proteinLocalCompartmentIndex) = 1;
            end

            this.localizationThreshold = sum(sum(this.localizationReactions,3),2);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %no substrate and byproducts required
            
            %at least 1 copy of the foundational proteins
            minEnzExp([
                this.enzymeIndexs_HMW1;
                this.enzymeIndexs_HMW2;
                this.enzymeIndexs_HMW3;
                this.enzymeIndexs_P32]) = 1;
        end       

        %initialization: proteins initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end

        %simulation
        function evolveState(this)
            % incorporate proteins into terminal organelle
            % cytoplasm -> terminal organelle cytoplasma
            % membrane  -> terminal organelle membrane
            while any(this.substrates(:, this.compartmentIndexs_unincorporated))
                activeReactions = ...
                    this.localizationReactions(:, :, 1) * (this.substrates(:, 1) > 0) + ...
                    this.localizationReactions(:, :, 2) * (this.substrates(:, 2) > 0) >= ...
                    this.localizationThreshold;
                
                localizingProteins = logical(...
                    this.substrates(:, this.compartmentIndexs_unincorporated) .* ...
                    (this.localizationSubstrates' * activeReactions));
                if ~any(localizingProteins)
                    break; 
                end

                this.substrates(localizingProteins, this.compartmentIndexs_incorporated) = ...
                    this.substrates(localizingProteins, this.compartmentIndexs_incorporated) + ...
                    this.substrates(localizingProteins, this.compartmentIndexs_unincorporated);
                this.substrates(localizingProteins, this.compartmentIndexs_unincorporated) = 0;
            end
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.substrates, 3) == 1
                value = this.substrateMolecularWeights' * sum(this.substrates, 2) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = permute(...
                    this.substrateMolecularWeights' * permute(sum(this.substrates,2), [1 3 2]),...
                    [1 3 2]) / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
