%Protein Folding
%
% @wholeCellModelID Process_ProteinFolding
% @name             Protein Folding
% @description
%   Biology
%   ==============
%   Protein folding is the process whereby proteins which are produced as linear
%   amino acid polymers relax to their most energetically favorable, and often
%   more compact and catalytically active configuration. Some proteins relax
%   quickly and spontaneously to their free energy minimum, while other proteins
%   require chaperones to help them achieve their folded configuration more
%   quickly. Additionally, some proteins coordinate prosthetic groups such as
%   metal ions and other small molecules during this protein rearrangement
%   process. Prosthetic groups may help stabilize the protein's catalytically
%   active configuration, or may perform folding process.
%
%   In addition to protein chaperones, membrane protein folding may be assisted
%   by phosphatidyl ethanolamine (PE) and phosphatidyl glycerol (PG) [PUB_0646].
%   This effect is not currently modeled by this process.
%
%   This process simulates
%   1. Prosthetic group complexation
%      - Inorganic ions: Fe2, Fe3, K, Mg, Mn, Na, Zn
%   2. Protein folding, where required mediated by ATP-dependent chaperone(s)
%      [PUB_0014, PUB_0644]:
%      a. Trigger factor (Tig, MG_238_MONOMER) – Highly expressed protein which
%         possibly interacts with the L23 subunit of every ribosome,
%         co-translationally assisting in the early folding every protein
%         [PUB_0014, PUB_0644]. Binds vacant ribosomes with half-life of 11 s
%         [PUB_0005]; binds ribosome-peptide complexes with half-life of 15-50 s
%         [PUB_0005]; remains associated with the peptide for a up to more
%         seconds after peptide release from the ribosome [PUB_0005]. Trigger
%         factor-50S ribosome affinity: KD=1 ?M [PUB_0005]. We model trigger as
%         required for the proper folding of all proteins.
%      b. Chaperone (DnaK, MG_305_MONOMER) – Folds 5-18% of proteins [PUB_0014,
%         PUB_0644]. Typically folds proteins > 30 kDa, typically in less than 2
%         min [PUB_0014, PUB_0644]. Binds backbone of short, linear, unfolded,
%         hydrophobic peptide segments [PUB_0014, PUB_0644]. ATP hydrolysis and
%         peptide release catalyzed by GrpE (MG_201_DIMER) [PUB_0014, PUB_0644].
%         Regulated by co-chaperone DnaJ (MG_019_DIMER) which interacts with the
%         peptide side chains [PUB_0014, PUB_0644].
%      c. Chaperonin (GroEL, GroES, MG_392_393_21MER) – Folds 10-15% proteins by
%         ATP-dependent mechanism [PUB_0014, PUB_0644]. Typically folds proteins
%         20-60 kDa in size [PUB_0014, PUB_0644] with half-life > 10 min
%         [PUB_0014, PUB_0644], 30-60 s [PUB_0389].
%
%   Knowledge Base
%   ==============
%   The chaperones (except Tig, which all protein monomers require) and
%   prosthetic groups required to fold each protein were reconstructed from the
%   literature (see following sections), organized into the knowledge base, and
%   are retrieved from the knowledge base in initializeConstants. Where known
%   the stoichiometry of each prosthetic group was stored in the knowledge base;
%   prosthetic groups with unreported stoichiometry are indicated with
%   stoichiometry values of -1.
%
%   Chaperones
%   ++++++++++++++
%   The chaperone requirements for the folding of each protein were compiled
%   for other bacteria (see sources below), and mapped to M. genitalium by
%   homology.
%
%       Tig    may interact with all nascent peptides at the ribosome
%              polypeptide exit site [PUB_0005, PUB_0009, PUB_0388] and assist
%              in early folding. Binds vacant ribosomes with half-life of 11 s
%              [PUB_0005]; binds ribosome-peptide complexes with half-life of
%              15-50 s [PUB_0005]; remains associated with the peptide for a
%              up to more seconds after peptide release from the ribosome
%              [PUB_0005]. Trigger factor-50S ribosome affinity: KD=1 ?M
%              [PUB_0005]. We model trigger as required for the proper folding
%              of all proteins.
%       DnaK   Deuerling et al performed a proteome-scale search for DnaK
%              substrates in E. coli [PUB_0388]
%       GroEL  Kerner et al performed a proteome-scale search for GroEL
%              substrates in E. coli [PUB_0389] and Endo and Kurusu performed
%              a proteome-scale search for GroEL substrates in B. subtilis
%              [PUB_0391]
%       FtsH   may act as molecule chaperone for membrane proteins
%              [PUB_0014]. Several other functions have also been associated
%              with FtsH. To date no proteome-scale studies of FtsH activity
%              has been performed, and FtsH's chaperone substrates are not
%              well characterized. Conequently we chose not to model FtsH as
%              a molecular chaperone, but rather a protease.
%
%   Substrates of SecB, a molecular chaperone not present in M. genitalium, have
%   also been identified on proteome-scale in E. coli [PUB_0390].
%
%   The reconstruction found that complexes typically require no chaperones to
%   fold, with the notable exception of MG_392_393_21MER (which itself is a
%   chaperone), which requires four.
%
%   Prosthetic groups
%   ++++++++++++++
%   Reaction coenzymes and prosthetic groups were reconstructed from several
%   sources including the databases BioCyc [PUB_0006], BRENDA [PUB_0570],
%   GenoBase [PUB_0386], Kinetikon [PUB_0571], Metal-MaCiE [PUB_0387], and
%   UniProt [PUB_0096] and the primary literature [PUB_0131]. These sources
%   frequently report protein prosthetic groups as a list of chemical species
%   (eg. metal ions) which can each fill the same prosthetic group role. We
%   simplified these ambiguous lists of possible prosthetic groups to a single
%   chemical species by choosing the most common chemical species according to
%   cell composition (Na+, K+ > Mg2+, Cl- > Fe2+, Fe3+ > Ca2+ > Mn2+ > Cu2+ >
%   Mo6+, Zn2, Co2, Ni2 [PUB_0393, PUB_0394, PUB_0395]), or the species with the
%   highest protein affinity (Irving-Williams Series, Mn2+ < Fe2+, Fe3+ < Co2+ <
%   Ni2+ < Cu2+ > Zn2+ [PUB_0404]).
%
%   Roughly 20% of the monomers and 1% of the complexes bind at least one
%   prosthetic group.
%
%   Representation
%   ==============
%   substrates, enzymes, unfoldedMonomers, unfoldedComplexs, foldedMonomers, and
%   foldedComplexs represent the counts of free metabolites, chaperones,
%   unfolded protein monomers, unfolded protein complexes, folded protein
%   monomers, and folded protein complexes. Each of these component types have
%   been mapped from the several compartments of the simulation class to one
%   pseudo compartment in this process. That is this process only "sees" the
%   counts of metabolites, chaperones, and protein monomers and complexes that
%   are relevant to protein folding. The process contains no intermediate
%   representation of protein folding, that is protein folding is simulated as a
%   all-or-nothing process that either does not occur or proceeds to completion
%   within a single time step.
%
%   The proteinChaperoneMatrix and proteinProstheticGroupMatrix are adjacency
%   matrices which represent the chaperone and prosthetic groups each protein
%   monomer and complex requires to achieve its catalytically active
%   configuration. Where the stoichiometry of a prosthetic group is unknown and
%   indicated in the knowledge base with the value -1, here we assume a
%   stoichiometry of 1 and set the element of proteinProstheticGroupMatrix
%   accordingly.
%
%   Initialization
%   ==============
%   All protein monomers and complexes are initialized to the mature state. This
%   is accomplished by the simulation class initializeState method.
%
%   Simulation
%   ==============
%   As with the other protein maturation processes, protein folding is modeled as
%   as all-or-nothing process at the single time-step level. A protein only
%   proceeds through this phase if its chaperones (if any) and any prosthetic
%   group ions are all present in the same time step. Chaperone kinetics are not
%   well characterized, and are not currrently modeled. Proteins that require
%   neither prosthetic groups nor chaperones are considered folded after one
%   time step.
%
%   Algorithm
%   ++++++++++++++
%   While(true)
%     1. Calculate numbers of proteins that can fold based on prosthetic group,
%        chaperone, and unfolded protein availability.
%     2. Randomly select proteins to fold weighted by limits calculated in step
%        (1).
%     3. Update counts of prosthetic groups, chaperones, and unfolded and folded
%        proteins.
%     4. Repeat until insufficient resources to further fold proteins
%   End
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef ProteinFolding < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'proteinProstheticGroupMatrix';
            'proteinChaperoneMatrix';
            'monomerCompartments';
            'complexCompartments';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unfoldedMonomers';
            'unfoldedComplexs';
            'foldedMonomers';
            'foldedComplexs'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ATP';'ADP';'PI';'H2O';'H'};
        substrateIndexs_adp       %index within substrates of ADP
        substrateIndexs_atp       %index within substrates of ATP
        substrateIndexs_fe2       %index within substrates of Fe+2 ion
        substrateIndexs_hydrogen  %index within substrates of hydrogen
        substrateIndexs_water     %index within substrates of water
        substrateIndexs_k         %index within substrates of K+ ion
        substrateIndexs_mg        %index within substrates of Mg2+ ion
        substrateIndexs_mn        %index within substrates of Mn2+ ion
        substrateIndexs_sodium    %index within substrates of Na+ ion
        substrateIndexs_phosphate %index within substrates of inorganic phosphate
        substrateIndexs_zinc      %index within substrates of Zn2+ ion

        enzymeWholeCellModelIDs = {     %enzyme whole cell model ids
            'MG_238_MONOMER';           %trigger factor
            'MG_305_MONOMER';           %chaperone
            'MG_019_DIMER';             %chaperone protein DnaJ
            'MG_201_DIMER';             %co-chaperone GrpE
            'MG_392_393_21MER'};        %GroEL-GroES chaperonin complex
        enzymeIndexs_triggerFactor      %index of trigger factor within enzymeWholeCellModelIDs
        enzymeIndexs_dnaK               %index of DnaK within enzymeWholeCellModelIDs
        enzymeIndexs_dnaJ               %index of DnaJ within enzymeWholeCellModelIDs
        enzymeIndexs_grpE               %index of GrpE within enzymeWholeCellModelIDs
        enzymeIndexs_groELES            %index of GroELES within enzymeWholeCellModelIDs

        unfoldedMonomerWholeCellModelIDs %whole cell model IDs of unfolded monomers
        unfoldedComplexWholeCellModelIDs %whole cell model IDs of unfolded complexes
        foldedMonomerWholeCellModelIDs   %whole cell model IDs of folded monomers
        foldedComplexWholeCellModelIDs   %whole cell model IDs of folded complexes
        
        unfoldedMonomerCompartmentIndexs %indices of unfolded monomers within this.monomer.counts
        unfoldedComplexCompartmentIndexs %indices of unfolded complexs within this.complex.counts
        foldedMonomerCompartmentIndexs   %indices of folded monomers within this.monomer.counts
        foldedComplexCompartmentIndexs   %indices of folded complexs within this.complex.counts
                
        monomerIndexs_folding            %indices of monomers that actively fold within unfoldedMonomers
        complexIndexs_folding            %indices of complexs that actively fold within unfoldedComplexs
        monomerIndexs_notFolding         %indices of monomers that don't actively fold within unfoldedMonomers
        complexIndexs_notFolding         %indices of complexs that don't actively fold within unfoldedComplexs
        
        monomerComplexIndexs_folded      %indices of monomers and complexs that actively fold within proteinProstheticGroupMatrix
        
        speciesIndexs_monomers           %indices of monomers that actively fold within speciesReactantByproductMatrix
        speciesIndexs_complexs           %indices of complexs that actively fold within speciesReactantByproductMatrix
    end
    
    %fixed biological constants
    properties
        proteinProstheticGroupMatrix     %proteins X metabolites X compartments
        proteinChaperoneMatrix           %proteins X proteins (chaperone) X compartments

        monomerCompartments              %compartment indices with simulation.compartmentWholeCellModelIDs of monomers
        complexCompartments              %compartment indices with simulation.compartmentWholeCellModelIDs of complexes
                            
        speciesReactantByproductMatrix   %stoichiometry of susbtrates, enzymes, monomers and complexs in reactions necessary to fold each protein: [Protein foldings] X [subtrates; enzymes; unfolded monomers; unfolded complexs]
        speciesReactantMatrix            %reactant stoichiometry of susbtrates, enzymes, monomers and complexs in reactions necessary to fold each protein: [Protein foldings] X [subtrates; enzymes; unfolded monomers; unfolded complexs]
    end

    %global state (stored locally for convenience)
    properties
        unfoldedMonomers                 %amount of unfolded monomers
        unfoldedComplexs                 %amount of unfolded complexes
        foldedMonomers                   %amount of folded monomers
        foldedComplexs                   %amount of folded complexes
    end

    %constructor
    methods
        function this = ProteinFolding(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            %prosthetic groups
            this.proteinProstheticGroupMatrix = knowledgeBase.proteinProstheticGroupMatrix(:, :, this.compartment.cytosolIndexs);
            this.proteinProstheticGroupMatrix(this.proteinProstheticGroupMatrix==-1)=1;

            this.substrateWholeCellModelIDs = unique([...
                this.substrateWholeCellModelIDs;
                this.metabolite.wholeCellModelIDs(sum(this.proteinProstheticGroupMatrix)>0)]);

            this.substrateIndexs_adp       = this.substrateIndexs({'ADP'});
            this.substrateIndexs_atp       = this.substrateIndexs({'ATP'});
            this.substrateIndexs_fe2       = this.substrateIndexs({'FE2'});
            this.substrateIndexs_hydrogen  = this.substrateIndexs({'H'});
            this.substrateIndexs_water     = this.substrateIndexs({'H2O'});
            this.substrateIndexs_k         = this.substrateIndexs({'K'});
            this.substrateIndexs_mg        = this.substrateIndexs({'MG'});
            this.substrateIndexs_mn        = this.substrateIndexs({'MN'});
            this.substrateIndexs_sodium    = this.substrateIndexs({'NA'});
            this.substrateIndexs_phosphate = this.substrateIndexs({'PI'});
            this.substrateIndexs_zinc      = this.substrateIndexs({'ZN'});

            %chaperones
            numMonomers = length(this.monomer.matureIndexs);
            proteinChaperone = knowledgeBase.proteinChaperoneMatrix;
            proteinMonomerChaperoneIndexs = logical(sum(sum(proteinChaperone(1:numMonomers,:,:),3),2));
            proteinComplexChaperoneIndexs = logical(sum(sum(proteinChaperone(numMonomers+1:end,:,:),3),2));
            this.enzymeWholeCellModelIDs = unique([
                this.enzymeWholeCellModelIDs
                this.monomer.wholeCellModelIDs{this.monomer.matureIndexs(proteinMonomerChaperoneIndexs)}
                this.complex.wholeCellModelIDs{this.complex.matureIndexs(proteinComplexChaperoneIndexs)}
                ]);

            this.enzymeIndexs_triggerFactor = this.enzymeIndexs({'MG_238_MONOMER'});
            this.enzymeIndexs_dnaK          = this.enzymeIndexs({'MG_305_MONOMER'});
            this.enzymeIndexs_dnaJ          = this.enzymeIndexs({'MG_019_DIMER'});
            this.enzymeIndexs_grpE          = this.enzymeIndexs({'MG_201_DIMER'});
            this.enzymeIndexs_groELES       = this.enzymeIndexs({'MG_392_393_21MER'});

            %call super class method to get metabolite, enzyme indexs
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});

            %prosthetic groups
            this.proteinProstheticGroupMatrix = ...
                this.proteinProstheticGroupMatrix(:,this.substrateMetaboliteGlobalIndexs);

            %chaperones
            numEnzymes    = length(this.enzymeWholeCellModelIDs);
            numSubstrates = size(proteinChaperone,2);            

            enzymeCompartmentIndexs = zeros(size(this.enzymeWholeCellModelIDs));
            enzymeCompartmentIndexs(this.enzymeMonomerLocalIndexs) = this.enzymeMonomerCompartmentIndexs;
            enzymeCompartmentIndexs(this.enzymeComplexLocalIndexs) = this.enzymeComplexCompartmentIndexs;
            
            enzymeGlobalIndexs = this.enzymeGlobalIndexs;
            enzymeGlobalIndexs(this.enzymeComplexLocalIndexs) = ...
                enzymeGlobalIndexs(this.enzymeComplexLocalIndexs) + ...
                numMonomers;
            this.proteinChaperoneMatrix = proteinChaperone(sub2ind(...
                size(proteinChaperone),...
                repmat(enzymeGlobalIndexs',numSubstrates,1), ...
                repmat((1:numSubstrates)',[1 numEnzymes]),...
                repmat(enzymeCompartmentIndexs',numSubstrates,1)));

            this.proteinChaperoneMatrix(1:numMonomers,this.enzymeIndexs_triggerFactor) = 1;
            this.proteinChaperoneMatrix(:,this.enzymeIndexs_dnaJ) = this.proteinChaperoneMatrix(:,this.enzymeIndexs_dnaK);
            this.proteinChaperoneMatrix(:,this.enzymeIndexs_grpE) = this.proteinChaperoneMatrix(:,this.enzymeIndexs_dnaK);

            %monomers
            this.unfoldedMonomerWholeCellModelIDs = this.monomer.wholeCellModelIDs(this.monomer.processedIIIndexs);
            this.foldedMonomerWholeCellModelIDs   = this.monomer.wholeCellModelIDs(this.monomer.foldedIndexs);
            this.monomerCompartments = double(knowledgeBase.proteinMonomerCompartments);
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleCytosolIndexs)  = this.compartment.cytosolIndexs;
            this.monomerCompartments(this.monomerCompartments == this.compartment.terminalOrganelleMembraneIndexs) = this.compartment.membraneIndexs;
            this.unfoldedMonomerCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count], ...
                this.monomer.processedIIIndexs, ...
                this.monomerCompartments);
            this.foldedMonomerCompartmentIndexs = sub2ind(...
                [numel(this.monomer.wholeCellModelIDs) this.compartment.count], ...
                this.monomer.foldedIndexs, ...
                this.monomerCompartments);

            %complexs
            this.unfoldedComplexWholeCellModelIDs = this.complex.wholeCellModelIDs(this.complex.nascentIndexs);
            this.foldedComplexWholeCellModelIDs   = this.complex.wholeCellModelIDs(this.complex.matureIndexs);
            this.complexCompartments = double(knowledgeBase.proteinComplexCompartments);
            this.unfoldedComplexCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count], ...
                this.complex.nascentIndexs, ...
                this.complexCompartments);
            this.foldedComplexCompartmentIndexs = sub2ind(...
                [numel(this.complex.wholeCellModelIDs) this.compartment.count], ...
                this.complex.matureIndexs, ...
                this.complexCompartments);
            
            this.initializeSpeciesNetwork();
        end
        
        function initializeSpeciesNetwork(this)
            this.monomerComplexIndexs_folded = find(sum(this.proteinProstheticGroupMatrix,2)> 0 | sum(this.proteinChaperoneMatrix,2)> 0);
            
            numMonomers   = length(this.unfoldedMonomerWholeCellModelIDs);
            numComplexs   = length(this.unfoldedComplexWholeCellModelIDs);            
            
            [~, this.monomerIndexs_folding, this.speciesIndexs_monomers] = intersect(1:numMonomers, this.monomerComplexIndexs_folded); %use of intersect is safe because the inputs are already sorted
            [~, this.complexIndexs_folding, this.speciesIndexs_complexs] = intersect((1:numComplexs) + numMonomers, this.monomerComplexIndexs_folded);
            [~, this.monomerIndexs_notFolding] = setdiff(1:numMonomers,                 this.monomerComplexIndexs_folded);
            [~, this.complexIndexs_notFolding] = setdiff((1:numComplexs) + numMonomers, this.monomerComplexIndexs_folded);
            
            this.speciesReactantByproductMatrix = [
                this.proteinProstheticGroupMatrix(this.monomerComplexIndexs_folded,:)...
                this.proteinChaperoneMatrix(this.monomerComplexIndexs_folded,:)...
                eye(length(this.monomerComplexIndexs_folded))];
            this.speciesReactantMatrix = max(0, this.speciesReactantByproductMatrix);
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();

            numTimePoints = size(this.unfoldedMonomers, 3);
            
            if numTimePoints == 1
                this.unfoldedMonomers = this.monomer.counts(this.unfoldedMonomerCompartmentIndexs);
                this.unfoldedComplexs = this.complex.counts(this.unfoldedComplexCompartmentIndexs);
                this.foldedMonomers = this.monomer.counts(this.foldedMonomerCompartmentIndexs);
                this.foldedComplexs = this.complex.counts(this.foldedComplexCompartmentIndexs);
            else
                this.unfoldedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedMonomers,1) 1]),[2 3 1])));
                this.unfoldedComplexs = this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complex.nascentIndexs, [1 1 numTimePoints]), ...
                    repmat(this.complexCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedComplexs,1) 1]),[2 3 1])));
                this.foldedMonomers = this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.foldedIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedMonomers,1) 1]),[2 3 1])));
                this.foldedComplexs = this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complex.matureIndexs, [1 1 numTimePoints]), ...
                    repmat(this.complexCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedComplexs,1) 1]),[2 3 1])));
            end
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();
            
            numTimePoints = size(this.unfoldedMonomers, 3);
            
            if numTimePoints == 1
                this.monomer.counts(this.unfoldedMonomerCompartmentIndexs) = this.unfoldedMonomers;
                this.complex.counts(this.unfoldedComplexCompartmentIndexs) = this.unfoldedComplexs;
                this.monomer.counts(this.foldedMonomerCompartmentIndexs) = this.foldedMonomers;
                this.complex.counts(this.foldedComplexCompartmentIndexs) = this.foldedComplexs;
            else                
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.processedIIIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedMonomers,1) 1]),[2 3 1]))) = ...
                    this.unfoldedMonomers;
                this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complex.nascentIndexs, [1 1 numTimePoints]), ...
                    repmat(this.complexCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedComplexs,1) 1]),[2 3 1]))) = ...
                    this.unfoldedComplexs;
                this.monomer.counts(sub2ind(...
                    size(this.monomer.counts), ...
                    repmat(this.monomer.foldedIndexs, [1 1 numTimePoints]), ...
                    repmat(this.monomerCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedMonomers,1) 1]),[2 3 1]))) = ...
                    this.foldedMonomers;
                this.complex.counts(sub2ind(...
                    size(this.complex.counts), ...
                    repmat(this.complex.matureIndexs, [1 1 numTimePoints]), ...
                    repmat(this.complexCompartments,[1 1 numTimePoints]),...
                    permute(repmat((1:numTimePoints)',[1 size(this.foldedComplexs,1) 1]),[2 3 1]))) = ...
                    this.foldedComplexs;
            end
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.unfoldedMonomers = zeros(length(this.unfoldedMonomerWholeCellModelIDs), 1, numTimePoints);
            this.unfoldedComplexs = zeros(length(this.unfoldedComplexWholeCellModelIDs), 1, numTimePoints);
            this.foldedMonomers   = zeros(length(this.foldedMonomerWholeCellModelIDs),   1, numTimePoints);
            this.foldedComplexs   = zeros(length(this.foldedComplexWholeCellModelIDs),   1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% substrate and byproducts
            %prosthetic groups that complex with proteins
            bmProd = ...
                this.proteinProstheticGroupMatrix' * ...
                [states.monomerProductions; states.complexProductions];
            
            %no byproducts
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            
            %% enzymes: at least one of each enzyme
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            minEnzExp([
                this.enzymeIndexs_triggerFactor;
                this.enzymeIndexs_dnaK;
                this.enzymeIndexs_dnaJ;
                this.enzymeIndexs_grpE;
                this.enzymeIndexs_groELES]) = 2;
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization
        %- monomers initialized to mature/bound/inactivated state
        %  by simulation initializeState method
        %- complexs initialized by
        %  1. Macromolecular complexation
        %  2. Ribosome assembly
        %  3. Protein folding
        %  4. Protein activation
        %  Here we push complexs to folded state
        function initializeState(this)
            this.foldedComplexs = this.foldedComplexs + this.unfoldedComplexs;
            this.unfoldedComplexs(:) = 0;
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            if all(this.enzymes)
                result = this.proteinProstheticGroupMatrix' * ...
                    [this.unfoldedMonomers; this.unfoldedComplexs];
            else
                result = this.proteinProstheticGroupMatrix' * ...
                    ((min(this.enzymes(:, ones(size(this.proteinChaperoneMatrix, 1), 1))' ./ this.proteinChaperoneMatrix, [], 2) > 0) .* ...
                    [this.unfoldedMonomers; this.unfoldedComplexs]);
            end
        end

        %simulation
        function evolveState(this)
            %update monomers that don't require chaperones or prosthetic groups
            this.foldedComplexs(this.complexIndexs_notFolding) = ...
                this.foldedComplexs(this.complexIndexs_notFolding) + ...
                this.unfoldedComplexs(this.complexIndexs_notFolding);
            this.unfoldedComplexs(this.complexIndexs_notFolding) = 0;
            
            %stop early if no proteins to fold
            if ~any(this.unfoldedMonomers) && ~any(this.unfoldedComplexs)
                return;
            end
                        
            numSubstrates = length(this.monomerComplexIndexs_folded);
            
            species = max(0, [
                this.substrates;
                this.enzymes * Inf;
                this.unfoldedMonomers;
                this.unfoldedComplexs(this.complexIndexs_folding)]');

            fluxes = zeros(size(this.speciesReactantByproductMatrix, 1),1);
            anyFlux = false;

            % simulate modification
            while true
                %compute maximum number of each species that can be modified
                substrateLimits = species(ones(numSubstrates,1), :) ./ this.speciesReactantMatrix;
                substrateLimits(:, this.substrateIndexs_water) = NaN;
                substrateLimits(:, this.substrateIndexs_hydrogen) = NaN;
                substrateLimits = min(substrateLimits, [], 2)';
                substrateLimits(isinf(substrateLimits) | isnan(substrateLimits) | substrateLimits<1)=0;

                %stop if no more substrates can be modified
                if ~any(substrateLimits); break; end;
                
                anyFlux = true;

                %pick substrate, reaction, enzymes
                selectedSubstrate = this.randStream.randsample(numel(substrateLimits), 1, true, substrateLimits);
                fluxes(selectedSubstrate) = fluxes(selectedSubstrate) + 1;

                %decrement metabolites, enzymes, unmodified substrates; increment modified substrates
                species = species - this.speciesReactantByproductMatrix(selectedSubstrate, :);
            end
            
            %stop if no reactions can proceed
            if ~anyFlux
                return;
            end

            % compute results
            this.substrates = this.substrates - this.proteinProstheticGroupMatrix(this.monomerComplexIndexs_folded, :)' * fluxes;

            this.unfoldedMonomers = this.unfoldedMonomers - fluxes(this.speciesIndexs_monomers);
            this.foldedMonomers   = this.foldedMonomers   + fluxes(this.speciesIndexs_monomers);

            this.unfoldedComplexs(this.complexIndexs_folding) = ...
                this.unfoldedComplexs(this.complexIndexs_folding) - ...
                fluxes(this.speciesIndexs_complexs);
            this.foldedComplexs(this.complexIndexs_folding) = ...
                this.foldedComplexs(this.complexIndexs_folding) + ...
                fluxes(this.speciesIndexs_complexs);
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unfoldedComplexs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.monomer.molecularWeights(this.monomer.processedIIIndexs)' * this.unfoldedMonomers + ...
                    this.monomer.molecularWeights(this.monomer.foldedIndexs)'      * this.foldedMonomers + ...                
                    this.complex.molecularWeights(this.complex.nascentIndexs)'     * this.unfoldedComplexs + ...
                    this.complex.molecularWeights(this.complex.matureIndexs)'      * this.foldedComplexs) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.monomer.molecularWeights(this.monomer.processedIIIndexs)' * permute(this.unfoldedMonomers,[1 3 2]),[1 3 2]) + ...
                    permute(this.monomer.molecularWeights(this.monomer.foldedIndexs)'      * permute(this.foldedMonomers,  [1 3 2]),[1 3 2]) + ...      
                    permute(this.complex.molecularWeights(this.complex.nascentIndexs)'     * permute(this.unfoldedComplexs,[1 3 2]),[1 3 2]) + ...
                    permute(this.complex.molecularWeights(this.complex.matureIndexs)'      * permute(this.foldedComplexs  ,[1 3 2]),[1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end