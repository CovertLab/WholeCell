%Translation
%
% @wholeCellModelID Process_Translation
% @name             Translation
% @description
%   Biology
%   ===============
%   Translation is the first step in the synthesis of function proteins whereby
%   the ribosome, accessory enzymes, and tRNAs transcode mRNAs and produce
%   amino acid polymers. Following translation proteins must be matured to
%   become fully functional:
%    1. Protein Processing I: N-terminal methionine deformylation and cleavage
%    2. Protein Translocation: Translocation of integral membrane proteins,
%       lipoproteins, and secreted proteins into and through the cell membrane
%    3. Protein Processing II: Diacylglyceryl addition to lipoproteins. Signal
%       sequence cleavage of lipo- and secreted proteins.
%    4. Protein Folding: Folding and prosthetic group coordination of proteins.
%    5. Protein Modification: Addition of covalently attached chemical groups to
%       proteins.
%    6. Macromolecular Complexation: Incorporation of protein monomers into
%       larger macromolecular assemblies.
%    7. Ribosome Assembly: Special case of macromolecular complexation of the
%       ribosome which requires GTPases.
%    8. Terminal Organelle Assembly: Special case of translocation which
%       sequesters several protein species in a membrane-bound bleb referred to
%       as the terminal organelle.
%    9. Protein Activation: Functional activation and inactivation of proteins
%       in response to their surrounding chemical environment.
%
%   Translation begins with the recruitment of the 30S and 50S ribosomal
%   particles and initiation factor 3 (IF3) to an mRNA molecule. Next the
%   ribosomal particles assembly into a 70S ribosome on the mRNA molecule with
%   the help of initiation factors. Third the ribosome polymerizes amino
%   acids with the help of elongation factors. Finally release factor
%   (MG_258) recognizes the stop codon UAG, release factor hydolyzes the
%   peptidyl tRNA bond and dissociates, ribosome recycling factor
%   dissociates the E-site tRna, elongation factor G release the release
%   factor and 50S ribosome, and initiation factor 3 dissociates the 30S
%   ribosome, P-site tRNA, and mRNA.
%
%   When a ribosome stalls its lastest tRNA is expelled and replaced by the
%   tRNA-like domain of a tmRNA molecule. Next the mRNA-like domain of the tmRNA
%   expels the bound mRNA. Third the ribosome resumes polymerization, now
%   using the tmRNA's mRNA-like domain as its template. This results in the
%   production of an amino acid polymer containing a C-terminal proteolysis tag.
%   Finally, the proteolysis tag will be recognized by the protein degradation
%   machinery, and the amino acid polymer will be degraded into its individual
%   component amino acids.
%
%   This process simulates protein translation by ribosomes and acessory
%   initiation, elongation, and termination factors. The process also simulates
%   the identification of stalled ribosomes by the tmRNA, the replacement of the
%   tRNA and mRNA with the tmRNA, and the synthesis of the proteolysis tag
%   encoded by the tmRNA's mRNA-like domain.
%
%   Knowledge Base
%   ===============
%   The knowledge base contains the region of the genome each mRNA species
%   corresponds to, and the location of the tmRNA mRNA-like domain on the
%   genome. This information is converted into a sequence of tRNAs which are
%   required to polymerize each protein and proteolysis tag by the knowledge
%   base class' computeTRNASequences method.
%
%   Representation
%   ===============
%   The properties substrates, enzymes, monomers, mRNAs, freeTRNAs,
%   aminoacylatedTRNAs, and aminoacylatedTMRNA represents the counts of
%   metabolites, translation enzymes, nascent peptides, mature mRNAs, free
%   tRNAs, aminoacylated tRNAs, and aminoacylated tmRNAs.
%
%   monomerTRNASequences and tmRNAProteolysisTagTRNASequence represent the
%   sequences of tRNAs required to polymerize each protein and proteolysis tag.
%   These tRNA sequences are computed from the DNA sequence of the mRNA gene /
%   tmRNA gene by the knowledge  base classs' computeTRNASequences method.
%
%   ribosomeStates represents the state/pseudostate (actively translating, free,
%   non-existent) of each ribosome. boundMRNAs indicates the mRNA
%   species each actively translating ribosome is bound to. The non-existent
%   pseudostate is used to keep track of elements in ribosomeStates
%   and boundMRNAs which have been allocated.
%   nascentMonomerLengths and proteolysisTagLengths represeent
%   the position of actively  translating ribosomes on mRNAs, and if stalled on
%   a tmRNA.
%
%   Initialization
%   ===============
%   All protein monomers are initialized to the mature state, and in their
%   correct localization. This is achieved by the simulation class
%   initializeState method.
%
%   Ribosomes are initialized to their steady state:
%   - Each ribosome is randomly assigned (without replacement) to mRNA species,
%     weighted by the current expression of the mRNAs,.
%   - Each ribosome is randomly assigned to positions within the assigned mRNA
%     with uniform probably.
%   - No ribosomes are initialized to the stalled state, since the expected
%     occupancy is negligible. No tmRNAs are initialized to the bound state.
%
%   tRNAs and tmRNAs are initialized to the aminoacyated state by
%   simulation.initializeState.
%
%   Simulation
%   ===============
%   Evolves the states of ribosomes among two states:
%   - actively translating
%   - free
%
%   Transition to from the free state to the actively translating state is
%   allowed when
%   - ribosome binding factor A is present
%   - initiation factors (IF-1, IF-2, and IF-3) are present
%   - one unit of energy (GTP) is available.
%   At this time the ribosome binds an mRNA randomly with uniform probability.
%
%   At the next iteration the initiation factors are released when the first
%   amino acid, f-methionine binds and elongation begins, assuming that
%   elongation factors (EF-tu, TS, and G) and energy (GTP) is available.
%
%   ASSUMPTION made here is that one of each is sufficient for each ribosome,
%   but need a separate set for each ribosome. So each ribosome can only
%   translate if a full set exists for it. Also, for cases in which
%   translation of a peptide finishes partway through the time step, the EFs
%   are free, but haven't added the complexity of another ribosome being able
%   to take up EFs partway through the timestep.
%
%   Finally, once all amino acids of a protein have been translated,
%   termination occurs if
%   - at least one terminator (RF-1) available
%   - at least one recycling factor available
%   - at least one elongation factor G available
%   - at least one trigger factor available
%   - energy is available.
%   The ribosome will be available to bind mRNA at the following iteration.
%
%   Ribosomes are created in the free state.
%
%   Algorithm
%   +++++++++++++++
%   1. Up to limit of elongation factors, randomly select actively translating
%      ribosomes to elongate.
%   2. Up to limit of initiation factors and energy, randomly select ribosomal
%      particles to initiatate. These ribosome will be able to start elongating
%      at the next time step. Randomly select mRNA species for each initiating
%      mRNA to bind to, weighted by the counts of each mRNA species. Update
%      ribosomeStates, boundMRNAs. Update substrates.
%   3. Allocate available amino acids and energy among actively translating
%      ribosomes. Update ribosome nascentMonomerLengths, or
%      proteolysisTagLengths for stalled ribosomes with the number of
%      polymerized bases. Update substrates.
%   4. If ribosome has reached end of (t)mRNA and termination factors
%      available, increase count of protein and dissolve (t)mRNA-ribosome
%      complex. Update substrates.
%   5. If ribosome hasn't advanced, then with a small probability transition
%      ribosome to stalled state. Expel mRNA and replace with tmRNA. Update
%      ribosomeStates. Update substrates.
%
%   References
%   ===============
%   1. Petry S, Weixlbaumer A, Ramakrishnan V (2008). The termination of
%      translation. Curr Opin Struct Biol. 18(1):70-7. [PUB_0226].
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanfod.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010

%TODO: Bind specific mRNAs rather than kinds of mRNAs.
%
%TODO: Prevent ribosomes from colliding/overlapping on an mRNA. Currently this
%      is only somewhat ensured by allowing only one ribosome to bind each mRNA
%      per time step.
%
%TODO: Associate elongation factors to mRNA for just 1 amino acid, not
%      elongationRate amino acids. With this change, all mRNAs will
%      translate at each second but with a slower rate. Related to the
%      assumption above. That is we should translate total of min(energy,
%      min(EFs)*rate*this.simulationTime) amino acids, where rate is a
%      rate of EF recycling per second, that is how many times an
%      individual EF can be used per second.
classdef Translation < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'ribosomeElongationRate'
            'tmRNABindingProbability'
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'mRNAs'
            'freeTRNAs'
            'aminoacylatedTRNAs'
            'aminoacylatedTMRNA'
            'monomers'
            };
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE';
            'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'FMET';
            'GTP'; 'GDP'; 'PI'; 'H2O'; 'H'};
        substrateIndexs_aminoAcids  = (1:20)'; %indices within substrates of amino acids
        substrateIndexs_fMet        = 21;      %index within substrates of formylmethionine
        substrateIndexs_gtp         = 22;      %index within substrates of GTP
        substrateIndexs_gdp         = 23;      %index within substrates of GDP
        substrateIndexs_phosphate   = 24;      %index within substrates of phosphate
        substrateIndexs_water       = 25;      %index within substrates of water
        substrateIndexs_hydrogen    = 26;      %index within substrates of hydrogen

        enzymeWholeCellModelIDs = { %whole cell model ids of enzymes
            'MG_173_MONOMER';   %translation initiation factor IF-1
            'MG_142_MONOMER';   %translation initiation factor IF-2
            'MG_196_MONOMER';   %translation initiation factor IF-3
            'MG_089_DIMER';     %translation elongation factor G
            'MG_026_MONOMER';   %translation elongation factor P
            'MG_451_DIMER';     %translation elongation factor Tu
            'MG_433_DIMER';     %translation elongation factor Ts
            'MG_258_MONOMER';   %peptide chain release factor 1
            'MG_435_MONOMER';   %ribosome recycling factor
            'RIBOSOME_30S';     %30S ribosomal subunit
            'RIBOSOME_30S_IF3'; %30S ribosomal subunit - translation initiation factor IF-3 complex
            'RIBOSOME_50S';     %50S ribosomal subunit
            'RIBOSOME_70S';     %70S ribosome
            'MG_0004';          %transfer-messenger RNA, tmRNA, MCS6, 10Sa RNA
            'MG_059_MONOMER';   %SsrA-binding protein
            'MG_083_MONOMER'};  %peptidyl-tRNA hydrolase
        enzymeIndexs_translationFactors          = (1:9)';    %indices within enzymes of translation factors
        enzymeIndexs_initiationFactors           = (1:3)';    %indices within enzymes of initiation factors
        enzymeIndexs_initiationFactor1           = 1;         %index within enzymes of initiation factor 1
        enzymeIndexs_initiationFactor2           = 2;         %index within enzymes of initiation factor 2
        enzymeIndexs_initiationFactor3           = 3;         %index within enzymes of initiation factor 3
        enzymeIndexs_elongationFactors           = (4:7)';    %indices within enzymes of elongation factors
        enzymeIndexs_elongationGFactor           = 4;         %index within enzymes of elongation factor G
        enzymeIndexs_terminationFactor           = 8;         %index within enzymes of termination factor
        enzymeIndexs_recyclingFactor             = 9;         %index within enzymes of recylcing factor
        enzymeIndexs_ribosome30S                 = 10;        %index within enzymes of 30S ribosomal particle
        enzymeIndexs_ribosome30SIF3              = 11;        %index within enzymes of 30S ribosomal particle-initiation factor 3 complex
        enzymeIndexs_ribosome50S                 = 12;        %index within enzymes of 50S ribosomal particle
        enzymeIndexs_ribosome70S                 = 13;        %index within enzymes of 70S ribosome
        enzymeIndexs_abortiveTranslationEnzymes  = (14:16)';  %indices within enzymes of abortive translation enzymes
        enzymeIndexs_tmRNA                       = 14;        %index within enzymes of tmRNA
        enzymeIndexs_tmRNABindingProtein         = 15;        %index within enzymes of tmRNA binding protein
        enzymeIndexs_peptidylTRNAHydrolase       = 16;        %index within enzymes of peptidyl-tRNA hydrolase

        freeTRNAWholeCellModelIDs          %IDs of free (non-aminoacylated) tRNAs
        aminoacylatedTRNAWholeCellModelIDs %IDs of aminoacylated tRNAs
        aminoacylatedTMRNAWholeCellModelID %ID of aminoacylated tmRNA
        
        tmRNAGlobalIndex %index of tmRNA within simulation.matureRNAIndexs
    end

    %fixed biological constants
    properties
        ribosomeElongationRate             %ribosomse elongation rate (amino acids per second per ribosome) [PUB_0564]
        tmRNABindingProbability            %probability that tmRNA is recruited to non-elongating, active ribosome        
    end

    %global state (stored locally for convenience)
    properties
        monomers                       %numbers of protein monomers
        mRNAs                          %numbers of mRNAs
        freeTRNAs                      %numbers of free (non-aminoacylated) tRNAs
        freeTMRNA                      %numbers of free tmRNA
        aminoacylatedTRNAs             %numbers of aminoacylated tRNAs
        aminoacylatedTMRNA             %numbers of aminoacylated tmRNA
        boundTMRNA                     %numbers of bound tmRNA
    end
  
    %global state (referenced locally for convenience)
    properties
        polypeptide                    %New polypeptides state class
        ribosome                       %Ribosome state class
    end
    
    %constructor
    methods
        function this = Translation(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            
            this.polypeptide = simulation.state('Polypeptide');
            this.ribosome = simulation.state('Ribosome');
            this.states = [this.states; {this.polypeptide; this.ribosome}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)            
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});           
            
            %tmRNA
            this.tmRNAGlobalIndex = this.enzymeRNAGlobalIndexs(this.enzymeRNALocalIndexs==this.enzymeIndexs_tmRNA);
            
            %RNA
            this.freeTRNAWholeCellModelIDs          = this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureTRNAIndexs));
            this.aminoacylatedTRNAWholeCellModelIDs = this.rna.wholeCellModelIDs(this.rna.aminoacylatedIndexs(this.rna.matureTRNAIndexs));
            this.aminoacylatedTMRNAWholeCellModelID = this.rna.wholeCellModelIDs(this.rna.aminoacylatedIndexs(this.tmRNAGlobalIndex));
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();
            
            r = this.rna;
            m = this.monomer;
            
            if size(r.counts, 3) == 1
                this.mRNAs = ...
                    r.matureRNAGeneComposition(this.gene.mRNAIndexs, r.matureMRNAIndexs) * ...
                    r.counts(r.matureIndexs(r.matureMRNAIndexs), this.compartment.cytosolIndexs);
            else
                this.mRNAs = multiprod(...
                    r.matureRNAGeneComposition(this.gene.mRNAIndexs, r.matureMRNAIndexs), ...
                    r.counts(r.matureIndexs(r.matureMRNAIndexs), this.compartment.cytosolIndexs, :), [1 2], 1);
            end
            this.monomers           = m.counts(m.nascentIndexs, this.compartment.cytosolIndexs, :);
            this.freeTRNAs          = r.counts(r.matureIndexs(this.rna.matureTRNAIndexs), this.compartment.cytosolIndexs, :);
            this.freeTMRNA          = r.counts(r.matureIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :);
            this.aminoacylatedTRNAs = r.counts(r.aminoacylatedIndexs(this.rna.matureTRNAIndexs), this.compartment.cytosolIndexs, :);
            this.aminoacylatedTMRNA = r.counts(r.aminoacylatedIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :);
            this.boundTMRNA         = r.counts(r.boundIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :);
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();
            
            r = this.rna;
            m = this.monomer;
            
            m.counts(m.nascentIndexs, this.compartment.cytosolIndexs, :)  = this.monomers;
            r.counts(r.matureIndexs(this.rna.matureTRNAIndexs), this.compartment.cytosolIndexs, :) = this.freeTRNAs;
            r.counts(r.matureIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :) = this.freeTMRNA;
            r.counts(r.aminoacylatedIndexs(this.rna.matureTRNAIndexs), this.compartment.cytosolIndexs, :) = this.aminoacylatedTRNAs;
            r.counts(r.aminoacylatedIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :) = this.aminoacylatedTMRNA;
            r.counts(r.boundIndexs(this.tmRNAGlobalIndex), this.compartment.cytosolIndexs, :) = this.boundTMRNA;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.monomers           = zeros(length(this.polypeptide.monomerLengths), 1, numTimePoints);
            this.mRNAs              = zeros(length(this.gene.mRNAIndexs), 1, numTimePoints);
            this.freeTRNAs          = zeros(length(this.rna.matureTRNAIndexs), 1, numTimePoints);
            this.freeTMRNA          = zeros(1, 1, numTimePoints);
            this.aminoacylatedTRNAs = zeros(length(this.rna.matureTRNAIndexs), 1, numTimePoints);
            this.aminoacylatedTMRNA = zeros(1, 1, numTimePoints);
            this.boundTMRNA         = zeros(1, 1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        %        
        %Contributions, by process, to protein synthesis:
        %                            AAs --(translation)----> nascent monomer + water
        %          nascent monomer + H2O --(processing)-----> processed protein + AAs
        %processed monomer + metabolites --(folding)--------> folded monomer
        %   folded monomer + metabolites --(modification)---> mature monomer + metabolites
        %              modified monomers --(complexation)---> nascent complexs
        % nascent complexs + metabolites --(folding)--------> mature complexes
        %          monomer/complex + H2O --(decay)----------> AAs + modified AAs
        %        phosphorylated AA + H2O --(catabolism)-----> AA + Pi
        %            lipoyl lysine + H2O --(catabolism)-----> lysine + lipoate
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            
            %breaking tRNA-peptidyl bond 
            bmProd(this.substrateIndexs_water) = ...
                bmProd(this.substrateIndexs_water) + ...
                this.polypeptide.monomerLengths' * states.monomerProductions;
            byProd(this.substrateIndexs_hydrogen) = ...
                byProd(this.substrateIndexs_hydrogen) + ...
                this.polypeptide.monomerLengths' * states.monomerProductions;
            
            %peptide bond formation
            byProd(this.substrateIndexs_water) = ...
                byProd(this.substrateIndexs_water) + ...
                max(0, this.polypeptide.monomerLengths - 1)' * states.monomerProductions;
            
            %energy
            translationEnergy = ...
                2 * this.polypeptide.monomerLengths' * states.monomerProductions + ...  %elongation
                3 * sum(states.monomerProductions);                                     %initiation, termination
            bmProd(this.substrateIndexs_gtp)       = bmProd(this.substrateIndexs_gtp)       + translationEnergy;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + translationEnergy;
            byProd(this.substrateIndexs_gdp)       = byProd(this.substrateIndexs_gdp)       + translationEnergy;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + translationEnergy;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + translationEnergy;
            
            %% enzymes

            %Ribosome
            fractionInitiating = this.ribosomeElongationRate *  ...
                states.monomerProductions0' * (1 ./ this.polypeptide.monomerLengths) / ...
                sum(states.monomerProductions0); %fraction of ribosomes that are initiating
            aaPolymerized0 = (states.monomerProductions0' * this.polypeptide.monomerLengths);
            minEnzExp(this.enzymeIndexs_ribosome30S) = 1.65 * aaPolymerized0 / this.ribosomeElongationRate / (1 - fractionInitiating);
            minEnzExp(this.enzymeIndexs_ribosome50S) = 1.65 * aaPolymerized0 / this.ribosomeElongationRate / (1 - fractionInitiating);

            %initiation, elongation, termination factors
            minEnzExp(this.enzymeIndexs_initiationFactors)     = 3 * minEnzExp(this.enzymeIndexs_ribosome30S) *      fractionInitiating;
            minEnzExp(this.enzymeIndexs_elongationFactors)     = 3 * minEnzExp(this.enzymeIndexs_ribosome30S) * (1 - fractionInitiating);
            minEnzExp(this.enzymeIndexs_terminationFactor)     = 3;
            minEnzExp(this.enzymeIndexs_recyclingFactor)       = 3;
            minEnzExp(this.enzymeIndexs_elongationGFactor)     = 3 + minEnzExp(this.enzymeIndexs_elongationGFactor);

            %abortive translation enzymes
            minEnzExp(this.enzymeIndexs_tmRNA)                 = 3;
            minEnzExp(this.enzymeIndexs_tmRNABindingProtein)   = 3;
            minEnzExp(this.enzymeIndexs_peptidylTRNAHydrolase) = 3;
        end
        
        %Initialize ribosome, bound mRNAs, translation factor states
        %Assumptions:
        %- Metabolic cost of the transcription of these amino acids is negligible
        %- Probability that so many ribosomes bind a given mRNA that they can't
        %  pack along the mRNA without steric interactions is negligible. Thus we
        %  don't try to handle this case separately.
        function initializeState(this)
            %references to states
            rib = this.ribosome;
            pol = this.polypeptide;
            
            %Check ribosomes and polypeptides are synchronized
            if ...
                    ~isequal(rib.boundMRNAs, pol.boundMRNAs) || ...
                    ~isequal(rib.mRNAPositions, pol.nascentMonomerLengths) || ...
                    ~isequal(rib.tmRNAPositions, pol.proteolysisTagLengths)
                throw(MException('Translation:error', 'Ribosomes and polypeptides must be synchronized'));
            end
            
            %enzymes
            this.enzymes = this.enzymes + this.boundEnzymes;
            this.boundEnzymes(:) = 0;
            this.enzymes(this.enzymeIndexs_ribosome30S) = this.enzymes(this.enzymeIndexs_ribosome30S) + this.enzymes(this.enzymeIndexs_ribosome70S);
            this.enzymes(this.enzymeIndexs_ribosome50S) = this.enzymes(this.enzymeIndexs_ribosome50S) + this.enzymes(this.enzymeIndexs_ribosome70S);
            this.enzymes(this.enzymeIndexs_ribosome70S) = 0;
            this.enzymes(this.enzymeIndexs_ribosome30S) = this.enzymes(this.enzymeIndexs_ribosome30S) + this.enzymes(this.enzymeIndexs_ribosome30SIF3);
            this.enzymes(this.enzymeIndexs_initiationFactor3) = this.enzymes(this.enzymeIndexs_initiationFactor3) + this.enzymes(this.enzymeIndexs_ribosome30SIF3);
            this.enzymes(this.enzymeIndexs_ribosome30SIF3) = 0;
            
            freeTranslationFactors  = this.enzymes(this.enzymeIndexs_translationFactors);
            boundTranslationFactors = this.boundEnzymes(this.enzymeIndexs_translationFactors);
            ribosome30S             = this.enzymes(this.enzymeIndexs_ribosome30S);
            ribosome30SIF3          = this.enzymes(this.enzymeIndexs_ribosome30SIF3);
            ribosome50S             = this.enzymes(this.enzymeIndexs_ribosome50S);
            ribosome70S             = this.enzymes(this.enzymeIndexs_ribosome70S);

            boundRibosome70S = min([
                ribosome30S;
                ribosome50S;
                freeTranslationFactors(this.enzymeIndexs_elongationFactors);
                sum(this.mRNAs)]);
            ribosome30S = ribosome30S - boundRibosome70S;
            ribosome50S = ribosome50S - boundRibosome70S;
            freeTranslationFactors(this.enzymeIndexs_elongationFactors)  = freeTranslationFactors(this.enzymeIndexs_elongationFactors)  - boundRibosome70S;
            boundTranslationFactors(this.enzymeIndexs_elongationFactors) = boundTranslationFactors(this.enzymeIndexs_elongationFactors) + boundRibosome70S;

            %allocate
            rib.states                = repmat(rib.notExistValue, [2*boundRibosome70S, 1, 1]);
            pol.boundMRNAs            = zeros(2 * boundRibosome70S, 1, 1);
            pol.nascentMonomerLengths = zeros(2 * boundRibosome70S, 1, 1);
            pol.proteolysisTagLengths = zeros(2 * boundRibosome70S, 1, 1);

            %bind random mRNA
            for i = 1:boundRibosome70S
                pol.boundMRNAs(i) = this.randStream.randsample(numel(this.mRNAs), 1, true, this.mRNAs .* pol.monomerLengths);                
            end

            %actively bound
            %- partition sequence among RNA polymerases bound to transcription unit
            %- randomly select position of RNA polymerases within partition
            for i = 1:length(this.mRNAs)
                boundRibosomes = find(pol.boundMRNAs == i);
                partitions = round(...
                    (0:length(boundRibosomes)) * ...
                    (pol.monomerLengths(i) + this.ribosomeElongationRate - 1) / length(boundRibosomes) -...
                    1/2*this.ribosomeElongationRate);

                for j = 1:length(boundRibosomes)
                    minState = partitions(j)   + 1/2 * this.ribosomeElongationRate;
                    maxState = partitions(j+1) - 1/2 * this.ribosomeElongationRate;
                    rib.states(boundRibosomes(j)) = rib.activeValue;
                    if maxState < minState
                        pol.nascentMonomerLengths(boundRibosomes(j)) = round(max(1, min(pol.monomerLengths(i), mean(partitions(j:j+1)))));
                    else
                        pol.nascentMonomerLengths(boundRibosomes(j)) = round(minState + this.randStream.rand * (maxState - minState));
                    end
                end
            end
            
            %synchronize ribosome and polypeptide states
            rib.boundMRNAs = pol.boundMRNAs;
            rib.mRNAPositions = pol.nascentMonomerLengths;
            rib.tmRNAPositions = pol.proteolysisTagLengths;

            %store enzymes
            this.enzymes(this.enzymeIndexs_translationFactors)      = freeTranslationFactors;
            this.boundEnzymes(this.enzymeIndexs_translationFactors) = boundTranslationFactors;
            this.enzymes(this.enzymeIndexs_ribosome30S)             = ribosome30S;
            this.enzymes(this.enzymeIndexs_ribosome30SIF3)          = ribosome30SIF3;
            this.enzymes(this.enzymeIndexs_ribosome50S)             = ribosome50S;
            this.enzymes(this.enzymeIndexs_ribosome70S)             = ribosome70S;
            this.boundEnzymes(this.enzymeIndexs_ribosome70S)        = boundRibosome70S;
            
            %decrement counts of protein
            this.copyToState();
            
            g = this.gene;
            rna = this.rna;
            mon = this.monomer;
            cpx = this.complex;
            pcComp = sum(cpx.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3);
            
            initMonCnts = mon.counts;
            initCpxCnts = cpx.counts;
            matureProteinWt = max(0, sum(mon.dryWeight + cpx.dryWeight) - pol.dryWeight);
            while sum(mon.dryWeight + cpx.dryWeight) > matureProteinWt
                idx = this.randStream.randsample(numel(mon.counts) + numel(cpx.counts), 1, false, [mon.counts(:); cpx.counts(:)]);
                if idx <= numel(mon.counts)
                    mon.counts(idx) = mon.counts(idx) - 1;
                else
                    idx = idx - numel(mon.counts);
                    [i, ~] = ind2sub(size(cpx.counts), idx);
                    i = mod(i - 1, numel(cpx.matureIndexs)) + 1;
                    if any(pcComp(:, i))
                        rna.counts(rna.matureIndexs(setdiff(1:end, rna.matureMRNAIndexs))) = ...
                            + rna.counts(rna.matureIndexs(setdiff(1:end, rna.matureMRNAIndexs))) ...
                            + pcComp(:, i);
                    end
                    cpx.counts(idx) = cpx.counts(idx) - 1;
                end
            end
            mon.counts = ...
                + mon.counts ...
                + mon.updateExternalState(-(initMonCnts - mon.counts), true);
            cpx.counts = ...
                + cpx.counts ...
                + cpx.updateExternalState(-(initCpxCnts - cpx.counts), true);
            
            this.copyFromState();
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            initiationFactor3 = this.enzymes(this.enzymeIndexs_initiationFactor3);
            ribosome30S       = this.enzymes(this.enzymeIndexs_ribosome30S);
            ribosome30SIF3    = this.enzymes(this.enzymeIndexs_ribosome30SIF3);
            ribosome50S       = this.enzymes(this.enzymeIndexs_ribosome50S);
            nRibosomes = ...
                this.enzymes(this.enzymeIndexs_ribosome70S) + ...
                this.boundEnzymes(this.enzymeIndexs_ribosome70S) + ...
                min(ribosome30SIF3 + min(ribosome30S, initiationFactor3), ribosome50S);
                
            nInitiations = nRibosomes;
            nElongations = min([
                sum(this.aminoacylatedTRNAs) + this.aminoacylatedTMRNA
                nRibosomes * this.ribosomeElongationRate]);
            nTerminations = nRibosomes;

            %energy for initiation, elongation, termination
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_gtp)   = nInitiations + 2 * nElongations + 3 * nTerminations;
            result(this.substrateIndexs_water) = nInitiations + 2 * nElongations + 3 * nTerminations;
            
            %water to cleave peptidyl-tRNA bound
            result(this.substrateIndexs_water) = - nInitiations - nElongations;
            
            %water from peptide bond formation
            result(this.substrateIndexs_water) = + nElongations;
        end
        
        %simulation
        function evolveState(this)
            %% references to states
            rib = this.ribosome;
            pol = this.polypeptide;
            
            %% Check ribosomes and polypeptides are synchronized
            if ...
                    ~isequal(rib.boundMRNAs, pol.boundMRNAs) || ...
                    ~isequal(rib.mRNAPositions, pol.nascentMonomerLengths) || ...
                    ~isequal(rib.tmRNAPositions, pol.proteolysisTagLengths)
                throw(MException('Translation:error', 'Ribosomes and polypeptides must be synchronized'));
            end
            
            %% allocate memory and initialize variables
            
            % Numbers of enzymes
            freeTranslationFactors     = this.enzymes(this.enzymeIndexs_translationFactors);
            boundTranslationFactors    = this.boundEnzymes(this.enzymeIndexs_translationFactors);
            ribosome30S                = this.enzymes(this.enzymeIndexs_ribosome30S);
            ribosome30SIF3             = this.enzymes(this.enzymeIndexs_ribosome30SIF3);
            ribosome50S                = this.enzymes(this.enzymeIndexs_ribosome50S);
            ribosome70S                = this.enzymes(this.enzymeIndexs_ribosome70S);
            boundRibosome70S           = this.boundEnzymes(this.enzymeIndexs_ribosome70S);
            freeTMRNA                  = this.freeTMRNA; %#ok<*PROP>
            freeTMRNABindingProtein    = this.enzymes(this.enzymeIndexs_tmRNABindingProtein);
            freePeptidylTRNAHydrolase  = this.enzymes(this.enzymeIndexs_peptidylTRNAHydrolase);
            boundTMRNA                 = this.boundTMRNA;
            aminoacylatedTmRNA         = this.aminoacylatedTMRNA;
            
            % form ribosomes
            new_ribosome30SIF3 = min(ribosome30S, freeTranslationFactors(this.enzymeIndexs_initiationFactor3));
            ribosome30SIF3 = ribosome30SIF3 + new_ribosome30SIF3;
            ribosome30S    = ribosome30S - new_ribosome30SIF3;
            freeTranslationFactors(this.enzymeIndexs_initiationFactor3) = freeTranslationFactors(this.enzymeIndexs_initiationFactor3) - new_ribosome30SIF3;
            
            %allocate space to store states of new ribosomes
            if boundRibosome70S + ribosome30SIF3 > numel(rib.states)
                rib.states = [
                    rib.states;
                    repmat(rib.notExistValue, [boundRibosome70S + ribosome30SIF3 1])];
                pol.boundMRNAs = [
                    pol.boundMRNAs;
                    zeros([boundRibosome70S + ribosome30SIF3 1])];
                pol.nascentMonomerLengths = [
                    pol.nascentMonomerLengths;
                    zeros([boundRibosome70S + ribosome30SIF3 1])];
                pol.proteolysisTagLengths = [
                    pol.proteolysisTagLengths;
                    zeros([boundRibosome70S + ribosome30SIF3 1])];
            end

            %ribosomes
            elongationRate = this.ribosomeElongationRate;               %ribosome elongation rate (aa/s)
            actRibs = find(...                                          %indices of active ribosomes
                rib.states == rib.activeValue | ...
                rib.states == rib.stalledValue);
            notExistRibs = find(rib.states == rib.notExistValue);       %indices of free ribosomes 30SIF3
            elngRibs = false(size(rib.states, 1), 1);                   %indicate of ribosome actively translating
            elngSeqs = zeros(size(rib.states, 1), elongationRate);      %tRNA sequences of mRNAs segments that should be translated
            termRibs = false(size(rib.states, 1), 1);
            termRibs(rib.states == rib.activeValue, 1) = ...
                pol.monomerLengths(pol.boundMRNAs(rib.states == rib.activeValue, 1), 1) == ...
                pol.nascentMonomerLengths(rib.states == rib.activeValue, 1);
            termRibs = termRibs | pol.proteolysisTagLengths == pol.proteolysisTagLength;
            
            %mRNAs
            bndProbs = this.mRNAs;                                      %unnormalized probability of binding mRNAs
            
            %tRNAs
            nTRNAs = length(this.freeTRNAs);                            %number of tRNA species
            tRNAs  = this.freeTRNAs + this.aminoacylatedTRNAs;          %total numbers of tRNAs (free and aminoacylated)
            
            %energy
            energyAllocation = max(0, this.substrates(this.substrateIndexs_gtp));
            availableEnergy  = energyAllocation;                        %energy (GTP) allocated to translation
            
            %water
            availableWater = this.substrates(this.substrateIndexs_water);
     
            %% recycle elongation factors
            freeTranslationFactors(this.enzymeIndexs_elongationFactors) = ...
                freeTranslationFactors(this.enzymeIndexs_elongationFactors) + ...
                boundTranslationFactors(this.enzymeIndexs_elongationFactors);
            boundTranslationFactors(this.enzymeIndexs_elongationFactors) = 0;

            %% select bound mRNAs to elongate
            %elongate if all elongation factors available
            nElongatingRibosomes = min([
                numel(actRibs)
                boundRibosome70S
                floor(min(availableEnergy, availableWater) / 2)]);
            
            randomOrder = this.randStream.randperm(boundRibosome70S);
            for idx = 1:nElongatingRibosomes
                i = actRibs(randomOrder(idx));

                %Use elongation factors to elongate
                if ~termRibs(i)
                    if all(freeTranslationFactors(this.enzymeIndexs_elongationFactors))
                        freeTranslationFactors(this.enzymeIndexs_elongationFactors)  = freeTranslationFactors(this.enzymeIndexs_elongationFactors)  - 1;
                        boundTranslationFactors(this.enzymeIndexs_elongationFactors) = boundTranslationFactors(this.enzymeIndexs_elongationFactors) + 1;
                    else
                        continue;
                    end
                end

                %Extract portion of sequence to be translated
                switch rib.states(i)
                    case rib.activeValue
                        idxs = pol.nascentMonomerLengths(i) + 1:min([
                            pol.monomerLengths(pol.boundMRNAs(i));
                            pol.nascentMonomerLengths(i) + elongationRate]);
                        elngSeqs(i, :) = padarray(pol.monomerTRNASequences{pol.boundMRNAs(i)}(idxs)', [0 elongationRate-length(idxs)], 0, 'post');
                        elngRibs(i) = true;
                    case rib.stalledValue
                        idxs = pol.proteolysisTagLengths(i) + 1:min([
                            pol.proteolysisTagLength;
                            pol.proteolysisTagLengths(i) + elongationRate]);
                        elngSeqs(i, :) = padarray(pol.proteolysisTagTRNASequence(idxs)', [0 elongationRate-length(idxs)], 0, 'post');
                        elngRibs(i) = true;
                end
            end
            nElongatingRibosomes = sum(elngRibs);

            %% make/update free ribosomes, initiate translation
            %Ribosome 70S assembly proceeds as:
            %30S-IF3 + IF2 + IF1 --> 30S-IF1-IF2-IF3-mRNA --> 30S-IF1-IF2-IF3-mRNA-tRNA
            %--> 30S-50S-IF1-IF2-mRNA-tRNA --> 30S-50S-mRNA-tRNA
            %All of these steps occur quite rapidly, so in our model if enough factors,
            %subunits, RNAs exist, then a 70S ribosome can be formed and
            %translation initiated.
            %(tRNA is bound in elongation, not initiation.)

            %Initiate if
            %- 30S-IF3 and 50S available
            %- initiation factors IF1 IF2 available
            %- energy available
            %- ribosome binding factor available
            %- mRNAs available
            %- (already verified by knowledge base computeTRNASequences
            %   method that first codon is formylmethionine)
            nInitiatingRibosomes = min([
                ribosome30SIF3;
                ribosome50S;
                freeTranslationFactors(this.enzymeIndexs_initiationFactor1);
                freeTranslationFactors(this.enzymeIndexs_initiationFactor2);
                max(0, min(availableEnergy, availableWater) - 2 * nElongatingRibosomes);
                sum(this.mRNAs)]);
            if nInitiatingRibosomes > 0
                ribosome30SIF3   = ribosome30SIF3   - nInitiatingRibosomes;
                ribosome50S      = ribosome50S      - nInitiatingRibosomes;
                boundRibosome70S = boundRibosome70S + nInitiatingRibosomes;
                freeTranslationFactors(this.enzymeIndexs_initiationFactor3) = ...
                    freeTranslationFactors(this.enzymeIndexs_initiationFactor3) + ...
                    nInitiatingRibosomes;
                availableEnergy = availableEnergy - nInitiatingRibosomes;
                availableWater = availableWater - nInitiatingRibosomes;
                for idx = 1:nInitiatingRibosomes
                    i = notExistRibs(idx);
                    
                    %randomly choose mRNA to bind
                    pol.boundMRNAs(i) = this.randStream.randsample(numel(bndProbs), 1, true, bndProbs);
                    
                    %update state of ribosome
                    rib.states(i) = rib.activeValue;
                    
                    %update binding probabilities to prevent additional ribosomes
                    %from initiating translation with same mRNA at same time
                    bndProbs(pol.boundMRNAs(i)) = max(0, bndProbs(pol.boundMRNAs(i)) - 1);
                end
            end

            %% Translate active sequences
            %compute progress of translation
            nInit = sum((pol.nascentMonomerLengths(elngRibs) + pol.proteolysisTagLengths(elngRibs)) == 0);
            availableEnergyWater = floor(min(availableEnergy, availableWater - nInit));
            [translationProgress, this.aminoacylatedTRNAs, translationCost, availableEnergyWater2] = ...
                edu.stanford.covert.cell.sim.util.polymerize(...
                elngSeqs(elngRibs, :), ...
                this.aminoacylatedTRNAs, 1:nTRNAs, 0, availableEnergyWater, 2, this.randStream);            
            if any(translationProgress)
                %update water, gtp
                availableWater = ...
                    + availableWater ...
                    - (availableEnergyWater - availableEnergyWater2) ...
                    - sum((pol.nascentMonomerLengths(elngRibs) + pol.proteolysisTagLengths(elngRibs)) == 0 & translationProgress > 0);
                availableEnergy = availableEnergy - (availableEnergyWater - availableEnergyWater2);
                
                %advance positions of ribosomes
                pol.nascentMonomerLengths(elngRibs) = ...
                    pol.nascentMonomerLengths(elngRibs) + ...
                    translationProgress .* (rib.states(elngRibs) == rib.activeValue);
                if any(rib.states(elngRibs) == rib.stalledValue)
                    pol.proteolysisTagLengths(elngRibs) = ...
                        pol.proteolysisTagLengths(elngRibs) + ...
                        translationProgress .* (rib.states(elngRibs) == rib.stalledValue);
                end
            end

            %% Terminate, remaining in elongating state, or convert to stalled state
            %If these conditions are true
            %(1) ribosome has translated entire sequence,
            %(2) termination factor is available (pfrA)--dissociates
            %    peptide from complex needs 1GTP,
            %(3) ribosome recycling factor (frr) and EF-G (fusA) available
            %    (use 1GTP) to split off the 50S subunit,
            %(4) IF-3 (infC) available to split off mRNA and tRNA and
            %(5) energy is available (total 2GTPs)
            %
            %Then terminate translation:
            %- use energy to terminate
            %- release ribosomes which have completed translation
            %- increment expression of the corresponding protein
            %- IF3 remains bound to 30S

            randomOrder = this.randStream.randperm(length(elngRibs));
            tmpProg = zeros(size(elngRibs));
            tmpProg(elngRibs) = translationProgress;
            for idx = 1:length(elngRibs)
                i = randomOrder(idx);
                if elngRibs(i) && (tmpProg(i) || pol.nascentMonomerLengths(i) >= pol.monomerLengths(pol.boundMRNAs(i)))
                    if (pol.nascentMonomerLengths(i) >= pol.monomerLengths(pol.boundMRNAs(i)) || ...
                            pol.proteolysisTagLengths(i) >= pol.proteolysisTagLength) && ...
                            freeTranslationFactors(this.enzymeIndexs_terminationFactor) && ...
                            freeTranslationFactors(this.enzymeIndexs_recyclingFactor) && ...
                            (freeTranslationFactors(this.enzymeIndexs_elongationGFactor) || tmpProg(i)) && ...
                            freeTranslationFactors(this.enzymeIndexs_initiationFactor3) && ...
                            availableEnergy >= 2 && ...
                            availableWater >= 2
                        
                        %use energy, (tRNAs have already been released in elongation)
                        availableEnergy = availableEnergy - 2;
                        availableWater = availableWater - 2;
                        
                        switch rib.states(i)
                            case rib.activeValue
                                %update protein count
                                this.monomers(pol.boundMRNAs(i)) = this.monomers(pol.boundMRNAs(i)) + 1;
                            case rib.stalledValue
                                %update stalled polypeptide count
                                pol.abortedPolypeptides = [
                                    pol.abortedPolypeptides;
                                    pol.boundMRNAs(i) pol.nascentMonomerLengths(i) pol.proteolysisTagLength];
                                
                                %release tmRNA
                                boundTMRNA = boundTMRNA - 1;
                                freeTMRNA  = freeTMRNA  + 1;
                        end
                        
                        %release ribosome
                        rib.states(i) = rib.notExistValue;
                        
                        %release mRNA
                        pol.nascentMonomerLengths(i) = 0;
                        pol.proteolysisTagLengths(i) = 0;
                        pol.boundMRNAs(i)            = 0;

                        %release ribosome, split into 30S and 50S
                        boundRibosome70S = boundRibosome70S - 1;
                        ribosome30S      = ribosome30S      + 1;
                        ribosome50S      = ribosome50S      + 1;
                    end
                elseif rib.states(i) == rib.activeValue && ...
                        pol.nascentMonomerLengths(i) > 0 && ...
                        aminoacylatedTmRNA > 0 && ...
                        freeTMRNABindingProtein > 0  && ...
                        freePeptidylTRNAHydrolase > 0 && ...
                        this.randStream.rand < this.tmRNABindingProbability

                    %set ribosome to stalled state
                    rib.states(i) = rib.stalledValue;

                    %increment bound tmRNAs
                    boundTMRNA = boundTMRNA + 1;

                    %decrement free aminoacylated tmRNAs
                    aminoacylatedTmRNA = aminoacylatedTmRNA - 1;

                    %increment proteolysis tag length by 1 for alanine
                    %bound to tRNA-like portion of tmRNA
                    pol.proteolysisTagLengths(i) = 1;
                end
            end
            
            %% synchronize ribosome and polypeptide states
            rib.boundMRNAs = pol.boundMRNAs;
            rib.mRNAPositions = pol.nascentMonomerLengths;
            rib.tmRNAPositions = pol.proteolysisTagLengths;

            %% store enzymes
            this.enzymes(this.enzymeIndexs_translationFactors)         = freeTranslationFactors;
            this.boundEnzymes(this.enzymeIndexs_translationFactors)    = boundTranslationFactors;
            this.enzymes(this.enzymeIndexs_ribosome30S)                = ribosome30S;
            this.enzymes(this.enzymeIndexs_ribosome30SIF3)             = ribosome30SIF3;
            this.enzymes(this.enzymeIndexs_ribosome50S)                = ribosome50S;
            this.enzymes(this.enzymeIndexs_ribosome70S)                = ribosome70S;
            this.boundEnzymes(this.enzymeIndexs_ribosome70S)           = boundRibosome70S;
            this.freeTMRNA                                             = freeTMRNA;
            this.enzymes(this.enzymeIndexs_tmRNABindingProtein)        = freeTMRNABindingProtein;
            this.enzymes(this.enzymeIndexs_peptidylTRNAHydrolase)      = freePeptidylTRNAHydrolase;
            this.boundTMRNA                                            = boundTMRNA;
            this.aminoacylatedTMRNA                                    = aminoacylatedTmRNA;

            %% acount for used substrates
            %aminoacylated tRNAs
            this.freeTRNAs = tRNAs - this.aminoacylatedTRNAs;

            %energy (energy-availableEnergy), water
            this.substrates(this.substrateIndexs_gtp)       = this.substrates(this.substrateIndexs_gtp)       - (energyAllocation - availableEnergy);
            this.substrates(this.substrateIndexs_water)     = availableWater;
            this.substrates(this.substrateIndexs_gdp)       = this.substrates(this.substrateIndexs_gdp)       + (energyAllocation - availableEnergy); 
            this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + (energyAllocation - availableEnergy); 
            this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + (energyAllocation - availableEnergy) + sum(translationCost);
        end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.monomers, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.monomer.molecularWeights(this.monomer.nascentIndexs)' * this.monomers           + ...
                    this.rna.molecularWeights(this.rna.matureIndexs(this.rna.matureTRNAIndexs))'          * this.freeTRNAs          + ...
                    this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.rna.matureTRNAIndexs))'   * this.aminoacylatedTRNAs + ...
                    this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.tmRNAGlobalIndex))'       * this.aminoacylatedTMRNA) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.monomer.molecularWeights(this.monomer.nascentIndexs)' * permute(this.monomers,           [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.matureIndexs(this.rna.matureTRNAIndexs))'        * permute(this.freeTRNAs,          [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.rna.matureTRNAIndexs))' * permute(this.aminoacylatedTRNAs, [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.tmRNAGlobalIndex))'     * permute(this.aminoacylatedTMRNA, [1 3 2]),[1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
