%RNAProcessing
%
% @wholeCellModelID Process_RNAProcessing
% @name             RNA Processing
% @description
%   Biology
%   ==================
%   Transcription produces mono- as well as polycistronic RNAs.
%   Following transcription polycistronic rRNA, tRNA, and tmRNA transcripts are
%   cleaved, modified, and aminoacylated. First, 45S rRNA transcripts are
%   cleaved by ribonuclease III (MG_367_DIMER) into 5S, 16S, and 23S rRNA
%   precursors [PUB_0038, PUB_0039]. Next these precursors are cleaved at their
%   3' and 5' ends by ribonuclease J and RsgA [PUB_0038]. Similarly the 5' ends
%   of each tRNA and tmRNA species are cleaved by ribonuclease P [PUB_0039,
%   PUB_0649]. Ribonuclease III also cleaves scRNA at two sites [PUB_0039].
%   Second, 13 enzymes formylate, lysidinate, methylate, pseudouridylate, and
%   thiolate 86 specific rRNA and tRNA bases. Third 19 tRNA synthetases
%   conjugate 19 amino acids to 37 tRNA and tmRNA species, and 2 tRNA
%   transferases complete the aminoacylation of the formylmethionine and
%   glutamine tRNAs.
%
%   rRNA, tRNA, and tmRNA cleavage, modification, and aminoacylation are modeled
%   motivated by mass-action kinetics. First, the maximum rate of maturation of
%   each RNA species is calculated based on (1) the abundance of immature
%   transcripts, (2) the abundance of substrates and enzymes required to cleave
%   or modify each RNA species, and (3) the experimentally measured kinetic rate
%   of each enzyme required to cleave or modify each RNA species. Second, RNA
%   maturations are randomly selected according to the calculated maximum rates.
%   Finally, steps (1) and (2) are repeated until insufficient resources are
%   available to mature additional RNAs.
%
%   This process simulates the cleavage of polycistronic r/s/t RNA transcripts
%   into individual genes:
%   - 30S RNA        -> pre 16/17S RNA, p23S RNA, 9S RNA    RNAseIII (MG_367)
%   - pre 16/17S RNA -> 16S RNA                             5' end:rnjA (MG_139); 3' end:rsgA (MG_110)
%   - p23S RNA       -> 23S RNA                             5' end:rnjA (MG_139), deaD (MG_425)
%   - 9S RNA         ->  5S RNA
%   - 3' end tRNA precursors                                RNAseIII (MG_367)
%   - 5' end tRNA precursors                                RNaseP (MG_0003, MG_465)
%   - scRNA (MG_0001) precusors                             RNAseIII (MG_367)
%   - 5' end tmRNA precursors                               RNaseP (MG_0003, MG_465)
%
%   Reactions
%   ++++++++++++++++++
%   DeaD           phosphorolytic                      3.6.4.13
%   RsgA           phosphorolytic                      3.6.5.2, 3.6.1.-
%   RNAaseIII      hydrolytic       endoribonuclease   3.1.26.3
%   RNAseP         hydrolytic       endoribonuclease   3.1.26.5
%   RNAseJ         hydrolytic       endoribonuclease   3.1.27.2
%
%   Cofactors
%   ++++++++++++++++++
%   DeaD           Mg2+                [PUB_0071]
%   RsgA       (1) Zn2+                [PUB_0096]
%   RNAaseIII      Mg2+                [PUB_0039]
%   RNAseJ     (2) Zn2+                [PUB_0096]
%   RNAseP     (3) Mg2+ or (3) Mn2+    [PUB_0012, PUB_0013, PUB_0043]
%
%   Kinetics:
%   ++++++++++++++++++
%   DeaD       1.48   1/s   B. subtilis  [PUB_0688]
%   RsgA       0.2917 1/s   E. coli      [PUB_0103]
%   RNAseP     0.027  1/s   B. subtilis  [PUB_0013]
%   RNAseP     6      1/s   B. subtilis  [PUB_0690]
%   RNAseIII   7.7    1/s   E. coli      [PUB_0691]
%   RNAseJ     0.37   1/s   B. subtilis  [PUB_0692]
%
%   Energy-Dependence
%   ++++++++++++++++++
%   DeaD       2.37  ATP  1/rxn   B. subtilis  [PUB_0069] (210 ATP/min B. subtilis [PUB_0069])/(1.48 rxn/s B. subtilis  [PUB_0688])
%   RsgA       1     GTP  1/rxn   E. coli      [PUB_0103]
%
%   Knowledge Base
%   ==================
%   The transcription unit organization was predicted by mapping the
%   experimentally determined genome organization of M. pnemoniae by Güell et al
%   [Table S5 "Suboperons", PUB_0418] onto the M. genitalium genome by homology
%   with 4 modifications:
%   - r/s/tRNAs were organized into transcription units according to their
%     "reference operons" (table s4)
%   - All mRNAs Güell Serrano et al did not assign to a suboperon because of
%     insufficient evidence (low quality expression data), were assigned to
%     their own transcription units. None of these mRNAs are located within
%     suboperons.
%   - All 5 mRNA genes which don't have homologs in M. pneumoniae were assigned
%     to their own transcription units. These occur outside other transcription
%     units.
%   - Genes which have undergone rearrangements between M. pnuemoniae and M.
%     genitalium were assigned to their own transcription units
%
%   The promoter (-35 and -10 boxes and TSS) for each transcription unit was set
%   predicted for its 3' gene by Weiner et al [PUB_0411].
%
%   The composition of each transcription unit, and the location of its -35 and
%   -10 boxes and TSS relative to the start coordinate of the 3' gene was stored
%   in the knowledge base, and it used by the TranscriptionUnit and Gene classes
%   to compute the sequences of transcription unit, gene, and intergenic
%   segment.
%
%   The kinetics and energy requirement of each RNA processing enzyme are
%   organized as parameters in the knowledge base. These values were curated
%   several papers [PUB_0069, PUB_0103, PUB_0688, PUB_0690, PUB_0691, PUB_0692].
%
%   Representation
%   ==================
%   substrates and enzymes represent the counts of metabolites and RNA
%   processing enzymes. unprocessedRNAs, processedRNAs, and intergenicRNAs
%   represent counts of RNAs. unprocessedRNAs represents the counts of nascent
%   RNAs produced by the RNA polymerase. processedRNA represents the counts of
%   r/s/tRNA precursors produced from cleavage of nascent transcripts, and mRNAs
%   that remain joined as transcription units. intergenic RNAs represents the
%   counts of RNA segments between r/s/tRNA genes in transcripts that are
%   released during the processing of these transcripts.
%
%   enzymeSpecificRate_* are the kcats of the five RNA processing enzymes. We
%   assume that each enzyme has the same kcat for all reactions it catalyzes.
%   enzymeEnergyCost_* are the amount of ATP required by the phosphorolytic RNA
%   processing enzymes (DeaD and Rsga) per cleavage reaction.
%
%   rna.nascentRNAMatureRNAComposition, rna.intergenicRNAMatrix,
%   reactantByproductMatrix, and catalysisMatrix are adjacency matrices.
%   nascentRNAMatureRNAComposition represents the the processed RNAs that arise
%   from each unprocessed RNA (eg. the r/s/t RNAs that arise from transcripts
%   containing multiple genes). intergenicRNAMatrix represents the intergenic
%   RNA segments that arise from the cleavage of r/s/tRNA transcripts into their
%   individual genes. reactantByproductMatrix represents the amounts of
%   metabolites required to cleave each polycistronic transcript, and the
%   metabolic byproducts of their cleavage. catalysisMatrix represents the
%   amount of enzyme required to cleave each polycistronic transcript.
%
%   Initialization
%   ==================
%   All RNAs are initialized to the mature state. This is accomplished by the
%   simulation class initializeState method.
%
%   Simulation
%   ==================
%   For each RNA class (mRNA, rRNA, sRNA, tRNA) in a random order:
%   1. Determine maximum number of RNAs that can be processed based on available
%      - RNA processing enzymes and their kinetics
%      - free metabolites
%      - RNAs that need processing
%   2. Randomly select among unprocessed RNAs proportional to count
%   3. Update state:
%      a. Decrement unprocessed RNAs
%      b. Increment processed RNAs, intergenic RNAs
%      c. Decrement metabolic reactants
%      d. Increment metabolic byproducts%
%
%   References
%   ==================
%   1. Güell M, van Noort V, Yus E, Chen WH, Leigh-Bell J, Michalodimitrakis K,
%      Yamada T, Arumugam M, Doerks T, Kühner S, Rode M, Suyama M, Schmidt S,
%      Gavin AC, Bork P, Serrano L (2009). Transcriptome complexity in a
%      genome-reduced bacterium. Science. 326(5957): 1268-71. [PUB_0418]
%   2. Weiner J 3rd, Herrmann R, Browning GF (2000). Transcription in Mycoplasma
%      pneumoniae. Nucleic Acids Res. 28: 4488-96. [PUB_0411].
%   3. Cordin O, Banroques J, Tanner NK, Linder P (2006). The DEAD-box protein
%      family of RNA helicases. Gene. 367: 17-37. [PUB_0069]
%   4. Himeno H, Hanawa-Suetsugu K, Kimura T, Takagi K, Sugiyama W, Shirata S,
%      Mikami T, Odagiri F, Osanai Y, Watanabe D, Goto S, Kalachnyuk L, Ushida
%      C, Muto A (2004). A  novel GTPase activated by the small subunit of
%      ribosome. Nucleic Acids Res. 32 (17): 5303-9. [PUB_0103]
%   5. Theissen B, Karow AR, Köhler J, Gubaev A, Klostermeier D (2008).
%      Cooperative binding of ATP and RNA induces a closed conformation in a
%      DEAD box RNA helicase. Proc Natl Acad Sci U S A. 105(2): 548-53.
%      [PUB_0688]
%   6. Xiao S, Scott F, Fierke CA, Engelke DR (2002). Eukaryotic ribonuclease P:
%      a plurality of ribonucleoprotein enzymes. Annu Rev Biochem. 71: 165-89.
%      [PUB_0690]
%   7. Amarasinghe AK, Calin-Jageman I, Harmouch A, Sun W, Nicholson AW (2001).
%      Escherichia coli ribonuclease III: affinity purification of
%      hexahistidine-tagged enzyme and assays for substrate binding and
%      cleavage. Methods Enzymol. 342: 143-58. [PUB_0691]
%   8. Niranjanakumari S, Day-Storms JJ, Ahmed M, Hsieh J, Zahler NH, Venters
%      RA, Fierke CA (2007). Probing the architecture of the B. subtilis RNase P
%      holoenzyme active site by cross-linking and affinity cleavage. RNA.
%      13(4): 521-35. [PUB_0692]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
classdef RNAProcessing < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'enzymeSpecificRate_DeaD';
            'enzymeSpecificRate_RsgA';
            'enzymeSpecificRate_RNAseP';
            'enzymeSpecificRate_RNAseIII';
            'enzymeSpecificRate_RNAseJ';
            'enzymeEnergyCost_DeaD';
            'enzymeEnergyCost_RsgA';            
            'reactantByproductMatrix';
            'catalysisMatrix';
            };			
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'unprocessedRNAs';
            'processedRNAs';
            'intergenicRNAs'};
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs   = {}; %whole cell model IDs of stimuli
        substrateWholeCellModelIDs = {  %whole cell model IDs of substrates
            'ATP'; 'GTP';
            'ADP'; 'GDP';
            'PI'; 'H2O'; 'H'};
        substrateIndexs_NTPs = (1:2)';
        substrateIndexs_ATP  = 1;
        substrateIndexs_GTP  = 2;
        substrateIndexs_NDPs = (3:4)';
        substrateIndexs_ADP  = 3;
        substrateIndexs_GDP  = 4;
        substrateIndexs_PI   = 5;
        substrateIndexs_H2O  = 6;
        substrateIndexs_H    = 7;
        
        enzymeWholeCellModelIDs = {     %whole cell model IDs of enzymes
            'MG_367_DIMER';             %ribonuclease III
            'MG_139_DIMER';             %ribonuclease J
            'MG_0003_465';              %ribonuclease P
            'MG_110_MONOMER';           %ribosome small subunit-dependent GTPase A
            'MG_425_DIMER'};            %ATP-dependent RNA helicase
        enzymeIndexs_RNAseIII = 1;      %index within enzymes of ribonuclease III
        enzymeIndexs_RNAseJ   = 2;      %index within enzymes of ribonuclease J
        enzymeIndexs_RNAseP   = 3;      %index within enzymes of ribonuclease P
        enzymeIndexs_RsgA     = 4;      %index within enzymes of ribosome small subunit-dependent GTPase A
        enzymeIndexs_DeaD     = 5;      %index within enzymes of ATP-dependent RNA helicase
        
        unprocessedRNAWholeCellModelIDs %whole cell model IDs of unprocessed RNAs
        processedRNAWholeCellModelIDs   %whole cell model IDs of processed RNAs
        intergenicRNAWholeCellModelIDs  %whole cell model IDs of intergenic RNAs
        
        unprocessedRNAIndexs_mRNA       %indices of unprocessed mRNAs within unprocessedRNAs
        unprocessedRNAIndexs_rRNA       %indices of unprocessed rRNAs within unprocessedRNAs
        unprocessedRNAIndexs_sRNA       %indices of unprocessed sRNAs within unprocessedRNAs
        unprocessedRNAIndexs_tRNA       %indices of unprocessed tRNAs within unprocessedRNAs
        unprocessedRNAIndexs_scRNA      %indices of unprocessed scRNAs within unprocessedRNAs
        unprocessedRNAIndexs_tmRNA      %indices of unprocessed tmRNAs within unprocessedRNAs
        
        processedRNAIndexs_mRNA         %indices of unprocessed mRNAs within processedRNAs
        processedRNAIndexs_rRNA         %indices of unprocessed rRNAs within processedRNAs
        processedRNAIndexs_sRNA         %indices of unprocessed sRNAs within processedRNAs
        processedRNAIndexs_tRNA         %indices of unprocessed tRNAs within processedRNAs
        processedRNAIndexs_scRNA        %indices of unprocessed scRNAs within processedRNAs
        processedRNAIndexs_tmRNA        %indices of unprocessed tmRNAs within processedRNAs
    end
    
    %fixed biological constants
    properties
        %enzyme kcats (reactions catalyzed per second per enzyme)
        enzymeSpecificRate_DeaD         % 1.48   1/s   B. subtilis  [PUB_0688]
        enzymeSpecificRate_RsgA         % 0.2917 1/s   E. coli      [PUB_0103]
        enzymeSpecificRate_RNAseP       % 6      1/s   B. subtilis  [PUB_0690]
        enzymeSpecificRate_RNAseIII     % 7.7    1/s   E. coli      [PUB_0691]
        enzymeSpecificRate_RNAseJ       % 0.37   1/s   B. subtilis  [PUB_0692]
        
        %NTP hydrolyzed per reaction
        enzymeEnergyCost_DeaD           % 2.37 ATP/reaction   (210 ATP/min B. subtilis [PUB_0069])/(1.48 rxn/s B. subtilis  [PUB_0688])
        enzymeEnergyCost_RsgA           % 1    GTP/reaction   [PUB_0103]
        
        %adjacency matrices
        reactantByproductMatrix         %metabolites required to cleavage each transcription unit, and metabolic byproducts of cleavage [substrates X unprocessed RNAs]
        catalysisMatrix                 %enzymes required to cleave each transcription unit, including the effect of their kcat [enzymes X unprocessed RNAs]
    end
    
    %global state (copied locally for convenience)
    properties
        unprocessedRNAs                 %counts of unprocessed RNAs
        processedRNAs                   %counts of processed RNAs
        intergenicRNAs                  %counts of intergenic RNAs
    end
    
    %constructor
    methods
        function this = RNAProcessing(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation, varargin{:});

            %% RNAs
            s = this.rna;
            this.unprocessedRNAWholeCellModelIDs = s.wholeCellModelIDs(s.nascentIndexs);
            this.processedRNAWholeCellModelIDs   = s.wholeCellModelIDs(s.processedIndexs);
            this.intergenicRNAWholeCellModelIDs  = s.wholeCellModelIDs(s.intergenicIndexs);
            
            this.processedRNAIndexs_mRNA = s.matureMRNAIndexs;
            this.processedRNAIndexs_rRNA = s.matureRRNAIndexs;
            this.processedRNAIndexs_sRNA = s.matureSRNAIndexs;
            this.processedRNAIndexs_tRNA = s.matureTRNAIndexs;
            this.processedRNAIndexs_scRNA = find(strcmp('MG_0001', this.processedRNAWholeCellModelIDs));
            this.processedRNAIndexs_tmRNA = find(strcmp('MG_0004', this.processedRNAWholeCellModelIDs));                    
            
            %% processing pathways
            %unprocessed RNAs that are precursors to processed RNAs           
            this.unprocessedRNAIndexs_mRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_mRNA, :),1));
            this.unprocessedRNAIndexs_rRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_rRNA, :),1));
            this.unprocessedRNAIndexs_sRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_sRNA, :),1));
            this.unprocessedRNAIndexs_tRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_tRNA, :),1));
            this.unprocessedRNAIndexs_scRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_scRNA, :),1));
            this.unprocessedRNAIndexs_tmRNA = find(sum(this.rna.nascentRNAMatureRNAComposition(this.processedRNAIndexs_tmRNA, :),1));
            
            %metabolites required to process each RNA, and released as
            %byproducts
            this.reactantByproductMatrix = zeros(numel(this.substrateWholeCellModelIDs), numel(this.unprocessedRNAWholeCellModelIDs));
            
            this.reactantByproductMatrix(this.substrateIndexs_H2O, :) = -(sum(this.rna.nascentRNAMatureRNAComposition, 1) + sum(this.rna.intergenicRNAMatrix, 1) - 1);
            this.reactantByproductMatrix(this.substrateIndexs_H,   :) =  (sum(this.rna.nascentRNAMatureRNAComposition, 1) + sum(this.rna.intergenicRNAMatrix, 1) - 1);
            
            energyIdxs_DeaD = [this.substrateIndexs_ATP; this.substrateIndexs_ADP; this.substrateIndexs_PI; this.substrateIndexs_H2O; this.substrateIndexs_H];
            this.reactantByproductMatrix(energyIdxs_DeaD, this.unprocessedRNAIndexs_rRNA) = ...
                this.reactantByproductMatrix(energyIdxs_DeaD, this.unprocessedRNAIndexs_rRNA) + ...
                round(this.enzymeEnergyCost_RsgA) * repmat([-1; 1; 1; -1; 1], 1, numel(this.unprocessedRNAIndexs_rRNA));
            
            energyIdxs_RsgA = [this.substrateIndexs_GTP; this.substrateIndexs_GDP; this.substrateIndexs_PI; this.substrateIndexs_H2O; this.substrateIndexs_H];
            this.reactantByproductMatrix(energyIdxs_RsgA, this.unprocessedRNAIndexs_rRNA) = ...
                this.reactantByproductMatrix(energyIdxs_RsgA, this.unprocessedRNAIndexs_rRNA) + ...
                round(this.enzymeEnergyCost_DeaD) * repmat([-1; 1; 1; -1; 1], 1, numel(this.unprocessedRNAIndexs_rRNA));
            
            %enzymes required to process each RNA (including affect of kcats)
            catalysis = zeros(numel(this.enzymeWholeCellModelIDs), numel(this.unprocessedRNAWholeCellModelIDs));
            catalysis(this.enzymeIndexs_DeaD,     this.unprocessedRNAIndexs_rRNA) = 1;
            catalysis(this.enzymeIndexs_RNAseJ,   this.unprocessedRNAIndexs_rRNA) = 2;
            catalysis(this.enzymeIndexs_RNAseIII, this.unprocessedRNAIndexs_rRNA) = 2;
            catalysis(this.enzymeIndexs_RsgA,     this.unprocessedRNAIndexs_rRNA) = 1;
            
            catalysis(this.enzymeIndexs_RNAseIII, this.unprocessedRNAIndexs_scRNA) = 2;
            catalysis(this.enzymeIndexs_RNAseP,   this.unprocessedRNAIndexs_tmRNA) = 1;
            
            catalysis(this.enzymeIndexs_RNAseIII, this.unprocessedRNAIndexs_tRNA) = 1;
            catalysis(this.enzymeIndexs_RNAseP,   this.unprocessedRNAIndexs_tRNA) = 1;
            
            enzymeSpecificRates = zeros(5,1);
            enzymeSpecificRates(this.enzymeIndexs_DeaD)     = this.enzymeSpecificRate_DeaD;
            enzymeSpecificRates(this.enzymeIndexs_RsgA)     = this.enzymeSpecificRate_RsgA;
            enzymeSpecificRates(this.enzymeIndexs_RNAseP)   = this.enzymeSpecificRate_RNAseP;
            enzymeSpecificRates(this.enzymeIndexs_RNAseIII) = this.enzymeSpecificRate_RNAseIII;
            enzymeSpecificRates(this.enzymeIndexs_RNAseJ)   = this.enzymeSpecificRate_RNAseJ;
            
            this.catalysisMatrix = catalysis ./ repmat(enzymeSpecificRates, 1, numel(this.unprocessedRNAWholeCellModelIDs));
        end
        
        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();
            
            this.unprocessedRNAs = this.rna.counts(this.rna.nascentIndexs,    this.compartment.cytosolIndexs, :);
            this.processedRNAs   = this.rna.counts(this.rna.processedIndexs,  this.compartment.cytosolIndexs, :);
            this.intergenicRNAs  = this.rna.counts(this.rna.intergenicIndexs, this.compartment.cytosolIndexs, :);
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.rna.counts(this.rna.nascentIndexs,    this.compartment.cytosolIndexs, :) = this.unprocessedRNAs;
            this.rna.counts(this.rna.processedIndexs,  this.compartment.cytosolIndexs, :) = this.processedRNAs;
            this.rna.counts(this.rna.intergenicIndexs, this.compartment.cytosolIndexs, :) = this.intergenicRNAs;
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            this.unprocessedRNAs = zeros(length(this.unprocessedRNAWholeCellModelIDs), 1, numTimePoints);
            this.processedRNAs   = zeros(length(this.processedRNAWholeCellModelIDs),   1, numTimePoints);
            this.intergenicRNAs  = zeros(length(this.intergenicRNAWholeCellModelIDs),  1, numTimePoints);
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            import edu.stanford.covert.util.ComputationUtil;
            invMat = ComputationUtil.invertCompositionMatrix(this.rna.nascentRNAMatureRNAComposition);
            
            %substrate and byproducts
            bmProd = max(0, -this.reactantByproductMatrix) * invMat * states.rnaProductions;
            byProd = max(0,  this.reactantByproductMatrix) * invMat * states.rnaProductions;
            
            %enzymes
            minEnzExp = 2 * this.catalysisMatrix * invMat * states.rnaProductions0;
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization: RNAs initialized to the mature/aminoacylated state by
        %simulation initializeState method
        function initializeState(~)
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = max(0, -this.reactantByproductMatrix) * min( ...
                ceil(min(this.enzymes(:, ones(size(this.catalysisMatrix, 2), 1)) ./ this.catalysisMatrix * this.stepSizeSec, [], 1))', ...
                this.unprocessedRNAs); %#ok<UDIM>
        end

        %simulation
        function evolveState(this)
            rnaTypes = {'mRNA'; 'rRNA'; 'sRNA'; 'tRNA'};
            
            order = this.randStream.randperm(4);
            for i = 1:4
                this.evolveState_Helper(rnaTypes{order(i)});
            end
        end
        
        function evolveState_Helper(this, rnaType)
            %indices of unprocessed RNAs of selected type
            unprocessedRNAIndexs = this.(['unprocessedRNAIndexs_' rnaType]);
            
            %limit RNA processing by availability of
            %- unprocessed RNAs
            %- enzymes
            %- free reactants
            %
            %Assumes that each RNA of this type uses same amount of substrates
            %(except water and hydrogen) and enzymes. There is an assertion in
            %the test class to check that this is true. For water and hydrogen
            %this upper bound calculation uses that most amount of water and
            %hydrogen required across all RNAs of the current type (Note: the
            %correct stoichiometry of water and hydrogen are used to update
            %state).
            totRNAs = sum(this.unprocessedRNAs(unprocessedRNAIndexs));
            if totRNAs == 0
                return;
            end
            
            substrates = this.substrates ./ max(0, -min(this.reactantByproductMatrix(:, unprocessedRNAIndexs), [], 2));
            numReactions = floor(min([
                totRNAs;
                this.randStream.stochasticRound((this.enzymes * this.stepSizeSec) ./ this.catalysisMatrix(:, unprocessedRNAIndexs(1)));
                substrates ...
                ]));
            if numReactions <= 0
                return;
            end
            
            %randomly select RNAs to process
            processingRNAs = this.randStream.randCounts(this.unprocessedRNAs(unprocessedRNAIndexs), numReactions);
            
            %increment processed RNAs
            this.processedRNAs = this.processedRNAs + ...
                this.rna.nascentRNAMatureRNAComposition(:, unprocessedRNAIndexs) * processingRNAs;
            
            %increment intergenic RNAs
            this.intergenicRNAs = this.intergenicRNAs + ...
                this.rna.intergenicRNAMatrix(:, unprocessedRNAIndexs) * processingRNAs;
            
            %decrement unprocessed RNAs
            this.unprocessedRNAs(unprocessedRNAIndexs) = ...
                this.unprocessedRNAs(unprocessedRNAIndexs) - processingRNAs;
            
            %decrement reactants, increment byproducts
            this.substrates = this.substrates + ...
                this.reactantByproductMatrix(:, unprocessedRNAIndexs) * processingRNAs;
        end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.unprocessedRNAs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    this.rna.molecularWeights(this.rna.nascentIndexs)'    * this.unprocessedRNAs + ...
                    this.rna.molecularWeights(this.rna.processedIndexs)'  * this.processedRNAs   + ...
                    this.rna.molecularWeights(this.rna.intergenicIndexs)' * this.intergenicRNAs) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(this.rna.molecularWeights(this.rna.nascentIndexs)'    * permute(this.unprocessedRNAs,[1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.processedIndexs)'  * permute(this.processedRNAs,  [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.intergenicIndexs)' * permute(this.intergenicRNAs, [1 3 2]),[1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
