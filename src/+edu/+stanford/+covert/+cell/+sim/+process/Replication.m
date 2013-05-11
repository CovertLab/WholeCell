%Replication
%
% @wholeCellModelID Process_Replication
% @name             Replication
% @description
%   Biology
%   =================
%   DNA replication by replisomes is initiated by DnaA complex formation near
%   the oriC and proceeds bidirectionally from oriC to terC in the 5'->3'
%   direction along leading strand as well as in Okazaki fragments in the 5'->3'
%   direction along lagging strands. This process models the locations, kinetics,
%   and biochemistry of all the replication proteins:
%   - replicative DNA helicase
%   - DNA primase
%   - DNA polymerase core
%   - Beta-clamp loading complex / gamma complex
%   - Beta-clamp
%   - DNA ligase
%   - single stranded binding proteins
%
%   The exact mechanism replication initiation in M. genitalium is unknown.
%   Furthermore, because M. genitalium does not contain a DnaC homolog it is
%   difficult to infer M. genitalium repliciation initiation from studies of
%   other bacterial species.
%
%   Knowledge Base
%   =================
%   The knowledge base contains data curated from the literature and databases:
%   - the chromosome DNA sequence
%   - the DNA footprints of all replication proteins
%   - DNA binding protein displacement reactions (that is which proteins can
%     displace which other proteins from the chromosome)
%   - values of various structural and kinetic parameters
%   - subunit composition of the replication proteins
%
%   Representation
%   =================
%   The properties substrates, enzymes, and boundEnzymes represent the counts of
%   free metabolites, free replication proteins, and chromosally-bound
%   replication proteins. The chromosomes property polymerizedRegions represents
%   the regions of the chromosomes which have been polymerized (and the
%   base-pairing of strands). The chromosomes property strandBreaks represents
%   strand breaks 5' to each base. The chromosomes property complexBoundSites
%   represents represents the specific chromosomal location of all chromosomally
%   bound proteins.
%
%   Spatial Model
%   +++++++++++++++++
%   - Helicase is centered on the boundary between ssDNA and dsDNA. The position
%     over which it is centered is the next position to be melted.
%   - Polymerase core is centered on the boundary between ssDNA and dsDNA. The
%     position over which it is centered is the next position to be polymerized.
%   - There is no gap between the helicase and polymerase core or between the
%     polymerase core and the beta clamp.
%   - Backup beta clamps bind slightly upstream of the start site of Okazaki
%     fragments such that there will be no gap between the polymerase and beta
%     clamp, and the polymerase core will be centered on the Okazaki fragment
%     start site
%   - At replication initiation a the mother strands are separated such that the
%     leading polymerase cores are centered at oriC+-1 and the helicases are
%     slightly (11 nt) ahead
%   - During replication initiation (and the final step of replication after the
%     last Okazaki fragment has completed) the lagging polymerase and primase
%     are accounted for as part of a complex on the leading strand (containing
%     also the helicase, leading polymerase, gamma complex, and leading beta
%     clamp). At all other times, the lagging polymerase, lagging beta clamp,
%     and primase are accounted for as a complex on a different strand. This
%     allows us to separately keep track of the leading and lagging polymerase
%     positions.
%   - Backup beta clamps can't binding until half of the previous Okazaki
%     fragment has been polymerized
%   - In general there will be approximately 1 Okazaki fragment length gap
%     between the progess of DNA polymerase on the leading and lagging strands
%
%   Initialization
%   =================
%   The mother chromosome is initialized by the Chromosome class and
%   several processes:
%   - completely synthesized
%   - mother strands base pairing
%   - undamged, except methylated at restriction/modification sites (DNA repair
%     process)
%   - bound by various proteins
%     - SMCs (DNA condensation)
%     - RNA polymerase (Transcription)
%     - Transcription factors (Transcriptional regulation)
%     - Topoisomerases (supercoiling)
%     - DnaA (replication initiation)
%
%   Simulation
%   =================
%   The simulation consists of 8 subfunctions executed in a random order:
%   - Initiate replication (initiateReplication)
%     If DnaA complex assembled at OriC and sufficient protein and metabolites,
%     unwind small segment of DNA and binding helicase, primase, polymerase,
%     gamma complex, and beta clamp to chromosomes. Associate all proteins with
%     the leading strand. Chromosome will take care of dissassembly of the
%     DnaA complex. ReplicationInitiation will take care of dissociating the
%     released DnaA-ATP polymers.
%
%   - Advance replisomes: unwind and polymerize DNA, release SSBs (unwindAndPolymerizeDNA)
%     1. If first Okazaki fragment starting, associate primase, lagging polymerase
%        with lagging strand.
%     2. Advance leading and lagging polymerases and helicases up to limits
%        1. Polymerase and primase kinetics
%        2. Available dNTP for polymerization and energy for unwinding
%        3. Prevent leading strand progressing if no SSBs bound to lagging
%           strand
%        4. Accessibility of upstream regions to helicase and polymerases
%        5. Don't polymerize leading strand past terC or lagging strands past
%           ends of Okazaki fragments
%        6. Don't allow the leading strand and helicase to go way beyond the
%           lagging strand progress.
%        7. Pause progress is RNA polymerase is encountered.
%           If the replication loop (helicase) hits an RNAP, polymerization
%           pauses, the RNAP falls off, and its transcript is degraded. If it's
%           a head-on collision, the replication loop will not proceed for some
%           unknown amount of time (a fittable parameter). If it's a codirec-
%           tional collision, polymerization will continue at full speed the
%           following time step. (Mirkin 2004, Mirkin 2006, Mirkin 2007).
%           If occupied by RNA-polymerase, then calculate "occupied DNA" by base
%           polymerase is currently on -9 through +2. (Neidhardt 1990).
%
%   - Bind SSBs (freeAndBindSSBs)
%     Bind free SSB 4mers to stochastically to deterministically selected
%     positions within single-stranded regions as SSB 8mers
%
%   - Dissociate free SSB 8mers into 2 SSB 4mers (dissociateFreeSSBComplexes)
%     Dissociate released SSB 8mers into 2 SSB 4mers
%
%   - Initiate Okazaki fragment by bind beta-clamp (initiateOkazakiFragment)
%     Binding beta-clamp just downstream of the start position of the next
%     Okazaki fragment. Requires that
%     - There is beta-clamp monomer to form new beta-clamp dimer on chromosome
%     - There is energy for the gamma-complex to form the new beta-clamp
%     - The position is accessible to the beta-clamp
%     - The leading helicase has already passed the position
%     - The lagging strand is at least 1/2 done with the current Okazaki
%       fragment
%
%   - Terminate Okazaki fragment by releasing beta-clamp (terminateOkazakiFragment)
%     Release beta-clamp and associate lagging primase and polymerase with
%     backup beta-clamp (or if terminating the last Okazaki fragment,
%     associate the lagging primase and with the leading strand machinery).
%     Mark end of Okazaki fragment as having a single strand break to be ligated
%     by ligase.
%     Occurs if:
%     - Okazaki fragment finished polymerizing
%     - SSBs bound to lagging strand
%     - Lagging backup beta-clamp bound
%     - Leading machinery has advanced beyond the Okazaki fragment start site
%
%   - Terminate replication (terminateReplication)
%     Release bound replisomes from leading strands. Mark terC as having single
%     strand breaks to be ligated by ligase. Occurs if:
%     - Polymerization completed for both leading and lagging strands
%     - Lagging strand ligated (except for last ligation at terC which will
%       occur after the execution of this subfunction)
%
%   - Ligate DNA (ligateDNA)
%     Stochastically ligate single strand breaks up to
%     - Ligase kinetics
%     - Ligase availability
%     - NAD availability
%
%   References
%   =================
%   1. Mirkin, E.V., Mirkin, S.M. (2007). Replication fork stalling at natural
%      impediments. Microbiology and molecular biology reviews 71: 13-35.
%   2. Mirkin, E.V., Mirkin, S.M. (2005). Mechanisms of
%      Transcription-replication collisions in bacteria. molecular and cellular
%      biology 25: 888-895.
%   3. Mirkin, E.V., Roa, D.C., Nudler, E., Mirkin, S.M. (2006). Transcription
%      regulatory elements are punctuation marks for dna replication. PNAS 103:
%      7276-7281.
%   4. Miyata, M. "Cell Division". Molecular biology and pathogenicity of
%      mycoplasmas. Razin, S., Herrmann, R. Kluwer Academic/Plenum Publishers,
%      2002. 117-130.
%   5. Rudolph, C.J., Dhillon, P., Moore, T., Lloyd, R.G. (2007). Avoiding and
%      resolving conflicts between DNA replication and transcription. DNA Repair
%      6: 981-993.
%   6. McGlynn, P., Guy, C.P. (2008). Replication forks blocked by protein-DNA
%      complexes have limited stability in vitro. J. Mol. Biol. 381: 249-255.
%   7. Kozlov AG and Lohman TM. (2002). Kinetic Mechanism
%      of Direct Transfer of Escherichia coli SSB Tetramers between
%      Single-Stranded DNA Molecules. Biochemistry. 41 (39): 11611-11627.
%      [PUB_0856]
%   8. Kunzelmann S, Morris C, Chavda AP, Eccleston JF, Webb MR (2010).
%      Mechanism of Interaction between Single-Stranded DNA Binding Protein
%      and DNA. Biochemistry 49(5): 843-52. [PUB_0857]
%   9. Roy R, Kozlov AG, Lohman TM, Ha T (2009). SSB protein diffusion on
%      single-stranded DNA stimulates RecA filament formation. Nature.
%      461(7267):1092-7. [PUB_0858]
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/19/2010

%TODO: consider moving requirement for SSBs in okazaki fragment termination to
%      initiation of next fragment
classdef Replication < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'oriCPosition'
            'terCPosition'
            'enzymeComposition'
            'primaseBindingLocations'
            'primerLength'
            'dnaPolymeraseElongationRate'
            'ssbComplexSpacing'
            'okazakiFragmentMeanLength'
            'ligaseRate'
            'laggingBackupClampReloadingLength'
            'startingOkazakiLoopLength'
            'rnaPolymeraseCollisionMeanDwellTime'
            'ssbDissociationRate'
            'dnaAFunctionalBoxStartPositions'
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
        
        leadingStrandIndexs = [1 4];
        laggingStrandIndexs = [3 2];
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli
        
        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'DATP'; 'DCTP'; 'DGTP'; 'DTTP';
            'ATP'; 'CTP'; 'GTP'; 'UTP';
            'PPI'; 'H2O'; 'H'; 'NAD'; 'NMN'; 'ADP'; 'AMP'; 'PI';
            };
        substrateIndexs_dntp        = (1:4)';
        substrateIndexs_ntp         = (5:8)';
        substrateIndexs_diphosphate = 9;
        substrateIndexs_water       = 10;
        substrateIndexs_hydrogen    = 11;
        substrateIndexs_nad         = 12;
        substrateIndexs_nmn         = 13;
        substrateIndexs_atp         = 5;
        substrateIndexs_adp         = 14;
        substrateIndexs_amp         = 15;
        substrateIndexs_phosphate   = 16;
        
        enzymeWholeCellModelIDs = { %enzyme whole cell model ids
            'REPLISOME';                                            %replisome
            'DNA_POLYMERASE_2CORE_BETA_CLAMP_GAMMA_COMPLEX_PRIMASE';%DNA-directed DNA polymerase (2) core + beta-clamp + gamma-complex + primase
            'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX';         %DNA-directed DNA polymerase core + beta-clamp + gamma-complex
            'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE';               %DNA-directed DNA polymerase core + beta-clamp + primase
            'DNA_POLYMERASE_CORE';                                  %DNA-directed DNA polymerase core
            'DNA_POLYMERASE_GAMMA_COMPLEX';                         %DNA-directed DNA polymerase gamma complex
            'MG_001_MONOMER';                                       %DNA polymerase III, beta subunit
            'MG_001_DIMER';                                         %DNA polymerase III, beta clamp
            'MG_094_HEXAMER';                                       %replicative DNA helicase
            'MG_254_MONOMER';                                       %DNA ligase, NAD-dependent
            'MG_250_MONOMER';                                       %DNA primase
            'MG_091_TETRAMER';                                      %single-strand binding protein family, tetramer
            'MG_091_OCTAMER';                                       %single-strand binding protein family, octamer
            };
        enzymeIndexs_replisome                         = 1;
        enzymeIndexs_2coreBetaClampGammaComplexPrimase = 2;
        enzymeIndexs_coreBetaClampGammaComplex         = 3;
        enzymeIndexs_coreBetaClampPrimase              = 4;
        enzymeIndexs_core                              = 5;
        enzymeIndexs_gammaComplex                      = 6;
        enzymeIndexs_betaClampMonomer                  = 7;
        enzymeIndexs_betaClamp                         = 8;
        enzymeIndexs_helicase                          = 9;
        enzymeIndexs_ligase                            = 10;
        enzymeIndexs_primase                           = 11;
        enzymeIndexs_ssb4mer                           = 12;
        enzymeIndexs_ssb8mer                           = 13;
        
        complexIndexs_DnaA_1mer_ATP       %index of DnaA-ATP 1mer within simulation.complexs
        complexIndexs_DnaA_7mer_ATP       %index of DnaA-ATP 7mer within simulation.complexs

        dnaAFunctionalBoxIndexs_R1234     %index of R1-4 functional DnaA boxes within dnaAFunctionalBoxStartPositions
        dnaAFunctionalBoxIndexs_R5        %index of R5 functional DnaA box within dnaAFunctionalBoxStartPositions
    end
    
    %fixed biological constants
    properties
        oriCPosition                       %oriC position
        terCPosition                       %terC position
        
        enzymeComposition                  %replisome subunit composition       
        
        primaseBindingLocations            %locations (randomly determined) where primers will bind the lagging strand to make okazaki fragments
        
        primerLength                       %length of primer that needs to form before an okazaki fragment can be replicated by DNA polymerase (nt) (11) [PUB_0606]
        dnaPolymeraseElongationRate        %DNA polymerase elongation rate (nt/s) (100) [Miyata 2002]
        ssbComplexSpacing                  %spacing between SSB octamers (nt) [PUB_0637]
        okazakiFragmentMeanLength          %Average length of an Okazaki fragment (nt) [PUB_0623, Kornberg 1980, Xie 2008]
        ligaseRate                         %Ligase Vmax (1/s) (0.04) [PUB_0488]
        laggingBackupClampReloadingLength  %Length at which a second backup clamp binds the loading complex in terms of how many bases have been polymerized of the Okazaki fragment (750)
        startingOkazakiLoopLength          %Length of DNA (bases) in the okazaki fragment loop right at the beginning of the fragment polymerization (50)
        rnaPolymeraseCollisionMeanDwellTime%mean time of DNA polymerase stalling upon interaction with anti-directional RNA polymerase (s) [1.7 PUB_0799] 
        ssbDissociationRate                %dissociation rate of SSBs (1/s) [0.3 PUB_0857]
        
        dnaAFunctionalBoxStartPositions    %sorted start coordinates of functional DnaA boxes
    end
    
    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent = true)
        isDnaAORIComplexAssembled      %boolean indicating whether or not a DnaA-ATP complex has formed at the oriC, if true replication can initiate
        isAnyHelicaseBound             %boolean indicating whether or not any helicase is bound
        isAnyPolymeraseBound           %boolean indicating whether or not any polymerase is bound
        leadingStrandElongating        %whether replication elongation is in progress
        laggingStrandElongating        %whether there is a backup clamp bound to each loading complex
        leadingStrandPolymerized       %whether leading strands have been polymerized
        laggingStrandPolymerized       %whether lagging strands have been polymerized
        strandPolymerized              %whether leading and lagging strands have been polymerized
        numLigations                   %number of ligations that have occurred (in replication order along lagging strand)
        strandLigated                  %whether or not strands are completely ligated
        strandDuplicated               %whether replication is complete and chromosome segregation can start
        
        helicasePosition               %helicase start coordinates
        leadingPolymerasePosition      %leading polymerase start coordinates
        laggingPolymerasePosition      %lagging polymerase start coordinates
        leadingPosition                %next position to be polymerized on leading strands
        laggingPosition                %next position to be polymerized on lagging strands
        laggingBackupBetaClampPosition %lagging backup beta clamp start coordinates
        
        okazakiFragmentIndex           %index of the current/most recent Okazaki fragments
        okazakiFragmentPosition        %5' position of the current/most recent Okazaki fragments, including the primer
        okazakiFragmentLength          %length of the current/most recent Okazaki fragments, including the primer
        okazakiFragmentProgress        %number of nucleotides polymerized of the current/most recent Okazaki fragments
        
        leadingStrandBoundSSBs         %sparse vector indicating where SSB complexes are bound
        laggingStrandBoundSSBs         %sparse vector indicating where SSB complexes are bound
        numLeadingTemplateBoundSSBs    %number of SSB 8mers bound to the leading strand's template
        numLaggingTemplateBoundSSBs    %number of SSB 8mers bound to the lagging strand's template
        areLaggingStrandSSBSitesBound  %true/false whether or not the lagging strand's template is sufficiently covered in SSB 8mers for replication to progress
    end
    
    %constructor
    methods
        function this = Replication(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(knowledgeBase, simulation, varargin{:});
            
            c = this.chromosome;
            
            %oriC position
            oriC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'oriC');
            this.oriCPosition = oriC.startCoordinate;
            
            %terC position
            terC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'terC');
            this.terCPosition = terC.startCoordinate;
            
            %enzyme subunit composition
            this.enzymeComposition = zeros(numel(this.enzymeWholeCellModelIDs));
            this.enzymeComposition(this.enzymeMonomerLocalIndexs, this.enzymeComplexLocalIndexs) = sum(knowledgeBase.proteinComplexMonomerComposition(this.gene.mRNAIndexs(this.enzymeMonomerGlobalIndexs), this.enzymeComplexGlobalIndexs, :), 3);
            this.enzymeComposition(this.enzymeComplexLocalIndexs, this.enzymeComplexLocalIndexs) = sum(knowledgeBase.proteinComplexComplexComposition(this.enzymeComplexGlobalIndexs, this.enzymeComplexGlobalIndexs, :), 3);
                         
            %primase binding locations
            this.primaseBindingLocations = {
                c.sequenceLen + 1 - this.calculatePrimaseBindingLocations(c.sequenceLen - this.terCPosition)
                this.calculatePrimaseBindingLocations(this.terCPosition)
                }';
            
            %DnaA and its binding sites
            this.complexIndexs_DnaA_1mer_ATP = this.complex.getIndexs({'MG_469_1MER_ATP'});
            this.complexIndexs_DnaA_7mer_ATP = this.complex.getIndexs({'MG_469_7MER_ATP'});
            
            dnaABoxes = findobj(knowledgeBase.genomeFeatures, 'type', 'DnaA box');            
            [~, order] = sort([dnaABoxes.startCoordinate]);
            dnaABoxes = dnaABoxes(order);
            
            dnaABoxNames = {dnaABoxes.name}';
            [tfs, idxs] = ismember(dnaABoxNames, {
                'Functional box R1';
                'Functional box R2';
                'Functional box R3';
                'Functional box R4';
                'Functional box R5';
                });
            
            dnaABoxes = dnaABoxes(tfs);
            this.dnaAFunctionalBoxStartPositions = ...
                ceil([dnaABoxes.startCoordinate]' + [dnaABoxes.sequenceLength]'/2 - c.complexDNAFootprints(this.complexIndexs_DnaA_1mer_ATP)/2);
            
            this.dnaAFunctionalBoxIndexs_R1234 = [
                find(idxs(tfs) == 1);
                find(idxs(tfs) == 2);
                find(idxs(tfs) == 3);
                find(idxs(tfs) == 4)];
            this.dnaAFunctionalBoxIndexs_R5 = find(idxs(tfs) == 5);
        end
        
        %Choose Okazaki fragment lengths randomly according to poisson distribution
        %with mean = okazakiFragmentMeanLength.
        function pos = calculatePrimaseBindingLocations(this, len)
            pos = this.randStream.random('poiss', this.okazakiFragmentMeanLength, ceil(1.2 * len/this.okazakiFragmentMeanLength), 1);
            while sum(pos) < len
                pos = [pos; this.randStream.random('poiss', this.okazakiFragmentMeanLength, 10, 1)]; %#ok<AGROW>
            end
            
            pos = [0; cumsum(pos)];
            pos = len - pos(end:-1:1);
            idx = find(pos <= 0, 1, 'last');
            
            pos(1:idx) = [];
            if pos(1) < 0.6 * this.okazakiFragmentMeanLength
                pos(1) = [];
            end
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, constants, ~)
            %% import classes
            import edu.stanford.covert.util.ConstantUtil;
            
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            c = this.chromosome;
            m = this.metabolite;
            mass = constants.states.Mass;
            nOkazakiFragments = sum(cellfun(@numel, this.primaseBindingLocations));
            
            baseComp = getBaseCounts(c.sequence);
            dnmpComp = [baseComp; nnz(c.damagedBases == m.m6ADIndexs)];
            dnmpComp(1) = dnmpComp(1) - nnz(c.damagedBases == m.m6ADIndexs);
            chrWt = dnmpComp' * (m.molecularWeights([m.dnmpIndexs; m.getIndexs({'m6dAMP'})]) - ConstantUtil.elements.O - ConstantUtil.elements.H) / ConstantUtil.nAvogadro;
            dntpProd = ...
                baseComp * min(1, mass.dryWeightFractionDNA * mass.cellInitialDryWeight / chrWt) + ...       
                baseComp * max(0, mass.dryWeightFractionDNA * mass.cellInitialDryWeight - chrWt) / ...
                (baseComp' * m.molecularWeights(m.dntpIndexs)) * ConstantUtil.nAvogadro;
            basesPol = c.sequenceLen * min(1, mass.dryWeightFractionDNA * mass.cellInitialDryWeight / chrWt);
           
            %ATP for helicase to unwind mother chromosome
            bmProd(this.substrateIndexs_atp)         = bmProd(this.substrateIndexs_atp)         + basesPol;
            bmProd(this.substrateIndexs_water)       = bmProd(this.substrateIndexs_water)       + basesPol;
            byProd(this.substrateIndexs_adp)         = byProd(this.substrateIndexs_adp)         + basesPol;
            byProd(this.substrateIndexs_phosphate)   = byProd(this.substrateIndexs_phosphate)   + basesPol;
            byProd(this.substrateIndexs_hydrogen)    = byProd(this.substrateIndexs_hydrogen)    + basesPol;
            
            %ATP used by gamma-complex to form beta-clamps
            nBetaClamps = sum(cellfun(@(x) ~isempty(x), this.primaseBindingLocations)) + nOkazakiFragments;
            bmProd(this.substrateIndexs_atp)         = bmProd(this.substrateIndexs_atp)         + nBetaClamps;
            bmProd(this.substrateIndexs_water)       = bmProd(this.substrateIndexs_water)       + nBetaClamps;
            byProd(this.substrateIndexs_adp)         = byProd(this.substrateIndexs_adp)         + nBetaClamps;
            byProd(this.substrateIndexs_phosphate)   = byProd(this.substrateIndexs_phosphate)   + nBetaClamps;
            byProd(this.substrateIndexs_hydrogen)    = byProd(this.substrateIndexs_hydrogen)    + nBetaClamps;            
            
            %dNTPs to polymerize new strands by end of replication phase
            bmProd(this.substrateIndexs_dntp)        = bmProd(this.substrateIndexs_dntp)        + dntpProd;
            byProd(this.substrateIndexs_diphosphate) = byProd(this.substrateIndexs_diphosphate) + 2 * basesPol;
            
            %Ligation: (2) DR5P + NAD ==> AMP + dRibose5P_dRibose5P + H + NMN
            nLigations = sum(cellfun(@(x) ~isempty(x), this.primaseBindingLocations)) + nOkazakiFragments;
            bmProd(this.substrateIndexs_nad)         = bmProd(this.substrateIndexs_nad)          + nLigations;
            byProd(this.substrateIndexs_nmn)         = byProd(this.substrateIndexs_nmn)          + nLigations;
            byProd(this.substrateIndexs_amp)         = byProd(this.substrateIndexs_amp)          + nLigations;
            byProd(this.substrateIndexs_hydrogen)    = byProd(this.substrateIndexs_hydrogen)     + nLigations;
            
            %% enzymes            
            %replisome: 
            %- helicase
            %- polymerase core (2)
            %- gamma-complex
            %- beta-clamp (2 dimers = 4 monomers)
            %- primase
            minEnzExp(this.enzymeIndexs_replisome) = 2;
            
            %lagging stand backup beta-clamp
            minEnzExp(this.enzymeIndexs_betaClamp) = 2;
            
            %ligase
            nLigations = sum(cellfun(@(x) ~isempty(x), this.primaseBindingLocations)) + sum(cellfun(@numel, this.primaseBindingLocations));
            minEnzExp(this.enzymeIndexs_ligase) = nLigations / this.ligaseRate / constants.states.Time.replicationDuration;

            %single-stranded binding protein
            ssbFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_ssb8mer);
            ssDNALen = 2 * this.okazakiFragmentMeanLength;
            minEnzExp(this.enzymeIndexs_ssb8mer) = 2 * ceil(ssDNALen / ssbFtpt);
            
            %account for timing of replication
            minEnzExp = minEnzExp / exp(log(2) * constants.states.Time.replicationInitiationDuration / constants.states.Time.cellCycleLength);
        end
        
        %initialization: simulation initialize state method initializes 1
        %chromosome
        function initializeState(~)
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            c = this.chromosome;
            
            %initialize
            result = zeros(size(this.substrates));
            
            %ATP, water for initiation
            if ~this.isAnyPolymeraseBound && this.isDnaAORIComplexAssembled && all(this.enzymes >= 2 * this.enzymeComposition(:, this.enzymeIndexs_2coreBetaClampGammaComplexPrimase))
                helFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_helicase);
                corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
                uwdLen = helFtpt3 + corFtpt5 + 1;
                result(this.substrateIndexs_atp)   = result(this.substrateIndexs_atp)   + 2 * (1 + uwdLen);
                result(this.substrateIndexs_water) = result(this.substrateIndexs_water) + 2 * (1 + uwdLen);
            end
            
            %ATP for helicase to unwind mother chromosome
            if this.isAnyHelicaseBound && all(this.leadingStrandElongating)
                nLeadingStrandElongating = 2;
                nLaggingStrandElongating = nnz(this.laggingStrandElongating);
                isOriCCplx = this.isDnaAORIComplexAssembled;
                
                unwindLen = min(2, 2 * isOriCCplx + nLeadingStrandElongating) * this.dnaPolymeraseElongationRate;
                result(this.substrateIndexs_atp)   = result(this.substrateIndexs_atp)   + unwindLen;
                result(this.substrateIndexs_water) = result(this.substrateIndexs_water) + unwindLen;
                
                %ATP used by gamma-complex to form beta-clamps
                nBetaClamps = min(2, 2 * isOriCCplx + nLaggingStrandElongating) - nnz(this.laggingBackupBetaClampPosition);
                result(this.substrateIndexs_atp)   = result(this.substrateIndexs_atp)   + nBetaClamps;
                result(this.substrateIndexs_water) = result(this.substrateIndexs_water) + nBetaClamps;
                
                %dNTPs to polymerize new strands
                dnmpComp = [
                    1 - c.sequenceGCContent;
                    c.sequenceGCContent;
                    c.sequenceGCContent;
                    1-c.sequenceGCContent] / 2;
                result(this.substrateIndexs_dntp) = ...
                    + result(this.substrateIndexs_dntp) ...
                    + min(4, 2 * isOriCCplx + nLeadingStrandElongating + nLaggingStrandElongating) * ...
                    this.dnaPolymeraseElongationRate * dnmpComp;
            end
            
            %Ligation: (2) DR5P + NAD ==> AMP + dRibose5P_dRibose5P + H + NMN
            result(this.substrateIndexs_nad) = ...
                + result(this.substrateIndexs_nad) ...
                + min(nnz(c.strandBreaks), ceil(this.enzymes(this.enzymeIndexs_ligase) * this.stepSizeSec * this.ligaseRate));
        end
        
        %simulation
        function evolveState(this)
            %% check that replisome state is in sync
            cnts = this.boundEnzymes([
                this.enzymeIndexs_2coreBetaClampGammaComplexPrimase
                this.enzymeIndexs_coreBetaClampGammaComplex
                this.enzymeIndexs_coreBetaClampPrimase
                this.enzymeIndexs_helicase
                ]);
            totPolCnts = [2 1 1] * cnts(1:3);
            if ~((totPolCnts == 0 || totPolCnts == 4) && ((cnts(4) == 0 && totPolCnts == 0) || cnts(4) == 2))
                msg = [];
                msg = [msg sprintf('Replisome state out of sync\n')];
                msg = [msg sprintf('- 2coreBetaClampGammaComplexPrimase %d\n', cnts(1))];
                msg = [msg sprintf('- coreBetaClampGammaComplex %d\n', cnts(2))];
                msg = [msg sprintf('- coreBetaClampPrimase %d\n', cnts(3))];
                msg = [msg sprintf('- helicase %d', cnts(4))];
                throw(MException('Replication:error', msg));
            end
            
            %% simulate
            subfunctions = {
                @this.initiateReplication;        %Initiate replication
                @this.ligateDNA;                  %Ligate DNA
                @this.dissociateFreeSSBComplexes; %Dissociate free SSB 8mers into 2 SSB 4mers
                @this.freeAndBindSSBs;            %Free, Bind SSBs
                };
            
            %Only execute if active replisomes
            if this.isAnyHelicaseBound && all(this.leadingStrandElongating)
                subfunctions = [
                    subfunctions;{
                    @this.unwindAndPolymerizeDNA;     %Advance replisomes: unwind and polymerize DNA, release SSBs
                    @this.initiateOkazakiFragment;    %Initiate Okazaki fragment by bind beta-clamp
                    @this.terminateOkazakiFragment;   %Terminate Okazaki fragment by releasing beta-clamp
                    @this.terminateReplication;       %Terminate replication
                    }];
            end
            
            order = this.randStream.randperm(numel(subfunctions));
            for i = 1:numel(subfunctions)
                subfunctions{order(i)}();
            end
        end
        
        %Initiate replication if DnaA complex at oriC and sufficient proteins
        %to form two replisomes:
        %- form two replisomes at origin
        %- indirectly (via Chromosome) disassemble DnaA complex
        function initiateReplication(this)
            chrLen = this.chromosome.sequenceLen;
            holFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            helFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_helicase);
            helFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_helicase);
            corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
            helCorSpacing = helFtpt3 + corFtpt5 + 1;
            
            %check if DnaA complex at oriC and sufficient protein to form 2
            %replisomes
            if ~this.isDnaAORIComplexAssembled || any(this.leadingPolymerasePosition)
                return;
            end
            
            %warn if this would be a secondary initiation
            if collapse(this.chromosome.polymerizedRegions) > 2 * chrLen
                warning('WholeCell:warning', 'Unable to represent multiple replication initiation events per cell cycle');
                return;
            end
            
            helPos = this.helicasePosition;
            maxSteps = floor(min(this.substrates([this.substrateIndexs_atp; this.substrateIndexs_water])) / 2);
            if all(helPos)
                uwdOld = mod(-(1 - helFtpt3 - helPos(2)) - 1, chrLen) + 1;
                uwdLen = min(helCorSpacing - uwdOld, maxSteps);
                atpCost = 2 * uwdLen;
                if uwdLen <= 0
                    return;
                end
                this.releaseProteinFromSites([helPos; this.leadingStrandIndexs]', false);
            elseif ~any(helPos)
                uwdLen = min(helCorSpacing, maxSteps - 1);
                uwdOld = 0;
                if uwdLen < 0 || this.enzymes(this.enzymeIndexs_helicase) < 2
                    return;
                end
                atpCost = 2 * (1 + uwdLen);
            else
                throw(MException('Replication:error', 'All or no helicases should be bound'));
            end
            
            %bind helicase and unwind chromosome
            this.chromosome.setRegionUnwound(chrLen-uwdOld, -uwdLen);
            this.chromosome.setRegionUnwound(1+uwdOld, uwdLen);
            if ~all(this.bindProteinToChromosome([chrLen-helFtpt5-uwdOld 1-helFtpt3+uwdOld; this.leadingStrandIndexs]', this.enzymeIndexs_helicase, ...
                    [], [], true, false, [-uwdLen-1; uwdLen+1], true, []))
                throw(MException('Replication:error', 'Both helicase must bind'));
            end
            
            %hydrolyze ATP used by gamma-complex to assemble beta-clamp, and by
            %helicase to melt DNA
            this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - atpCost;
            this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - atpCost;
            this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + atpCost;
            this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + atpCost;
            this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + atpCost;
            
            %Form two replisomes (each with 2 cores, 1 gamma-complex, and 1 beta-clamp (MG_001 dimer)) at oriC
            if uwdOld + uwdLen == helCorSpacing && all(this.enzymes >= 2 * this.enzymeComposition(:, this.enzymeIndexs_2coreBetaClampGammaComplexPrimase))
                this.enzymes = this.enzymes - 2 * this.enzymeComposition(:, this.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
                this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) + 2;
                
                %bind DNA polymerase core and beta-clamp to initiate leading strand
                %synthesis
                if ~all(this.bindProteinToChromosome([chrLen-corFtpt5 1+corFtpt5-holFtpt+1; this.leadingStrandIndexs]', ...
                        this.enzymeIndexs_2coreBetaClampGammaComplexPrimase, [], [], true, false, 1, true, []))
                    throw(MException('Replication:error', 'Both polymerases must bind'));
                end
            end
        end
        
        %Advance replisomes:
        %- unwind DNA
        %- polymerize DNA
        %
        %TODO: RNA primer
        function unwindAndPolymerizeDNA(this)
            %% Terminate early if no active replisomes
            if ~this.isAnyHelicaseBound || ~all(this.leadingStrandElongating)
                return;
            end
            
            %% indices, parameters
            helGblIdx = this.enzymeGlobalIndexs(this.enzymeIndexs_helicase);
            leadPolLclIdx = [this.enzymeIndexs_2coreBetaClampGammaComplexPrimase; this.enzymeIndexs_coreBetaClampGammaComplex];
            leadPolGblIdx = this.enzymeGlobalIndexs(leadPolLclIdx);
            lagPolGblIdx = this.enzymeGlobalIndexs(this.enzymeIndexs_coreBetaClampPrimase);
            helFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_helicase);
            helFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_helicase);
            helFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_helicase);
            holFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_core);
            corFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_core);
            corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
            bClmpFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_betaClamp);
            
            c = this.chromosome;
            chrLen = c.sequenceLen;
            firstBetaClampPos = [
                this.primaseBindingLocations{1}(1)-corFtpt3-bClmpFtpt
                this.primaseBindingLocations{2}(1)+corFtpt3+1]';
            
            %% bind lagging polymerase to lagging strand
            leadingPolPos  = this.leadingPolymerasePosition;
            laggingPolPos  = this.laggingPolymerasePosition;
            
            tfs = ...
                laggingPolPos == 0 & ...
                this.laggingBackupBetaClampPosition == firstBetaClampPos & ...
                c.complexBoundSites([leadingPolPos; this.leadingStrandIndexs]')' == leadPolGblIdx(1);
            
            if tfs(1); laggingPolPos(1) = firstBetaClampPos(1);           end;
            if tfs(2); laggingPolPos(2) = firstBetaClampPos(2) - corFtpt; end;
            
            n = sum(tfs);
            
            this.releaseProteinFromSites([leadingPolPos(:, tfs);     this.leadingStrandIndexs(:,tfs)]', false);
            this.releaseProteinFromSites([firstBetaClampPos(:, tfs); this.laggingStrandIndexs(:,tfs)]', false);
            
            this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) - n;
            this.enzymes(this.enzymeIndexs_coreBetaClampGammaComplex)         = this.enzymes(this.enzymeIndexs_coreBetaClampGammaComplex)  + n;
            this.enzymes(this.enzymeIndexs_core)                              = this.enzymes(this.enzymeIndexs_core)                       + n;
            this.enzymes(this.enzymeIndexs_primase)                           = this.enzymes(this.enzymeIndexs_primase)                    + n;
            
            this.enzymes(this.enzymeIndexs_betaClamp)                         = this.enzymes(this.enzymeIndexs_betaClamp)                  - n;
            this.enzymes(this.enzymeIndexs_betaClampMonomer)                  = this.enzymes(this.enzymeIndexs_betaClampMonomer)           + 2*n;
            
            this.enzymes(this.enzymeIndexs_betaClampMonomer)                  = this.enzymes(this.enzymeIndexs_betaClampMonomer)           - 2*n;
            this.enzymes(this.enzymeIndexs_core)                              = this.enzymes(this.enzymeIndexs_core)                       - n;
            this.enzymes(this.enzymeIndexs_primase)                           = this.enzymes(this.enzymeIndexs_primase)                    - n;
            this.enzymes(this.enzymeIndexs_coreBetaClampPrimase)              = this.enzymes(this.enzymeIndexs_coreBetaClampPrimase)       + n;
            
            if ~all(this.bindProteinToChromosome([leadingPolPos(:, tfs); this.leadingStrandIndexs(:,tfs)]', this.enzymeIndexs_coreBetaClampGammaComplex, [], [], true, false, 1, true, []))
                throw(MException('Replication:error', 'All polymerases should bind.'));
            end
            if ~all(this.bindProteinToChromosome([laggingPolPos(:, tfs); this.laggingStrandIndexs(:,tfs)]', this.enzymeIndexs_coreBetaClampPrimase, [], [], true, false, 1, true, []))
                throw(MException('Replication:error', 'All polymerases should bind.'));
            end
            
            %% Calculate maximum unwinding and polymerization extent
            helicasePos = this.helicasePosition;
            leadingPolPos  = this.leadingPolymerasePosition;
            laggingPolPos  = this.laggingPolymerasePosition;
            leadingPos  = this.leadingPosition;
            laggingPos  = this.laggingPosition;
            fPos = this.okazakiFragmentPosition;
            fPrg = this.okazakiFragmentProgress;
            fLen = this.okazakiFragmentLength;
            
            %initialize limits
            limits = zeros(2, 2); %[leading 1  leading 2; lagging 1  lagging 2]
            
            %primase kinetics
            limits(1, 1) = this.primerLength - (chrLen - leadingPos(1));
            limits(1, 2) = this.primerLength - (leadingPos(2) - 1);
            limits(2, 1) = this.primerLength - (laggingPos(1) - fPos(1));
            limits(2, 2) = this.primerLength - (fPos(2) - laggingPos(2));
            
            %DNA polymerase kinetics
            limits(limits <= 0) = this.dnaPolymeraseElongationRate;
            
            %polymerase must already be bound
            limits(1, :) = limits(1, :) .* (leadingPos ~= 0);
            limits(2, :) = limits(2, :) .* (laggingPos ~= 0);
            
            %prevent leading strand from progressing if no SSBs protecting lagging strand
            limits(1, :) = limits(1, :) .* this.areLaggingStrandSSBSitesBound;
            
            %helicase
            this.releaseProteinFromSites([helicasePos; this.leadingStrandIndexs]', false);
            [~, ~, ~, extents] = c.isRegionAccessible([helicasePos + [0 helFtpt-1]; this.leadingStrandIndexs]', [-1; 1] .* (limits(1, :) + helFtpt + 1)', ...
                [], helGblIdx, false, [], true, true);
            limits(1, :) = max(0, abs(extents)-helFtpt-1);
            
            %polymerase
            bndLeadPolIdx = ismembc2(c.complexBoundSites([leadingPolPos; this.leadingStrandIndexs]'), leadPolGblIdx);
            this.releaseProteinFromSites([leadingPolPos; this.leadingStrandIndexs]', false);
            [~, ~, ~, leadingExtents] = c.isRegionAccessible([leadingPolPos + [0 holFtpt-1]; this.leadingStrandIndexs]', [-1; 1] .* (limits(1, :) + holFtpt + 1)', ...
                [], leadPolGblIdx(1), false, [], true, true);
            limits(1, :) = max(0, abs(leadingExtents)-holFtpt-1);
            
            if laggingPos(1) ~= 0
                this.releaseProteinFromSites([laggingPolPos(1) this.laggingStrandIndexs(1)], false);
                [~, ~, ~, laggingExtents] = c.isRegionAccessible([laggingPolPos(1) this.laggingStrandIndexs(1)], limits(2, 1) + holFtpt + 1, ...
                    [], lagPolGblIdx, false, [], true, true);
                limits(2, 1) = max(0, abs(laggingExtents)-holFtpt-1);
            end
            if laggingPos(2) ~= 0
                this.releaseProteinFromSites([laggingPolPos(2) this.laggingStrandIndexs(2)], false);
                [~, ~, ~, laggingExtents] = c.isRegionAccessible([laggingPolPos(2)+holFtpt-1 this.laggingStrandIndexs(2)], -(limits(2, 2) + holFtpt + 1), ...
                    [], lagPolGblIdx, false, [], true, true);
                limits(2, 2) = max(0, abs(laggingExtents)-holFtpt-1);
            end
            
            %no polymerized of lagging strands past end of Okazaki fragment
            limits(2, :) = min(limits(2, :), fLen - fPrg);
            
            %prevent leading strand from getting too far ahead of lagging strand
            tmp = fPos;
            if fPos(1) == 0; tmp(1) = chrLen; end;
            if fPos(2) == 0; tmp(2) = 1; end;
            limits(1, 1) = limits(1, 1) * ((tmp(1) - helicasePos(1) - helFtpt) < 2 * this.okazakiFragmentMeanLength);
            limits(1, 2) = limits(1, 2) * ((helicasePos(2) - tmp(2)) < 2 * this.okazakiFragmentMeanLength);
            
            %no helicase progress if the double-stranded region about the terC
            %has non-zero linking number and would be annilated by unwinding
            if ...
                    helicasePos(1) + helFtpt5 >= c.terCPosition + 1 && ...
                    helicasePos(1) + helFtpt5 - limits(1, 1) < c.terCPosition+1 && ...
                    helicasePos(2) + helFtpt5 - 1 <= c.terCPosition && ...
                    helicasePos(2) + helFtpt5 - 1 + limits(1, 2) > c.terCPosition && ...
                    abs(c.linkingNumbers([min(c.terCPosition + 1, helicasePos(2) + helFtpt5 - 1) 1])) > 1e-6
                limits(1, ceil(2 * this.randStream.rand())) = 0;
            elseif ...
                    helicasePos(1) + helFtpt5 >= c.terCPosition + 1 && ...
                    helicasePos(1) + helFtpt5 - limits(1, 1) < c.terCPosition+1 && ...
                    (helicasePos(2) + helFtpt5 - 1 > c.terCPosition || helicasePos(2) + helFtpt5 - 1 + limits(1, 2) <= c.terCPosition) && ...
                    abs(c.linkingNumbers([min(c.terCPosition + 1, helicasePos(2) + helFtpt5 - 1) 1])) > 1e-6 && ...
                    ~c.isRegionPolymerized([c.terCPosition 2], 1, false)
                limits(1, 1) = 0;
            elseif ...
                    helicasePos(2) + helFtpt5 - 1 <= c.terCPosition && ...
                    helicasePos(2) + helFtpt5 - 1 + limits(1, 2) > c.terCPosition && ...
                    (helicasePos(1) + helFtpt5 < c.terCPosition + 1 || helicasePos(1) + helFtpt5 - limits(1, 1) >= c.terCPosition+1) && ...
                    abs(c.linkingNumbers([min(c.terCPosition + 1, helicasePos(2) + helFtpt5 - 1) 1])) > 1e-6 && ...
                    ~c.isRegionPolymerized([c.terCPosition + 1 2], 1, false)
                limits(1, 2) = 0;
            end
            
            %pause if collision with RNA polymerase
            [posStrands, vals] = find(c.complexBoundSites);
            tfs = ismembc(vals, this.complex.rnaPolymeraseIndexs);
            if any(tfs)
                pos = posStrands(tfs, 1);
                strnds = posStrands(tfs, 2);
                rnaPolFtpt = c.getDNAFootprint([], this.complex.rnaPolymeraseIndexs(1));
                
                %calculate limits based on RNA polymerase position
                tmpLimits = zeros(2, 2);
                tfs = pos < helicasePos(1)   & mod(strnds, 2) == 1; if any(tfs), tmpLimits(1, 1) = min(helicasePos(1) - (pos(tfs) + rnaPolFtpt));   end
                tfs = pos > helicasePos(2)   & mod(strnds, 2) == 0; if any(tfs), tmpLimits(1, 2) = min(pos(tfs) - (helicasePos(2) + helFtpt));      end
                tfs = pos > laggingPolPos(1) & mod(strnds, 2) == 0; if any(tfs), tmpLimits(2, 1) = min(pos(tfs) - (laggingPolPos(1) + holFtpt));    end
                tfs = pos < laggingPolPos(2) & mod(strnds, 2) == 1; if any(tfs), tmpLimits(2, 2) = min(laggingPolPos(2) - (pos(tfs) + rnaPolFtpt)); end
                 
                %apply limits to DNA polymerase which are stalled
                tfs = this.randStream.random('poisson', 1 / this.rnaPolymeraseCollisionMeanDwellTime, [2 2]) == 0;
                limits(tfs) = min(limits(tfs), max(0, tmpLimits(tfs)));
            end
            
            %no polymerization of leading strands past terC
            %(move leading polymerization slightly past terC to make room for
            %lagging polymerase)
            polLimits = limits;
            polLimits(1, 1) = min(limits(1, 1), max(0, leadingPos(1) - this.terCPosition));
            polLimits(1, 2) = min(limits(1, 2), max(0, this.terCPosition - leadingPos(2) + 1));

            limits(1, 1) = min(limits(1, 1), max(0, leadingPos(1) - (this.terCPosition-(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3))));
            limits(1, 2) = min(limits(1, 2), max(0, (this.terCPosition+1+(corFtpt5+bClmpFtpt+corFtpt5+helFtpt3)) - leadingPos(2)));
            
            %dNTPs
            sequences = char(zeros(4, max(polLimits(:))));
            sequences(1, 1:polLimits(1, 1)) = c.sequence(leadingPos(1) - (0:polLimits(1, 1)-1), 2)';
            sequences(2, 1:polLimits(1, 2)) = c.sequence(leadingPos(2) + (0:polLimits(1, 2)-1), 3)';
            sequences(3, 1:polLimits(2, 1)) = c.sequence(laggingPos(1) + (0:polLimits(2, 1)-1), 3)';
            sequences(4, 1:polLimits(2, 2)) = c.sequence(laggingPos(2) - (0:polLimits(2, 2)-1), 2)';
            
            pols = min(polLimits, reshape(edu.stanford.covert.cell.sim.util.polymerize(...
                sequences, this.substrates(this.substrateIndexs_dntp), 'ACGT', ' ', 0, 0, this.randStream), [], 2)');
            limits(pols ~= polLimits) = pols(pols ~= polLimits);
            polLimits = pols;
            
            %energy to unwind
            unwinds = min(polLimits(1, :), max(0, [
                (helicasePos(1) + helFtpt5) - this.terCPosition;
                this.terCPosition + 1 - (helicasePos(2)+helFtpt-helFtpt5-1)]'));
            unwindLimits = min(unwinds, floor(unwinds / sum(unwinds) * min(this.substrates([this.substrateIndexs_atp; this.substrateIndexs_water]))));
            
            polLimits(1, unwinds ~= unwindLimits) = unwindLimits(unwinds ~= unwindLimits);
            limits(1, unwinds ~= unwindLimits) = unwindLimits(unwinds ~= unwindLimits);
            
            %% Unwind DNA
            %advance helicase
            if ~all(this.bindProteinToChromosome([helicasePos; this.leadingStrandIndexs]' , this.enzymeIndexs_helicase, ...
                    [], [], true, false, [-1; 1] .* (limits(1,:) + 1)', true, []))
                throw(MException('Replication:error', 'helicases must bind'));
            end
            
            %unwind chromosome
            c.setRegionUnwound(helicasePos(1)+helFtpt5,           -unwindLimits(1, 1));
            c.setRegionUnwound(helicasePos(2)+helFtpt-helFtpt5-1,  unwindLimits(1, 2));
            
            %update metabolites
            n = sum(unwindLimits(1, :));
            this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - n;
            this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - n;
            this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + n;
            this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + n;
            this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + n;
            
            %% Polymerize DNA
            %advance DNA polymerase
            if ~this.bindProteinToChromosome([leadingPolPos(1) this.leadingStrandIndexs(1)], leadPolLclIdx(bndLeadPolIdx(1)), [], [], true, false, -(limits(1, 1)+1), true, [])
                throw(MException('Replication:error', 'polymerases must bind'));
            end
            if ~this.bindProteinToChromosome([leadingPolPos(2) this.leadingStrandIndexs(2)], leadPolLclIdx(bndLeadPolIdx(2)), [], [], true, false,  (limits(1, 2)+1), true, [])
                throw(MException('Replication:error', 'polymerases must bind'));
            end
            if laggingPos(1) ~= 0
                if ~this.bindProteinToChromosome([laggingPolPos(1) this.laggingStrandIndexs(1)], this.enzymeIndexs_coreBetaClampPrimase, [], [], true, false, limits(2, 1)+1, true, [])
                    throw(MException('Replication:error', 'polymerases must bind'));
                end
            end
            if laggingPos(2) ~= 0
                if ~this.bindProteinToChromosome([laggingPolPos(2) this.laggingStrandIndexs(2)], this.enzymeIndexs_coreBetaClampPrimase, [], [], true, false, -(limits(2, 2)+1), true, [])
                    throw(MException('Replication:error', 'polymerases must bind'));
                end
            end
            
            %set regions polymerized
            c.setRegionPolymerized([leadingPos; 1 2]', [-1; 1] .* polLimits(1, :)');
            c.setRegionPolymerized([laggingPos; 2 1]', [1; -1] .* polLimits(2, :)');
            
            %update metabolites
            usedDNTPs = ...
                c.sequence.subsequenceBaseCounts(leadingPos(1) - (0:polLimits(1, 1)-1), 2) + ...
                c.sequence.subsequenceBaseCounts(leadingPos(2) + (0:polLimits(1, 2)-1), 3) + ...
                c.sequence.subsequenceBaseCounts(laggingPos(1) + (0:polLimits(2, 1)-1), 3) + ...
                c.sequence.subsequenceBaseCounts(laggingPos(2) - (0:polLimits(2, 2)-1), 2);
            
            this.substrates(this.substrateIndexs_dntp)        = this.substrates(this.substrateIndexs_dntp)        - usedDNTPs;
            this.substrates(this.substrateIndexs_diphosphate) = this.substrates(this.substrateIndexs_diphosphate) + sum(usedDNTPs);
        end
        
        %Dissociate free SSB 8mers into 2 SSB 4mers
        function dissociateFreeSSBComplexes(this)
            this.enzymes(this.enzymeIndexs_ssb4mer) = ...
                this.enzymes(this.enzymeIndexs_ssb4mer) + ...
                2 * this.enzymes(this.enzymeIndexs_ssb8mer);
            this.enzymes(this.enzymeIndexs_ssb8mer) = 0;
        end
        
        %Free, Bind SSBs 8mers to single-stranded DNA
        function freeAndBindSSBs(this)
            %% free SSBs
            positionsStrands = this.releaseProteinFromChromosome(this.enzymeIndexs_ssb8mer, this.ssbDissociationRate, [], []);
            
            this.enzymes(this.enzymeIndexs_ssb4mer) = ...
                + this.enzymes(this.enzymeIndexs_ssb4mer) ...
                + 2 * size(positionsStrands, 1);
            this.enzymes(this.enzymeIndexs_ssb8mer) = ...
                + this.enzymes(this.enzymeIndexs_ssb8mer) ...
                - size(positionsStrands, 1);
            
            %% bind SSBs
            nPossibleSSB8mers = this.enzymes(this.enzymeIndexs_ssb4mer) / 2;
            if nPossibleSSB8mers < 1
                return;
            end
            
            c = this.chromosome;
            
            %find accessible regions
            [rgnPosStrnds, rgnLens] = c.getAccessibleRegions([], this.enzymeGlobalIndexs(this.enzymeIndexs_ssb8mer));
            if isempty(rgnPosStrnds)
                return;
            end
            
            %SSB 8mer footprint and spacing
            ftpt = this.enzymeDNAFootprints(this.enzymeIndexs_ssb8mer);
            spcg = this.ssbComplexSpacing;
            
            %exclude regions near other SSB 8mers
            ssbPosStrnds = find(c.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_ssb8mer));
            [rgnPosStrnds, rgnLens] = c.excludeRegions(rgnPosStrnds, rgnLens, ...
                [ssbPosStrnds(:, 1) - this.ssbComplexSpacing ssbPosStrnds(:, 2)], ftpt + 2*spcg);
            if isempty(rgnPosStrnds)
                return;
            end
            
            %split regions over terC
            idx = reshape(find(rgnPosStrnds(:, 1) <= this.terCPosition & rgnPosStrnds(:, 1) + rgnLens - 1 > this.terCPosition), [], 1);
            rgnPosStrnds = [rgnPosStrnds; repmat(this.terCPosition + 1, size(idx)) rgnPosStrnds(idx, 2)];
            rgnLens = [rgnLens; rgnLens(idx) - (this.terCPosition - rgnPosStrnds(idx, 1) + 1)];
            rgnLens(idx) = this.terCPosition - rgnPosStrnds(idx, 1) + 1;
            
            %switch directions of region on far side of terC
            tfs = rgnPosStrnds(:, 1) > this.terCPosition;
            rgnPosStrnds(tfs, 1) = rgnPosStrnds(tfs, 1) + rgnLens(tfs, 1) - ftpt;
            rgnLens(tfs) = -rgnLens(tfs);
            
            %split regions into pieces large enough to bind just 1 SSB 8mer
            tfs1 = abs(rgnLens) == ftpt;
            tfs2 = abs(rgnLens) > ftpt;
            rgnPosStrnds = [
                c.splitRegions(rgnPosStrnds(tfs2, :), rgnLens(tfs2), ftpt + spcg)
                rgnPosStrnds(tfs1, :)];
            if isempty(rgnPosStrnds)
                return;
            end
            
            %bind SSB 8mers to regions
            nBindings = min([
                size(rgnPosStrnds, 1);
                floor(nPossibleSSB8mers)]);
            
            this.enzymes(this.enzymeIndexs_ssb4mer) = this.enzymes(this.enzymeIndexs_ssb4mer) - 2*nBindings;
            this.enzymes(this.enzymeIndexs_ssb8mer) = this.enzymes(this.enzymeIndexs_ssb8mer) +   nBindings;
            
            if nBindings ~= sum(this.bindProteinToChromosome(rgnPosStrnds, this.enzymeIndexs_ssb8mer, nBindings, [], [], false, 1, false, []))
                throw(MException('Replication:error', 'SSBs must bind'));
            end
        end
        
        %Bind beta-clamp to start of next Okazaki fragment
        function initiateOkazakiFragment(this)
            %Terminate early if no active replisomes
            if ~this.isAnyHelicaseBound || ~all(this.leadingStrandElongating)
                return;
            end
            
            helFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_helicase);
            helFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_helicase);
            corFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_core);
            bClmpFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_betaClamp);
            
            fIdx = this.okazakiFragmentIndex;
            fPrg = this.okazakiFragmentProgress;
            helPos = this.helicasePosition;
            bClmpPos = this.laggingBackupBetaClampPosition;
            strndPol = this.strandPolymerized;
            
            %Compute list of starting points of next Okazaki fragments that
            %aren't yet bound by core/beta-clamps
            posStrnds = zeros(0, 2);
            
            if ...
                    fIdx(1) < numel(this.primaseBindingLocations{1}) && ...
                    (fIdx(1) == 0  || (fIdx(1) > 0 && fPrg(1) >= this.laggingBackupClampReloadingLength)) && ...
                    helPos(1) + helFtpt5 + 1 < this.primaseBindingLocations{1}(fIdx(1)+1)-corFtpt3-bClmpFtpt && ...
                    bClmpPos(1) ~= this.primaseBindingLocations{1}(fIdx(1)+1)-corFtpt3-bClmpFtpt && ...
                    ~strndPol(1)
                
                posStrnds = [posStrnds; this.primaseBindingLocations{1}(fIdx(1)+1)-corFtpt3-bClmpFtpt this.laggingStrandIndexs(1)];
            end
            
            if ...
                    fIdx(2) < numel(this.primaseBindingLocations{2}) && ...
                    (fIdx(2) == 0 || (fIdx(2) > 0 && fPrg(2) >= this.laggingBackupClampReloadingLength)) && ...
                    helPos(2)+helFtpt3 > this.primaseBindingLocations{2}(fIdx(2)+1)+corFtpt3+1+bClmpFtpt && ...
                    bClmpPos(2) ~= this.primaseBindingLocations{2}(fIdx(2)+1)+corFtpt3+1 && ...
                    ~strndPol(2)
                
                posStrnds = [posStrnds; this.primaseBindingLocations{2}(fIdx(2)+1)+corFtpt3+1 this.laggingStrandIndexs(2)];
            end
            
            %Bind beta-clamp to starting point of next Okazaki fragments
            nBinding = floor(min([
                size(posStrnds, 1)
                this.substrates(this.substrateIndexs_atp)
                this.substrates(this.substrateIndexs_water)
                this.enzymes(this.enzymeIndexs_betaClampMonomer)/2]));
            if nBinding ~= sum(this.bindProteinToChromosome(posStrnds, this.enzymeIndexs_betaClamp, nBinding, [], true, false, 1, false, []))
                throw(MException('Replication:error', 'Beta clamps must bind'));
            end
            
            this.enzymes(this.enzymeIndexs_betaClampMonomer) = this.enzymes(this.enzymeIndexs_betaClampMonomer) - 2 * nBinding;
            this.enzymes(this.enzymeIndexs_betaClamp)        = this.enzymes(this.enzymeIndexs_betaClamp)        +     nBinding;
            
            %hydrolyze ATP used by gamma-complex to form beta-clamp
            this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - nBinding;
            this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - nBinding;
            this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + nBinding;
            this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + nBinding;
            this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + nBinding;
        end
        
        %Terminate Okazaki fragments which have been completely synthesized
        %- release beta-clamp
        %- binding core to beta-clamp of next Okazaki fragment
        function terminateOkazakiFragment(this)
            %Terminate early if no active replisomes
            if ~this.isAnyHelicaseBound || ~all(this.leadingStrandElongating)
                return;
            end
            
            c = this.chromosome;
            chrLen = c.sequenceLen;
            fIdx = this.okazakiFragmentIndex;
            fPrg = this.okazakiFragmentProgress;
            fLen = this.okazakiFragmentLength;
            helPos = this.helicasePosition;
            helFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_helicase);
            holFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampPrimase);
            corFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_core);
            corFtpt3 = this.enzymeDNAFootprints3Prime(this.enzymeIndexs_core);
            corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
            areLaggingStrandSSBSitesBound = this.areLaggingStrandSSBSitesBound;
            laggingBackupBetaClampPosition = this.laggingBackupBetaClampPosition;
            leadingPolymerasePosition = this.leadingPolymerasePosition;
            laggingPolymerasePosition = this.laggingPolymerasePosition;
            
            %calculate positions of core/beta-clamps that have reached end of
            %Okazaki fragments and for which a core/beta-clamp is bound to the
            %next Okazaki fragment
            releasePolPosStrnds = zeros(0, 2);
            releaseBClmpPosStrnds = zeros(0, 2);
            bindPolPosStrnds = zeros(0, 2);
            terminatedPosStrnds = zeros(0, 2);
            if ...
                    fIdx(1) > 0 && ...
                    fPrg(1) == fLen(1) && ...
                    areLaggingStrandSSBSitesBound(1) && ...
                    (fIdx(1) == numel(this.primaseBindingLocations{1}) || (...
                    (fIdx(1) == numel(this.primaseBindingLocations{1}) - 1 || ...
                    this.primaseBindingLocations{1}(fIdx(1)+1)-(helPos(1)+helFtpt) > this.startingOkazakiLoopLength) && ...
                    laggingBackupBetaClampPosition(1) == this.primaseBindingLocations{1}(fIdx(1)+1)-(holFtpt-corFtpt5)+1))
                
                releasePolPosStrnds = [releasePolPosStrnds;
                    laggingPolymerasePosition(1) this.laggingStrandIndexs(1)];
                if fIdx(1) < numel(this.primaseBindingLocations{1})
                    releaseBClmpPosStrnds = [releaseBClmpPosStrnds;
                        laggingBackupBetaClampPosition(1) this.laggingStrandIndexs(1)];
                    bindPolPosStrnds = [bindPolPosStrnds;
                        laggingBackupBetaClampPosition(1) this.laggingStrandIndexs(1)];
                else
                    terminatedPosStrnds = [terminatedPosStrnds;
                        leadingPolymerasePosition(1) this.leadingStrandIndexs(1)];
                end
                
                if fIdx(1) == 1
                    c.strandBreaks(chrLen, 3) = 1;
                else
                    c.strandBreaks(this.primaseBindingLocations{1}(fIdx(1)-1)-1, 3) = 1;
                end
            end
            if ...
                    fIdx(2) > 0 && ...
                    fPrg(2) == fLen(2) && ...
                    areLaggingStrandSSBSitesBound(2) && ...
                    (fIdx(2) == numel(this.primaseBindingLocations{2}) || (...
                    (fIdx(2) == numel(this.primaseBindingLocations{2})-1 || ...
                    helPos(2) - this.primaseBindingLocations{2}(fIdx(2)+1) > this.startingOkazakiLoopLength) && ...
                    laggingBackupBetaClampPosition(2) == this.primaseBindingLocations{2}(fIdx(2)+1)+corFtpt3+1))
                
                releasePolPosStrnds = [releasePolPosStrnds;
                    laggingPolymerasePosition(2) this.laggingStrandIndexs(2)];
                if fIdx(2) < numel(this.primaseBindingLocations{2})
                    releaseBClmpPosStrnds = [releaseBClmpPosStrnds;
                        laggingBackupBetaClampPosition(2) this.laggingStrandIndexs(2)];
                    bindPolPosStrnds = [bindPolPosStrnds;
                        laggingBackupBetaClampPosition(2)-corFtpt this.laggingStrandIndexs(2)];
                else
                    terminatedPosStrnds = [terminatedPosStrnds;
                        leadingPolymerasePosition(2) this.leadingStrandIndexs(2)];
                end
                
                if fIdx(2) == 1
                    c.strandBreaks(chrLen, 2) = 1;
                else
                    c.strandBreaks(this.primaseBindingLocations{2}(fIdx(2)-1), 2) = 1;
                end
            end
            
            %release the core/beta-clamp that have reached the end of the Okazki
            %fragment
            this.releaseProteinFromSites(releasePolPosStrnds, false);
            this.enzymes(this.enzymeIndexs_coreBetaClampPrimase) = this.enzymes(this.enzymeIndexs_coreBetaClampPrimase)    -   size(releasePolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_primase)          = this.enzymes(this.enzymeIndexs_primase)          +   size(releasePolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_core)             = this.enzymes(this.enzymeIndexs_core)             +   size(releasePolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_betaClampMonomer) = this.enzymes(this.enzymeIndexs_betaClampMonomer) + 2*size(releasePolPosStrnds, 1);
            
            this.releaseProteinFromSites(releaseBClmpPosStrnds, false);
            this.enzymes(this.enzymeIndexs_betaClamp)        = this.enzymes(this.enzymeIndexs_betaClamp)        -   size(releaseBClmpPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_betaClampMonomer) = this.enzymes(this.enzymeIndexs_betaClampMonomer) + 2*size(releaseBClmpPosStrnds, 1);
            
            if ~all(this.bindProteinToChromosome(bindPolPosStrnds, this.enzymeIndexs_coreBetaClampPrimase, [], [], true, false, 1, false, []))
                throw(MException('Replication:error', 'lagging strand polymerases must bind'));
            end
            this.enzymes(this.enzymeIndexs_primase)              = this.enzymes(this.enzymeIndexs_primase)              -   size(bindPolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_core)                 = this.enzymes(this.enzymeIndexs_core)                 -   size(bindPolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_betaClampMonomer)     = this.enzymes(this.enzymeIndexs_betaClampMonomer)     - 2*size(bindPolPosStrnds, 1);
            this.enzymes(this.enzymeIndexs_coreBetaClampPrimase) = this.enzymes(this.enzymeIndexs_coreBetaClampPrimase) +   size(bindPolPosStrnds, 1);
            
            if ~isempty(terminatedPosStrnds)
                n = size(terminatedPosStrnds, 1);
                this.modifyProteinOnChromosome(terminatedPosStrnds, this.enzymeIndexs_2coreBetaClampGammaComplexPrimase(ones(n, 1)));
                
                this.enzymes(this.enzymeIndexs_coreBetaClampGammaComplex)    = this.enzymes(this.enzymeIndexs_coreBetaClampGammaComplex)  - n;
                this.enzymes(this.enzymeIndexs_core)                         = this.enzymes(this.enzymeIndexs_core)                       + n;
                this.enzymes(this.enzymeIndexs_betaClampMonomer)             = this.enzymes(this.enzymeIndexs_betaClampMonomer)           + 2*n;
                this.enzymes(this.enzymeIndexs_gammaComplex)                 = this.enzymes(this.enzymeIndexs_gammaComplex)               + n;
                
                this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase)   = this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) + n;
                this.enzymes(this.enzymeIndexs_primase)                             = this.enzymes(this.enzymeIndexs_primase)                    - n;
                this.enzymes(this.enzymeIndexs_core)                                = this.enzymes(this.enzymeIndexs_core)                       - 2*n;
                this.enzymes(this.enzymeIndexs_betaClampMonomer)                    = this.enzymes(this.enzymeIndexs_betaClampMonomer)           - 2*n;
                this.enzymes(this.enzymeIndexs_gammaComplex)                        = this.enzymes(this.enzymeIndexs_gammaComplex)               - n;
            end
        end
        
        %Ligate DNA
        %(2) DR5P + NAD ==> AMP + dRibose5P_dRibose5P + H + NMN
        %
        %Note: this function is redundant with evolveState_Ligate in DNARepair
        function ligateDNA(this)
            c = this.chromosome;
            
            positionsStrands = find(c.singleStrandBreaks);
            if isempty(positionsStrands)
                return;
            end
                
            numReactions = max(0, floor(min([
                size(positionsStrands, 1);
                this.randStream.stochasticRound(this.enzymes(this.enzymeIndexs_ligase) * this.stepSizeSec * this.ligaseRate);
                this.substrates(this.substrateIndexs_nad)])));
            if numReactions == 0
                return;
            end
            
            [~, ~, positionsStrands, numReactions] = this.bindProteinToChromosome(positionsStrands, ...
                this.enzymeIndexs_ligase, numReactions, [], false, true, 1, false, [1 2]);
            if numReactions == 0
                return;
            end
            
            c.strandBreaks(positionsStrands) = 0;
            
            this.substrates(this.substrateIndexs_nad)      = this.substrates(this.substrateIndexs_nad)      - numReactions;
            this.substrates(this.substrateIndexs_nmn)      = this.substrates(this.substrateIndexs_nmn)      + numReactions;
            this.substrates(this.substrateIndexs_amp)      = this.substrates(this.substrateIndexs_amp)      + numReactions;
            this.substrates(this.substrateIndexs_hydrogen) = this.substrates(this.substrateIndexs_hydrogen) + numReactions;
        end
        
        %Terminate replication if replisomes have reached terC (meaning the
        %chromosome has been duplicated)
        %- dissolve replisomes
        function terminateReplication(this)
            if ...
                    ~all(this.leadingStrandElongating) || ...  %Terminate early if no active replisomes
                    ~all(this.strandPolymerized)       || ...  %check if replisomes have reached terC
                    any(this.laggingStrandElongating)          %check if last Okazaki fragment already terminated
                return;
            end
            
            %set strand break
            this.chromosome.strandBreaks(this.terCPosition, 2:3) = 1;
            
            %dissolve replisomes
            if size(this.releaseProteinFromChromosome(this.enzymeIndexs_helicase, Inf, [], []), 1) ~= 2
                throw(MException('Replication:error', 'Two helicases must be released'));
            end
            if size(this.releaseProteinFromChromosome(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase, Inf, [], []), 1) ~= 2
                throw(MException('Replication:error', 'Two polymerase complexes must be released'));
            end
            
            this.enzymes = this.enzymes + 2 * this.enzymeComposition(:, this.enzymeIndexs_2coreBetaClampGammaComplexPrimase);
            this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) = this.enzymes(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase) - 2;
        end
    end
    
    %get methods of dependent local state
    methods
        function value = get.isDnaAORIComplexAssembled(this)
            c = this.chromosome;
            [posStrnds, complexs] = find(c.complexBoundSites);
            
            %check if any DnaA-ATP 7mers bound
            idxs7 = find(complexs == this.complexIndexs_DnaA_7mer_ATP);
            if length(idxs7) < 4
                value = false;
                return;
            end
            
            %check if DnaA-ATP 1mer bound at R5 box on chromosome 1
            idxs1 = find(posStrnds(:, 1) == this.dnaAFunctionalBoxStartPositions(this.dnaAFunctionalBoxIndexs_R5));
            if ~any(complexs(idxs1) == this.complexIndexs_DnaA_1mer_ATP & posStrnds(idxs1, 2) == 1)
                value = false;
                return;
            end
            
            %check if DnaA-ATP 7mers bound at R1-4 boxes on chromosome 1
            value = all(ismembc(this.dnaAFunctionalBoxStartPositions(this.dnaAFunctionalBoxIndexs_R1234, 1), posStrnds(idxs7(posStrnds(idxs7, 2) == 1), 1)));
        end
        
        function result = get.isAnyHelicaseBound(this)
            [~, complexs] = find(this.chromosome.complexBoundSites);
            result = any(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_helicase));
        end
        
        function result = get.isAnyPolymeraseBound(this)
            [~, complexs] = find(this.chromosome.complexBoundSites);
            result = ...
                any(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_2coreBetaClampGammaComplexPrimase)) || ...
                any(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_coreBetaClampGammaComplex)) || ...
                any(complexs == this.enzymeGlobalIndexs(this.enzymeIndexs_coreBetaClampPrimase));
        end
        
        function result = get.leadingStrandElongating(this)
            result = ...
                this.helicasePosition & ...
                this.leadingPolymerasePosition;
        end
        
        function result = get.laggingStrandElongating(this)
            result = this.laggingPolymerasePosition ~= 0;
        end
        
        function result = get.leadingStrandPolymerized(this)
            c = this.chromosome;
            result = false(1, 2);
            
            [pos2, lens2] = find(c.polymerizedRegions(:, 2));
            [pos4, lens4] = find(c.polymerizedRegions(:, 4));
            result(1) = ...
                any(pos2(:,1) <= this.terCPosition + 1 & pos2(:,1) + lens2 - 1 == c.sequenceLen) && ...
                any(pos4(:,1) <= this.terCPosition + 1 & pos4(:,1) + lens4 - 1 == c.sequenceLen);
            
            result(2) = ...
                all(c.polymerizedRegions([1 1; 1 3]) >= this.terCPosition);
        end
        
        function result = get.laggingStrandPolymerized(this)
            c = this.chromosome;
            result = false(1, 2);
            
            [pos1, lens1] = find(c.polymerizedRegions(:, 1));
            [pos3, lens3] = find(c.polymerizedRegions(:, 3));
            result(1) = ...
                any(pos1(:,1) <= this.terCPosition + 1 & pos1(:,1) + lens1 - 1 == c.sequenceLen) && ...
                any(pos3(:,1) <= this.terCPosition + 1 & pos3(:,1) + lens3 - 1 == c.sequenceLen);
            
            result(2) = ...
                all(c.polymerizedRegions([1 2; 1 4]) >= this.terCPosition);
        end
        
        function result = get.strandPolymerized(this)
            result = ...
                this.leadingStrandPolymerized & ...
                this.laggingStrandPolymerized;
        end
        
        function result = get.numLigations(this)
            fIdx = this.okazakiFragmentIndex;
            fPrg = this.okazakiFragmentProgress;
            
            pos1 = [];
            pos2 = [];
            
            if fIdx(1) > 0; pos1 = this.primaseBindingLocations{1}(1:fIdx(1)-2)-1; end;
            if fIdx(2) > 0; pos2 = this.primaseBindingLocations{2}(1:fIdx(2)-2)  ; end;
            
            if fIdx(1) > 1
                pos1 = [pos1; 0];
            end
            if fIdx(2) > 1
                pos2 = [pos2; 0];
            end
            if fIdx(1) == numel(this.primaseBindingLocations{1}) && fPrg(1) > 0
                pos1 = [pos1; this.terCPosition];
            end
            if fIdx(2) == numel(this.primaseBindingLocations{2}) && fPrg(2) > 0
                pos2 = [pos2; this.terCPosition];
            end
            
            strandPolymerized = this.strandPolymerized;
            leadingPosition = this.leadingPosition;
            if strandPolymerized(1)
                pos1 = [0; this.primaseBindingLocations{1}(1:end)-1];
            end
            if strandPolymerized(2) && (strandPolymerized(1) || any(leadingPosition))
                pos2 = [0; this.primaseBindingLocations{2}(1:end)];
            end
            
            result = zeros(1, 2);
            result(1) = numel(pos1) - nnz(this.chromosome.strandBreaks(pos1, 3));
            result(2) = numel(pos2) - nnz(this.chromosome.strandBreaks(pos2, 2));
        end
        
        %Every Okazaki fragment needs to have been ligated on one end, and
        %the last one on each strand needs to have been ligated on both
        %ends.
        function result = get.strandLigated(this)
            result = this.numLigations == ...
                [numel(this.primaseBindingLocations{1}) + 1 ...
                numel(this.primaseBindingLocations{2}) + 1];
        end
        
        function result = get.strandDuplicated(this)
            result = ...
                this.strandPolymerized & ...
                this.strandLigated & ...
                this.helicasePosition == 0 & ...
                this.leadingPosition == 0 & ...
                this.laggingPosition == 0 &  ...
                this.laggingBackupBetaClampPosition == 0;
        end
        
        %helicase position
        function result = get.helicasePosition(this)
            positionsStrands = find(this.chromosome.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_helicase));
            
            idxs1 = find(positionsStrands(:, 2) == this.leadingStrandIndexs(1));
            idxs2 = find(positionsStrands(:, 2) == this.leadingStrandIndexs(2));
            
            if numel(idxs1) > 1 || numel(idxs2) > 1
                throw(MException('Replication:error', 'More than 1 active replisome'));
            end
            
            result = zeros(1, 2);
            if ~isempty(idxs1)
                result(1) = positionsStrands(idxs1, 1);
            end
            if ~isempty(idxs2)
                result(2) = positionsStrands(idxs2, 1);
            end
        end
        
        function result = get.leadingPolymerasePosition(this)
            polLclIdxs = [this.enzymeIndexs_coreBetaClampGammaComplex; this.enzymeIndexs_2coreBetaClampGammaComplexPrimase];
            [positionsStrands, complexs] = find(this.chromosome.complexBoundSites);
            tfs = complexs == this.enzymeGlobalIndexs(polLclIdxs(1)) | complexs == this.enzymeGlobalIndexs(polLclIdxs(2));
            positionsStrands = positionsStrands(tfs, :);
            
            result = zeros(1, 2);
            
            idx = find(positionsStrands(:, 2) == this.leadingStrandIndexs(1));
            if numel(idx) > 1
                throw(MException('Replication:error', 'More than 1 active replisome'));
            elseif numel(idx) == 1
                result(1) = positionsStrands(idx, 1);
            end
            
            idx = find(positionsStrands(:, 2) == this.leadingStrandIndexs(2));
            if numel(idx) > 1
                throw(MException('Replication:error', 'More than 1 active replisome'));
            elseif numel(idx) == 1
                result(2) = positionsStrands(idx, 1);
            end
        end
        
        function result = get.laggingPolymerasePosition(this)
            positionsStrands = find(this.chromosome.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_coreBetaClampPrimase));
            idxs1 = find(positionsStrands(:, 2) == this.laggingStrandIndexs(1));
            idxs2 = find(positionsStrands(:, 2) == this.laggingStrandIndexs(2));
            
            result = zeros(1, 2);
            if ~isempty(idxs1); result(1) = max(positionsStrands(idxs1, 1)); end;
            if ~isempty(idxs2); result(2) = min(positionsStrands(idxs2, 1)); end
            
            if numel(idxs1) > 2 || numel(idxs2) > 2
                throw(MException('Replication:error', 'More than 1 active replisome'));
            end
        end
        
        function result = get.leadingPosition(this)
            holFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
            
            result  = this.leadingPolymerasePosition;
            if result(1) ~= 0; result(1) = result(1) + corFtpt5;               end;
            if result(2) ~= 0; result(2) = result(2) + holFtpt - corFtpt5 - 1; end;
            
            result(result ~= 0) = mod(result(result ~= 0) - 1, this.chromosome.sequenceLen) + 1;
        end
        
        function result = get.laggingPosition(this)
            holFtpt  = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampGammaComplex);
            corFtpt5 = this.enzymeDNAFootprints5Prime(this.enzymeIndexs_core);
            
            result = this.laggingPolymerasePosition;
            if result(1) ~= 0; result(1) = result(1) + holFtpt - corFtpt5 - 1; end;
            if result(2) ~= 0; result(2) = result(2) + corFtpt5;               end;
            
            result(result ~= 0) = mod(result(result ~= 0) - 1, this.chromosome.sequenceLen) + 1;
        end
        
        function result = get.laggingBackupBetaClampPosition(this)
            positionsStrands = find(this.chromosome.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_betaClamp));
            idxs1 = find(positionsStrands(:, 2) == this.laggingStrandIndexs(1));
            idxs2 = find(positionsStrands(:, 2) == this.laggingStrandIndexs(2));
            
            result = zeros(1, 2);
            if ~isempty(idxs1); result(1) = max(positionsStrands(idxs1, 1)); end;
            if ~isempty(idxs2); result(2) = min(positionsStrands(idxs2, 1)); end
            
            if numel(idxs1) > 1 || numel(idxs2) > 1
                throw(MException('Replication:error', 'More than 1 active replisome'));
            end
        end
        
        function result = get.okazakiFragmentIndex(this)
            c = this.chromosome;
            laggingPos = this.laggingPosition;
            
            adjLaggingPos = laggingPos;
            if laggingPos(1) ~= 0 && laggingPos(1) < 1/2 * this.terCPosition
                adjLaggingPos(1) = laggingPos(1) + this.chromosome.sequenceLen;
            end
            if laggingPos(2) ~= 0 && laggingPos(2) > 3/2 * this.terCPosition
                adjLaggingPos(2) = laggingPos(2) - this.chromosome.sequenceLen;
            end
            
            result = zeros(1, 2);
            
            idx1 = find(this.primaseBindingLocations{1} <= adjLaggingPos(1) - (c.isRegionPolymerized([adjLaggingPos(1)-1 3], 1, false) && adjLaggingPos(1) > this.terCPosition+1), 1, 'first');
            idx2 = find(this.primaseBindingLocations{2} >= adjLaggingPos(2) + (c.isRegionPolymerized([adjLaggingPos(2)+1 2], 1, false) && adjLaggingPos(2) < this.terCPosition),   1, 'first');
            if laggingPos(1) ~= 0 && ~isempty(idx1); result(1) = idx1; end;
            if laggingPos(2) ~= 0 && ~isempty(idx2); result(2) = idx2; end;
        end
        
        function result = get.okazakiFragmentPosition(this)
            okazakiFragmentIndex = this.okazakiFragmentIndex;
            result = zeros(1, 2);
            if okazakiFragmentIndex(1) > 0
                result(1) = this.primaseBindingLocations{1}(okazakiFragmentIndex(1));
            end
            if okazakiFragmentIndex(2) > 0
                result(2) = this.primaseBindingLocations{2}(okazakiFragmentIndex(2));
            end
        end
        
        function result = get.okazakiFragmentLength(this)
            okazakiFragmentIndex = this.okazakiFragmentIndex;
            
            starts = zeros(1, 2);
            if okazakiFragmentIndex(1) > 0
                starts(1) = this.primaseBindingLocations{1}(okazakiFragmentIndex(1));
            end
            if okazakiFragmentIndex(2) > 0
                starts(2) = this.primaseBindingLocations{2}(okazakiFragmentIndex(2));
            end
            
            ends = zeros(1, 2);
            if okazakiFragmentIndex(1) > 1
                ends(1) = this.primaseBindingLocations{1}(okazakiFragmentIndex(1)-1) - 1;
            elseif okazakiFragmentIndex(1) == 1
                ends(1) = this.chromosome.sequenceLen;
            end
            if okazakiFragmentIndex(2) > 1
                ends(2) = this.primaseBindingLocations{2}(okazakiFragmentIndex(2)-1) + 1;
            elseif okazakiFragmentIndex(2) == 1
                ends(2) = 1;
            end
            
            result = abs(ends-starts)+1;
            result(okazakiFragmentIndex == 0) = 0;
        end
        
        function result = get.okazakiFragmentProgress(this)
            laggingPos  = this.laggingPosition;
            fIdx = this.okazakiFragmentIndex;
            
            adjLaggingPos = laggingPos;
            if laggingPos(1) ~= 0 && laggingPos(1) < 1/2* this.terCPosition
                adjLaggingPos(1) = laggingPos(1) + this.chromosome.sequenceLen;
            end
            if laggingPos(2) ~= 0 && laggingPos(2) > 3/2 * this.terCPosition
                adjLaggingPos(2) = laggingPos(2) - this.chromosome.sequenceLen;
            end
            
            result = zeros(1, 2);
            if laggingPos(1) ~= 0 && fIdx(1) ~= 0; result(1) = adjLaggingPos(1) - this.primaseBindingLocations{1}(fIdx(1)); end;
            if laggingPos(2) ~= 0 && fIdx(2) ~= 0; result(2) = this.primaseBindingLocations{2}(fIdx(2)) - adjLaggingPos(2); end;
        end
        
        function result = get.leadingStrandBoundSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            helPos = this.helicasePosition;
            polPos = this.leadingPolymerasePosition;
            
            posStrnds = find(this.chromosome.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_ssb8mer));
            pos = posStrnds(:, 1);
            str = posStrnds(:, 2);
            
            idx1 = [];
            idx2 = [];
            
            if polPos(1) ~= 0
                idx1 = find(...
                    str == 1 & ...
                    pos < polPos(1) & pos > helPos(1));
            end
            if polPos(2) ~= 0
                idx2 = find(...
                    str == 4 & ...
                    pos > polPos(2) & pos < helPos(2));
            end
            
            result = CircularSparseMat([pos(idx1) ones(size(idx1)); pos(idx2) repmat(2, size(idx2))], ones(numel(idx1) + numel(idx2), 1), [this.chromosome.sequenceLen, 2], 1);
        end
        
        function result = get.laggingStrandBoundSSBs(this)
            import edu.stanford.covert.util.CircularSparseMat;
            
            helPos = this.helicasePosition;
            
            posStrnds = find(this.chromosome.complexBoundSites == this.enzymeGlobalIndexs(this.enzymeIndexs_ssb8mer));
            pos = posStrnds(:, 1);
            str = posStrnds(:, 2);
            
            fIdx = this.okazakiFragmentIndex;
            
            starts = zeros(1, 2);
            if fIdx(1) == 0
                starts(1) = this.chromosome.sequenceLen;
            else
                starts(1) = this.primaseBindingLocations{1}(fIdx(1));
            end
            if fIdx(2) == 0
                starts(2) = 1;
            else
                starts(2) = this.primaseBindingLocations{2}(fIdx(2));
            end
            
            idx1 = find(...
                str == 4 & ...
                pos < starts(1) & ...
                pos > helPos(1));
            idx2 = find(...
                str == 1 & ...
                pos > starts(2) & ...
                pos < helPos(2));
            
            result = CircularSparseMat([pos(idx1) ones(size(idx1)); pos(idx2) repmat(2, size(idx2))], ones(numel(idx1) + numel(idx2), 1), [this.chromosome.sequenceLen, 2], 1);
        end
        
        function result = get.numLeadingTemplateBoundSSBs(this)
            result = full(sum(this.leadingStrandBoundSSBs, 1));
        end
        
        function result = get.numLaggingTemplateBoundSSBs(this)
            result = full(sum(this.laggingStrandBoundSSBs, 1));
        end
        
        function result = get.areLaggingStrandSSBSitesBound(this)
            nBndSSBs = this.numLaggingTemplateBoundSSBs;
            ssbFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_ssb8mer);
            ssbSpcg = this.ssbComplexSpacing;
            helPos = this.helicasePosition;
            fIdx = this.okazakiFragmentIndex;
            leadFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampGammaComplex);
            lagFtpt = this.enzymeDNAFootprints(this.enzymeIndexs_coreBetaClampPrimase);
            
            starts = zeros(1, 2);
            if fIdx(1) == 0
                starts(1) = this.chromosome.sequenceLen;
            else
                starts(1) = this.primaseBindingLocations{1}(fIdx(1));
            end
            if fIdx(2) == 0
                starts(2) = 1;
            else
                starts(2) = this.primaseBindingLocations{2}(fIdx(2));
            end
            
            result = false(1, 2);
            if helPos(1) ~= 0; result(1) = nBndSSBs(1) >= floor((starts(1) - helPos(1) - leadFtpt - lagFtpt) / (ssbFtpt + ssbSpcg) - 2); end;
            if helPos(2) ~= 0; result(2) = nBndSSBs(2) >= floor((helPos(2) - starts(2) - leadFtpt - lagFtpt) / (ssbFtpt + ssbSpcg) - 2); end;
        end
    end
end
