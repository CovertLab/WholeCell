%Protein Decay
%
% @wholeCellModelID Process_ProteinDecay
% @name             Protein Decay
% @description
%   Biology
%   ===========
%   Poisson simulation of protein damage, repair, and decay mediated by
%   - chaperone: clpB chaperone refolds/disaggregates proteins by an
%     energy-dependent mechanism
%   - proteases: lon, ftsH; each is believed to processively cleave
%     peptides into smaller peptides of approximately 20 amino acids
%     - ftsH: cleaves peptides with tmRNA stalled translation proteolysis
%       tag; parameterized by
%       - fragment length (aa)
%       - kinetic rate (cleaves/s)
%       - energy per cleavage (ATP/cleavage)
%     - lon: cleaves all other peptides; parameterized by
%       - fragment length (aa)
%       - kinetic rate (cleaves/s)
%       - energy per cleavage (ATP/cleavage)
%   - several peptidases: the specific functions and kinetics of the
%     individual peptidases are unknown and not modelled; decay proceeds as
%     long as at least one of each peptidase is present
%
% 1. Misfold proteins
%    - Protein misfolding is modeled as a poisson process
%      with a small rate constant
% 2. Refold cytosolic proteins
%    - Occurs if clpB is present
% 3. Decay macromolecular complexes
%    - Model decay is poisson process with rate parameter given by the
%      inverse weighted average half life of the complex's subunits
%    - Salvage bound prosthetic groups
%    - Mark subunits as damaged to be degraded by either the
%      protease/peptidase or ribonuclease machinery
% 4. Decay protein monomers and proteolysis tagged polypeptides
%    - Model selection of monomers that decay as poisson process with rate
%      parameters equal to the decay rates (ln(2)/half life as computed by
%      the N-end rule) of the monomers.
%    - Execute decay reactions as long as the following are available
%      - Lon/FtsH protease for normal / ribosome-stalled peptides
%      - Complete peptidase complement
%      - Water for hydrolysis of peptide bonds
%
%   Knowledge Base
%   ===============
%   The knowledge base contains implements the N-end rule which predicts the
%   half life of every protein monomer based on the average experimentally
%   measured (in E. coli) protein half lives for each possible N-terminal amino
%   acid. Exceptions:
%   - signal sequences are assumed to have 0 half life
%   - proteolysis tagged polypeptides are assumed to have 0 half life
%
%   The knowledge base predicts the half lives of complexes as the weighted mean
%   of that of their constituent protein monomers and RNAs.
%
%   The knowledge base ProteinMonomer and ProteinComplex classes also computes
%   hydrolytic degradation reactions of every protein monomer and complex:
%   - the (modified) amino acids released by hydrolytic cleavage of every
%     protein monomer
%   - water required for hydrolytic peptide bond cleavage
%   - prosthetic groups released by cleavage
%
%   The ProteinMonomer class is also used to compute hydrolytic degradation
%   reactions for proteolysis tagged polypeptides.
%
%   Representation
%   ===============
%   The substrates, enzymes, monomers, RNAs, and complexs properties represent
%   the counts of metabolites, proteases and peptidases, protein monomers,
%   damaged RNAs, protein complexes. The substrates and enzymes properties have
%   only 1 compartment (cytosol). The monomers, RNAs, and complexs properties
%   have 5 compartments. The monomers and complexs properties represent all
%   forms of proteins (nascent, mature, damaged, folded, misfolded, etc.)
%
%   abortedSequences represents the sequences of every proteolysis tagged
%   monomer in the cytosol.
%
%   The process contains intermediate representation of protein complex
%   degradation, but not of protein monomer degradation, misfolding, or
%   refolding. Degrading complexes are broken up into the constituent RNAs and
%   monomers, and these components are marked as "damaged" and recognized by RNA
%   and protein monomer degradation as molecules with 0 half life. Protein
%   monomer degradation is treated as an all-or-nothing event that either
%   proceeds to complete with a time step or doesn't progress at all.
%
%   proteinMisfoldingRate represents the rate at which every protein misfolds in
%   seconds.
%
%   monomerDecayRates and complexDecayRates represent the decay rate of each
%   protein monomer and complex species in seconds. misfolded and damaged
%   proteins, signal sequences, and proteolysis tagged polypeptides have 0 half
%   lives. All other proteins have positive half lives.
%
%   monomerDecayReactions and complexDecayReactions represent the metabolites
%   required and released by hydrolytic cleavage of protein monomers and the
%   break down of macromolecular complexes into their constituent subunits and
%   sequestered prosthetic groups. monomerDecayReactions and
%   complexDecayReactions are computed by the knowledge base ProteinMonomer and
%   ProteinComplex classes.
%
%   Initialization
%   ===============
%   All protein monomers and complexs are initialized to the mature state. This
%   is accomplished by the simulation class initializeState method.
%
%   Simulation
%   ===============
%   This process models misfolding, refolding, and protein monomer and complex
%   degradation as an enyzme-dependent (excepct misfolding) poisson processes
%   with rate parameter:
%      lambda = proteins .* decayRates * stepSizeSec
%
%   Algorithm
%   +++++++++++++++
%   Each of misfolding, refolding, protein complex degradation, proteolysis
%   tagged monomer degradation, and protein monomer degradation use a variant of
%   the algorithm:
%
%   1. Stochastically select proteins to misfold/refolding/degrade based on
%      poission distribution with
%         lambda = proteins .* rates * stepSizeSec
%   2. Limit refolding/degradation by availability of water and energy
%   3. Limit protein refolding/degradation by available enzyme activity
%      a. Refolding requires ClpB
%      b. Protein monomer degradation requires Lon protease and 6 peptidases
%      c. Proteolysis tagged polypeptide degradation requires FtsH protease and
%         6 peptidases
%   4. Update counts of proteins
%   5. Update counts of metabolic reactants and byproducts of protein
%      refolding/degradation
%
%   Compartments
%   +++++++++++++++
%   - Misfolding: occurs in all compartments
%   - Refolding: occurs only in compartments where ClpB is present, that is the
%     cytosol (and terminal organelle cytosol)
%   - Complex degradation: occurs only in compartments where the protein decay
%     machinery (Lon protease and peptidases) are presents, that is the cytosol
%     (and terminal organelle cytosol)
%   - Proteolysis tagged polypeptide degradation: only exist in cytosol, and are
%     only degraded their (assumes they are accessible from within the cytosol
%     to the integral membrane FtsH protease)
%   - Protein monomer degradation:  occurs only in compartments where the Lon
%     protease and peptidases are present, that is the cytosol.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/4/2010

classdef ProteinDecay < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'lonProteaseSpecificRate';
            'lonProteaseEnergyCost';
            'lonProteaseFragmentLength';
            'ftsHProteaseSpecificRate';
            'ftsHProteaseEnergyCost';
            'ftsHProteaseFragmentLength';
            'oligoendopeptidaseFSpecificRate';
            'proteinMisfoldingRate';
            'proteinComplexRNAComposition';
            'proteinComplexMonomerComposition';
            'monomerLonProteaseCleavages';
            'monomerDecayReactions';
            'complexDecayReactions';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'monomers';
            'complexs';
            'RNAs'};
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};   %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = {}; %whole cell model IDs of substrates
        substrateIndexs_atp              %index within substrates of ATP
        substrateIndexs_adp              %index within substrates of ADP
        substrateIndexs_phosphate        %index within substrates of inorganic phosphate
        substrateIndexs_hydrogen         %index within substrates of hydrogen
        substrateIndexs_water            %index within substrates of water
        substrateIndexs_aminoAcids       %indices within substrates of amino acids
        substrateIndexs_methionine       %index within substrates of methionine
        substrateIndexs_fmethionine      %index within substrates of formylmethionine
        substrateIndexs_glutamate        %index within substrates of glutamate
        substrateIndexs_glutamine        %index within substrates of glutamine
        substrateIndexs_formate          %index within substrates of formate
        substrateIndexs_ammonia          %index within substrates of ammonia

        enzymeWholeCellModelIDs = { %enzyme whole cell model ids
            'MG_355_HEXAMER';    %protease clpB            cytoplasmic
            'MG_239_HEXAMER';    %protease La              cytoplasmic
            'MG_457_HEXAMER';    %metalloprotease FtsH     integral membrane
            'MG_183_MONOMER';    %oligoendopeptidase F     cytoplasmic
            'MG_324_MONOMER';    %aminopeptidase           cytoplasmic
            'MG_391_HEXAMER';    %cytosol aminopeptidase   cytoplasmic
            'MG_208_DIMER';      %glycoprotease            cytoplasmic
            'MG_046_DIMER';      %metalloendopeptidase     cytoplasmic
            'MG_020_MONOMER'};   %proline iminopeptidase   cytoplasmic
        enzymeIndexs_clpBProtease          = 1;      %index within enzymes of protease clpB
        enzymeIndexs_lonProtease           = 2;      %index within enzymes of protease La
        enzymeIndexs_ftsHProtease          = 3;      %index within enzymes of metalloprotease FtsH
        enzymeIndexs_oligoendopeptidaseF   = 4;      %index within enzymes of oligoendopeptidase F
        enzymeIndexs_aminopeptidase        = 5;      %index within enzymes of aminopeptidase
        enzymeIndexs_cytosolAminopeptidase = 6;      %index within enzymes of cytosol aminopeptidase
        enzymeIndexs_glycoprotease         = 7;      %index within enzymes of glycoprotease
        enzymeIndexs_metalloendopeptidase  = 8;      %index within enzymes of metalloendopeptidase
        enzymeIndexs_prolineIminopeptidase = 9;      %index within enzymes of proline iminopeptidase
        enzymeIndexs_peptidases            = (4:9)'; %index within enzymes of peptidases
        
        rnaWholeCellModelIDs                        %whole cell model IDs of RNA subunits
        damagedRNAGlobalIndexs                      %indices of damaged rRNA, sRNA within simulation.RNAs
        
        complexIndexs_modified                      %indices of modified protein complexes within this.complexs
        complexIndexs_unmodified                    %indices of unmodified protein complexes within this.complexs
    end

    %fixed biological constants
    properties
        lonProteaseSpecificRate           %lon protease kinetic rate (1.667) [PUB_0029]
        lonProteaseEnergyCost    	      %lon protease energy requirement (6) [PUB_0029]
        lonProteaseFragmentLength         %length of peptide fragments after cleavage by lon protease (20) [PUB_0029]
        ftsHProteaseSpecificRate      	  %ftsH protease kinetic rate (0.03) [PUB_0031]
        ftsHProteaseEnergyCost            %ftsH protease energy requirement (8.3) [PUB_0031]
        ftsHProteaseFragmentLength        %length of peptide fragments after cleavage by ftsH protease (15) [PUB_0030]
        oligoendopeptidaseFSpecificRate   %oligoendopeptidase F kinetic rate (27.474) [PUB_0035]
        proteinMisfoldingRate             %expected number of misfolds per protein per time step (1e-6)

        proteinComplexRNAComposition      %RNA subunit composition of protein complexes (RNAs X protein complexes X compartments)
        proteinComplexMonomerComposition  %protein monomer subunit composition of protein complexes (protein monomers X protein complexes X compartments)

        monomerLonProteaseCleavages       %number of lon protease cleavages required to degrade each protein monomer
        monomerDecayReactions             %substrates consumed and produced by decay of each protein monomer
        complexDecayReactions             %substrates consumed and produced by decay of each protein complex
    end

    %global state (stored locally for convenience)
    properties
        RNAs                              %counts of RNAs
        monomers                          %counts of protein monomers
        complexs                          %counts of protein complexes
    end

    %global state (referenced locally for convenience)
    properties
        polypeptide                       %New polypeptides state class
    end
    
    %constructor
    methods
        function this = ProteinDecay(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            
            this.polypeptide = simulation.state('Polypeptide');
            this.states = [this.states; {this.polypeptide}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)            
            %include all necessary metabolites
            decayReaction_nascentMonomers     = reshape([knowledgeBase.proteinMonomers.decayReaction],                    [], knowledgeBase.numProteinMonomers);
            decayReaction_processedIMonomers  = reshape([knowledgeBase.proteinMonomers.processedISequenceDecayReaction],  [], knowledgeBase.numProteinMonomers);
            decayReaction_processedIIMonomers = reshape([knowledgeBase.proteinMonomers.processedIISequenceDecayReaction], [], knowledgeBase.numProteinMonomers);
            decayReaction_signalSequence      = reshape([knowledgeBase.proteinMonomers.signalSequenceDecayReaction],      [], knowledgeBase.numProteinMonomers);
            decayReaction_foldedMonomers      = reshape([knowledgeBase.proteinMonomers.foldedSequenceDecayReaction],      [], knowledgeBase.numProteinMonomers);
            decayReaction_matureMonomers      = reshape([knowledgeBase.proteinMonomers.matureSequenceDecayReaction],      [], knowledgeBase.numProteinMonomers);
            decayReaction_nascentComplexs     = reshape([knowledgeBase.proteinComplexs.decayReaction],                    [], knowledgeBase.numProteinComplexs);
            decayReaction_matureComplexs      = reshape([knowledgeBase.proteinComplexs.matureDecayReaction],              [], knowledgeBase.numProteinComplexs);
            
            m = this.metabolite;
            gtpIdxs = m.getIndexs({'GTP'; 'H2O'; 'PI'; 'H'});
            gdpIdxs = m.getIndexs('GDP');
            for i = 1:size(decayReaction_nascentComplexs, 2)
                tmp = decayReaction_nascentComplexs(gtpIdxs, i) ./ [1; 1; -1; -1];
                if tmp(1) && all(tmp == tmp(1))
                    decayReaction_nascentComplexs(gdpIdxs, i) = tmp(1);
                    decayReaction_nascentComplexs(gtpIdxs, i) = 0;
                end
                
                tmp = decayReaction_matureComplexs(gtpIdxs, i) ./ [1; 1 ;-1; -1];
                if tmp(1) && all(tmp == tmp(1))
                    decayReaction_matureComplexs(gdpIdxs, i) = tmp(1);
                    decayReaction_matureComplexs(gtpIdxs, i) = 0;
                end
            end
            
            this.substrateWholeCellModelIDs = unique([
                this.metabolite.wholeCellModelIDs(...
                any(decayReaction_nascentMonomers,     2) | ...
                any(decayReaction_processedIMonomers,  2) | ...
                any(decayReaction_processedIIMonomers, 2) | ...
                any(decayReaction_signalSequence,      2) | ...
                any(decayReaction_foldedMonomers,      2) | ...
                any(decayReaction_matureMonomers,      2) | ...
                any(decayReaction_nascentComplexs,     2) | ...
                any(decayReaction_matureComplexs,      2));
                'ATP'; 'ADP'; 'PI'; 'H'; 'H2O'; 'NH3'; 'FOR']);
            this.substrateIndexs_atp         = this.substrateIndexs({'ATP'});
            this.substrateIndexs_adp         = this.substrateIndexs({'ADP'});
            this.substrateIndexs_phosphate   = this.substrateIndexs({'PI'});
            this.substrateIndexs_hydrogen    = this.substrateIndexs({'H'});
            this.substrateIndexs_water       = this.substrateIndexs({'H2O'});
            this.substrateIndexs_aminoAcids  = this.substrateIndexs({'ALA'; 'ARG'; 'ASN'; 'ASP'; 'CYS'; 'GLN'; 'GLU'; 'GLY'; 'HIS'; 'ILE'; 'LEU'; 'LYS'; 'MET'; 'PHE'; 'PRO'; 'SER'; 'THR'; 'TRP'; 'TYR'; 'VAL'; 'FMET'});
            this.substrateIndexs_methionine  = this.substrateIndexs({'MET'});
            this.substrateIndexs_fmethionine = this.substrateIndexs({'FMET'});
            this.substrateIndexs_glutamate   = this.substrateIndexs({'GLU'});
            this.substrateIndexs_glutamine   = this.substrateIndexs({'GLN'});
            this.substrateIndexs_ammonia     = this.substrateIndexs({'NH3'});
            this.substrateIndexs_formate     = this.substrateIndexs({'FOR'});

            %call super class method
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});

            %RNA, monomer, complex whole cell model IDs
            this.rnaWholeCellModelIDs = this.rna.wholeCellModelIDs(this.rna.damagedIndexs([this.rna.matureRRNAIndexs; this.rna.matureSRNAIndexs]));
            
            %RNA indices
            this.damagedRNAGlobalIndexs = this.rna.damagedIndexs([this.rna.matureRRNAIndexs; this.rna.matureSRNAIndexs]);

            %complex composition
            cplxRNAComp     = knowledgeBase.proteinComplexAllRNAComposition([this.gene.rRNAIndexs; this.gene.sRNAIndexs], :, this.compartment.cytosolIndexs);
            cplxMonomerComp = knowledgeBase.proteinComplexAllMonomerComposition(:, :, this.compartment.cytosolIndexs);

            this.proteinComplexRNAComposition = zeros(size(cplxRNAComp,1), numel(this.complex.wholeCellModelIDs));
            this.proteinComplexRNAComposition(:, this.complex.nascentIndexs)     = cplxRNAComp;
            this.proteinComplexRNAComposition(:, this.complex.matureIndexs)      = cplxRNAComp;
            this.proteinComplexRNAComposition(:, this.complex.inactivatedIndexs) = cplxRNAComp;
            this.proteinComplexRNAComposition(:, this.complex.boundIndexs)       = cplxRNAComp;
            this.proteinComplexRNAComposition(:, this.complex.misfoldedIndexs)   = cplxRNAComp;
            this.proteinComplexRNAComposition(:, this.complex.damagedIndexs)     = cplxRNAComp;

            this.proteinComplexMonomerComposition = zeros(size(cplxMonomerComp,1), numel(this.complex.wholeCellModelIDs));
            this.proteinComplexMonomerComposition(:, this.complex.nascentIndexs)     = cplxMonomerComp;
            this.proteinComplexMonomerComposition(:, this.complex.matureIndexs)      = cplxMonomerComp;
            this.proteinComplexMonomerComposition(:, this.complex.inactivatedIndexs) = cplxMonomerComp;
            this.proteinComplexMonomerComposition(:, this.complex.boundIndexs)       = cplxMonomerComp;
            this.proteinComplexMonomerComposition(:, this.complex.misfoldedIndexs)   = cplxMonomerComp;
            this.proteinComplexMonomerComposition(:, this.complex.damagedIndexs)     = cplxMonomerComp;

            %decay reactions           
            this.monomerDecayReactions = zeros(numel(this.substrateWholeCellModelIDs), numel(this.monomer.wholeCellModelIDs));
            this.complexDecayReactions = zeros(numel(this.substrateWholeCellModelIDs), numel(this.complex.wholeCellModelIDs));
            
            this.monomerDecayReactions(:, this.monomer.nascentIndexs)        = decayReaction_nascentMonomers(    this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.processedIIndexs)     = decayReaction_processedIMonomers( this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.processedIIIndexs)    = decayReaction_processedIIMonomers(this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.signalSequenceIndexs) = decayReaction_signalSequence(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.foldedIndexs)         = decayReaction_foldedMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.matureIndexs)         = decayReaction_matureMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.inactivatedIndexs)    = decayReaction_matureMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.boundIndexs)          = decayReaction_matureMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.misfoldedIndexs)      = decayReaction_matureMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            this.monomerDecayReactions(:, this.monomer.damagedIndexs)        = decayReaction_matureMonomers(     this.substrateMetaboliteGlobalIndexs, :);
            
            this.complexDecayReactions(:, this.complex.nascentIndexs)        = decayReaction_nascentComplexs(    this.substrateMetaboliteGlobalIndexs, :);
            this.complexDecayReactions(:, this.complex.matureIndexs)         = decayReaction_matureComplexs(     this.substrateMetaboliteGlobalIndexs, :);
            this.complexDecayReactions(:, this.complex.inactivatedIndexs)    = decayReaction_matureComplexs(     this.substrateMetaboliteGlobalIndexs, :);
            this.complexDecayReactions(:, this.complex.boundIndexs)          = decayReaction_matureComplexs(     this.substrateMetaboliteGlobalIndexs, :);
            this.complexDecayReactions(:, this.complex.misfoldedIndexs)      = decayReaction_matureComplexs(     this.substrateMetaboliteGlobalIndexs, :);
            this.complexDecayReactions(:, this.complex.damagedIndexs)        = decayReaction_matureComplexs(     this.substrateMetaboliteGlobalIndexs, :);
            
            this.complexIndexs_modified = find(any(this.complexDecayReactions > 0, 1)');
            this.complexIndexs_unmodified = setdiff((1:numel(this.complex.wholeCellModelIDs))', this.complexIndexs_modified);
            
            this.monomerLonProteaseCleavages = max(0, ceil(this.monomer.lengths/this.lonProteaseFragmentLength)-1);
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.stanford.covert.cell.sim.Process();
            
            this.monomers = this.monomer.counts;
            this.complexs = this.complex.counts;
            this.RNAs = this.rna.counts(this.damagedRNAGlobalIndexs, :, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.stanford.covert.cell.sim.Process();

            this.monomer.counts = this.monomers;
            this.complex.counts = this.complexs;
            this.rna.counts(this.damagedRNAGlobalIndexs, :, :) = this.RNAs;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.stanford.covert.cell.sim.Process(numTimePoints);

            numCompartments = this.compartment.count;

            this.RNAs     = zeros(size(this.rnaWholeCellModelIDs,1),      numCompartments, numTimePoints);
            this.monomers = zeros(size(this.monomer.wholeCellModelIDs,1), numCompartments, numTimePoints);
            this.complexs = zeros(size(this.complex.wholeCellModelIDs,1), numCompartments, numTimePoints);
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
            %refolding
            nRefolding = (sum(states.monomers) + sum(states.complexs)) * this.proteinMisfoldingRate;
            
            bmProd(this.substrateIndexs_atp)       = bmProd(this.substrateIndexs_atp)       + nRefolding;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + nRefolding;
            byProd(this.substrateIndexs_adp)       = byProd(this.substrateIndexs_adp)       + nRefolding;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + nRefolding;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + nRefolding;
            
            %complex degradation
            bmProd = bmProd + max(0, -this.complexDecayReactions(:, this.complex.matureIndexs)) * states.complexDecays;
            byProd = byProd + max(0,  this.complexDecayReactions(:, this.complex.matureIndexs)) * states.complexDecays;
            
            %mature monomer degradation
            bmProd = bmProd + max(0, -this.monomerDecayReactions(:, this.monomer.matureIndexs)) * states.monomerDecays;
            byProd = byProd + max(0,  this.monomerDecayReactions(:, this.monomer.matureIndexs)) * states.monomerDecays;
            
            lonCleavages = this.lonProteaseEnergyCost * states.monomerDecays' * this.monomerLonProteaseCleavages(this.monomer.matureIndexs);
            bmProd(this.substrateIndexs_atp)       = bmProd(this.substrateIndexs_atp)       + lonCleavages;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + lonCleavages;
            byProd(this.substrateIndexs_adp)       = byProd(this.substrateIndexs_adp)       + lonCleavages;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + lonCleavages;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + lonCleavages;
            
            %signal sequence degradation
            bmProd = bmProd + max(0, -this.monomerDecayReactions(:, this.monomer.signalSequenceIndexs)) * states.monomerProductions;
            byProd = byProd + max(0,  this.monomerDecayReactions(:, this.monomer.signalSequenceIndexs)) * states.monomerProductions;
            
            lonCleavages = this.lonProteaseEnergyCost * states.monomerProductions' * ...
                (this.monomerLonProteaseCleavages(this.monomer.signalSequenceIndexs) .* ...
                any(this.monomerDecayReactions(:, this.monomer.signalSequenceIndexs), 1)');
            bmProd(this.substrateIndexs_atp)       = bmProd(this.substrateIndexs_atp)       + lonCleavages;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + lonCleavages;
            byProd(this.substrateIndexs_adp)       = byProd(this.substrateIndexs_adp)       + lonCleavages;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + lonCleavages;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + lonCleavages;
            
            %stalled polypeptide degradation
            ftsHCleavages = this.ftsHProteaseEnergyCost * 0;
            bmProd(this.substrateIndexs_atp)       = bmProd(this.substrateIndexs_atp)       + ftsHCleavages;
            bmProd(this.substrateIndexs_water)     = bmProd(this.substrateIndexs_water)     + ftsHCleavages;
            byProd(this.substrateIndexs_adp)       = byProd(this.substrateIndexs_adp)       + ftsHCleavages;
            byProd(this.substrateIndexs_phosphate) = byProd(this.substrateIndexs_phosphate) + ftsHCleavages;
            byProd(this.substrateIndexs_hydrogen)  = byProd(this.substrateIndexs_hydrogen)  + ftsHCleavages;
            
            %% enzymes
            minEnzExp(this.enzymeIndexs_lonProtease) = ...
                2 * states.monomerDecays0' * this.monomerLonProteaseCleavages(this.monomer.matureIndexs) / ...
                this.lonProteaseSpecificRate;
            minEnzExp(this.enzymeIndexs_clpBProtease) = 2;
            minEnzExp(this.enzymeIndexs_ftsHProtease) = 70; %to keep pace with degrading aborted polypeptides
            minEnzExp(this.enzymeIndexs_peptidases)   = 2;
        end

        %initialization: proteins initialized to mature/bound/inactivated state
        %by simulation initializeState method
        function initializeState(~)
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            
            mnmers = sum(this.monomers(:, [this.compartment.cytosolIndexs this.compartment.terminalOrganelleCytosolIndexs]), 2);
            cmplxs = sum(this.complexs(:, [this.compartment.cytosolIndexs this.compartment.terminalOrganelleCytosolIndexs]), 2);
            
            %refolding
            nRefolding = (this.enzymes(this.enzymeIndexs_clpBProtease) > 0) * ( ...
                sum(mnmers(this.monomer.misfoldedIndexs)) + ...
                sum(cmplxs(this.complex.misfoldedIndexs)));
            
            if nRefolding > 0 && this.enzymes(this.enzymeIndexs_clpBProtease)
                result(this.substrateIndexs_atp)   = result(this.substrateIndexs_atp)   + nRefolding;
                result(this.substrateIndexs_water) = result(this.substrateIndexs_water) + nRefolding;
            end
            
            %complex degradation
            cmplxDecays = cmplxs .* min(1, this.complex.decayRates);
            result = result + max(0, -this.complexDecayReactions) * cmplxDecays;
            
            %monomer degradation
            mnmerDecays = mnmers .* min(1, this.monomer.decayRates);
            result = result + max(0, -this.monomerDecayReactions) * mnmerDecays;
            result([this.substrateIndexs_atp; this.substrateIndexs_water]) = ...
                + result([this.substrateIndexs_atp; this.substrateIndexs_water]) ...
                + all(this.enzymes(this.enzymeIndexs_peptidases)) * this.lonProteaseEnergyCost * min(...
                    this.enzymes(this.enzymeIndexs_lonProtease) * this.lonProteaseSpecificRate, ...
                    mnmerDecays' * this.monomerLonProteaseCleavages);
            
            %stalled polypeptide degradation
            if ~isempty(this.polypeptide.abortedPolypeptides)
                result([this.substrateIndexs_atp; this.substrateIndexs_water]) = ...
                    result([this.substrateIndexs_atp; this.substrateIndexs_water]) + ...
                    ceil(this.ftsHProteaseEnergyCost * min(...
                        ceil(this.enzymes(this.enzymeIndexs_ftsHProtease) * this.ftsHProteaseSpecificRate), ...
                        sum(ceil(this.polypeptide.abortedSequenceLengths / this.ftsHProteaseFragmentLength) - 1)));
                        
                result(this.substrateIndexs_water) = ...
                    result(this.substrateIndexs_water) + ...
                    sum(this.polypeptide.abortedSequenceLengths - 1);
            end
        end

        %simulation
        function evolveState(this)
            this.evolveState_MisfoldProteins();   %misfold proteins
            this.evolveState_RefoldProteins();    %refold proteins
            this.evolveState_DegradeComplexes();  %degrade complexes to damaged monomers and RNAs
            this.evolveState_DegradeAbortedPolypeptides(); %degrade proteolysis tagged monomers to amino acids
            this.evolveState_DegradeMonomers();   %degrade monomers to amino acids (including monomers marked as damaged by evolveState_DegradeComplexes)
        end

        %Model protein misfolding (of mature; inactive; and inactivated forms)
        %as a passive poisson process with rate proteinMisfoldingRate.
        %
        %Protein misfolding does not require metabolites or enzymes.
        function evolveState_MisfoldProteins(this)
            %% monomers
            misfolding = this.randStream.stochasticRound(...
                this.monomers([this.monomer.matureIndexs; this.monomer.boundIndexs; this.monomer.inactivatedIndexs], :) ...
                * min(1, this.proteinMisfoldingRate * this.stepSizeSec));
            misfolding(:, this.compartment.membraneIndexs) = 0;
            misfolding(:, this.compartment.terminalOrganelleMembraneIndexs) = 0;
            
            if any(misfolding(:))
                tmp = zeros(size(this.monomers));
                tmp([this.monomer.matureIndexs; this.monomer.boundIndexs; this.monomer.inactivatedIndexs], :) = -misfolding;
                tmp(this.monomer.misfoldedIndexs, :) = sum(permute(reshape(misfolding', [size(this.monomers, 2) numel(this.monomer.misfoldedIndexs) 3]), [2 1 3]), 3);
                notUpdatingProteins = this.monomer.updateExternalState(tmp, true);
                misfolding = misfolding - notUpdatingProteins([this.monomer.matureIndexs; this.monomer.boundIndexs; this.monomer.inactivatedIndexs], :);
                
                this.monomers([this.monomer.matureIndexs; this.monomer.boundIndexs; this.monomer.inactivatedIndexs], :) = ...
                    this.monomers([this.monomer.matureIndexs; this.monomer.boundIndexs; this.monomer.inactivatedIndexs], :) - ...
                    misfolding;
                this.monomers(this.monomer.misfoldedIndexs, :) = ...
                    this.monomers(this.monomer.misfoldedIndexs, :) + ...
                    sum(permute(reshape(misfolding', [size(this.monomers, 2) numel(this.monomer.misfoldedIndexs) 3]), [2 1 3]), 3);
            end
            
            %% complexs
            misfolding = this.randStream.stochasticRound(...
                this.complexs([this.complex.matureIndexs; this.complex.boundIndexs; this.complex.inactivatedIndexs], :) * ...
                min(1, this.proteinMisfoldingRate * this.stepSizeSec));
            misfolding(:, this.compartment.membraneIndexs) = 0;
            misfolding(:, this.compartment.terminalOrganelleMembraneIndexs) = 0;
            
            if any(misfolding(:))
                tmp = zeros(size(this.complex.counts));
                tmp([this.complex.matureIndexs; this.complex.boundIndexs; this.complex.inactivatedIndexs], :) = -misfolding;
                tmp(this.complex.misfoldedIndexs, :) = sum(permute(reshape(misfolding', [size(this.complexs, 2) numel(this.complex.misfoldedIndexs) 3]), [2 1 3]), 3);
                notUpdatingProteins = this.complex.updateExternalState(tmp, true);
                misfolding = misfolding - notUpdatingProteins([this.complex.matureIndexs; this.complex.boundIndexs; this.complex.inactivatedIndexs], :);
                
                this.complexs([this.complex.matureIndexs; this.complex.boundIndexs; this.complex.inactivatedIndexs], :) = ...
                    this.complexs([this.complex.matureIndexs;this.complex.boundIndexs; this.complex.inactivatedIndexs], :) - ...
                    misfolding;
                this.complexs(this.complex.misfoldedIndexs, :) = ...
                    this.complexs(this.complex.misfoldedIndexs, :) + ...
                    sum(permute(reshape(misfolding', [size(this.complexs, 2) numel(this.complex.misfoldedIndexs) 3]), [2 1 3]), 3);
            end
        end

        %Model protein refolding as a binary metabolite-limited process:
        %If active ClpB is present in the cytosol and ATP is present, refold
        %all misfolded cytoplasmic proteins up to the amount of available ATP.
        %
        %Because the ATP cost of ClpB-mediated protein refolding/disagggregation
        %is unknown we assume that 1 ATP is required per refolding.
        %
        %Refolding proteins transition to the active state.
        function evolveState_RefoldProteins(this)
            %cumulative sum of misfolded monomers and complexes
            if ...
                    ~any(this.monomers(this.monomer.misfoldedIndexs, this.compartment.cytosolIndexs))                  && ...
                    ~any(this.monomers(this.monomer.misfoldedIndexs, this.compartment.terminalOrganelleCytosolIndexs)) && ...
                    ~any(this.complexs(this.complex.misfoldedIndexs, this.compartment.cytosolIndexs))                  && ...
                    ~any(this.complexs(this.complex.misfoldedIndexs, this.compartment.terminalOrganelleCytosolIndexs))
                return;
            end
            
            %indices
            cIdxs = [this.compartment.cytosolIndexs; this.compartment.terminalOrganelleCytosolIndexs];
            reactantIdxs = [this.substrateIndexs_atp; this.substrateIndexs_water];
            productIdxs  = [this.substrateIndexs_adp; this.substrateIndexs_phosphate; this.substrateIndexs_hydrogen];
            
            %stochastically select proteins to refold, up to limit of available
            %reactants
            cumsumMisfoldedProteins = cumsum([0;
                this.monomers(this.monomer.misfoldedIndexs, this.compartment.cytosolIndexs);
                this.monomers(this.monomer.misfoldedIndexs, this.compartment.terminalOrganelleCytosolIndexs);
                this.complexs(this.complex.misfoldedIndexs, this.compartment.cytosolIndexs);
                this.complexs(this.complex.misfoldedIndexs, this.compartment.terminalOrganelleCytosolIndexs)]);
            numFoldingProteins = (this.enzymes(this.enzymeIndexs_clpBProtease) > 0) * min([
                this.substrates(reactantIdxs)
                cumsumMisfoldedProteins(end)]);
            if numFoldingProteins == 0
                return;
            end
            
            foldingProteinIdxs = this.randStream.randperm(cumsumMisfoldedProteins(end));
            
            foldingMonomers = zeros(size(this.monomer.matureIndexs,1), 2);
            foldingComplexs = zeros(size(this.complex.matureIndexs,1), 2);
            
            for i = 1:numFoldingProteins
                j = find(foldingProteinIdxs(i) <= cumsumMisfoldedProteins, 1, 'first') - 1;
                if j <= numel(this.monomer.matureIndexs)
                    foldingMonomers(j,1) = foldingMonomers(j,1)+1;
                elseif j <= 2*numel(this.monomer.matureIndexs)
                    foldingMonomers(j - numel(this.monomer.matureIndexs), 2) = foldingMonomers(j - numel(this.monomer.matureIndexs), 2) + 1;
                elseif j <= 2*numel(this.monomer.matureIndexs)+numel(this.complex.matureIndexs)
                    foldingComplexs(j - 2*numel(this.monomer.matureIndexs), 1) = foldingComplexs(j - 2*numel(this.monomer.matureIndexs), 1) + 1;
                else
                    foldingComplexs(j - 2*numel(this.monomer.matureIndexs) - numel(this.complex.matureIndexs), 2) = foldingComplexs(j - 2*numel(this.monomer.matureIndexs)-numel(this.complex.matureIndexs), 2) + 1;
                end
            end
            
            %update metabolites
            this.substrates(reactantIdxs) = this.substrates(reactantIdxs) - numFoldingProteins;
            this.substrates(productIdxs)  = this.substrates(productIdxs)  + numFoldingProteins;
            
            %update monomers
            this.monomers(this.monomer.matureIndexs, cIdxs) = ...
                this.monomers(this.monomer.matureIndexs, cIdxs) + ...
                foldingMonomers;
            this.monomers(this.monomer.misfoldedIndexs, cIdxs) = ...
                this.monomers(this.monomer.misfoldedIndexs, cIdxs) - ...
                foldingMonomers;
            
            %update complexs
            this.complexs(this.complex.matureIndexs, cIdxs) = ...
                this.complexs(this.complex.matureIndexs, cIdxs) + ...
                foldingComplexs;
            this.complexs(this.complex.misfoldedIndexs, cIdxs) = ...
                this.complexs(this.complex.misfoldedIndexs, cIdxs) - ...
                foldingComplexs;
        end

        %Model complex degradation as a stochastic enzyme-dependent poisson
        %process. To simplify the modeling of complex degradation, this method
        %first stochastically selects cytosolic complexes to degrade, and then
        %breaks up the selected complexes into their constituents parts:
        %prosthetic groups, protein monomers, and RNAs, and finally relies on
        %the protein  monomer and RNA degradation models
        %(evolveState_DegradeMonomers and RNA Decay process, respectively) to
        %degrade the individual components. The method marks the component
        %protein monomers and RNAs and selected complexes as "damaged" (that is
        %in a state with 0 half life), so that the the protein monomer and RNA
        %machinery know to degrade these components.
        %
        %Note: only cytosolic complexes are degraded because the M. genitalium
        %protein degradation machinery (protease La) is a cytosolic protein, except that
		%all "knocked" out proteins are degraded.
        function evolveState_DegradeComplexes(this)
            %compartment indices
            iC = this.compartment.cytosolIndexs;
            iTC = this.compartment.terminalOrganelleCytosolIndexs;
            
            %stochastically choose complexes to decay
            decayingRates = this.complexs .* ...
                min(1e6, this.complex.decayRates(:, ones(1, this.compartment.count))) * this.stepSizeSec;
			decayingRates(this.complex.decayRates < 1e6, [1:iC-1  iC+1:iTC-1  iTC+1:end]) = 0;
            decayingProteins = min(this.complexs, this.randStream.stochasticRound(decayingRates));
            decayingRates(decayingProteins == 0) = 0;
            if ~any(decayingProteins(:))
                return;
            end
            decayedProteins = zeros(size(decayingProteins));
            
            %% decay unmodified complexes
            decayingSimpleProteins = decayingProteins(this.complexIndexs_unmodified, :);
            
            %decrement unmodified complexes
            if any(decayingSimpleProteins(:))
                decayedProteins(this.complexIndexs_unmodified, :) = decayingSimpleProteins;
            end
            
            %% decay modified complexes up to limit of metabolites
            %(complex decay requires metabolites to reverse the affects of any
            %modifications)
            decayingModifiedRates = decayingRates(this.complexIndexs_modified, :);
            decayingModifiedProteins = decayingProteins(this.complexIndexs_modified, :);
            decayedModifiedProteins = zeros(size(decayingModifiedProteins));
            substrates = this.substrates;
            while any(decayingModifiedProteins(:))
                %randomly select complex, weighted by copy number * decay rate
                idx = this.randStream.randsample(numel(decayingModifiedRates), 1, true, decayingModifiedRates(:));
                iProtein = this.complexIndexs_modified(mod(idx - 1, numel(this.complexIndexs_modified)) + 1);
                
                %check if sufficient resources available
                if any(substrates < max(0, -this.complexDecayReactions(:, iProtein)))
                    break;
                end
                
                %decrement complexes
                decayingModifiedProteins(idx) = decayingModifiedProteins(idx) - 1;
                decayedModifiedProteins(idx) = decayedModifiedProteins(idx) + 1;
                if decayingModifiedProteins(idx) == 0
                    decayingModifiedRates(idx) = 0;
                else
                    decayingModifiedRates(idx) = max(0, decayingModifiedRates(idx) - min(1e6, this.complex.decayRates(iProtein)) * this.stepSizeSec);
                end
                
                %release prosthetic groups,
                %invert any modifications
                substrates = substrates + this.complexDecayReactions(:, iProtein);
            end
            decayedProteins(this.complexIndexs_modified, :) = decayedModifiedProteins;
            
            tmp = -decayedProteins;
            notUpdatingProteins = this.complex.updateExternalState(tmp, true);
            decayedProteins = decayedProteins - notUpdatingProteins;
            
            %decrease complexes. mark protein and RNA monomers for decay
            %(that is put them in a state with 0 half life) by
            %evolveState_DegradeMonomers and RNADecay
            this.complexs = this.complexs - decayedProteins;
            this.monomers(this.monomer.damagedIndexs, iC) = ...
                this.monomers(this.monomer.damagedIndexs, iC) + ...
                this.proteinComplexMonomerComposition * sum(decayedProteins, 2);
            this.RNAs(:, iC) = ...
                this.RNAs(:, iC) + ...
                this.proteinComplexRNAComposition * sum(decayedProteins, 2);
            this.substrates = ...
                this.substrates + ...
                this.complexDecayReactions * sum(decayedProteins, 2);
        end

        %Stochastically degrade proteolysis tagged monomers up to limit of
        %available enzymes (metalloprotease FtsH and 6 peptidases) and their
        %kinetics (only of FtsH; peptidase kinetics haven't been characterized;
        %in this case we model that degradation can occur as only as at least 1
        %of each peptidase species is present) and metabolites (hydrolytic
        %cleavage requires water, and FtsH requires ATP). Because the specific
        %function and mechanism of each peptidase hasn't been well characterized
        %we assume that degradation requires all peptidases, and they don't
        %require ATP. Proteolysis tagged polypeptides are assumed to be degraded
        %immediately, that is they have a half life of 0.
        function evolveState_DegradeAbortedPolypeptides(this)
            %terminate early if no proteolysis tagged monomers
            if isempty(this.polypeptide.abortedPolypeptides)
                return;
            end
            
            % import classes
            import edu.stanford.covert.cell.kb.ProteinMonomer;
            
            %indices
            aaWaterIdxs = [
                this.substrateIndexs_aminoAcids;
                this.substrateIndexs_water];
            energyIdxs = [
                this.substrateIndexs_atp;
                this.substrateIndexs_adp;
                this.substrateIndexs_phosphate;
                this.substrateIndexs_water;
                this.substrateIndexs_hydrogen];

            %randomly order proteolysis tagged monomer
            if numel(this.polypeptide.abortedPolypeptides) > 1
                this.polypeptide.abortedPolypeptides = ...
                    this.polypeptide.abortedPolypeptides(this.randStream.randperm(size(this.polypeptide.abortedPolypeptides, 1)), :);
            end
            
            %determine number of monomers that can be cleaved
            cleavages = ceil(this.polypeptide.abortedSequenceLengths / this.ftsHProteaseFragmentLength) - 1;
            peptidaseCosts = (1:numel(cleavages))';
            if any(cleavages)
                ftsHCosts = this.randStream.stochasticRound((this.enzymes(this.enzymeIndexs_ftsHProtease) * this.stepSizeSec) ./ cumsum(cleavages / this.ftsHProteaseSpecificRate));
                ftsHCosts(find(ftsHCosts == 0, 1, 'first'):end) = 0;
            else
                ftsHCosts = ones(size(cleavages));
            end
            energyCosts = cumsum(this.randStream.stochasticRound(cleavages * this.ftsHProteaseEnergyCost));
            waterCosts = energyCosts + cumsum(this.polypeptide.abortedSequenceLengths - 1);
            numReactions = find(...
                [0; peptidaseCosts] <= min(this.enzymes(this.enzymeIndexs_peptidases)) * this.stepSizeSec & ...
                [1; ftsHCosts]      >= 1 & ...
                [0; energyCosts]    <= this.substrates(this.substrateIndexs_atp) & ...
                [0; waterCosts]     <= this.substrates(this.substrateIndexs_water), ...
                1, 'last') - 1;
            
            %stop if no degradation
            if isempty(numReactions) || numReactions == 0
                return;
            end
            
            %loop over tagged monomers
            decayReactions = zeros(22, 1);
            abortedSeqs = this.polypeptide.abortedSequences;
            for i = 1:numReactions
                decayReactions = ...
                    + decayReactions ...
                    + ProteinMonomer.computeDecayReaction(...
                    ProteinMonomer.computeBaseCount(abortedSeqs{i}, 22, 1:21, true),...
                    length(abortedSeqs{i}), 22)';
            end

            %% update state
            %metabolites
            this.substrates(aaWaterIdxs) = this.substrates(aaWaterIdxs) + decayReactions;
            this.substrates(energyIdxs)  = this.substrates(energyIdxs) + [-1; 1; 1; -1; 1] * energyCosts(numReactions);

            %proteolysis tagged monomers
            this.polypeptide.abortedPolypeptides = this.polypeptide.abortedPolypeptides(numReactions+1:end, :);
        end

        %Stochastically degrade monomers weighted by copy number and decay rate
        %up to limit of available enzymes (monomers require protease La and 6
        %peptidases) and their kinetics (except for the peptidases where
        %their kinetics haven't been characterized; in this case we model that
        %degradation can occur as only as at least 1 of each peptidase species
        %is present) and metabolites (hydrolytic cleavage requires water, and
        %Lon requires ATP). Because the specific function and  mechanism
        %of each peptidase hasn't been well characterized we assume that
        %degradation requires all peptidases, and they don't require ATP.
        function evolveState_DegradeMonomers(this)
            %indices
            iC = this.compartment.cytosolIndexs;
            iTC = this.compartment.terminalOrganelleCytosolIndexs;
            
            %decaying monomers
            decayingRates = this.monomers .* ...
                min(1e6, this.monomer.decayRates(:, ones(1, this.compartment.count))) * this.stepSizeSec;
			decayingRates(this.monomer.decayRates < 1e6, [1:iC-1  iC+1:iTC-1  iTC+1:end]) = 0;
            decayingProteins = min(this.monomers, this.randStream.stochasticRound(decayingRates));
            decayingRates(decayingProteins == 0) = 0;
            if ~any(decayingProteins(:))
                return;
            end
            decayedProteins = zeros(size(decayingProteins));
            
            iEnergy = [
                this.substrateIndexs_atp
                this.substrateIndexs_adp
                this.substrateIndexs_phosphate
                this.substrateIndexs_water
                this.substrateIndexs_hydrogen
                ];
            
            protease  = this.enzymes(this.enzymeIndexs_lonProtease) * this.stepSizeSec * this.lonProteaseSpecificRate;
            peptidase = this.randStream.stochasticRound(min(this.enzymes(this.enzymeIndexs_peptidases)) * this.stepSizeSec);
            
            %decay monomers up to limit of metabolites and enzymes
            substrates = this.substrates;
            while any(decayingProteins(:))
                %randomly select monomer, weighted by copy number * decay rate
                idx = this.randStream.randsample(numel(decayingRates), 1, true, decayingRates(:));
                iProtein = mod(idx - 1, size(this.monomers, 1)) + 1;
                
                %check if sufficient resources available
                substrateCost = -this.monomerDecayReactions(:, iProtein);
                substrateCost(iEnergy) = substrateCost(iEnergy) + ...
                    [1; -1; -1; 1; -1] * this.lonProteaseEnergyCost * this.monomerLonProteaseCleavages(iProtein);
                proteaseCost = this.monomerLonProteaseCleavages(iProtein);
                
                if ...
                        any(substrates(iEnergy) < max(0, substrateCost(iEnergy))) || ...
                        this.randStream.stochasticRound(protease / proteaseCost) < 1 || ...
                        peptidase < 1
                    break;
                end
                
                %update state
                decayingProteins(idx) = decayingProteins(idx) - 1;
                decayedProteins(idx) = decayedProteins(idx) + 1;
                if decayingProteins(idx) <= 0
                    decayingRates(idx) = 0;
                else
                    decayingRates(idx) = max(0, decayingRates(idx) - min(1e6, this.monomer.decayRates(iProtein)) * this.stepSizeSec);
                end
                
                substrates = substrates - substrateCost;
                protease = max(0, protease - proteaseCost);
            end
            
            tmp = -decayedProteins;
            notUpdatingProteins = this.monomer.updateExternalState(tmp, true);
            decayedProteins = decayedProteins - notUpdatingProteins;
            
            this.monomers = this.monomers - decayedProteins;
            this.substrates = ...
                this.substrates + ...
                this.monomerDecayReactions * sum(decayedProteins, 2);
            this.substrates(iEnergy) = ...
                this.substrates(iEnergy) - ...
                [1; -1; -1; 1; -1] * this.lonProteaseEnergyCost * (this.monomerLonProteaseCleavages' * sum(decayedProteins, 2));
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            rnaMolecularWeights = this.rna.molecularWeights(this.rna.damagedIndexs([this.rna.matureRRNAIndexs; this.rna.matureSRNAIndexs]));
            
            if size(this.RNAs, 3) == 1
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    rnaMolecularWeights'           * sum(this.RNAs,     2) + ...
                    this.monomer.molecularWeights' * sum(this.monomers, 2) + ...
                    this.complex.molecularWeights' * sum(this.complexs, 2)) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.stanford.covert.cell.sim.Process() + (...
                    permute(rnaMolecularWeights'           * permute(sum(this.RNAs,     2), [1 3 2]), [1 3 2]) + ...
                    permute(this.monomer.molecularWeights' * permute(sum(this.monomers, 2), [1 3 2]), [1 3 2]) + ...
                    permute(this.complex.molecularWeights' * permute(sum(this.complexs, 2), [1 3 2]), [1 3 2])) ...
                    / edu.stanford.covert.util.ConstantUtil.nAvogadro;
            end
        end
    end
end