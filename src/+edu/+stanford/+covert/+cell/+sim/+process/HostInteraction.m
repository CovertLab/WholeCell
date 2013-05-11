%Host Interaction
%
% @wholeCellModelID Process_HostInteraction
% @name             Host Interaction
% @description
%   Introduction
%   ======================
%   Medical relevance
%   ++++++++++++++++++++++
%   14 Mycoplasmas and 2 ureaplasmas have been detected in humans
%   [PUB_0840]. 6 of these species, including M. genitalium have been found
%   in the urogenital tract. M. genitalium is a common urogenital pathogen,
%   and one of the most common causes of non-gnoccocal urethritis. M.
%   genitalium is estimated to colonize the urogenital tract of 3% of all
%   women [PUB_0840], and up to 9-10% in some populations [PUB_0841]. Most
%   M. genitalium infections are asymptomatic, and many researchers have
%   written that M. genitalium is close to "ideal parasite" [PUB_0313],
%   typically living in  harmony with its host. A small fraction of M.
%   genitalium infections manifest as urethritis, cervicitis, and pelvic
%   inflammatory disease.
%
%   Pathophysiology
%   ++++++++++++++++++++++
%   Little is known about the pathophysiology of M. genitalium. The
%   literature on M. genitalium pathophysiology has focused on a few
%   topics:
%   1. Adhesion
%      Several studies have suggested that the terminal organelle
%      lipoproteins MG_191, MG_192, MG_217, MG_318, and MG_386 are involved
%      in adhesion [PUB_0088, PUB_0313]
%   2. Host immune response activation
%      Several studies have suggested that the gene products of MG_149,
%      MG_309, and MG_412 activate the host NF-kB response through
%      interaction with host TLR receptors 1, 2, and 6 [PUB_0194, PUB_0842,
%      PUB_0844].
%   3. Antigenicity and antigenic variation
%      Several studies have focused on antigenticity and antigenic
%      variation. Razin reported that MgPa (MG_191) is one of the main
%      antigenic components, and that antibodies against it may be found in
%      infected patients [PUB_0088]. Shimizu et al have observed two
%      differentially expressed forms of MG_149 -- P18 and P20 [PUB_0844].
%   4. Toxic byproducts
%      M. genitalium byproducts including hydrogen peroxide and superoxide
%      radicals have been suggested to cause oxidative damage to host
%      membranes [PUB_0313].
%
%   The following table summarizes genes that have been implicated
%   interactions with mammalian hosts.
%
%     Gene    Symbols                  Localization         Function
%     ======  =======================  ============         ===========
%     MG_075                           lipoprotein          antigen
%     MG_149  N-ALP2, P18, P20         lipoprotein          activates NF-kB as TLR 1/2 ligand
%     MG_191  MgPa, P140               Trm Org lipoprotein  adhesion
%     MG_192  P110, P40, P90, pre-B/C  Trm Org lipoprotein  adhesion
%     MG_200                           cytosol              motility
%     MG_217  P65                      Trm Org lipoprotein  adhesion
%     MG_218  HMW2                     Trm Org core         motility, attachment
%     MG_288                           cytosol              homologous to protein L, but annotation as protein L is inconsistent with cytosolic localization; not currently modeled as participating in interaction with host  
%     MG_309                           lipoprotein          activates NF-kB as TLR 2/6 ligand
%     MG_312  HMW1, P65                Trm Org core         motility, attachment
%     MG_317  HMW3, P200, P69          Trm Org core         motility, attachment
%     MG_318  P30, P32                 Trm Org lipoprotein  adhesion
%     MG_386  P200                     Trm Org lipoprotein  adhesion
%     MG_412  N-ALP1                   lipoprotein          activates NF-kB as TLR 1/2 ligand
%
%   Diagnosis
%   ++++++++++++++++++++++
%   M. genitalium is typically not directly diagnosed, but rather is
%   diagnosed by negative results for other urogenital pathogen. Nucleic
%   acid amplification testing (NAAT) is commonly used to rule out
%   chlamydia infection, and suggest M. genitalium infection. Less commonly
%   NAAT is used directly to diagnose M. genitalium infection.
%
%   Treatment
%   ++++++++++++++++++++++
%   M. genitalium is most commonly treated using macrolides such as
%   azithromycin and tetracycline [PUB_0840]. Azithromycin-resistant M.
%   genitalium may be treated with fluoroquinolones [PUB_0840].
%   Taylor-Robinson and Bebear have catalogued the susceptibility of M.
%   genitalium to 18 antibiotics [PUB_0846].
%
%   Tetracycline and aminoglycoside antibiotics inhibit the 30S ribosomal
%   particle. Macrolides and lincosamide antibiotics inhibit the 50S
%   ribosomal particle. Fluoroquinolone antibiotics inhibit topoisomerase
%   II/DNA gyrase and topoisomerase IV.
%
%   Knowledge Base
%   ======================
%   The antibiotic susceptibilities of M. genitalium proteins are encoded
%   as protein activation rules in the knowledge base, and are simulated by
%   the ProteinActivation class. 32 antibiotics are represented as
%   metabolites in the knowledge base.
%
%   Representation
%   ======================
%   The state of the M. genitalium proteins involved in interaction with
%   the mammalian urogenital epithelium is represented by the enzymes
%   property. This property contains the current copy number of each
%   protein.
%
%   The state of the host mammalian urogenital epithelium is represented by
%   4 boolean properties in the host class:
%   - isBacteriumAdherent: indicates whether or not cell is adherent to the epithelial surface
%   - isTLRActivated: indicates the state of host TLRs 1, 2, and 6
%   - isNFkBActivated: indicates the state of host NF-kB
%   - isInflammatoryResponseActivated: indicates the state of the host inflammatory response
%
%   Initialization
%   ======================
%   The state of the M. genitalium host interaction proteins is initialized
%   by the initializeState method of the Simulation class. All proteins are
%   initialized to the mature state and localization.
%
%   The state of the host urogenital epithelium is initialized using the
%   same Boolean rules as in the evolveState method. That is we assume that
%   the urogenital epithelium is in a steady-state with respect to its
%   interaction with M. genitalium.
%
%   Simulation
%   ======================
%   Very little quantitative is available for the interaction between M.
%   genitalium and its mammalian host. The little data available is
%   patchwork and qualitative. Thus the model here for the interaction of
%   M. genitalium with its host is necessarily qualitative.
%
%   The model consists of several Boolean rules which calculate the value
%   of the 4 properties of the Host class:
%   1. M. genitalium is adherent if the terminal organelle is properly
%      formed and all the adhesion proteins are expressed.
%   2. TLRs 1, 2, and 6 are activated if M. genitalium is adherent (see
%      rule 1) and the M. genitalium lipoprotein ligands are expressed.
%   3. NF-kB is activated if TLRs 1 and 2 or TLRs 2 and 6 are activated.
%   4. The host inflammatory response is activated if NF-kB is activated or
%      if MG_075 is present.
%
%   The model assumes that the mammalian urogenital epithelium reaches a
%   steady-state quickly with respect to its interaction with M.
%   genitalium.
%
%   References
%   ======================
%   1. Taylor-Robinson D, Lamont RF (2011). Mycoplasmas in pregnancy. BJOG.
%      118(2):164-74. [PUB_0840]
%   2. Larsen B, Hwang J (2010). Mycoplasma, Ureaplasma, and adverse
%      pregnancy outcomes: a fresh look. Infect Dis Obstet Gynecol. 521921.
%      [PUB_0841]
%   3. Razin S, Yogev D, Naot Y (1998). Molecular biology and pathogenicity
%      of mycoplasmas. Microbiol Mol Biol Rev. 62(4):1094-156. [PUB_0313]
%   4. Razin S, Jacobs E (1992). Mycoplasma adhesion. J Gen Microbiol.
%      138(3): 407-22. [PUB_0088]
%   5. McGowin CL, Ma L, Martin DH, Pyles RB (2009). Mycoplasma
%      genitalium-encoded MG309 activates NF-kappaB via Toll-like receptors
%      2 and 6 to elicit proinflammatory cytokine secretion from human
%      genital epithelial cells. Infect Immun. 77(3): 1175-81. [PUB_0194]
%   6. Shimizu T, Kida Y, Kuwano K (2008). A triacylated lipoprotein from
%      Mycoplasma genitalium activates NF-kappaB through Toll-like receptor
%      1 (TLR1) and TLR2. Infect Immun. 76(8): 3672-8. [PUB_0842]
%   7. Shimizu T, Kida Y, Kuwano K (2007). Triacylated lipoproteins derived
%      from Mycoplasma pneumoniae activate nuclear factor-kappaB through
%      toll-like receptors 1 and 2. Immology. 121(4): 473-83. [PUB_0844]
%   8. Duffy MF, Walker ID, Browning GF (1997). The immunoreactive 116 kDa
%      surface protein of Mycoplasma pneumoniae is encoded in an operon.
%      Microbiology. 143: 3391-402. [PUB_0176]
%   9. Pich OQ, Burgos R, Ferrer-Navarro M, Querol E, Pinol J (2006).
%      Mycoplasma genitalium mg200 and mg386 genes are involved in gliding
%      motility but not in cytadherence. Mol Microbiol. 60(6):1509-19.
%      [PUB_0813]
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/27/2011
classdef  HostInteraction < edu.stanford.covert.cell.sim.Process
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {}; %names of fixed constant properties
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs   = {}; %whole cell model IDs of stimuli
        substrateWholeCellModelIDs = {}; %whole cell model IDs of substrates
        
        enzymeWholeCellModelIDs = {      %whole cell model IDs of enzymes
            'MG_075_MONOMER'
            'MG_149_MONOMER' %N-ALP2, P18, P20
            'MG_191_MONOMER' %MgPa, P140
            'MG_192_MONOMER' %P110, P40, P90, pre-B/C
            'MG_200_MONOMER'
            'MG_217_MONOMER' %P65
            'MG_218_MONOMER' %HMW2
            'MG_288_MONOMER'
            'MG_309_MONOMER'
            'MG_312_MONOMER' %HMW1, P65
            'MG_317_MONOMER' %HMW3, P200, P69
            'MG_318_MONOMER' %P30, P32
            'MG_386_MONOMER' %P200
            'MG_412_MONOMER' %N-ALP1
            'MG_410_411_412_PENTAMER' %phosphonate ABC transporter
            };
        enzymeIndexs_terminalOrganelle = [3; 4; 6; 7; 10; 11; 12; 13];
        enzymeIndexs_terminalOrganelleCore = [7; 10; 11];
        enzymeIndexs_terminalOrganelleLiprotein = [3; 4; 5; 12; 13];
        enzymeIndexs_liprotein = [1; 2; 3; 4; 5; 9; 12; 13; 14; 15];
        enzymeIndexs_adhesin = [3; 4; 6; 12];
        enzymeIndexs_antigen = [1; 2; 9; 14; 15];
        enzymeIndexs_tlrLigand = [2; 9; 14; 15];
        enzymeIndexs_tlr12Ligand = [2; 14; 15];
        enzymeIndexs_tlr26Ligand = 9;
    end
    
    %state references
    properties
        host
    end
    
    %constructor
    methods
        function this = HostInteraction(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            
            this.host = simulation.state('Host');
            
            this.states = [
                this.states;
                {this.host}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end
        
        %initialization: chromosome state initializes chromosome to state with 1
        %chromosome (which trivially is not segregated)
        function initializeState(this)
            this.evolveState();
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end
        
        %simulation
        function evolveState(this)
            h = this.host;
            
            %adherence requires
            %- properly assembled terminal organelle
            %- all adhesins
            h.isBacteriumAdherent = ...
                all(this.enzymes(this.enzymeIndexs_terminalOrganelle)) && ...
                all(this.enzymes(this.enzymeIndexs_adhesin));
            
            %TLR activation requires ligands
            h.isTLRActivated(h.tlrIndexs_1) = ...
                h.isBacteriumAdherent && ...
                any(this.enzymes(this.enzymeIndexs_tlr12Ligand));
            h.isTLRActivated(h.tlrIndexs_2) = ...
                h.isBacteriumAdherent && (...
                any(this.enzymes(this.enzymeIndexs_tlr12Ligand)) || ...
                any(this.enzymes(this.enzymeIndexs_tlr26Ligand)));
            h.isTLRActivated(h.tlrIndexs_6) = ...
                h.isBacteriumAdherent && ...
                any(this.enzymes(this.enzymeIndexs_tlr26Ligand));
            
            %NFkB activation requires TLRs to be activated by ligands
            h.isNFkBActivated = ...
                (h.isTLRActivated(h.tlrIndexs_2) && h.isTLRActivated(h.tlrIndexs_1)) || ...
                (h.isTLRActivated(h.tlrIndexs_2) && h.isTLRActivated(h.tlrIndexs_6));
            
            %inflammatory response requires NF-kB activation or antigen
            %with unknock mechanism of action
            h.isInflammatoryResponseActivated = ...
                h.isNFkBActivated || ...
                (h.isBacteriumAdherent && any(this.enzymes(this.enzymeIndexs_antigen)));
        end
    end
end
