%Transcriptional Regulation
%
% @wholeCellModelID Process_TranscriptionalRegulation
% @name             Transcriptional Regulation
% @description
%   Biology
%   ==============
%   The the rate of transcription of each transcription unit is known to be
%   regulated by proteins referred to as transcription factors which modulate
%   the affinity of the RNA polymerase for each promoter. Transcription factors
%   can both positively and negatively modulate RNA polymerase - promoter
%   affinity. Transcription factors can stabilize RNA polymerase - promoter
%   complexes by contributing additional negative free energy to the complex by
%   providing an additional surface for the RNA polymerase binding.
%   Transcription factors can destabilize RNA polymerase - promoter complexes
%   for example by sterically blocking promoters and preventing the RNA
%   polymerase from binding promoters.
%
%   This process simulates the binding of transcription factors to the promoters
%   of each transcription unit (affinity), as well as the affect of the
%   transcription factors on the recruiment of the RNA polymerase to each
%   promoter (activity). The process was built using experimentally observed fold
%   change expression effect of each transcription factor on each promoter.
%
%   Affinity
%   ++++++++++++++
%   Binds enzymes to promoters assuming:
%   1) Transcription factors have high affinity for promoters
%   2) Transcription factors bind promoters rapidly
%   3) Transcription factors bind promoters stably over time
%   4) Only 1 copy of a transcription factor can bind each promoter
%   5) Transcription factors only compete within their species for
%      promoters. That is for each transcription unit we assume each
%      transcription factor species binds a distinct promoters region.
%
%   Consequently, at each time step we simulate that each free
%   transcription factor binds randomly binds unoccupied promoters (no copy
%   of that transcription factor is already bound to the promoter). Random
%   transcription factor-promoter binding is weighted by the affinity of
%   each transcription factor for each promoter.
%
%   Because transcription factor-promoter affinities are generally not
%   experimentally observed, we base them on the transcription factor fold
%   change activities.
%
%   Activity
%   ++++++++++++++
%   The effect of bound transcription factors on the recruitment of RNA
%   polymerase and the expression of transcription units is simulated here
%   and incorporated into the calculation of the RNA polymerase
%   transcription unit promoter binding probabilities in the transcription
%   process. Specifically, the wild-type average RNA polymerase
%   transcription unit promoter binding probabilities are multiplied by the
%   binding probability fold change effects simulated in this process.
%
%   The RNA polymerase binding probability fold change is simulated for
%   each promoter as the product of the observed expression fold change
%   effects of each bound transcription factor. When a promoter is bound by
%   a single transcription factor, the net RNA polymerase binding
%   probability fold change is the observed expression fold change of that
%   transcription factor. When a promoter is bound by multiple
%   transcription units, the net RNA polymerase binding probability fold
%   change is given by the product of the individual fold change effects of
%   the bound transcription factors.
%
%   Knowledge Base
%   ==============
%   The list of transcriptional regulatory relationships is maintained in the
%   the knowledge base. As of 8/10/2010, it contained 31 such relationships
%   between 5 transcription factors and 29 transcription units containing 37
%   genes. The knowledge base was built from a variety of literature sources
%   and databases including:
%   - PUB_0096
%   - PUB_0110
%   - PUB_0112
%   - PUB_0196
%   - PUB_0418-20
%   - PUB_0433-8
%   - PUB_0505
%
%   Initialization
%   =================
%   Because we assume that transcription factors have high affinity for DNA and
%   bind DNA stably, we initialize as many transcription factors as possible to
%   the promoter-bound state. We randomly assign transcription factors to
%   promoters using their relative affinities.
%
%   Simulation
%   ==============
%   For each kind of transcription factor, bind any free transcription factors
%   to promoters that aren't already occupied by this kind of transcription
%   factor. Choose the binding sites randomly, weighted by the transcription
%   factor's affinity to them.
%
%   References
%   ==============
%   1) Lacramioara Bintu and Nicolas E Buchler and Hernan G Garcia and
%      Ulrich Gerland and Terence Hwa, Jane Kondev and Rob Phillips (2005).
%      Transcriptional regulation by the numbers: models. Curr Opin Genet
%      Dev. 15:116-24.
%   2) Lacramioara Bintu and Nicolas E Buchler and Hernan G Garcia and
%      Ulrich Gerland and Terence Hwa, Jane Kondev and Rob Phillips (2005).
%      Transcriptional regulation by the numbers: applications. Curr Opin
%      Genet Dev. 15:125-35.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/19/2010
classdef TranscriptionalRegulation < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'tfIndexs';
            'tuIndexs';
            'tfAffinities';
            'tfActivities';
            'tfPositionStrands';
            'otherActivities';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};
        substrateWholeCellModelIDs = {};

        enzymeWholeCellModelIDs
        enzymeIndexs_fur         %MG_236_MONOMER  ferric uptake repressor
        enzymeIndexs_gntR        %MG_101_MONOMER  Uncharacterized HTH-type transcriptional regulator
        enzymeIndexs_hrcA        %MG_205_DIMER    heat-inducible transcription repressor HrcA, putative
        enzymeIndexs_luxR        %MG_428_DIMER    LuxR bacterial regulatory protein, putative
        enzymeIndexs_spx         %MG_127_MONOMER  Spx subfamily protein

        transcriptionUnitWholeCellModelIDs
    end

    %fixed biological constants
    properties
        tfIndexs           %identifying index of transcription factor that binds each site [sites X 2]
        tuIndexs           %identifying index of transcription unit associated with each TF binding site [sites X 2]
        tfAffinities       %affinity of each transcription factor for each promoter it binds [sites X 2]
        tfActivities       %fold change expression of each transcription unit when its promoter is bound by a particular transcription factor [sites X 2]
        tfPositionStrands  %start position and strand of each transcription factor DNA binding site [2*sites X 2]
        otherActivities    %fold change each non-DNA-binding transcription factor has on each transcription unit's expression [nTFs X nTUs]
    end

    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent, SetAccess = protected)
        tfBoundPromoters       %number of each transcription factor bound to promoter of each transcription unit [transcription factors X promoters X 2]
        boundTFs               %count of how many of each transcription factor are bound [transcription factors X 1]
    end
    
    %references to cell state
    properties
        rnaPolymerase             %reference to RNA polymerase state
    end

    %constructor
    methods
        function this = TranscriptionalRegulation(wholeCellModelID, name)
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
        function initializeConstants(this, knowledgeBase, varargin)
            %build object arrays of transcription factor monomers and complexes
            monomers = [];
            complexs = [];
            for i = 1:knowledgeBase.numTranscriptionUnits
                m = knowledgeBase.transcriptionUnits(i).transcriptionFactorProteinMonomers;
                c = knowledgeBase.transcriptionUnits(i).transcriptionFactorProteinComplexs;
                if ~isempty(m)
                    monomers = [monomers; m]; %#ok<AGROW>
                end
                if ~isempty(c)
                    complexs = [complexs; c]; %#ok<AGROW>
                end
            end
            
            %whole cell model ids of all transcription factors
            this.enzymeWholeCellModelIDs = unique({
                monomers.wholeCellModelID...
                complexs.wholeCellModelID}');

            %transcription factors indices
            this.enzymeIndexs_fur  = this.enzymeIndexs({'MG_236_MONOMER'});
            this.enzymeIndexs_gntR = this.enzymeIndexs({'MG_101_MONOMER'});
            this.enzymeIndexs_hrcA = this.enzymeIndexs({'MG_205_DIMER'});
            this.enzymeIndexs_luxR = this.enzymeIndexs({'MG_428_DIMER'});
            this.enzymeIndexs_spx  = this.enzymeIndexs({'MG_127_MONOMER'});

            %call super class method to compute mapping between enzymes of
            %this process, and of protein monomers and complexes of the simulation
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(...
                knowledgeBase, varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, varargin{:});

            %whole cell model ids of transcription units
            this.transcriptionUnitWholeCellModelIDs = ...
                {knowledgeBase.transcriptionUnits.wholeCellModelID}';

            %build matrices of transcription factor positions, affinities, activities
            n = length(monomers) + length(complexs);
            this.tuIndexs = zeros(n, 2);
            this.tfIndexs = zeros(n, 2);
            this.tfActivities = zeros(n, 2);
            this.tfPositionStrands = zeros(n*2, 2);
            c = this.chromosome;
            j = 1;
            for i = 1:knowledgeBase.numTranscriptionUnits
                tu = knowledgeBase.transcriptionUnits(i);

                if ~isempty(tu.transcriptionFactorProteinMonomers)
                    [tf, proteinIndexs] = ismember(...
                        [tu.transcriptionFactorProteinMonomers.idx]',...
                        this.enzymeMonomerGlobalIndexs);
                    k = nnz(tf) - 1;
                    this.tuIndexs(j:j+k) = i;

                    proteinIndexs = proteinIndexs(tf);
                    this.tfIndexs(j:j+k) = this.enzymeMonomerLocalIndexs(proteinIndexs);

                    regulationIndexs = ...
                        this.enzymeMonomerCompartmentIndexs(proteinIndexs) == ...
                        [tu.transcriptionFactorProteinMonomerCompartments.idx]';
                    this.tfActivities(j:j+k) = ...
                        tu.transcriptionFactorProteinMonomerActivitys(regulationIndexs);

                    footprints = c.monomerDNAFootprints([tu.transcriptionFactorProteinMonomers.idx]');
                    startCoordinates = tu.transcriptionFactorProteinMonomerBindingSiteStartCoordinates;
                    lengths = tu.transcriptionFactorProteinMonomerBindingSiteLengths;
                    this.tfPositionStrands(j:j+k) = ...
                        mod(ceil(startCoordinates + lengths/2 - footprints/2) - 1, c.sequenceLen) + 1;
                    j = j + 1 + k;
                end

                if ~isempty(tu.transcriptionFactorProteinComplexs)
                    [tf, proteinIndexs] = ismember(...
                        [tu.transcriptionFactorProteinComplexs.idx]',...
                        this.enzymeComplexGlobalIndexs);
                    k = nnz(tf) - 1;
                    this.tuIndexs(j:j+k) = i;

                    proteinIndexs = proteinIndexs(tf);
                    this.tfIndexs(j:j+k) = this.enzymeComplexLocalIndexs(proteinIndexs);
                   
                    regulationIndexs = ...
                        this.enzymeComplexCompartmentIndexs(proteinIndexs) == ...
                        [tu.transcriptionFactorProteinComplexCompartments.idx]';
                    this.tfActivities(j:j+k) = ...
                        tu.transcriptionFactorProteinComplexActivitys(regulationIndexs);

                    footprints = c.complexDNAFootprints([tu.transcriptionFactorProteinComplexs.idx]');
                    startCoordinates = tu.transcriptionFactorProteinComplexBindingSiteStartCoordinates;
                    lengths = tu.transcriptionFactorProteinComplexBindingSiteLengths;
                    this.tfPositionStrands(j:j+k) = ...
                        mod(ceil(startCoordinates + lengths/2 - footprints/2) - 1, c.sequenceLen) + 1;
                    j = j + 1 + k;
                end
            end

            validateattributes(this.tfActivities, {'numeric'}, {'nonnegative'});
            
            this.tuIndexs(:,2) = this.tuIndexs(:,1);
            this.tfIndexs(:,2) = this.tfIndexs(:,1);
            this.tfActivities(:,2) = this.tfActivities(:,1);

            this.tfPositionStrands(1:n,2) = 1;
            this.tfPositionStrands(n+1:2*n,1) = this.tfPositionStrands(1:n,1);
            this.tfPositionStrands(n+1:2*n,2) = 3;

            this.tfAffinities = this.tfActivities;
            this.tfActivities(this.tfActivities<1 & this.tfAffinities>0) = ...
                this.tfActivities(this.tfActivities<1 & this.tfAffinities>0) .^ -1;
 
            bind = ~isnan(this.tfPositionStrands(1:n,1));
            this.otherActivities = ones(...
                length(this.transcriptionUnitWholeCellModelIDs), ...
                length(this.enzymeWholeCellModelIDs));
            this.otherActivities(sub2ind(...
                size(this.otherActivities), ...
                this.tuIndexs(~bind), ...
                this.tfIndexs(~bind))) = ...
                this.tfActivities(~bind);
            this.tuIndexs = this.tuIndexs(bind,:);
            this.tfIndexs = this.tfIndexs(bind,:);
            this.tfAffinities = this.tfAffinities(bind,:);
            this.tfActivities = this.tfActivities(bind,:);
            this.tfPositionStrands = this.tfPositionStrands(repmat(bind, [2 1]), :);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %no substrate and byproducts required
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            
            %no enzymes required
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end

        %initialization:
        %- Simulation initializeState method initializes 1 undamaged chromosome
        %- Various processes may bind proteins to chromosome
        %- Here we bind transcription factors to the initial chromosome
        function initializeState(this)
            this.evolveState();
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end

        %simulation
        function evolveState(this)
            %Bind transcription factors to accessible sites
            tfBoundPromoters = this.bindTranscriptionFactors();
            
            %Compute new fold changes
            this.rnaPolymerase.transcriptionFactorBindingProbFoldChange = ...
                this.calcBindingProbabilityFoldChange(tfBoundPromoters);
        end
        
        %Bind transcription factors to accessible sites.
        function tfBoundPromoters = bindTranscriptionFactors(this)
            %Estimate which binding sites exist and are accessible to the
            %transcription unit. (Note chromosome class which ensure that sites
            %are truly accessible, and only allow transcription factors to bind
            %to these sites).
            tfBoundPromoters = this.tfBoundPromoters;
            accessible = ~tfBoundPromoters & reshape(this.chromosome.isRegionPolymerized(this.tfPositionStrands, 1, false), [], 2);
            
            for i = 1:size(this.enzymes, 1)
                %For each transcription factor:
                %  Stochastically pick transcription unit promoters elements to
                %  bind from among free promoter elements weighted by the
                %  transcription factors' affinities for the promoters. (Note:
                %  each transcription factor is assumed to bind a distinct
                %  promoter region near each transcription unit such that
                %  transcription factors only compete within their species
                %  for binding to promoters.)
                
                if this.enzymes(i) <= 0
                    continue;
                end
                
                sites = find((this.tfIndexs == i) & accessible);
                if isempty(sites)
                    continue;
                end
                
                tf = this.bindProteinToChromosome(...
                    this.tfPositionStrands(sites, :), i, this.enzymes(i), this.tfAffinities(sites), [], false, 1, false, []);
                
                tfBoundPromoters(sites(tf)) = true;
            end
        end
        
        function foldChange = calcBindingProbabilityFoldChange(this, tfBoundPromoters)
            %reset fold change
            foldChange = ones(size(this.chromosome.transcriptionUnitStartCoordinates, 1), 2);
            
            %fold change effect of bound TFs
            for i = 1:size(this.tfIndexs, 1)
                for j = 1:size(this.tfIndexs, 2)
                    if ~tfBoundPromoters(i, j)
                        continue;
                    end
                
                    foldChange(this.tuIndexs(i, j), j) = foldChange(this.tuIndexs(i, j), j) * this.tfActivities(i, j);
                end
            end
            
            otherFoldChanges = prod(this.otherActivities .^ (this.enzymes(:, ones(1, size(foldChange, 1)))' > 0), 2);
            foldChange = foldChange .* [otherFoldChanges otherFoldChanges];
        end
    end
    
    %get methods of dependent local state
    methods
        function result = get.tfBoundPromoters(this)
            result = reshape(this.isDnaBound(this.tfPositionStrands, this.tfIndexs(:)), [], 2);
        end
        
        function result = get.boundTFs(this)
            result = histc(this.tfIndexs(logical(this.tfBoundPromoters(:,1))), 1:numel(this.enzymes));
        end
    end
end
