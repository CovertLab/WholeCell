%Cytokinesis
%
% @wholeCellModelID Process_Cytokinesis
% @name             Cytokinesis
% @description
%
%   Cytokinesis is the division of the cytoplasm, achieved by the pinching of
%   the cell membrane. At a high level, it occurs by the formation of a
%   contractile Z ring along the interior surface of the cell membrane. The ring
%   is comprised of GTP-activated ftsZ polymer filaments that bend when their
%   GTP is hydrolized. Cytokinesis is thus a cycle of filament binding, bending,
%   and dissociation.
%
%   This process is based on the model described in Li et al. 2007.
%     1. Filaments bind the membrane end-to-end to form a regular polygon
%        (inscribed in the pinched cell circumference). Two filaments bind at
%        each polygon edge.
%     2. If there was a previous cycle, then when at least one filament has
%        bound each polygon edge, the last remaining ring of bent filaments from
%        the previous cycle can begin to dissociate.
%     3. When all polygon edges have two filaments bound and all residual bent
%        filaments have dissociated, the bending of the edges of the newly
%        completed polygon may begin. For simplicity, all of the GTP in the pair
%        of filaments at a particular edge hydrolyze at the same time.
%     4. When all the filaments have been bent, the bent filaments can begin
%        dissociating -- but only one from each polygon edge. The other ring
%        must remain to maintain the new smaller pinched circumference.
%     5. When only one ring of bent filaments remains, the cycle repeats.
%     6. Cytokinesis concludes when the pinched diameter is smaller than the
%        length of one filament.
%
%   How bending works:
%    - Filaments are bound in the straight configuration, forming a polygon
%      inscribed in the cell circumference.
%    - When the filaments bend, their length does not change. They bend just
%      enough so that when they've all bent, a new circle is formed. Its
%      circumference equals the old polygon's perimeter. Each fragment is now
%      an arc.
%
%   In this version of the model, all filaments are of a fixed length. When they
%   are joined end-to-end to form a regular polygon, the polygon does not quite
%   fully inscribe the entire circumference. This remaining portion of the
%   circumference does not bent when the polygon does. It is preserved and
%   accounted for in the next iteration. For more details, see Li et al. 2007.
%
%   Note: Lluch-Senar et al have shown that M. genitalium cell division can occur 
%   in the absence of FtsZ. Lluch-Senar showed that division occurs in ftsZ knockouts
%   through motility.      
%
%   References
%   ==========
%   1. Lluch-Senar M, Querol E, Pinol J (2010). Cell division in a minimal bacterium 
%      in the absence of ftsZ. Mol Microbiol. 78(2): 278-89. [PUB_0796]
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/7/2010

classdef Cytokinesis < edu.stanford.covert.cell.sim.Process & edu.stanford.covert.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'rateFtsZGtpHydrolysis';    
            'rateFilamentBindingMembrane';
            'rateFilamentDissociation';
			};
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};           %whole cell model IDs of stimuli

        substrateWholeCellModelIDs = {           %whole cell model IDs of substrates
            'PI';'H2O';'H'};
        substrateIndexs_phosphate   = 1;
        substrateIndexs_water       = 2;
        substrateIndexs_hydrogen    = 3;

        enzymeWholeCellModelIDs = {              %whole cell model IDs of enzymes
            'MG_224_9MER_GTP';      %cell division protein ftsZ polymer activated
            'MG_224_9MER_GDP';      %cell division protein ftsZ polymer deactivated and bound to cell membrane
            'MG_224_MONOMER_GDP';   %cell division protein ftsZ deactivated (monomer)
            'MG_224_MONOMER_GTP'};  %cell division protein ftsZ activated (monomer)
        enzymeIndexs_ftsZ_GTP_polymer = 1;
        enzymeIndexs_ftsZ_GDP_polymer = 2;
        enzymeIndexs_ftsZ_GDP         = 3;
        enzymeIndexs_ftsZ_GTP         = 4;
    end

    %fixed biological constants
    properties        
        rateFtsZGtpHydrolysis           %FtsZ-GTP filament hydrolysis rate (1/s) [0.15 1/s; Huecas 2007]
        rateFilamentBindingMembrane     %FtsZ-GTP filament binding rate (1/s) [0.7; not measured, fit by sensitivity analysis]
        rateFilamentDissociation        %FtsZ-GDP filament dissociation rate (1/s) [0.7; not measured, fit by sensitivity analysis]
    end

    %global state (referenced locally for convenience)
    properties
        ftsZRing
    end

    %constructor
    methods
        function this = Cytokinesis(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    methods
        %retrieve references to state objects from simulation
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.stanford.covert.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.ftsZRing = simulation.state('FtsZRing');
            this.states = [this.states; {this.ftsZRing}];
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, varargin)
            this.initializeConstants@edu.stanford.covert.cell.sim.Process(varargin{:});
            this.initializeConstants@edu.stanford.covert.cell.sim.ChromosomeProcessAspect(varargin{:});
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
            
            %substrate and byproducts: FtsZ polymerization uses GTP to make
            %FtsZ-GTP filaments for cytokinesis. Hydrolysis of the bound GTP is
            %coupled to contraction in cytokinesis. Here we account for the free
            %metabolites involved in the hydrolysis of bound GTP to GDP. GTP and
            %GDP are accounted for in FtsZPolymerization.
            [~, gtpCost] = this.calcRequiredPinchingCycles(this.geometry.width, this.ftsZRing.numFtsZSubunitsPerFilament, this.ftsZRing.filamentLengthInNm);
            
            bmProd(this.substrateIndexs_water)     = gtpCost;
            byProd(this.substrateIndexs_phosphate) = gtpCost;
            byProd(this.substrateIndexs_hydrogen)  = gtpCost;
            
            %enzymes: three Z-ring filament equivalents are needed
            minEnzExp(this.enzymeIndexs_ftsZ_GTP_polymer) = ...
                4 * this.ftsZRing.calcNumEdges(this.geometry.width, this.ftsZRing.filamentLengthInNm);
        end

        %initialization
        function initializeState(~)
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_water) = ...
                this.ftsZRing.numFtsZSubunitsPerFilament * this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer);
        end

        %simulation
        function evolveState(this)
            if ~this.chromosome.segregated
                return
            end
            
            ring = this.ftsZRing;

            %% Bind two straight filaments at each polygon edge.
            if ring.numEdgesTwoBent == 0 && ~this.geometry.pinched
                for j = 1:2  %either filament can bind an edge first
                    for i = 1:(ring.numEdges-ring.numEdgesOneStraight-ring.numEdgesTwoStraight)
                        if this.randStream.rand() <= this.rateFilamentBindingMembrane && ...
                           this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer) >= 1
                            this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer)      = this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer)      - 1;
                            this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) = this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) + 1;
                            ring.numEdgesOneStraight = ring.numEdgesOneStraight + 1;
                        end
                    end
                end
                for i = 1:ring.numEdgesOneStraight
                    if this.randStream.rand() <= this.rateFilamentBindingMembrane && ...
                       this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer) >= 1
                        this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer)      = this.enzymes(this.enzymeIndexs_ftsZ_GTP_polymer)      - 1;
                        this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) = this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) + 1;
                        ring.numEdgesOneStraight = ring.numEdgesOneStraight - 1;
                        ring.numEdgesTwoStraight = ring.numEdgesTwoStraight + 1;
                    end
                end
            end

            %% Unbind any residual bent ftsZ filaments when new ring has formed.
            if ring.numEdgesOneStraight + ring.numEdgesTwoStraight == ring.numEdges || this.geometry.pinched
                for i = 1:ring.numResidualBent
                    if this.randStream.rand() <= this.rateFilamentDissociation
                        ring.numResidualBent = ring.numResidualBent - 1;
                        
                        %the filaments fall off as monomers
                        this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) = ...
                            this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) - 1;
                        this.enzymes(this.enzymeIndexs_ftsZ_GDP) = ...
                            this.enzymes(this.enzymeIndexs_ftsZ_GDP) + ...
                            ring.numFtsZSubunitsPerFilament;
                    end
                end
            end

            %% Hydrolize GTP to bend filaments. Update diameter when done.
            if ring.numEdgesTwoBent + ring.numEdgesTwoStraight == ring.numEdges && ...
               ring.numEdgesTwoStraight > 0 && ring.numResidualBent == 0
                for i = 1:ring.numEdgesTwoStraight
                    n = 2 * ring.numFtsZSubunitsPerFilament;
                    if this.randStream.rand() <= this.rateFtsZGtpHydrolysis && this.substrates(this.substrateIndexs_water) >= n
                        ring.numEdgesTwoStraight = ring.numEdgesTwoStraight - 1;
                        ring.numEdgesTwoBent = ring.numEdgesTwoBent + 1;
                        
                        %account for hydrolysis: H2O, PI, H
                        this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - n;
                        this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + n;
                        this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + n;
                        
                        %bound polymers are now ftsZ-GDP
                        this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) = this.boundEnzymes(this.enzymeIndexs_ftsZ_GTP_polymer) - 2;
                        this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) = this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) + 2;
                    end
                end
                if ring.numEdgesTwoStraight == 0 && ~this.geometry.pinched
                    this.geometry.pinchedDiameter = this.calcNextPinchedDiameter(this.geometry.pinchedDiameter, this.ftsZRing.filamentLengthInNm);
                end
            end

            %% The first of the two bent ftsZ rings falls off.
            if ring.numEdgesTwoStraight == 0
                for i = 1:ring.numEdgesTwoBent
                    if this.randStream.rand() <= this.rateFilamentDissociation
                        ring.numEdgesTwoBent = ring.numEdgesTwoBent - 1;
                        ring.numResidualBent = ring.numResidualBent + 1;
                        
                        %filaments fall off depolymerized
                        this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) = ...
                            this.boundEnzymes(this.enzymeIndexs_ftsZ_GDP_polymer) - 1;
                        this.enzymes(this.enzymeIndexs_ftsZ_GDP) = ...
                            this.enzymes(this.enzymeIndexs_ftsZ_GDP) + ...
                            ring.numFtsZSubunitsPerFilament;
                    end
                end
            end
        end
    end

    %model helpers
    methods (Static)
        function result = calcNextPinchedDiameter(pinchedDiameter, filamentLengthInNm)
            import edu.stanford.covert.cell.sim.state.FtsZRing;
            
            numEdges = FtsZRing.calcNumEdges(pinchedDiameter, filamentLengthInNm);
            
            %the diameter of the circle that circumscribes the regular polygon
            %having numEdges sides
            flooredDiameter = filamentLengthInNm * 1e-9 / sin(pi/numEdges);
            
            %the diameter of the circle whose circumference consists of the
            %bound, bent filaments
            newDiameter = filamentLengthInNm * 1e-9 * numEdges / pi;
            
            %the parenthesized expression makes up for floor error in numEdges
            result = newDiameter + (pinchedDiameter - flooredDiameter);
            
            %if no additional pinching necessary set pinched diameter to zero
            if result <= filamentLengthInNm * 1e-9
                result = 0;
            end
        end
        
        %calculate GTP used by FtsZ polymerization to make FtsZ-GTP
        %filaments for cytokinesis and account for hydrolysis of FtsZ-GTP
        %filaments during contraction in cytokinesis
        function [nCycles, nGTP] = calcRequiredPinchingCycles(diameter, numFtsZSubunitsPerFilament, filamentLengthInNm)
            import edu.stanford.covert.cell.sim.process.Cytokinesis;
            import edu.stanford.covert.cell.sim.state.FtsZRing;
            
            nCycles  = 0;
            nGTP     = 0;
            while diameter > 0
                nCycles  = nCycles + 1;
                nGTP     = nGTP + 2 * numFtsZSubunitsPerFilament * FtsZRing.calcNumEdges(diameter, filamentLengthInNm);
                diameter = Cytokinesis.calcNextPinchedDiameter(diameter, filamentLengthInNm);
            end
        end
    end
end
