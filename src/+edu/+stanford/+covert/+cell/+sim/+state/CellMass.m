%CellWeight
% Calculates weight of various fractions of the simulation:
% - cell/media
% - metabolite/RNA/protein
% - cytosol/membrane/terminal organelle/extracellular
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef CellMass < edu.stanford.covert.cell.sim.CellState
    properties (Constant)
        dryWeight = 0;        %dry weight of this class' state properties
    end
    
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'fractionWetWeight';
            'initialBiomassConcentration';
            'cellInitialDryWeight';
            'initialFractionNTPsInRNAs';
            'initialFractionAAsInMonomers';
            'timeAveragedCellWeight';
            'cellInitialMassVariation';
            };
        fittedConstantNames     = {   %names of process properties that are considered fitted constants, and should be stored with the simulation as such
            'dryWeightFractionDNA'
            'dryWeightFractionRNA'
            'dryWeightFractionProtein'
            'dryWeightFractionLipid'
            'dryWeightFractionPolyamine'
            'dryWeightFractionCarbohydrate'
            'dryWeightFractionVitamin'
            'dryWeightFractionIon'
            'dryWeightFractionNucleotide'
            };
        stateNames              = {}; %names of properties which are part of the simulation's state
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'total'
            'cell'
            'cellDry'
            'media'
            'waterWt'
            'metaboliteWt'
            'dnaWt'
            'rnaWt'
            'proteinWt'
            };
    end
        
    %constants
    properties
        fractionWetWeight                %percent cell mass that is wet (ie. water) [PUB_0554 pp. 303]
        initialBiomassConcentration      %initial concentration of biomass (g/L)
        cellInitialDryWeight             %cell dry mass -- DNA, RNA, protein, membrane, etc. (g) [PUB_0553, PUB_0554, PUB_0555]
        initialFractionNTPsInRNAs        %fraction of initial NTPs that should be incorporated into RNAs
        initialFractionAAsInMonomers     %fraction of initial AAs that should be incorporated into protein monomers
        timeAveragedCellWeight           %time average of normalized cell weight, (current cell weight / one cell weight) [PUB_0094, PUB_0553, PUB_0554, PUB_0555, PUB_0556, PUB_0559, PUB_0560, PUB_0561]
        cellInitialMassVariation         %ratio of standard deviation to mean initial cell mass
        
        dryWeightFractionDNA             %DNA fraction of cell dry weight
        dryWeightFractionRNA             %RNA fraction of cell dry weight
        dryWeightFractionProtein         %Protein fraction of cell dry weight
        dryWeightFractionLipid           %Lipid fraction of cell dry weight
        dryWeightFractionPolyamine       %Polyamine fraction of cell dry weight
        dryWeightFractionCarbohydrate    %Carbohydrate fraction of cell dry weight
        dryWeightFractionVitamin         %Vitamin fraction of cell dry weight
        dryWeightFractionIon             %Ion fraction of cell dry weight
        dryWeightFractionNucleotide      %Free nucleotide fraction of cell dry weight
    end
    
    %component weight
    properties (SetAccess = protected)
        waterWt               %water weight (g)
        metaboliteWt          %metabolite weight (g)
        dnaWt                 %DNA weight (g)
        rnaWt                 %RNA weight (g)
        proteinWt             %protein weight (g)
        
        total                 %total weight of cell and media (g)
        cell                  %total weight of cell (g)
        cellDry               %total dry weight of cell (g)
        media                 %total weight of media (g)
    end
    
    %references to other parts of cell state
    properties
        compartment
        
        states
        
        metabolite
        rna
        monomer
        complex
        chromosome
    end
    
    %constructor
    methods
        function this = CellMass(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.compartment = simulation.compartment;
            
            this.states = simulation.states;
            
            this.metabolite = simulation.state('Metabolite');
            this.rna        = simulation.state('Rna');
            this.monomer    = simulation.state('ProteinMonomer');
            this.complex    = simulation.state('ProteinComplex');
            this.chromosome = simulation.state('Chromosome');
        end
    end
    
    methods
        function allocateMemory(~, ~)
        end
        
        function initialize(this)
            this.calcMass();
        end
        
        function calcMass(this)
            %fractions
            this.waterWt = this.metabolite.wetWeight;
            
            %total
            this.proteinWt = zeros(size(this.waterWt));
            this.total = this.waterWt;
            for i = 1:numel(this.states)
                dryWt = this.states{i}.dryWeight;
                if size(dryWt, 2) == 1
                    this.total(:, this.compartment.cytosolIndexs, :) = ...
                        this.total(:, this.compartment.cytosolIndexs, :) + ...
                        dryWt;
                else
                    this.total = this.total + dryWt;
                end
                if isequal(this.states{i}, this.metabolite)
                    this.metaboliteWt = dryWt;
                elseif isequal(this.states{i}, this.chromosome)
                    this.dnaWt = zeros(size(this.waterWt));
                    this.dnaWt(:, this.compartment.cytosolIndexs, :) = dryWt;
                elseif isequal(this.states{i}, this.rna)
                    this.rnaWt = dryWt;
                elseif isequal(this.states{i}, this.monomer) || isequal(this.states{i}, this.complex)
                    this.proteinWt = this.proteinWt + dryWt;
                end
            end
            
            %cell
            this.cell = this.total;
            this.cell(:, this.compartment.extracellularIndexs, :) = 0;
            
            this.cellDry = zeros(size(this.cell));
            this.cellDry(:, this.compartment.cellularIndexs, :) = ...
                + this.cell(:, this.compartment.cellularIndexs, :) ...
                - this.waterWt(:, this.compartment.cellularIndexs, :);
            
            %media
            this.media = this.total;
            this.media(:, this.compartment.cellularIndexs, :) = 0;
        end
    end
end