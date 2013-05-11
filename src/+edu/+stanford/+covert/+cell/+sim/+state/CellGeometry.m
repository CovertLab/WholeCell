%CellGeometry
%
% @wholeCellModelID State_CellGeometry
% @name             Cell geometry
% @description
%
% The Cell Shape process calculates the length and surface area of the cell
% across its lifespan. It also keeps track of the cell volume. The cell is
% approximated to be of a rod shape similar to E. coli, even though
% M. genitalium is known to have a less uniform flask shape. Initially, the
% cell is modeled as a cylinder with two hemispherical caps. Once cell
% pinching commences at the midline of the cell, the shape and size of a
% "septum region" is also modeled.
%
% General geometric equations to represent the shape of a cell, and the
% idea that the density of the cell does not change across the cell cycle
% are borrowed from Domach et al. (1983). We add the assumption that the
% width of the cell remains constant during the lifespan of a cell. Thus,
% the density and cell width are inputs into our model. The cell density of
% E. coli is used 1100g/L (Baldwin et al., 1995). The cell width is
% calculated based on the initial cell mass and density and the assumption
% that the cell is a sphere. The initial cell mass is fit to result in a
% cell width of 200nm (Lind et al, 1984).
%
% Our model calculates the mass of the cell at all timesteps. The mass and
% constant density provide us with the volume at all time points.
% Volume = [2 hemispheres] + [2 cylinders] + [septum region]
% The volume of the septum region is calculated as a cylinder of length,
% 2*septumLength, and width of the cell. Then two cones (height=septum,
% radius=septum) are subtracted from this cylinder.
% (This approximation was used in Shuler et al. 1979).
%
% Since we know the volume of the cell from the cell's mass and density,
% and the septum length from our cytokinesis process, the volume formula
% gives us the length of the cell at each timestep.
%
% Similarly, we can also calculate the surface area of the cell.
% Surface Area = [2 hemispheres] + [2 cylinders] + [septum region]
% The surface area of the septum region is calculated as a cylinder of
% length, 2*septum, and width of the cell.
%
% References
% ==================
% 1. Shuler, M.L., Leung, S., Dick, C.C. (1979). A Mathematical Model for the
%    Growth of a Single Bacteria Cell. Annals of the New York Academy of Sciences
%    326: 35-52.
% 2. Domach, M.M., Leung, S.K., Cahn, R.E., Cocks, G.G., Shuler, M.L. (1983).
%    Computer model for glucose-limited growth of a single cell of Escherichia
%    coli B/r-A. Biotechnology and Bioengineering 26: 203-216.
% 3. Lind, K., Lindhardt, B., Schutten, H.J., Blom, J., Christiansen, C.
%    (1984). Serological Cross-Reactions Between Mycoplasma genitalium and
%    Mycoplasma pneumoniae. Journal of Clinical Microbiology 20: 1036-1043.
% 4. Baldwin WW, Myer R, Powell N, Anderson E, Koch AL. (1995). Buoyant
%    density of Escherichia coli is determined solely by the osmolarity of the
%    culture medium. Arch Microbiology 164: 155-157.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/3/2010
classdef CellGeometry < edu.stanford.covert.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
		    'density';
			};
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such  
        stateNames              = {   %names of properties which are part of the simulation's state
            'width';
            'pinchedDiameter'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'volume'
            'cylindricalLength'
            'surfaceArea'
            'totalLength'
            'pinchedCircumference'
            'pinched'
            'chamberVolume'
            };
    end
    
    properties (Constant)
        dryWeight = 0;      %dry weight of this class' state properties
    end

    %fixed biological constants
    properties
        density             %g/L [PUB_0553]
    end

    %state
    properties
        width               %m; width of the cell; stays constant across simulation; calculated during initialization        
        pinchedDiameter     %m; diameter where cell is pinched the most
    end   
    
    %dependent state
    properties
        volume               %L
    end
    
    properties (Dependent = true, SetAccess = protected)
        cylindricalLength    %m; length of just the cylindrical part of the cell
        surfaceArea          %m^2; surface area of the cell
        totalLength          %m; total cell length = width + length + (2*septumLength)
        pinchedCircumference %circumference where cell is pinched the most (m)
        pinched              %whether cytokinesis is finished
        
        chamberVolume        %volume of chamber in which cell is simulated (L)
    end
    
    %references to other parts of the cell state
    properties
        mass
    end

    %constructor
    methods
        function this = CellGeometry(wholeCellModelID, name)
            this = this@edu.stanford.covert.cell.sim.CellState(wholeCellModelID, name);            
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.mass = simulation.state('Mass');
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)            
            this.width             = zeros(1, 1, numTimePoints);
            this.pinchedDiameter   = zeros(1, 1, numTimePoints);
        end
    end

    %model
    methods
        %initialize to spherical shape
        function initialize(this)
            this.calculateVolume();
            
            w = this.calculateWidth(this.volume);
            this.width = w;
            this.pinchedDiameter = w;
        end
                
        %calculates width
        function w = calculateWidth(~, volume)
            w = 2 * (volume/1000 * 3/(4*pi))^(1/3);
        end
        
        %updates volume; called by metabolism
        function calculateVolume(this)
            this.volume = sum(this.mass.cell) / this.density;
        end
    end
       
    %getters
    methods
        function value = get.width(this)
            if this.pinchedDiameter > 0
                value = this.width;
            else
                value = -1;
            end
        end       
        
        function value = get.cylindricalLength(this)
            [value, ~, ~] = this.calcGeometry();
        end
        
        function value = get.surfaceArea(this)
            [~, value, ~] = this.calcGeometry();
        end
        
        function value = get.totalLength(this)
            [~, ~, value] = this.calcGeometry();
        end
        
        function result = get.pinchedCircumference(this)
            result = this.pinchedDiameter * pi;
        end
        
        function result = get.pinched(this)
            result = this.pinchedDiameter == 0;
        end
       
        function [cylindricalLength, surfaceArea, totalLength] = calcGeometry(this)
            import edu.stanford.covert.cell.sim.state.CellGeometry;
            [cylindricalLength, surfaceArea, totalLength] = ...
                CellGeometry.calculateGeometry(this.width,this.pinchedDiameter, this.volume);
        end
        
        function value = get.chamberVolume(this)
            value = this.mass.cellInitialDryWeight / (1 - this.mass.fractionWetWeight) / ...
                this.mass.initialBiomassConcentration;
        end 
    end
    methods (Static)
        function [cylindricalLength, surfaceArea, totalLength] = calculateGeometry(width, pinchedDiameter, volume)
            if pinchedDiameter > 0
                %calculate dimensions of cell whose shape is a cyclinder with
                %hemispherical caps
                
                w = width;
                s = (w - pinchedDiameter)/2; %m
                
                cylindricalLength = ...                                             %(m)
                    max(0, ...
                    + volume / 1000 ...                                             %total volume
                    - 1/6 * pi * width^3 ...                                        %spherical caps
                    - pi/2*(s*(8*s^2-4*s*w+w^2) + 2*s^2*(-2*s+w)*pi/2 - 4/3*s^3) ... %septum
                    ) * 4/(pi * width^2);
                surfaceArea = ...                            %(m^2)
                    pi * width^2 + ...                       %spherical caps
                    pi * width * cylindricalLength + ...     %cylinder
                    pi * 4*s*(w-2*s+s);                      %septum
                totalLength = ...
                    cylindricalLength + ...
                    width + ...
                    2 * s; %m
            else
                % cell has divided
                warning('WholeCell:warning', 'Cell has divided. Cell shape undefined');
                cylindricalLength = -1;
                surfaceArea = -1;
                totalLength = -1;
            end
        end
    end
end
