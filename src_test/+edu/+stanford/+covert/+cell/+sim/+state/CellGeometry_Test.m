%CellShape test cases
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef CellGeometry_Test < edu.stanford.covert.cell.sim.CellStateTestCase
    methods
        function this = CellGeometry_Test(name)
            this = this@edu.stanford.covert.cell.sim.CellStateTestCase(name);
        end
    end
    
    methods
        function setCellMass(this, mass)
            s = this.state;
            m = s.mass;
            
            m.metabolite.counts(:) = 0;
            m.rna.counts(:) = 0;
            m.monomer.counts(:) = 0;
            m.complex.counts(:) = 0;
            
            stateIDs = cellfun(@(state) state.wholeCellModelID(7:end), m.states, 'UniformOutput', false);
            rnaPol = m.states{strcmp(stateIDs, 'RNAPolymerase')};
            transcript = m.states{strcmp(stateIDs, 'Transcript')};
            rnaPol.states(:) = rnaPol.notExistValue;
            rnaPol.positionStrands(:) = 0;
            transcript.boundTranscriptionUnits(:) = 0;
            transcript.boundTranscriptProgress(:) = 0;
            transcript.boundTranscriptChromosome(:) = 0;
            transcript.abortedTranscripts = cell(0, 1);
            
            rib = m.states{strcmp(stateIDs, 'Ribosome')};
            pol = m.states{strcmp(stateIDs, 'Polypeptide')};
            rib.states(:) = rib.notExistValue;
            rib.boundMRNAs(:) = 0;
            rib.mRNAPositions(:) = 0;
            rib.tmRNAPositions(:) = 0;
            pol.boundMRNAs = rib.boundMRNAs;
            pol.nascentMonomerLengths = rib.mRNAPositions;
            pol.proteolysisTagLengths = rib.tmRNAPositions;
            pol.abortedPolypeptides = cell(0, 1);
                        
            m.chromosome.allocateMemory(1);
            
            m.calcMass();
            assertEqual(0, sum(m.cell));
            m.metabolite.counts(1, 1) = mass / m.metabolite.molecularWeights(1) * ...
                edu.stanford.covert.util.ConstantUtil.nAvogadro; %g
            m.calcMass();
            assertElementsAlmostEqual(mass, sum(m.cell));
            
            s.calculateVolume();
        end
    end
    
    methods
        % Verifies that initialize sets the cell shape properties using
        % only cellWeight and cellDensity as inputs. It currently assumes that
        % the cell begins as a sphere.  V = 4/3 * pi * r^3, SA = 4 * pi * r^2.
        function testInitialize_sphere(this)
            s = this.state;
            
            this.setCellMass(8/3 * pi);
            s.density = 2; %g/L
            s.width = 0.2;
            s.pinchedDiameter = 0.2;
            
            s.initialize();
            assertElementsAlmostEqual(4/3 * pi, s.volume); %L
            assertElementsAlmostEqual(0.2, s.width); %m
            assertElementsAlmostEqual(0, s.cylindricalLength); %m
            assertElementsAlmostEqual(0.2^2 * pi, s.surfaceArea); %m^2
            assertElementsAlmostEqual(0.2, s.totalLength); %m
        end
        
        % Verifies that Update updates the cell shape properties using
        % inputs: cellWeight, cellDensity, width, and pinchedDiameter.
        % This case covers the initial cell shape -- a sphere.
        function testUpdate_sphere(this)
            s = this.state;
            
            s.density = 2; %g/L
            this.setCellMass(8/3 * pi);
            
            s.width = 0.2; %m; just right to make a sphere
            s.pinchedDiameter = 0.2; %m; same as width -- no pinching
            
            assertElementsAlmostEqual(4/3 * pi, s.volume);
            assertElementsAlmostEqual(0.2, s.width);
            assertElementsAlmostEqual(0, s.cylindricalLength);
            assertElementsAlmostEqual(0.2^2 * pi, s.surfaceArea);
            assertElementsAlmostEqual(0.2, s.totalLength);
        end
        
        % Verifies that Update updates the cell shape properties using
        % inputs: cellWeight, cellDensity, width, and pinchedDiameter.
        % This case covers the cell after some growth -- a rod.
        function testUpdate_rod(this)
            s = this.state;
            
            s.density = 2; %g/L
            this.setCellMass((8/3 + 20) * pi);
            
            s.width = 0.2; %m
            s.pinchedDiameter = 0.2; %m; same as width -- no pinching
            
            assertElementsAlmostEqual((4/3 + 10) * pi, s.volume); %L
            assertElementsAlmostEqual(0.2, s.width); %m
            assertElementsAlmostEqual(1, s.cylindricalLength); %m
            assertElementsAlmostEqual((0.2^2 + 0.2) * pi, s.surfaceArea); %m^2
            assertElementsAlmostEqual(1.2, s.totalLength); %m
        end
        
        % Verifies that Update updates the cell shape properties using
        % inputs: cellWeight, cellDensity, width, and pinchedDiameter.
        % This case covers the dividing cell -- a rod with a septum.
        function testUpdate_rodWithSeptum(this)
            s = this.state;
            
            s.density = 2; %g/L
            this.setCellMass((8/3 + 20 - 1/6) * pi);
            
            s.width = 0.2; %m
            s.pinchedDiameter = 0.1; %m; half of cell width
            
            assertElementsAlmostEqual((4/3 + 10 - 1/12) * pi, s.volume); %L
            assertElementsAlmostEqual(0.2, s.width); %m
            assertElementsAlmostEqual(0.9107, s.cylindricalLength, 'relative', 1e-3); %m
            assertElementsAlmostEqual(0.7921, s.surfaceArea, 'relative', 1e-3); %m^2
            assertElementsAlmostEqual(s.width + s.cylindricalLength + 2*(s.width-s.pinchedDiameter)/2, s.totalLength); %m
        end
        
        function testUpdate_2spheres(this)
            s = this.state;
            
            s.density = 2; %g/L
            this.setCellMass(2 * (8/3) * pi);
            
            s.width = 0.2; %m
            s.pinchedDiameter = 0; %m; half of cell width
            
            assertElementsAlmostEqual(2 * (4/3) * pi, s.volume); %L
            assertElementsAlmostEqual(-1, s.width); %m
            assertElementsAlmostEqual(-1, s.cylindricalLength); %m
            assertElementsAlmostEqual(-1, s.surfaceArea); %m^2
            assertElementsAlmostEqual(-1, s.totalLength); %m
        end
    end
end
