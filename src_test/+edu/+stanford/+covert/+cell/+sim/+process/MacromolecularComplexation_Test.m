%Macromolecular complexation test cases
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef MacromolecularComplexation_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = MacromolecularComplexation_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function testNoComplexationWithoutAnySubunits(this)
            m = this.process;
            m.substrates(:) = 0;  %no subunits
            m.complexs(:) = 0;

            m.evolveState();
            assertAllEqual(0, m.complexs);
        end

        function testSubunitAndComplexQuantitiesIntegerValued(this)
            m = this.process;

            m.evolveState();
            assertAllEqual(0, rem(m.substrates, 1));
            assertAllEqual(0, rem(m.complexs, 1));
        end

        function testMultipleComplexsRequiringTheSameSubunit(this)
            this.loadSimpleTestFixture();
            m = this.process;
            nSubunits = size(m.complexComposition,1);
            nComplexs = size(m.complexComposition,2);
            nCompartments = size(m.complexComposition,3);
            m.substrates = repmat(20, [nSubunits nCompartments]);
            m.complexs = zeros(nComplexs, 1);

            m.evolveState();
            assertEqual([10 6 4 6]', m.complexs);  % some of each made
            assertEqual([10 0 0 8 12 2]', m.substrates);
        end

        function testComplexRequiresEachSubunit(this)
            [s1, s2] = this.substrateIdx('MG_0001', 'MG_048_MONOMER');
            c = this.complexIdx('MG_0001_048');
            m = this.process;

            % both subunits present
            m.substrates(:) = 0;
            m.substrates([s1 s2]) = 1;
            m.complexs(:) = 0;
            m.evolveState();
            assertEqual(1, m.complexs(c));
            assertEqual([0; 0], m.substrates([s1 s2], 1));

            % only MG_0001 present
            m.substrates(:) = 0;
            m.substrates(s1) = 1;
            m.complexs(:) = 0;
            m.evolveState();
            assertEqual(0, m.complexs(c));
            assertEqual([1; 0], m.substrates([s1 s2], 1));

            % only MG_048_MONOMER present
            m.substrates(:) = 0;
            m.substrates(s2) = 1;
            m.complexs(:) = 0;
            m.evolveState();
            assertEqual(0, m.complexs(c));
            assertEqual([0; 1], m.substrates([s1 s2], 1));
        end

        function testSubunitInCorrectCompartment(this)
            [s1, s2, s3] = this.substrateIdx(...
                'MG_467_MONOMER',...
                'MG_468_MONOMER',...
                'MG_526_MONOMER');
            c = this.complexIdx('MG_467_468_526_TETRAMER');
            m = this.process;

            % in correct compartments
            m.substrates(:) = 0;
            m.substrates(s1, 1) = 1;  % cmpt 1 is cytosol
            m.substrates(s2, 1) = 2;  % cmpt 4 is membrane
            m.substrates(s3, 1) = 1;
            m.complexs(:) = 0;
            m.evolveState();
            assertEqual(1, m.complexs(c));
            assertAllEqual(0, m.substrates);
        end

        function testGeneEssentiality(this)
            this.process.substrates(:) = 100;
            this.helpTestGeneEssentiality(...
                {}, @this.isProperlyFunctioning, struct('knockoutSubstrates',false));
        end
    end

    %helper functions
    methods (Access = private)
        function result = isProperlyFunctioning(~, m, i)
            initialMaxComplexs = ...
                double(m.complexComposition'>0) * double(i.substrates>0) ...
                == sum(m.complexComposition>0, 1)';
            result = all(initialMaxComplexs == 0 | m.complexs > i.complexs);
        end

        function varargout = substrateIdx(this, varargin)
            [~,i] = ismember(varargin, this.process.substrateWholeCellModelIDs);
            varargout = num2cell(i);
        end
        
        function varargout = complexIdx(this, varargin)
            [~,i] = ismember(varargin, this.process.complexWholeCellModelIDs);
            varargout = num2cell(i);
        end
        
        function loadSimpleTestFixture(this)
            m = this.process;

            %whole cell model IDs of subunits, macromolecular complexes
            m.substrateWholeCellModelIDs = {'substrate1','substrate2','substrate3','substrate4','substrate5','substrate6'}';
            m.complexWholeCellModelIDs   = {'complex1','complex2','complex3','complex4'}';

            %protein complex composition (monomers X complexes X compartments)
            %complex1 = substrate1[c] + (2) substrate2[m]
            %complex2 = (2) substrate3[c] + (2) substrate4[c]
            %complex3 = (2) substrate3[c] + (2) substrate5[c]
            %complex4 = (3) substrate6[c]
            nSubunits = length(m.substrateWholeCellModelIDs);
            nComplexs = length(m.complexWholeCellModelIDs);
            nCompartments = 1;

            m.complexComposition = zeros(nSubunits, nComplexs, nCompartments);
            m.complexComposition(1,1) = 1;
            m.complexComposition(2,1) = 2;
            m.complexComposition(3,2) = 2;
            m.complexComposition(4,2) = 2;
            m.complexComposition(3,3) = 2;
            m.complexComposition(5,3) = 2;
            m.complexComposition(6,4) = 3;

            [m.substrates2complexNetworks, m.complexs2complexNetworks, m.complexNetworks] = ...
                edu.stanford.covert.util.findNonInteractingRowsAndColumns(...
                reshape(permute(m.complexComposition,[2 1 3]),size(m.complexComposition,2),[])');

            %weights of substrates
            m.substrateMolecularWeights = [1 2 3 4 5 6]';
            m.complex.molecularWeights = sum(m.complexComposition,3)'*m.substrateMolecularWeights;

            m.substrates = repmat(20, nSubunits, nCompartments);
            m.complexs = zeros(nComplexs, 1);
        end
    end
end
