%Cytokinesis process test case
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 9/8/2010
classdef Cytokinesis_Test < edu.stanford.covert.cell.sim.ProcessTestCase
    methods
        function this = Cytokinesis_Test(name)
            this = this@edu.stanford.covert.cell.sim.ProcessTestCase(name);
        end

        function loadTestFixture(this)
            this.loadTestFixture@edu.stanford.covert.cell.sim.ProcessTestCase();
            m = this.process;
            
            c = m.chromosome;
            c.segregated = true;
            
            g = m.geometry;
            g.width = 200 * 1e-9;
            g.pinchedDiameter = 200 * 1e-9;
        end
    end

    %tests
    methods
        function testNumEdges(this)
            m = this.process;
            g = m.geometry;
            r = m.ftsZRing;
            assertTrue(r.numEdges > 0);
            assertTrue(r.numEdges * r.filamentLengthInNm*1e-9 < g.pinchedCircumference);
        end

        function testBindingUnlimited(this)
             m = this.process;
             r = m.ftsZRing;
             
             m.substrates(:) = 1e6;
             m.enzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = 1e3;
             m.rateFilamentBindingMembrane = 1;  %force binding
             m.rateFtsZGtpHydrolysis = 0;

             m.evolveState();
             assertEqual(0, r.numEdgesOneStraight);
             assertEqual(r.numEdges, r.numEdgesTwoStraight);
             assertEqual(1e3-2*r.numEdges, m.enzymes(m.enzymeIndexs_ftsZ_GTP_polymer));
             assertEqual(2*r.numEdges, m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer));
        end

        function testBindingPolymerLimited(this)
             m = this.process;
             r = m.ftsZRing;
             
             m.substrates(:) = 1e6;
             m.enzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = 1;
             m.rateFilamentBindingMembrane = 1;  %force binding

             m.evolveState();
             assertEqual(1, r.numEdgesOneStraight);
             assertEqual(0, m.enzymes(m.enzymeIndexs_ftsZ_GTP_polymer));
             assertEqual(1, m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer));
        end

        function testBendingUnlimited(this)
            m = this.process;
            g = m.geometry;
            r = m.ftsZRing;
            
            r.numEdgesTwoStraight = r.numEdges;
            m.substrates(:) = 1e6;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = 2*r.numEdges;
            m.rateFtsZGtpHydrolysis = 1;  %force bending
            m.rateFilamentDissociation = 0;
            
            initialDiameter = g.pinchedDiameter;

            m.evolveState();
            assertEqual(0, r.numEdgesTwoStraight);
            assertEqual(r.numEdges, r.numEdgesTwoBent);
            assertEqual(0, r.numResidualBent);
            assertEqual(2*r.numEdges, m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer));
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer));
            assertEqual(1e6-(2*r.numEdges*r.numFtsZSubunitsPerFilament), m.substrates(m.substrateIndexs_water));
            assertEqual(1e6+(2*r.numEdges*r.numFtsZSubunitsPerFilament), m.substrates(m.substrateIndexs_phosphate));
            assertEqual(1e6+(2*r.numEdges*r.numFtsZSubunitsPerFilament), m.substrates(m.substrateIndexs_hydrogen));
            assertTrue(g.pinchedDiameter < initialDiameter);
        end

        function testBendingRateLimited(this)
            m = this.process;
            g = m.geometry;
            r = m.ftsZRing;
            
            r.numEdgesTwoStraight = r.numEdges;
            m.substrates(:) = 1e6;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = 30;

            initialDiameter = g.pinchedDiameter;

            m.evolveState();
            assertTrue(r.numEdgesTwoBent > 0);
            assertTrue(r.numEdgesTwoBent < r.numEdges);
            assertEqual(initialDiameter, g.pinchedDiameter);

            m.rateFtsZGtpHydrolysis = 1;  %force bending
            m.rateFilamentDissociation = 0;

            m.evolveState();
            assertEqual(0, r.numEdgesTwoStraight);
            assertEqual(r.numEdges, r.numEdgesTwoBent);
            assertTrue(g.pinchedDiameter < initialDiameter);
        end

        function testFirstRingFallingOffUnlimited(this)
            m = this.process;
            r = m.ftsZRing;
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer) = 30;
            r.numEdgesTwoBent = r.numEdges;
            m.rateFilamentDissociation = 1;  %force falling off

            m.evolveState();
            assertEqual(0, r.numEdgesTwoBent);
            assertEqual(r.numEdges, r.numResidualBent);
            assertEqual(r.numEdges, m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer));
            assertEqual(...
                r.numEdges*r.numFtsZSubunitsPerFilament,...
                m.enzymes(m.enzymeIndexs_ftsZ_GDP));
        end

        function testFirstRingFallingOffRateLimited(this)
            m = this.process;
            r = m.ftsZRing;
            
            m.substrates(:) = 0;
            m.enzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer) = 30;
            r.numEdgesTwoBent = r.numEdges;

            m.evolveState();
            assertTrue(r.numEdgesTwoBent < r.numEdges);
            assertTrue(r.numEdgesTwoBent > 0);
            assertTrue(r.numResidualBent > 0);
            assertTrue(r.numResidualBent < r.numEdges);

            m.rateFilamentDissociation = 1;  %force falling off

            m.evolveState();
            assertEqual(0, r.numEdgesTwoBent);
            assertEqual(r.numEdges, r.numResidualBent);
        end
        
        function testSecondRingFallingOffUnlimited(this)
            m = this.process;
            r = m.ftsZRing;
            
            r.numEdgesOneStraight = r.numEdges;
            r.numResidualBent = 5;
            m.enzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = r.numEdgesOneStraight;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer) = r.numResidualBent;

            m.rateFilamentDissociation = 1;

            m.evolveState();
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer));
            assertEqual(0, r.numResidualBent);
            assertEqual(...
                5 * r.numFtsZSubunitsPerFilament,...
                m.enzymes(m.enzymeIndexs_ftsZ_GDP));
        end
        
        function testSecondRingFallingOffRateLimited(this)
            m = this.process;
            r = m.ftsZRing;
            
            r.numEdgesOneStraight = r.numEdges;
            r.numResidualBent = 12;
            m.enzymes(:) = 0;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GTP_polymer) = r.numEdgesOneStraight;
            m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer) = r.numResidualBent;

            m.evolveState();
            assertTrue(r.numResidualBent < 15);
            assertTrue(r.numResidualBent > 0);

            m.rateFilamentDissociation = 1;

            m.evolveState();
            assertEqual(0, m.boundEnzymes(m.enzymeIndexs_ftsZ_GDP_polymer));
            assertEqual(0, r.numResidualBent);
            assertEqual(...
                12 * r.numFtsZSubunitsPerFilament,...
                m.enzymes(m.enzymeIndexs_ftsZ_GDP));
        end

        function testPinchCell(this)
            m = this.process;
            g = m.geometry;
            r = m.ftsZRing;
            
            m.enzymes(:) = 1e6;
            m.substrates(:) = 1e6;
            m.rateFtsZGtpHydrolysis = 1;
            m.rateFilamentBindingMembrane = 1;
            m.rateFilamentDissociation = 1;
            
            for j = 1:75
                m.evolveState();
            end
            
            assertEqual(0, g.pinchedDiameter);
            assertTrue(g.pinched);
            assertEqual(0, r.numEdges);
            assertEqual(0, r.numEdgesOneStraight);
            assertEqual(0, r.numEdgesTwoStraight);
            assertEqual(0, r.numEdgesTwoBent);
            assertEqual(0, r.numResidualBent);
        end

        function testNoEnzymes(this)
            m = this.process;
            m.enzymes(:) = 0;
            m.substrates(:) = 0;
            m.rateFtsZGtpHydrolysis = 1;
            m.rateFilamentBindingMembrane = 1;
            m.rateFilamentDissociation = 1;
            for j = 1:10
                m.evolveState();
            end
            
            g = m.geometry;
            assertEqual(g.width, g.pinchedDiameter);
        end
        
        function testNoWater(this)
            m = this.process;
            m.enzymes(:) = 1e6;
            m.substrates(:) = 1e6;
            m.substrates(m.substrateIndexs_water) = 0;
            m.rateFtsZGtpHydrolysis = 1;
            m.rateFilamentBindingMembrane = 1;
            m.rateFilamentDissociation = 1;
            for j = 1:10
                m.evolveState();
            end
            
            g = m.geometry;
            assertEqual(g.width, g.pinchedDiameter);
        end

        function testGeneEssentiality(this)
            m = this.process;
            g = m.geometry;
            g.pinchedDiameter = g.width;
            m.enzymes(:) = 1e6;
            m.substrates(:) = 1e6;
            m.rateFtsZGtpHydrolysis = 1;
            m.rateFilamentBindingMembrane = 1;
            m.rateFilamentDissociation = 1;
            
            initialPinchedDiameter = g.pinchedDiameter;

            this.helpTestGeneEssentiality(...
                {'MG_224'},... %cell division protein ftsZ
                @(m, i) m.geometry.pinchedDiameter < initialPinchedDiameter,...
                struct('lengthSec',10));
        end
    end
end
