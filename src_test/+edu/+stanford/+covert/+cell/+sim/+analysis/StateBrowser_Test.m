%StateBrowser_Test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 5/9/2011
classdef StateBrowser_Test < TestCase
    properties
        browser
    end
    
    methods
        function this = StateBrowser_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            import edu.stanford.covert.cell.sim.analysis.StateBrowser;
            import edu.stanford.covert.cell.sim.SimulationFixture;
            import edu.stanford.covert.cell.sim.util.DiskLogger;
            
            simDir = edu.stanford.covert.cell.sim.util.SimulationDiskUtil.getLatestSimulation();
            metadata = DiskLogger.loadMetadata(simDir);
            simLength = metadata.lengthSec;
            segmentSizeStep = metadata.segmentSizeStep;
            
            plotConfiguration = [
                struct(...
                'state', 'Metabolite - counts', ...
                'stateSubset', []);
                struct(...
                'state', 'Geometry - width', ...
                'stateSubset', [])];
            plotConfiguration(1).stateSubset = {'ATP'; 'CTP'; 'GTP'; 'UTP'};
            this.browser = StateBrowser(simDir, plotConfiguration, 1, simLength, segmentSizeStep, 0);
        end
        
        function tearDown(this)
            this.browser.close();
        end
    end
    
    
    methods
        function testBrowser(this)
            this.browser.open();
            this.browser.plot();
            this.browser.save('output/runAnalysisTests/StateBrowser.mat');
            %this.browser.save('output/runAnalysisTests/StateBrowser.pdf');
        end
    end
end