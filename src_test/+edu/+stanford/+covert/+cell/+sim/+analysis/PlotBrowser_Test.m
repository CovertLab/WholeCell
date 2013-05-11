%PlotBrowser_Test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/6/2011
classdef PlotBrowser_Test < TestCase
    properties
        simulation
    end
    
    methods
        function this = PlotBrowser_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        function setUp(this)
            this.simulation = edu.stanford.covert.cell.sim.SimulationFixture.load([], true);
        end
        
        function tearDown(~)
            close('all');
        end
    end
    
    methods
        function test(this)
            import edu.stanford.covert.cell.sim.analysis.PlotBrowser;
            browser = PlotBrowser(this.simulation);
            
            %test all plots
            messages = {};
            for i = 1:numel(browser.plotNames)
                try
                    browser.selectPlots(browser.plotNames{i});
                catch exception
                    messages{end+1} = sprintf('%20s: %s', browser.plotFunctions{i}, exception.message); %#ok<AGROW>
                end
            end
            assert(isempty(messages), ...
                sprintf('The following plot methods failed:\n%s',strjoin(sprintf('\n'),messages{:})));
            
            %close browser
            browser.close();
        end
    end
end