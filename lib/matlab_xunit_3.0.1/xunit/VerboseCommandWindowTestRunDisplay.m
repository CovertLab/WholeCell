classdef VerboseCommandWindowTestRunDisplay < CommandWindowTestRunDisplay
%VerboseCommandWindowTestRunDisplay Print test suite execution results to Command Window.
%   VerboseCommandWindowTestRunDisplay is a subclass of
%   CommandWindowTestRunDisplay.  It supports the -verbose option of runtests.
%
%   Overriddent methods:
%       testComponentStarted  - Update Command Window display
%       testComponentFinished - Update Command Window display
%       testRunFinished       - Update Command Window display at end of run
%
%   See also CommandWindowTestRunDisplay, TestRunLogger, TestRunMonitor, TestSuite

%   Steven L. Eddins
%   Copyright 2010 The MathWorks, Inc.     
    
    properties (SetAccess = private, GetAccess = private)
        TicStack = uint64([])
    end
    
    methods
        
        function testComponentStarted(self, component)
            %testComponentStarted Update Command Window display
            
            self.pushTic();
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            end
            
            fprintf('%s%s', self.indentationSpaces(), component.Name);
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            else
                fprintf(' %s ', self.leaderDots(component.Name));
            end
        end    
            
        function testComponentFinished(self, component, did_pass)
            %testComponentFinished Update Command Window display

            if ~isa(component, 'TestCase')
                fprintf('%s%s %s ', self.indentationSpaces(), component.Name, ...
                    self.leaderDots(component.Name));
            end
            
            component_run_time = toc(self.popTic());
            
            if did_pass
                fprintf('passed in %12.6f seconds\n', component_run_time);
            else
                fprintf('FAILED in %12.6f seconds\n', component_run_time);
            end
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            end
            
            if isempty(self.TicStack)
                self.testRunFinished();
            end
                
        end
        
    end
    
    methods (Access = protected)
        function testRunFinished(self)
            %testRunFinished Update Command Window display
            %    obj.testRunFinished(component) displays information about the test
            %    run results, including any test failures, to the Command
            %    Window.
            
            self.displayFaults();
        end
    end
    
    methods (Access = private)
        function pushTic(self)
            self.TicStack(end+1) = tic;
        end
        
        function t1 = popTic(self)
            t1 = self.TicStack(end);
            self.TicStack(end) = [];
        end
        
        function str = indentationSpaces(self)
            str = repmat(' ', 1, self.numIndentationSpaces());
        end
        
        function n = numIndentationSpaces(self)
            indent_level = numel(self.TicStack) - 1;
            n = 3 * indent_level;
        end
        
        function str = leaderDots(self, name)
            num_dots = max(0, 60 - self.numIndentationSpaces() - numel(name));
            str = repmat('.', 1, num_dots);
        end
        
    end
    
end
