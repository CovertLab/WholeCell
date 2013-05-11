classdef XMLTestRunDisplay < CommandWindowTestRunDisplay
    properties (Access = protected)
        TicStack = uint64([])
        
        ReportName
        ReportLabel
        XMLFileName
        
        TimeStamp
        iTestCase
        TestCases
        
        NumFailures = 0
        NumErrors = 0
        NumSkipped = 0;
    end
    
    methods
        function this = XMLTestRunDisplay(ReportName, ReportLabel, XMLFileName)
            this = this@CommandWindowTestRunDisplay();
            this.ReportName = ReportName;
            this.ReportLabel = ReportLabel;
            this.XMLFileName = XMLFileName;
        end
    end
    
    methods (Access = protected)
        function testRunStarted(this, component)
            this.TimeStamp = datestr(now, 31);
            this.iTestCase = 0;
            this.TestCases = repmat(...
                struct('classname', [], 'name', [], 'location', [], 'time', [], 'skipped', [], 'failure', [], 'error', [], 'stack', [], 'exception', []), ...
                component.numTestCases(), 1);
        end
        
        function testRunFinished(this, totalRunTime)
            import edu.stanford.covert.util.escapeXML;
            
            this.displayFaults();
            
            [~, hostName] = edu.stanford.covert.util.computerInfo();
            try
                revision = edu.stanford.covert.util.revision();
            catch %#ok<CTCH>
                revision = [];
            end
            
            fid = fopen(this.XMLFileName, 'w');
            fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
            fprintf(fid, '<testsuite name="%s" label="%s" time="%f" tests="%d" failures="%d" errors="%d" skipped="%d" hostname="%s" timestamp="%s">\n', ...
                this.ReportName, this.ReportLabel, ...
                totalRunTime, ...
                numel(this.TestCases), this.NumFailures, this.NumErrors, this.NumSkipped, ...
                hostName, this.TimeStamp);
            fprintf(fid, '\t<properties>\n');
            fprintf(fid, '\t\t<property name="%s" value="%d"/>\n', 'revision', revision);
            fprintf(fid, '\t\t<property name="%s" value="%s"/>\n', 'matlab.version', version);
            fprintf(fid, '\t\t<property name="%s" value="%s"/>\n', 'java.version', version('-java'));
            fprintf(fid, '\t</properties>\n');
            for i = 1:numel(this.TestCases)
                testCase = this.TestCases(i);
                if isempty(testCase.skipped) && isempty(testCase.error) && isempty(testCase.failure)
                    fprintf(fid, '\t<testcase classname="%s" name="%s" time="%f"/>\n', ...
                        testCase.classname, testCase.name, testCase.time);
                elseif ~isempty(testCase.skipped)
                    fprintf(fid, '\t<testcase classname="%s" name="%s" time="%f" skipped="1">\n', ...
                        testCase.classname, testCase.name, testCase.time);
                    fprintf(fid, '\t\t<skipped message="%s"/>\n', escapeXML(testCase.skipped));
                    fprintf(fid, '\t</testcase>\n');
                elseif ~isempty(testCase.failure)
                    fprintf(fid, '\t<testcase classname="%s" name="%s" time="%f">\n', ...
                        testCase.classname, testCase.name, testCase.time);
                    fprintf(fid, '\t\t<failure type="%s" message="%s"><![CDATA[%s\n%s]]></failure>\n', ...
                        escapeXML(testCase.exception), escapeXML(testCase.failure), ...
                        this.formatStack(testCase.stack), escapeXML(testCase.failure, false));
                    fprintf(fid, '\t</testcase>\n');
                else
                    fprintf(fid, '\t<testcase classname="%s" name="%s" time="%f">\n', ...
                        testCase.classname, testCase.name, testCase.time);
                    fprintf(fid, '\t\t<error type="%s" message="%s"><![CDATA[%s\n%s]]></error>\n', ...
                        escapeXML(testCase.exception), escapeXML(testCase.error), ...
                        this.formatStack(testCase.stack), escapeXML(testCase.error, false));
                    fprintf(fid, '\t</testcase>\n');
                end
            end
            fprintf(fid, '\t<system-out><![CDATA[%s]]></system-out>\n', '');
            fprintf(fid, '\t<system-err><![CDATA[%s]]></system-err>\n', '');
            fprintf(fid, '</testsuite>\n');
            fclose(fid);
        end
    end
    
    methods
        function testComponentStarted(this, component)
            this.testComponentStarted@CommandWindowTestRunDisplay(component);            
            
            %testComponentStarted Update Command Window display
            
            this.pushTic();
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            end
            
            fprintf('%s%s', this.indentationSpaces(), component.Name);
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            else
                fprintf(' %s ', this.leaderDots(component.Name));
            end
            
            %disabled tests
            if ~isa(component, 'TestCase')
                metaClass = meta.class.fromName(component.Name);
                if ~isempty(metaClass)
                    for i = 1:numel(metaClass.Methods)
                        if numel(metaClass.Methods{i}.Name) >= 13 && strcmp(metaClass.Methods{i}.Name(1:13), 'disabled_test')                            
                            this.TestCases = [this.TestCases; struct(...
                                'classname', metaClass.Name, ...
                                'name', metaClass.Methods{i}.Name(10:end), ...
                                'location', component.Location, ...
                                'time', 0, ...
                                'skipped', 'disabled', ...
                                'failure', [], ...
                                'error', [], ...
                                'stack', [], ...
                                'exception', [])];
                            
                            this.NumSkipped = this.NumSkipped + 1;
                            fprintf('%s%s', repmat(' ', 1, 3*numel(this.TicStack)), metaClass.Methods{i}.Name(10:end));
                            fprintf(' %s ', repmat('.', 1, max(0, 60 - 3*numel(this.TicStack) - numel(metaClass.Methods{i}.Name(10:end)))));
                            fprintf('skiped in %12.6f seconds\n', 0);
                        end
                    end
                end
            end
        end
        
        function testComponentFinished(this, component, did_pass)
            %testComponentFinished Update Command Window display
            
            if ~isa(component, 'TestCase')
                fprintf('%s%s %s ', this.indentationSpaces(), component.Name, ...
                    this.leaderDots(component.Name));
            end
            
            component_run_time = toc(this.popTic());
            
            if did_pass
                fprintf('passed in %12.6f seconds\n', component_run_time);
            else
                fprintf('FAILED in %12.6f seconds\n', component_run_time);
            end
            
            if ~isa(component, 'TestCase')
                fprintf('\n');
            end
            
            if isa(component, 'TestCase')
                this.iTestCase = this.iTestCase + 1;
                this.TestCases(this.iTestCase).classname = class(component);
                this.TestCases(this.iTestCase).name = component.MethodName;
                this.TestCases(this.iTestCase).location = component.Location;
                this.TestCases(this.iTestCase).time = component_run_time;
            end
            
            if isempty(this.TicStack)
                this.testRunFinished(component_run_time);
            end
        end
        
        function testCaseFailure(this, test_case, failure_exception)
            this.NumFailures = this.NumFailures + 1;
            
            this.TestCases(this.iTestCase+1).failure = failure_exception.message;
            this.TestCases(this.iTestCase+1).exception = failure_exception.identifier;
            this.TestCases(this.iTestCase+1).stack = failure_exception.stack;
            
            this.testCaseFailure@CommandWindowTestRunDisplay(test_case, failure_exception);
        end
        
        function testCaseError(this, test_case, error_exception)
            this.NumErrors = this.NumErrors + 1;
            
            this.TestCases(this.iTestCase+1).error = error_exception.message;
            this.TestCases(this.iTestCase+1).exception = error_exception.identifier;
            this.TestCases(this.iTestCase+1).stack = error_exception.stack;
            
            this.testCaseError@CommandWindowTestRunDisplay(test_case, error_exception);
        end
    end
    
    methods (Access = protected)
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
        
        function output = formatStack(this, stack)
            output = [];
            stack = this.filterStack(stack);
            for k = 1:numel(stack)
                filename = stack(k).file;
                linenumber = stack(k).line;
                output = [output sprintf('%s at line %d\n', filename, linenumber)]; %#ok<AGROW>
            end
        end
        
        function new_stack = filterStack(~, stack)
            %filterStack Remove unmeaningful stack trace calls
            %    new_stack = filterStack(stack) removes from the input stack trace calls
            %    that are framework functions and methods that are not likely to be
            %    meaningful to the user.
            
            % Testing stack traces follow this common pattern:
            %
            % 1. The first function call in the trace is often one of the assert functions
            % in the framework directory.  This is useful to see.
            %
            % 2. The next function calls are in the user-written test functions/methods and
            % the user-written code under test.  These calls are useful to see.
            %
            % 3. The final set of function calls are methods in the various framework
            % classes.  There are usually several of these calls, which clutter up the
            % stack display without being that useful.
            %
            % The pattern above suggests the following stack filtering strategy: Once the
            % stack trace has left the framework directory, do not follow the stack trace back
            % into the framework directory.
            
            mtest_directory = fileparts(which('runtests'));
            last_keeper = numel(stack);
            have_left_mtest_directory = false;
            for k = 1:numel(stack)
                directory = fileparts(stack(k).file);
                if have_left_mtest_directory
                    if strcmp(directory, mtest_directory)
                        % Stack trace has reentered mtest directory.
                        last_keeper = k - 1;
                        break;
                    end
                else
                    if ~strcmp(directory, mtest_directory)
                        have_left_mtest_directory = true;
                    end
                end
            end
            
            new_stack = stack(1:last_keeper);
        end
    end
end