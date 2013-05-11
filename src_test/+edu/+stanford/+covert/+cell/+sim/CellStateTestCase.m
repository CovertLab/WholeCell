%CellStateTestCase
% Base test class for cell state classes
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef CellStateTestCase < TestCase
    properties
        state               %state object
        stateClass          %full class name of state
        stateClassName      %short class name of state
        stateBaseLocation   %base directory of state
    end
    
    %constructor
    methods
        function this = CellStateTestCase(name)
            %parent class constructor
            this = this@TestCase(name);
            
            %get class name of state
            this.stateClass     = edu.stanford.covert.util.substr(class(this), 1, -5);
            this.stateClassName = this.stateClass(find(this.stateClass=='.',1,'last')+1:end);
            
            %find parent folder of state -- used to store test fixture and
            %expected outputs
            locationParts = strfind(this.Location, filesep);
            if this.Location(locationParts(end-1)+1) == '@'
                this.stateBaseLocation = this.Location(1:locationParts(end-1));
            else
                this.stateBaseLocation = this.Location(1:locationParts(end));
            end
        end
    end
    
    %setup
    methods
        function setUp(this)
            constructor = str2func(this.stateClass);
            this.state = constructor([], this.stateClassName);
            this.loadTestFixture();
            this.state.verbosity = 0;
            this.state.seed = 1;
            this.state.seedRandStream();
        end
        
        %Called following the execution of each test case.
        function tearDown(this)
            this.state = [];
        end
        
        %Loads knowledge base test fixture into state. Called by setUp prior
        %to the execution of each test. You may wish to override this method to
        %set additional properties.
        function loadTestFixture(this)
            this.state = edu.stanford.covert.cell.sim.CellStateFixture.load(this.state);
        end
    end
    
    %test interface
    methods
        function testImplementsStateInterface(this)
            import edu.stanford.covert.util.InterfaceUtil;
            
            %% assert interface
            metadata = struct();
            
            %properties
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'wholeCellModelID';
            metadata.Properties{2}.Name = 'name';
            metadata.Properties{3}.Name = 'optionNames';
            metadata.Properties{4}.Name = 'fixedConstantNames';
            metadata.Properties{5}.Name = 'fittedConstantNames';
            metadata.Properties{6}.Name = 'parameterNames';
            metadata.Properties{7}.Name = 'parameterIndexs';
            metadata.Properties{8}.Name = 'stateNames';
            metadata.Properties{8}.Name = 'dependentStateNames';
            
            metadata.Properties{9}.Name = 'verbosity';
            metadata.Properties{10}.Name = 'seed';
            metadata.Properties{11}.Name = 'randStream';
            
            metadata.Properties{12}.Name = 'dryWeight';
            
            %methods
            metadata.Methods{1}.Name = 'seedRandStream';
            metadata.Methods{2}.Name = 'initializeConstants';
            metadata.Methods{3}.Name = 'allocateMemory';
            metadata.Methods{4}.Name = 'initialize';
            
            %assert
            InterfaceUtil.assertInterface(this.state, metadata);
            
            %% verbosity and seed are marked as options
            assertAllEqual(true, ismember({'verbosity'; 'seed'}, this.state.optionNames));
        end
    end
    
    %tests for run-time errors
    methods
        function testAllocateMemoryForState(this)
            s = this.state;
            s.allocateMemory(1);
        end
        
        function testInitialize(this)
            s = this.state;
            s.initialize();
        end
        
        function testWeightCalculations(this)
            s = this.state;
            
            s.dryWeight;
        end
        
        %check for run-time errors in calculations
        function testDependentProperties(this)
            s = this.state;
            metaClass = metaclass(s);
            for i = 1:numel(metaClass.Properties)
                if ~metaClass.Properties{i}.Dependent && isempty(metaClass.Properties{i}.GetMethod)
                    continue;
                end
                val = s.(metaClass.Properties{i}.Name); %#ok<NASGU>
            end
        end
        
        function testDisp(this)
            s = this.state;
            s.disp();
            s.display();
        end
    end
end