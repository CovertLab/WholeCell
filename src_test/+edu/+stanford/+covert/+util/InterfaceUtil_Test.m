%InterfaceUtil
% Checks if a class (or class instance) properly implements a specificed
% interface.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef InterfaceUtil_Test < TestCase
    methods 
        function this = InterfaceUtil_Test(name)
            this = this@TestCase(name);
        end
    end    
    
    %tests
    methods
        function testAssertInterface(~)
            import edu.stanford.covert.util.InterfaceUtil;
            
            %ex 1-4
            InterfaceUtil.assertInterface('double', 'double');
            InterfaceUtil.assertInterface('double', metaclass(zeros(0,0,'double')));
            InterfaceUtil.assertInterface('double', meta.class.fromName('double'));            
            InterfaceUtil.assertInterface(zeros(1,1,'double'), 'double');            
            
            %ex 5
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            InterfaceUtil.assertInterface('double', metadata);
            
            %ex 6
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Name2 = 'plus';
            assertExceptionThrown(@() InterfaceUtil.assertInterface('double', metadata), 'Interface:invalidInput');
            
            %ex 7
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'data';
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            assertExceptionThrown(@() InterfaceUtil.assertInterface('double', metadata), 'Interface:incompleteInterface'); 
            
            %ex 8
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Access = 'protected';
            assertExceptionThrown(@() InterfaceUtil.assertInterface('double', metadata), 'Interface:incompatibleInterface'); 
            
            %ex 9
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Access = 'protected';
            assertExceptionThrown(@() InterfaceUtil.assertInterface('double', metadata), 'Interface:incompatibleInterface'); 
            
            %ex 10
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Access = {'protected','public'};
            InterfaceUtil.assertInterface('double', metadata); 
            
            %ex 11
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Abstract = 1;
            assertExceptionThrown(@() InterfaceUtil.assertInterface('double', metadata), 'Interface:incompatibleInterface'); 
            
            %ex 12
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Abstract = [0 1];
            InterfaceUtil.assertInterface('double', metadata); 
            
            %ex 13
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'plus';
            metadata.Methods{1}.Abstract = [];
            InterfaceUtil.assertInterface('double', metadata); 
            
            %ex 14
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'ind2sub';
            metadata.Methods{1}.InputNames = 3;
            InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata); 
            
            %ex 15
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'ind2sub';
            metadata.Methods{1}.InputNames = 4;
            assertExceptionThrown(@() InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata), 'Interface:incompatibleInterface'); 
            
            %ex 16
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Methods = cell(0,1);
            metadata.Methods{1} = struct;
            metadata.Methods{1}.Name = 'ind2sub';
            metadata.Methods{1}.InputNames = {'this','inds','siz'};
            InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata); 
            
            %ex 17
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'inds';
            metadata.Properties{1}.GetMethod = 1;
            metadata.Methods = cell(0,1);
            InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata); 
            
            %ex 18
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'inds';
            metadata.Properties{1}.GetMethod = 0;
            metadata.Methods = cell(0,1);
            assertExceptionThrown(@() InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata), 'Interface:incompatibleInterface'); 
            
            %ex 19
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'inds';
            metadata.Properties{1}.GetMethod = [0 1];
            metadata.Methods = cell(0,1);
            InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata); 
            
            %ex 20
            metadata = struct;
            metadata.Properties = cell(0,1);
            metadata.Properties{1}.Name = 'inds';
            metadata.Properties{1}.GetMethod = [];
            metadata.Methods = cell(0,1);
            InterfaceUtil.assertInterface('edu.stanford.covert.util.SparseMat', metadata); 
        end
    end
end