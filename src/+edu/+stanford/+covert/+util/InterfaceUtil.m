%InterfaceUtil
% Checks if a class (or class instance) properly implements a specificed
% interface.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef InterfaceUtil
    methods (Static = true)
        %Throws an exception is the class/class instance doesn't properly
        %implement the specified interface.
        %
        %Inputs:
        %- instanceOrClassName is an instance of a class, or a char equal to the
        %  name of a class
        %- metaclass is a struct of duck type meta.class with properties
        %  - Properties: cell array of objects of duck type meta.property. Each
        %    element of Properties must have the field Name,
        %    corresponding to the name of the property whose attributes should
        %    be validated. The other field names of Properties must equal
        %    names of property attributes, and the field values represent
        %    acceptable values of the corresponding property attribute. A
        %    Properties field value of null means that any value of the
        %    property attribute is acceptable. 
        %  - Methods: cell array of objects of duck type meta.method. Each
        %    element of Methods must have the field Name, corresponding
        %    to the name of the method whose attributes should be validated. The
        %    other field names of Methods must equal names of method
        %    attributes, and the field values represent acceptable values of the
        %    corresponding method attribute. A Methods field value of
        %    null means that any value of the method attribute is acceptable.        
        function assertInterface(instanceOrClassName, instanceOrMetaclass)
            import edu.stanford.covert.util.InterfaceUtil;
            
            if isa(instanceOrMetaclass, 'meta.class') || isstruct(instanceOrMetaclass)
                metadata = instanceOrMetaclass;
            elseif ischar(instanceOrMetaclass)
                metadata = meta.class.fromName(instanceOrMetaclass);
            else
                metadata = metaclass(instanceOrMetaclass);
            end
            
            InterfaceUtil.assertPropertyInterface(instanceOrClassName, metadata.Properties);
            InterfaceUtil.assertMethodInterface(instanceOrClassName, metadata.Methods);
        end
    end
    
    methods (Static = true, Access = protected)
        function assertPropertyInterface(instanceOrClassName, interface)
            import edu.stanford.covert.util.InterfaceUtil;
            
            if isempty(interface)
                return;
            end
            
            %validate input
            if ~isvector(interface) || ~iscell(interface) || ~all(cellfun(@InterfaceUtil.hasDuckTypeMetaProperty, interface))
                throw(MException('Interface:invalidInput', 'interface must be a cell array of objects of duck type meta.property'))
            end
            
            %get metadata
            if ischar(instanceOrClassName)
                metadata = meta.class.fromName(instanceOrClassName);
            else
                metadata = metaclass(instanceOrClassName);
            end
            
            %check property interface
            propertyNames = cellfun(@(property) property.Name, metadata.Properties, 'UniformOutput', false);
            for i = 1:numel(interface)
                iProperty = find(strcmp(propertyNames, interface{i}.Name), 1, 'first');                
                interfaceFieldNames = setdiff(fieldnames(interface{i}), {'Name', 'Description', 'DetailedDescription', 'DefiningClass'});

                if isempty(iProperty)
                    throw(MException('Interface:incompleteInterface', ...
                        'The class %s must implement a %s property', ...
                        metadata.Name, interface{i}.Name));
                end

                property = metadata.Properties{iProperty};
                for j = 1:numel(interfaceFieldNames)
                    if isempty(interface{i}.(interfaceFieldNames{j}))
                        continue;
                    end
                    
                    interfaceValue = interface{i}.(interfaceFieldNames{j});
                    interfaceValueStr = interface{i}.(interfaceFieldNames{j});
                    if ischar(interfaceValue)                        
                        interfaceValueStr = {interfaceValue};
                        interfaceValue = {interfaceValue};
                    elseif ~iscell(interfaceValue)                                   
                        interfaceValueStr = cellfun(@num2str, num2cell(interfaceValue), 'UniformOutput', false);
                    end
                    
                    propertyValue = property.(interfaceFieldNames{j});
                    if ischar(propertyValue)
                        %GetAccess, SetAccess                        
                        propertyValue = {propertyValue};
                    elseif isa(propertyValue, 'function_handle')
                        %GetMethod, SetMethod
                        propertyValue = ~isempty(propertyValue);
                    end
                    
                    if ~ismember(propertyValue, interfaceValue)
                        throw(MException('Interface:incompatibleInterface', ...
                            'The class %s %s property must have attribute %s with value {%s}', ...
                            metadata.Name, interface{i}.Name, interfaceFieldNames{j}, ['''' strjoin(''', ''',interfaceValueStr{:}) '''']));
                    end
                end
            end
        end
        
        function assertMethodInterface(instanceOrClassName, interface)
            import edu.stanford.covert.util.InterfaceUtil;
            
            if isempty(interface)
                return;
            end
            
            %validate input
            if ~isvector(interface) || ~iscell(interface) || ~all(cellfun(@InterfaceUtil.hasDuckTypeMetaMethod, interface))
                throw(MException('Interface:invalidInput', 'interface must be a cell array of objects of duck type meta.method'))
            end
            
            %get metadata
            if ischar(instanceOrClassName)
                metadata = meta.class.fromName(instanceOrClassName);
            else
                metadata = metaclass(instanceOrClassName);
            end
            
            %check method interface
            methodNames = cellfun(@(method) method.Name, metadata.Methods, 'UniformOutput', false);
            for i = 1:numel(interface)
                iMethod = find(strcmp(methodNames, interface{i}.Name), 1, 'first');                
                interfaceFieldNames = setdiff(fieldnames(interface{i}), {'Name', 'Description', 'DetailedDescription', 'DefiningClass'});

                if isempty(iMethod)
                    throw(MException('Interface:incompleteInterface', ...
                        'The class %s must implement a %s method', ...
                        metadata.Name, interface{i}.Name));
                end

                method = metadata.Methods{iMethod};
                for j = 1:numel(interfaceFieldNames)
                    if isempty(interface{i}.(interfaceFieldNames{j}))
                        continue;
                    end
                    
                    interfaceValue = interface{i}.(interfaceFieldNames{j});
                    interfaceValueStr = interface{i}.(interfaceFieldNames{j});
                    if ischar(interfaceValue)                        
                        interfaceValueStr = {interfaceValue};
                        interfaceValue = {interfaceValue};
                    elseif ~iscell(interfaceValue)                                   
                        interfaceValueStr = cellfun(@num2str, num2cell(interfaceValue), 'UniformOutput', false);
                    end
                    
                    methodValue = method.(interfaceFieldNames{j});
                    if ischar(methodValue)
                        %Access
                        methodValue = {methodValue};
                        
                        if ~ismember(methodValue, interfaceValue)
                            throw(MException('Interface:incompatibleInterface', ...
                                'The class %s %s method must have attribute %s with value {%s}', ...
                                metadata.Name, interface{i}.Name, interfaceFieldNames{j}, ['''' strjoin(''', ''',interfaceValueStr{:}) '''']));
                        end
                    elseif iscell(methodValue)
                        %InputNames, OutputNames
                        methodValue = numel(methodValue);
                        if iscell(interfaceValue)
                            interfaceValue = numel(interfaceValue);
                        elseif numel(interfaceValue)~=1
                            throw(MException('Interface:invalidInput', 'interface must specify number or names of %s',interfaceFieldNames{j}));
                        end
                        
                        if methodValue < interfaceValue
                            throw(MException('Interface:incompatibleInterface', ...
                                'The class %s %s method must have at at least %d %s', ...
                                metadata.Name, interface{i}.Name, interfaceValue, interfaceFieldNames{j}));
                        end
                    else
                        if ~ismember(methodValue, interfaceValue)
                            throw(MException('Interface:incompatibleInterface', ...
                                'The class %s %s method must have attribute %s with value {%s}', ...
                                metadata.Name, interface{i}.Name, interfaceFieldNames{j}, ['''' strjoin(''', ''',interfaceValueStr{:}) '''']));
                        end
                    end
                end
            end
        end
        
        function tf = hasDuckTypeMetaProperty(prop)
            tf = false;
            
            if isa(prop, 'meta.property')
                tf = true;
                return;
            end
            
            metaPropertyFieldNames = {'Name'; 'Description'; 'DetailedDescription'; 
                'GetAccess'; 'SetAccess'; 'Dependent'; 'Constant'; 'Abstract'; 'Transient'; 
                'Hidden'; 'GetObservable'; 'SetObservable'; 'AbortSet'; 'GetMethod'; 
                'SetMethod'; 'HasDefault'; 'DefaultValue'; 'DefiningClass'};
            
            if isstruct(prop) && ismember('Name',fieldnames(prop)) && all(ismember(fieldnames(prop), metaPropertyFieldNames))
                tf = true;
                return
            end           
        end
        
        function tf = hasDuckTypeMetaMethod(method)
            tf = false;
            
            if isa(method, 'meta.method')
                tf = true;
                return;
            end
            
            metaMethodFieldNames = {'Name'; 'Description'; 'DetailedDescription';
                'Access'; 'Static'; 'Abstract'; 'Sealed'; 'Hidden'; 'InputNames';
                'OutputNames'; 'DefiningClass'};
            
            if isstruct(method) && ismember('Name',fieldnames(method)) && all(ismember(fieldnames(method), metaMethodFieldNames))
                tf = true;
                return
            end
        end
    end
end