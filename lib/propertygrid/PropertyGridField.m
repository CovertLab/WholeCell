% Metadata on a single property.
% This class is used to provide type and shape information, initial value,
% categorization, and context-sensitive help when a property is added as a
% field to a property grid.
%
% See also: PropertyGrid, PropertyType

% Copyright 2010 Levente Hunyadi
classdef PropertyGridField < hgsetget
    properties
        % The short name of the property.
        Name;
        % The type constraints the property value must conform to.
        Type;
        % The initial value of the property.
        Value;
        % The name of the category which the property should be put into.
        Category;
        % A longer, descriptive name for the property.
        DisplayName;
        % A more detailed description of the property for novice users.
        Description;
        % Whether the property is read-only.
        ReadOnly = false;
        % Whether the property is computed based on other properties.
        Dependent = false;
        % Direct descendants of this property.
        Children = PropertyGridField.empty(1,0);
    end
    methods
        function self = PropertyGridField(name, value, varargin)
            self.Name = name;
            self.Value = value;
            self = constructor(self, varargin{:});
        end

        function self = set.Name(self, name)
            validateattributes(name, {'char'}, {'nonempty','row'});
            self.Name = name;
        end
        
        function self = set.Type(self, type)
            validateattributes(type, {'PropertyType'}, {'scalar'});
            assert(type.CanAccept(self.Value), ...
                'PropertyGridField:InvalidOperation', ...
                'Setting type "%s" would invalidate current property value.', char(type)); %#ok<MCSUP>
            self.Type = type;
            self.Value = type.ConvertFromMatLab(self.Value);  %#ok<MCSUP> % convert value type to match assigned type
        end
        
        function self = set.Value(self, value)
            if isempty(self.Type) %#ok<MCSUP>
                self.Value = value;
                self.Type = PropertyType.AutoDiscover(value); %#ok<MCSUP>
            else
                assert(self.Type.CanAccept(value), ...
                    'PropertyGridField:InvalidArgumentValue', ...
                    'This value does not conform to property type restrictions.'); %#ok<MCSUP>
                self.Value = value;
            end
        end
        
        function self = set.Category(self, category)
            validateattributes(category, {'char'}, {'nonempty','row'});
            self.Category = category;
        end
        
        function self = set.DisplayName(self, name)
            if isempty(name)
                self.DisplayName = [];
            else
                validateattributes(name, {'char'}, {'nonempty','row'});
                self.DisplayName = name;
            end
        end

        function self = set.Description(self, description)
            if isempty(description)
                self.Description = [];  % clear description
            else
                validateattributes(description, {'char'}, {'nonempty','row'});
                self.Description = description;
            end
        end
        
        function self = set.ReadOnly(self, readonly)
            validateattributes(readonly, {'logical'}, {'scalar'});
            self.ReadOnly = readonly;
        end
        
        function tf = HasCategory(selfarray)
        % True if any object in the array has a category specification.
            for k = 1 : numel(selfarray)
                self = selfarray(k);
                if ~isempty(self.Category) || self.Children.HasCategory()
                    tf = true;
                    return;
                end
            end
            tf = false;
        end

        function tf = HasDescription(selfarray)
        % True if any object in the array has a description.
            for k = 1 : numel(selfarray)
                self = selfarray(k);
                if ~isempty(self.Description) || self.Children.HasDescription()
                    tf = true;
                    return;
                end
            end
            tf = false;
        end
        
        function root = GetHierarchy(selfarray)
        % Converts a flat property list into a hierarchical structure.
        % Parent-child relationships are discovered from hierarchical
        % qualified names where components are separated by dot (.).
            root = arrayfilter(@(self) isempty(strfind(self.Name, '.')) , selfarray);  % get simple names only
            root.WireHierarchy(selfarray);
            selfarray.UnqualifyNames();
        end
        
        function UnqualifyNames(selfarray)
        % Reduces qualified names to unqualified names (without dot).
            for k = 1 : numel(selfarray)
                self = selfarray(k);
                dot = strfind(self.Name, '.');
                if ~isempty(dot)
                    self.Name = self.Name(dot(end)+1:end);
                end
            end
        end
        
        function self = FindByName(selfarray, name)
        % Looks up a property field by name.
            for k = 1 : numel(selfarray)
                self = selfarray(k);
                if strcmp(name, self.Name)
                    return
                end
            end
            self = [];  % not found
        end
    end
    methods (Access = private)
        function WireHierarchy(selfarray, descendants)
        % Wires parent-child relationships between nested properties.
            for k = 1 : numel(selfarray)
                self = selfarray(k);
                self.Children = descendants.FilterChildren(self.Name);  % add direct descendants as children
                self.Children.WireHierarchy(descendants);               % recurse for children
            end
        end

        function filteredarray = FilterChildren(selfarray, filterprefix)
        % An array of direct children.
        %
        % Input arguments:
        % filterprefix:
        %    a cell array that specifies which members to select.
        %
        % Example:
        % nodes.FilterChildren('hu.bme') selects
        %    hu.bme.aut, hu.bme.cs, hu.bme.mit
        % nodes.FilterChildren('hu.bme') does not select
        %    hu.bme.aut.www (not a direct child)
            names = getclassfield(selfarray, 'Name');
            if iscell(filterprefix)
                prefix = [strjoin('.', filterprefix) '.'];
            else
                prefix = [filterprefix '.'];
            end
            ix = strmatch(prefix, names);  % get names that begin with prefix
            if isempty(ix)
                filteredarray = PropertyGridField.empty(1,0);
                return;
            end
            
            len = numel(prefix);
            names = names(ix);                                             % drop names that do not begin with prefix
            names = cellfun(@(name) name(len+1:end), names, ...            % drop leading prefix from names
                'UniformOutput', false);
            filter = cellfun(@(name) isempty(strfind(name, '.')), names);  % property names that become simple after prefix is chopped
            filteredarray = selfarray(ix(filter));
        end
    end
    methods (Static)
        function fields = GenerateFrom(obj)
        % Property fields for a structure, a value or a handle object.
            if isstruct(obj)
                fields = PropertyGridField.GenerateFromStruct(obj);
            elseif isobject(obj)
                fields = PropertyGridField.GenerateFromClass(obj);
            else
                fields = PropertyGridField.empty(1,0);
            end
        end
        
        function fields = GenerateFromStruct(obj)
        % Automatically generated property fields for a structure.
            names = fieldnames(obj);
            n = numel(names);
            
            k = 0;
            fields = PropertyGridField.empty(1,0);
            for i = 1 : n
                name = names{i};
                value = obj.(name);
                k = k + 1;
                fields(k) = PropertyGridField(name, value);
                fields(k).Children = PropertyGridField.GenerateFrom(value);
            end
        end
        
        function fields = GenerateFromClass(obj)
        % Automatically generated property fields for an object.
            assert(isobject(obj), 'PropertyGridField:ArgumentTypeMismatch', ...
                'New-style MatLab object (defined with the classdef keyword) is required.');
            try
                clazz = metaclass(obj);
            catch %#ok<CTCH>
                fields = GenerateFromStruct(obj);  % fallback to method for old-style classes (i.e. not defined with the classdef keyword)
                return;
            end

            k = 0;  % number of properties found
            n = numel(clazz.Properties);  % maximum number of properties
            fields = PropertyGridField.empty(0, 1);
            for i = 1 : n
                property = clazz.Properties{i};
                if property.Abstract || property.Hidden || ~strcmp(property.GetAccess, 'public');
                    continue;
                end
                try
                    value = obj.(property.Name);
                catch %#ok<CTCH>
                    continue;  % could not fetch property value
                end
                description = property.Description;  % not currently used in MatLab
                if isempty(description)
                    text = helptext([clazz.Name '.' property.Name]);  % fetch help text as description
                    if ~isempty(text)
                        description = text{1};  % use first line
                    end
                end
                readonly = property.Constant || ~strcmp(property.SetAccess, 'public') || property.Dependent && isempty(property.SetMethod);
                dependent = property.Dependent;
                
                k = k + 1;
                fields(k) = PropertyGridField(property.Name, value, ...
                    'DisplayName', property.Name, ...
                    'Description', description, ...
                    'ReadOnly', readonly, ...
                    'Dependent', dependent);
                fields(k).Children = PropertyGridField.GenerateFrom(value);
            end
        end
    end
end
