% A property grid based on the JIDE grid implementation.
% A graphical user control based on the integrated JIDE PropertyGrid that
% mimics MatLab's property inspector. Unlike the inspector, it supports
% structures, new-style MatLab objects, both with value and handle
% semantics.
% 
% The property grid displays a list of (object) properties with values
% editable in-place. Each property has an associated semantics (or type)
% that restricts the possible values the property can take and helps
% visualize the property value in a natural manner, in particular:
% * a character array is displayed as a string and can be edited in-place
% * a scalar logical is mapped to a checkbox
% * an integer value that has a limited range is manipulated with a spinner
% * a selection from a set of values is presented as a drop-down list
% * a cell array of strings (either row or column vector) can be edited as
%   multi-line text in a pop-up text box
% * a logical vector that is an indicator for a set (e.g. [false false
%   true] for 'C' from the universe {'A','B','C'}]) is visualized as a
%   checkbox list
% * numeric vectors and matrices can be edited elementwise in a pop-up
%   window
% * parent-child relationships are presented with the help of expandable
%   properties
%
% Supported types include all shapes (scalar, vector and matrix) of all
% primitive types (logicals, integers, real/complex double/single) as well
% as cell arrays of strings, structures, and both value and handle MatLab
% objects with arbitrary level of nesting.
%
% If a property is selected use F1 to get help (a dialog is displayed with
% the help text of the property) or F2 to edit a numeric matrix in a pop-up
% dialog.
%
% References:
% The com.jidesoft.grid package by JIDE Software,
%    http://www.jidesoft.com/products/JIDE_Grids_Developer_Guide.pdf
%
% See also: inspect

% Copyright 2010 Levente Hunyadi
classdef PropertyGrid < UIControl
    properties (Dependent)
        % The handle graphics control that wraps the property grid.
        Control;
        % Properties listed in the property grid.
        Properties;
        % The MatLab structure or object bound to the property grid.
        Item;
    end
    properties (Access = private)
        % A uipanel that wraps the property grid.
        Container;
        % A com.jidesoft.grid.PropertyPane instance.
        % Encapsulates a property table and decorates it with icons to
        % choose sorting order, expand and collapse categories, and a 
        % description panel.
        Pane;
        % A com.jidesoft.grid.PropertyTable instance.
        % Lets the user view and edit property name--value pairs.
        Table;
        % A com.jidesoft.grid.PropertyTableModel instance.
        % Contains the properties enlisted in the property grid.
        Model;
        Fields = JidePropertyGridField.empty(1,0);
        % The MatLab structure or object bound to the property grid.
        BoundItem = [];
    end
    methods
        function self = PropertyGrid(varargin)
            self = self@UIControl(varargin{:});
        end
        
        function self = Instantiate(self, parent)
            if nargin < 2
                parent = figure;
            end
            
            self.Container = uipanel(parent, ...
                'Units', 'normalized', ...
                'Position', [0 0 1 1], ...
                'Tag', '__PropertyGrid__', ...
                'UserData', self);

            % initialize JIDE
            com.mathworks.mwswing.MJUtilities.initJIDE;
            com.jidesoft.grid.CellEditorManager.registerEditor(javaclass('cellstr',1), com.jidesoft.grid.StringArrayCellEditor);
            com.jidesoft.grid.CellEditorManager.registerEditor(javaclass('char',1), com.jidesoft.grid.MultilineStringCellEditor, com.jidesoft.grid.MultilineStringCellEditor.CONTEXT);
            com.jidesoft.grid.CellRendererManager.registerRenderer(javaclass('char',1), com.jidesoft.grid.MultilineStringCellRenderer, com.jidesoft.grid.MultilineStringCellEditor.CONTEXT);

            % create JIDE property pane
            self.Table = handle(objectEDT('com.jidesoft.grid.PropertyTable'), 'CallbackProperties');  % property grid (without table model)
            self.Pane = objectEDT('com.jidesoft.grid.PropertyPane', self.Table);  % property pane (with icons at top and help panel at bottom)

            % control = jcontrol(parent, pane, 'Position', [0 0 1 1]);
            panel = self.Container;
            pixelpos = getpixelposition(panel);
            [control,container] = javacomponent(self.Pane, [0 0 pixelpos(3) pixelpos(4)], panel); %#ok<ASGLU>
            set(container, 'Units', 'normalized');
            set(self.Table, 'KeyPressedCallback', @PropertyGrid.OnKeyPressed);
        end
        
        function ctrl = get.Control(self)
            ctrl = self.Container;
        end
        
        function properties = get.Properties(self)
        % Retrieves properties displayed in the grid.
            properties = self.Fields.GetProperties();
        end
        
        function self = set.Properties(self, properties)
        % Explicitly sets properties displayed in the grid.
        % Setting this property clears any object bindings.
            validateattributes(properties, {'PropertyGridField'}, {'vector'});
            self.BoundItem = [];
            
            if ~isempty(self.Model)
                set(self.Model, 'PropertyChangeCallback', []);  % clear callback
            end
            
            % create JIDE properties
            toolbar = properties.HasCategory();
            description = properties.HasDescription();
            self.Fields = JidePropertyGridField.empty(0,1);
            for k = 1 : numel(properties)
                self.Fields(k) = JidePropertyGridField(properties(k));
            end

            % create JIDE table model
            list = self.Fields.GetTableModel();
            model = handle(com.jidesoft.grid.PropertyTableModel(list), 'CallbackProperties');
            model.setMiscCategoryName('Miscellaneous');  % caption for uncategorized properties
            model.expandAll();
            self.Model = model;

            % set JIDE table model to property table
            self.Table.setModel(model);
            self.Pane.setShowToolBar(toolbar);
            if toolbar
                self.Pane.setOrder(0);
            else
                self.Pane.setOrder(1);
            end
            self.Pane.setShowDescription(description);

            % wire property change event hook
            set(model, 'PropertyChangeCallback', @PropertyGrid.OnPropertyChange);
        end
        
        function item = get.Item(self)
        % Retrieves the object bound to the property grid.
            item = self.BoundItem;
        end
        
        function self = set.Item(self, item)
        % Binds an object to the property grid.
        % Any changes made in the property grid are automatically reflected
        % in the bound object. Only handle objects (i.e. those that derive
        % from the handle superclass) are supported for on-line binding
        % (i.e. changes are immediately available in the original object).
        % In order to simplify usage, this property may even be set with
        % value objects (i.e. those not derived from handle) and structures
        % but the Item property has to be queried explicitly to fetch
        % changes (off-line binding).
            if ~isempty(item)
                assert(isstruct(item) || isobject(item), 'PropertyGrid:ArgumentTypeMismatch', ...
                    'Only structures and MatLab objects are bindable.');
            end
            if ~isempty(item)
                properties = PropertyGridField.GenerateFrom(item);
            else
                properties = PropertyGridField.empty(1,0);
            end
            self.Bind(item, properties);
        end
        
        function self = Bind(self, item, properties)
        % Binds an object to the property grid with a custom property list.
            self.Properties = properties;
            self.BoundItem = item;
        end
        
        function s = GetPropertyValues(self)
        % Returns the set of property names and values in a structure.
            s = struct;
            for k = 1 : numel(self.Fields)
                field = self.Fields(k);
                s = nestedassign(s, field.PropertyData.Name, field.PropertyData.Value);
            end
        end
        
        function DisableToolBar(self)
            self.Pane.setShowToolBar(false);
        end
    end
    methods (Access = private)
        function EditMatrix(self, name)
        % Opens the matrix editor to change the entries of a matrix.
        %
        % Input arguments:
        % name:
        %    the name of the property (in dot notation) for which to open
        %    the matrix editor
            field = self.Fields.FindByName(name);
            if ~isnumeric(field.PropertyData.Type) && ~islogical(field.PropertyData.Type) || ~is2d(field.PropertyData.Type)
                return;
            end
            self.Table.editingStopped(javax.swing.event.ChangeEvent(self.Table));  % commit value to edited cell
            fig = figure( ...
                'DockControls', 'off', ...
                'MenuBar', 'none', ...
                'Name', sprintf('Edit matrix "%s"', name), ...
                'NumberTitle', 'off', ...
                'Toolbar', 'none');
            editor = MatrixEditor(fig, ...
                'ReadOnly', field.PropertyData.ReadOnly, ...
                'Type', field.PropertyData.Type, ...
                'Item', field.PropertyData.Value);
            uiwait(fig);
            field.Value = editor.Item;
            self.UpdateField(name);
        end
        
        function UpdateDependentProperties(self, field)
        % Propagates changes triggered by dependent properties.
        %
        % Input arguments:
        % field:
        %    the JidePropertyGridField that has changed
            if isempty(self.BoundItem)
                return;
            end
            % requery affected property values as needed
            if field.PropertyData.Dependent  % dependent property set; requery all properties as the dependent property might have changed the value of any of them
                for k = 1 : numel(self.Fields)
                    f = self.Fields(k);
                    if f ~= field
                        value = nestedfetch(self.BoundItem, f.PropertyData.Name);  % query dependent property value
                        f.Value = value;
                        f.PropertyData.Value = value;
                    end
                end
                self.Table.repaint();
            else  % requery dependent properties only
                dependent = getdependentproperties(self.BoundItem);  % a cell array of dependent property names
                if ~isempty(dependent)
                    for k = 1 : numel(dependent)
                        name = dependent{k};
                        value = nestedfetch(self.BoundItem, name);  % query dependent property value
                        field = self.Fields.FindByName(name);
                        field.Value = value;               % update value displayed in grid
                        field.PropertyData.Value = value;  % update value stored internally
                    end
                    self.Table.repaint();
                end
            end
        end
        
        function UpdateField(self, name)
        % Updates a property value or reverts changes if value is illegal.
            field = self.Fields.FindByName(name);
            value = field.Value;
            if field.CanAccept(value)
                try
                    if ~isempty(self.BoundItem)  % reflect changes in bound object
                        self.BoundItem = nestedassign(self.BoundItem, name, value);
                    end
                    field.PropertyData.Value = value;  % persist changes in property value
                    self.UpdateDependentProperties(field);
                catch me
                    field.Value = field.PropertyData.Value;  % revert changes
                    self.Table.repaint();
                    rethrow(me);
                end
            else
                field.Value = field.PropertyData.Value;  % revert changes
                self.Table.repaint();
            end
        end
    end
    methods (Access = private, Static)
        function self = FindPropertyGrid(obj, member)
        % Finds the object property grid that contains the given field.
        %
        % Input arguments:
        % obj:
        %    a com.jidesoft.grid.DefaultProperty instance
            validateattributes(member, {'char'}, {'nonempty','row'});
            % find which PropertyGrid contains the object for which the callback is executing
            h = findobjuser(@(userdata) userdata.(member) == obj, '__PropertyGrid__');
            self = get(h, 'UserData');
        end

        function name = GetSelectedProperty(obj)
        % The name of the currently selected property (if any).
        % Like JIDE, this function also uses a hierarchical naming scheme
        % (dot notation).
        %
        % Output arguments:
        % name:
        %    a selected property in dot notation
            selectedfield = obj.getSelectedProperty();
            if isempty(selectedfield)
                name = [];
            else
                name = char(selectedfield.getFullName());
            end
        end
    end
    methods (Static)  % methods (Access = private, Static) for MatLab 2010a and up
        function OnKeyPressed(obj, event)
        % Fired when a key is pressed when the property grid has the focus.
            key = char(event.getKeyText(event.getKeyCode()));
            switch key
                case 'F1'
                    name = PropertyGrid.GetSelectedProperty(obj);
                    self = PropertyGrid.FindPropertyGrid(obj, 'Table');
                    if ~isempty(name) && ~isempty(self.BoundItem)  % help
                        nameparts = strsplit(name, '.');
                        if numel(nameparts) > 1
                            helpobject = nestedfetch(self.BoundItem, strjoin('.', nameparts(1:end-1)));
                        else
                            helpobject = self.BoundItem;
                        end
                        helpdialog([class(helpobject) '.' nameparts{end}]);
                    end
                case 'F2'
                    name = PropertyGrid.GetSelectedProperty(obj);
                    if ~isempty(name)  % edit property value
                        self = PropertyGrid.FindPropertyGrid(obj, 'Table');
                        self.EditMatrix(name);
                    end
            end
        end
        
        function OnPropertyChange(obj, event) %#ok<INUSL>
        % Fired when a property value in a property grid has changed.
        % This function is declared static because object methods cannot be
        % directly used with the @ operator. Even though the anonymous
        % function construct @(obj,evt) self.OnPropertyChange(obj,evt);
        % could be feasible, it leads to a memory leak.
            self = PropertyGrid.FindPropertyGrid(obj, 'Model');
            name = get(event, 'PropertyName');  % JIDE automatically uses a hierarchical naming scheme
            self.UpdateField(name);

            if 0  % debug mode
                oldvalue = var2str(get(event, 'OldValue'));  % old value as a string
                newvalue = var2str(get(event, 'NewValue'));  % new value as a string
                fprintf('Property value of "%s" has changed.\n', name);
                fprintf('Old value: %s\n', oldvalue);
                fprintf('New value: %s\n', newvalue);
            end
        end
    end
end
