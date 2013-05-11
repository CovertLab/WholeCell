% A single JIDE property field in a JIDE property grid.
% This is an internal class and is not intended to be used directly.
%
% References:
% http://www.jidesoft.com/products/JIDE_Grids_Developer_Guide.pdf
%
% See also: PropertyGrid

% Copyright 2010 Levente Hunyadi
classdef JidePropertyGridField < handle
    properties
        % A PropertyGridField instance.
        PropertyData;
    end
    properties (Dependent)
        Value;
    end
    properties (Access = private)
        % Children the property field might have.
        Children = JidePropertyGridField.empty(0,1);
        % A JIDE DefaultProperty instance.
        Control;
        % JIDE editor context (if any).
        Context;
        % JIDE editor context type (if any).
        ContextType;
    end
    methods
        function self = JidePropertyGridField(data)
            if nargin > 0
                self.Initialize(data);
            end
        end
        
        function self = Initialize(self, data)
        % Initializes a JIDE property based on a property grid field.
        % May register editor contexts with the global CellEditorManager.
            validateattributes(self, {'JidePropertyGridField'}, {'scalar'});
        
            self.PropertyData = data;
            field = com.jidesoft.grid.DefaultProperty();
            field.setName(data.Name);  % JIDE automatically uses a hierarchical naming scheme
            field.setCategory(data.Category);
            field.setDisplayName(data.DisplayName);
            field.setDescription(data.Description);
            field.setEditable(~data.ReadOnly);
            switch data.Type.Shape
                case 'scalar'
                    switch data.Type.PrimitiveType
                        case {'denserealdouble','sparserealdouble'}  % add a double-precision floating point property
                            matlabtype = 'double';  % MatLab type that is marshalled to Java
                        case {'denserealsingle','sparserealsingle'}  % add a single-precision floating point property
                            matlabtype = 'single';
                        case {'int8','uint8','int16'}
                            matlabtype = 'int16';
                        case {'uint16','int32'}  % add an integer property
                            matlabtype = 'int32';
                        case {'uint32','int64'}
                            matlabtype = 'int64';
                        case 'logical'  % add a logical property
                            matlabtype = 'logical';
                            field.setEditorContext(com.jidesoft.grid.BooleanCheckBoxCellEditor.CONTEXT);
                        case {'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle'}
                            matlabtype = [];
                            field.setType(javaclass('char',1));  % edit as string and convert with eval
                        otherwise
                            matlabtype = [];
                            field.setEditable(false);
                            %error('PropertyGrid:ArgumentTypeMismatch', 'Scalar %s is not supported.', data.Type.PrimitiveType);
                    end
                    if ~isempty(matlabtype)
                        javatype = javaclass(matlabtype);
                        field.setType(javatype);
                        if ~isempty(data.Type.Domain)
                            if iscell(data.Type.Domain)  % explicit enumeration of domain elements
                                self.AddComboBoxEditor(field, javatype, data.Type.Domain);
                            elseif isnumeric(data.Type.Domain) && isinteger(data.Type)  % domain expressed as interval
                                self.AddSpinnerEditor(field, javatype, min(data.Type.Domain), max(data.Type.Domain));
                            end
                        end
                    end
                case {'row','column'}
                    InitializeVector(self, field);
                case 'matrix'
                    InitializeMatrix(self, field);
                case 'empty'
                    field.setEditable(false);
                otherwise
                    error('PropertyGrid:ArgumentTypeMismatch', 'Data shape %s is not supported.', data.Type.Shape);
            end
            self.Control = field;
            self.Value = data.Value;
            
            for k = 1 : numel(data.Children)
                self.Children(k) = JidePropertyGridField(data.Children(k));
                self.Control.addChild(self.Children(k).Control);
            end
        end
        
        function delete(self)
        % Deallocates registered editor contexts.
            if ~isempty(self.Context)
                com.jidesoft.grid.CellEditorManager.unregisterEditor(self.ContextType, self.Context);
            end
        end
        
        function value = get.Value(self)
        % Gets the native MatLab value of a property.
            javavalue = self.Control.getValue();
            [value,stringconversion] = self.PropertyData.Type.ConvertFromJava(javavalue);
            if stringconversion  % display pretty-formatted value
                javavalue = self.PropertyData.Type.ConvertToJava(value);
                self.Control.setValue(javavalue);
            end
        end
        
        function set.Value(self, value)
        % Sets the Java value of a property based on the native value.
            javavalue = self.PropertyData.Type.ConvertToJava(value);
            self.Control.setValue(javavalue);
        end

        function tf = CanAccept(self, value)
        % Tests if the property can be assigned the native value.
        %
        % Input arguments:
        % value:
        %    a native MatLab value
            validateattributes(self, {'JidePropertyGridField'}, {'scalar'});
            tf = self.PropertyData.Type.CanAccept(value);
        end

        function self = FindByName(selfarray, name)
        % Finds the property field with the given name in an array.
            self = selfarray.FindByNameRecurse(strsplit(name, '.'));
        end

        function properties = GetProperties(selfarray)
            properties = PropertyGridField.empty(0,1);
            for k = 1 : numel(selfarray)
                properties(k) = selfarray(k).PropertyData;
            end
        end

        function list = GetTableModel(selfarray)
            n = numel(selfarray);
            ctrls = cell(1,n);
            for k = 1 : n
                ctrls{k} = selfarray(k).Control;
            end
            list = javaArrayList(ctrls);
        end
    end
    methods (Access = private)
        function self = FindByNameRecurse(selfarray, nameparts)
            names = getclassfield(getclassfield(selfarray, 'PropertyData'), 'Name');
            ix = strmatch(nameparts{1}, names, 'exact');
            if isempty(ix)
                self = JidePropertyGridField.empty(0,1);
                return;
            end
            ix = ix(1);
            namerest = nameparts(2:end);
            if isempty(namerest)
                self = selfarray(ix);
            else
                self = selfarray(ix).Children.FindByNameRecurse(namerest);
            end
        end
        
        function InitializeVector(self, field)
            data = self.PropertyData;
            switch data.Type.PrimitiveType
                case 'char'  % add a string property
                    switch data.Type.Shape
                        case 'row'
                            field.setType(javaclass('char',1));
                            if ~isempty(data.Type.Domain)
                                self.AddComboBoxEditor(field, javaclass('char',1), javaStringArray(data.Type.Domain));
                            end
                        otherwise
                            field.setType(javaclass('char',1));  % edit as string and convert with eval
                    end
                case 'cellstr'
                    field.setType(javaclass('char',1));
                    field.setEditorContext(com.jidesoft.grid.MultilineStringCellEditor.CONTEXT);
                case 'logical'
                    if ~isempty(data.Type.Domain)
                        field.setType(javaclass('cellstr',1));  % java.lang.String array
                        self.AddCheckBoxListEditor(field, data.Type.Domain);
                    else
                        field.setType(javaclass('char',1));  % edit as string and convert with eval
                    end
                case {...
                        'denserealdouble','sparserealdouble','denserealsingle','sparserealsingle',...
                        'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle',...
                        'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
                    field.setType(javaclass('char',1));  % edit as string and convert with eval
                otherwise
                    field.setEditable(false);
                    %error('PropertyGrid:ArgumentTypeMismatch', 'Unsupported type: %s %s.', data.Type.Shape, data.Type.PrimitiveType);
            end
        end
        
        function InitializeMatrix(self, field)
            data = self.PropertyData;
            switch data.Type.PrimitiveType
                case {...
                        'denserealdouble','sparserealdouble','denserealsingle','sparserealsingle',...
                        'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle',...
                        'char','int8','uint8','int16','uint16','int32','uint32','int64','uint64','logical'}
                    field.setType(javaclass('char',1));  % edit as string and convert with eval
                otherwise
                    field.setEditable(false);
                    %error('PropertyGrid:ArgumentTypeMismatch', 'Matrix %s is not supported.', data.Type.PrimitiveType);
            end
        end
        
        function ApplyContext(self, field, javatype, editor, editortype)
        % Registers a context for an editor.
        %
        % Input arguments:
        % field:
        %    a com.jidesoft.grid.Property instance
        % javatype:
        %    the Java type to be associated with the context
        % editor:
        %    a cell editor instance
        % editortype:
        %    the key used for generating the name of the context
            self.Context = com.jidesoft.grid.EditorContext(['JIDE_' editortype '_' int2str(self.GetContextCounterValue())]);
            self.ContextType = javatype;
            com.jidesoft.grid.CellEditorManager.registerEditor(javatype, editor, self.Context);  % register new context
            field.setEditorContext(self.Context);  % apply context to property
        end
        
        function AddCheckBoxListEditor(self, field, labels)
        % Registers a new check box list context.
        %
        % Input arguments:
        % field:
        %    a com.jidesoft.grid.Property instance
        % labels:
        %    a cell array of strings to label elements in the set shown
            editor = com.jidesoft.grid.CheckBoxListComboBoxCellEditor(javaStringArray(labels), javaclass('cellstr',1));
            self.ApplyContext(field, javaclass('cellstr',1), editor, 'checkboxlist');
        end
        
        function AddComboBoxEditor(self, field, javatype, javadomain)
        % Registers a new list selection context.
        % A drop-down combo box allows the user to select a single value
        % from a predefined set of values.
        %
        % Input arguments:
        % field:
        %    a com.jidesoft.grid.Property instance
        % javatype:
        %    a java.lang.Class instance
        % javadomain:
        %    an array of Java objects whose type corresponds to javatype
            editor = com.jidesoft.grid.ListComboBoxCellEditor(javadomain);
            self.ApplyContext(field, javatype, editor, 'combobox');
        end
        
        function AddSliderEditor(self, field, javatype, lower, upper)
        % Registers a new slider context.
        % The slider has a limited range.
        %
        % field:
        %    a com.jidesoft.grid.Property instance
        % javatype:
        %    a java.lang.Class instance
        % lower:
        %    an integer representing the upper lower of the interval
        % upper:
        %    an integer representing the upper bound of the interval
            editor = com.jidesoft.grid.SliderCellEditor(int32(lower), int32(upper));
            self.ApplyContext(field, javatype, editor, 'slider');
        end
        
        function AddSpinnerEditor(self, field, javatype, lower, upper)
        % Registers a new spinner context.
        % The spinner has a limited range and a fixed step.
        %
        % field:
        %    a com.jidesoft.grid.Property instance
        % javatype:
        %    a java.lang.Class instance
        % lower:
        %    an integer representing the upper lower of the interval
        % upper:
        %    an integer representing the upper bound of the interval
            spinner = javax.swing.SpinnerNumberModel(int32(lower), int32(lower), int32(upper), int32(1));
            editor = com.jidesoft.grid.SpinnerCellEditor(spinner);
            self.ApplyContext(field, javatype, editor, 'spinner');
        end
    end
    methods (Static, Access = private)
        function counter = GetContextCounterValue()
        % The current counter value for the JIDE editor context.
            persistent jide_counter;
            if isempty(jide_counter)  % if not yet assigned
                jide_counter = 0;
            end
            jide_counter = jide_counter + 1;  % increment global field counter
            counter = jide_counter;           % return current value
        end
    end
end
