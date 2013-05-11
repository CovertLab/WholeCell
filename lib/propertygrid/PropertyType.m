% Type information on a property.
% This class encapsulates type, shape and domain constraints on property
% values.
%
% Examples:
% PropertyType('denserealsingle', 'scalar')
%    creates a non-sparse single real scalar with unrestricted value
% PropertyType('denserealdouble', 'scalar', [-1,1])
%    creates a non-sparse double real scalar with value in range -1 to 1
% PropertyType('char', 'row', {'spring','summer','fall','winter'})
%    creates a character array which may be either be 'spring', 'summer',
%    'fall' or 'winter'
% PropertyType('logical', 'column', {'A','B','C'})
%    creates a set whose elements are 'A', 'B' and 'C'; the logical vector
%    is an indexer into this set, e.g. [1 1 0] maps to the set {'A','B'}
%
% See also: PropertyGridField

% Copyright 2010 Levente Hunyadi
classdef PropertyType
    properties
        % The underlying MatLab type for the property.
        PrimitiveType = 'denserealdouble';
        % The expected dimensions for the property.
        Shape = 'scalar';
        % A cell array of possible property values (or {} if unrestricted).
        Domain = [];
    end
    properties (Access = protected)
        ObjectType = '';
    end
    methods
        function self = PropertyType(type, shape, domain)
            self.PrimitiveType = type;
            self.Shape = shape;
            if nargin > 2
                self.Domain = domain;
            end
        end
        
        function s = char(self)
        % Compact textual representation of property type information.
            if ~isempty(self.Shape)
                s = sprintf('%s/%s', self.PrimitiveType, self.Shape);
            else
                s = self.PrimitiveType;
            end
            if ~isempty(self.Domain)
                if iscellstr(self.Domain)
                    s = [s ' ' strjoin(',', self.Domain)];
                elseif isnumeric(self.Domain)
                    s = [s ' ' mat2str(self.Domain)];
                end
            end
        end
        
        function tf = islogical(self)
            tf = strcmp(self.PrimitiveType, 'logical');
        end
        
        function tf = isinteger(self)
            tf = any(strmatch(self.PrimitiveType, {'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}, 'exact'));
        end
        
        function tf = isfloat(self)
            tf = any(strmatch(self.PrimitiveType, {'denserealsingle','denserealdouble','densecomplexsingle','densecomplexdouble','sparserealsingle','sparserealdouble','sparsecomplexsingle','sparsecomplexdouble'}, 'exact'));
        end
        
        function tf = isnumeric(self)
            tf = isinteger(self) || isfloat(self);
        end
        
        function tf = isreal(self)
            tf = any(strmatch(self.PrimitiveType, {'denserealsingle','denserealdouble','sparserealsingle','sparserealdouble'}, 'exact'));
        end
        
        function tf = is2d(self)
            tf = any(strmatch(self.Shape, {'row','column','matrix'}, 'exact'));
        end
        
        function tf = isvector(self)
            tf = any(strmatch(self.Shape, {'row','column'}, 'exact'));
        end
        
        function tf = isstring(self)
            tf = strcmp(self.PrimitiveType, 'char') && strcmp(self.Shape, 'row');
        end
        
        function tf = isset(self)
            tf = islogical(self) && isvector(self) && ~isempty(self.Domain);
        end

        function self = set.PrimitiveType(self, type)
            validateattributes(type, {'char'}, {'nonempty','row'});
            self.PrimitiveType = validatestring(type, { ...
                'logical', ...
                'char', ...
                'int8','uint8','int16','uint16','int32','uint32','int64','uint64', ...
                'denserealsingle','denserealdouble', ...
                'densecomplexsingle','densecomplexdouble', ...
                'sparserealsingle','sparserealdouble', ...
                'sparsecomplexsingle','sparsecomplexdouble', ...
                'cellstr', ...
                'object'});
        end
        
        function self = set.Shape(self, shape)
            validateattributes(shape, {'char'}, {'nonempty','row'});
            self.Shape = validatestring(shape, {'scalar','row','column','matrix','empty'});  % 'empty' is a pseudo-shape for the 0-by-0 empty matrix
        end
        
        function self = set.Domain(self, domain)
            if isempty(domain)
                self.Domain = [];
                return;
            end
            if ~isscalar(self) && ~isstring(self) && ~(islogical(self) && isvector(self))
                error('PropertyType:InvalidArgumentValue', ...
                    'Domain can be set only for scalars, strings or logical vectors.');
            end
            if iscell(domain)
                validateattributes(domain, {'cell'}, {'vector'});
                if ~(islogical(self) && isvector(self))  % interpret as exhaustive enumeration of domain elements
                    for k = 1 : numel(domain)
                        assert(self.CanAccept(domain{k}), ...
                            'PropertyType:InvalidArgumentValue', ...
                            'Invalid value specified for property domain.');
                    end
                end
            elseif isnumeric(domain)  % interpret as numeric range for domain elements
                validateattributes(domain, {'numeric'}, {'vector'});
                assert(isnumeric(self), 'PropertyType:ArgumentTypeMismatch', ...
                    'Closed interval domain is meaningful only for numeric data.');
                assert(length(domain) == 2, 'PropertyType:InvalidArgumentValue', ...
                    'Closed interval data is expected in the form [a,b] where a and b are numbers.');
                domain = [min(domain) max(domain)];
            else
                error('PropertyType:InvalidArgumentValue', ...
                    'Closed interval specification [a,b] or explicit enumeration {a,b,c,...} is expected.');
            end
            self.Domain = domain;
        end

        function tf = CanAccept(self, value, precisionloss)
        % Whether the property can accept the given value.
        %
        % Output arguments:
        % tf:
        %    false if any of the property constaints (type, shape or
        %    domain) would be violated
            validateattributes(self, {'PropertyType'}, {'scalar'});
            if nargin < 3
                precisionloss = false;
            end
            if isobject(value)
                tf = ~isempty(self.ObjectType) && metaclass(value) <= meta.class.fromName(self.ObjectType);
            elseif (isnumeric(value) || islogical(value) || ischar(value)) && any(isnan(value(:)));
                tf = false;
            else
                type = PropertyType.AutoDiscoverType(value);
                shape = PropertyType.AutoDiscoverShape(value);
                tf = PropertyType.CanAcceptType(self.PrimitiveType, type, precisionloss) ...
                    && PropertyType.CanAcceptShape(self.Shape, shape) ...
                    && self.CanAcceptValue(value);
            end
        end
        
        function clazz = GetPrimitiveMatLabType(self)
            validateattributes(self, {'PropertyType'}, {'scalar'});
            switch (self.PrimitiveType)
                case {'logical','char','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
                    clazz = self.PrimitiveType;
                case {'denserealsingle','densecomplexsingle','sparserealsingle','sparsecomplexsingle'}
                    clazz = 'single';
                case {'denserealdouble','densecomplexdouble','sparserealdouble','sparsecomplexdouble'}
                    clazz = 'double';
                case 'cellstr'
                    clazz = 'cell';
            end
        end
        
        function javatype = GetPrimitiveJavaType(self)
            validateattributes(self, {'PropertyType'}, {'scalar'});
            switch self.PrimitiveType
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
                case {'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle'}
                    matlabtype = 'cellstr';
                otherwise
                    error('MatrixEditor:ArgumentTypeMismatch', 'Type %s is not supported.', self.PrimitiveType);
            end
            javatype = javaclass(matlabtype);
        end
        
        function javavalue = GetPrimitiveJavaValue(self, value)
            validateattributes(self, {'PropertyType'}, {'scalar'});
            switch self.PrimitiveType
                case {'denserealdouble','sparserealdouble'}  % add a double-precision floating point property
                    javavalue = java.lang.Double(value);
                case {'denserealsingle','sparserealsingle'}  % add a single-precision floating point property
                    javavalue = java.lang.Float(value);
                case {'int8','uint8','int16'}
                    javavalue = java.lang.Short(value);
                case {'uint16','int32'}  % add an integer property
                    javavalue = java.lang.Integer(value);
                case {'uint32','int64'}
                    javavalue = java.lang.Long(value);
                case 'logical'  % add a logical property
                    javavalue = java.lang.Boolean(value);
                case {'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle'}
                    javavalue = java.lang.String(mat2str(value));
                otherwise
                    if isnumeric(value) || islogical(value)
                        javavalue = mat2str(value);
                    else
                        javavalue = char(value);
                    end
                    %error('PropertyType:ArgumentTypeMismatch', 'Type %s is not supported.', self.PrimitiveType);
            end
        end
        
        function javamatrix = GetJavaVectorOfVectors(self, matrix)
            assert(isnumeric(self) || islogical(self), 'PropertyType:InvalidOperation', ...
                'Operation supported on numeric and logical matrices only.');
            javamatrix = java.util.Vector(int32(size(matrix,1)));
            for m = 1 : size(matrix,1)
                javarow = java.util.Vector(int32(size(matrix,2)));
                javamatrix.add(javarow);
                for n = 1 : size(matrix,2)
                    javarow.add(self.GetPrimitiveJavaValue(matrix(m,n)));
                end
            end
        end
        
        function javamatrix = GetJavaMatrix(self, matrix)
        % Get MatLab value as Java 2d array of proper type.
            assert(isnumeric(self), 'PropertyType:InvalidOperation', ...
                'Operation supported on numeric matrices only.');
            assert(numel(matrix) > 0, 'PropertyType:InvalidOperation', ...
                'Empty values are not supported.');
            javatype = self.GetPrimitiveJavaType();
            javamatrix = javaArray(char(javatype.getName()), size(matrix,1), size(matrix,2));
            for m = 1 : size(matrix,1)
                for n = 1 : size(matrix,2)
                    javamatrix(m,n) = self.GetPrimitiveJavaValue(matrix(m,n));
                end
            end
        end
        
        function javavector = GetJavaVector(self, vector)
        % Get MatLab value as a Java vector of proper type.
            assert(isnumeric(self), 'PropertyType:InvalidOperation', ...
                'Operation supported on numeric matrices only.');
            assert(numel(vector) > 0, 'PropertyType:InvalidOperation', ...
                'Empty values are not supported.');
            javatype = self.GetPrimitiveJavaType();
            javavector = javaArray(char(javatype.getName()), numel(vector));
            for k = 1 : numel(vector)
                javavector(k) = self.GetPrimitiveJavaValue(vector(k));
            end
        end

        function javavalue = ConvertToJava(self, value)
        % Converts a MatLab value into the appropriate Java object.
            switch self.Shape
                case 'scalar'
                    javavalue = self.GetPrimitiveJavaValue(value);
                case {'row','column'}
                    switch self.PrimitiveType
                        case 'char'  % add a string property
                            if strcmp(self.Shape, 'row')
                                javavalue = value;
                            else
                                javavalue = mat2str(value);
                            end
                        case 'cellstr'
                            javavalue = java.lang.String(strjoin(sprintf('\n'), value));
                        case 'logical'
                            if ~isempty(self.Domain)
                                javavalue = javaStringArray(self.Domain(value));  % value is an indicator vector
                            else
                                javavalue = mat2str(value);
                            end
                        otherwise
                            if isnumeric(value) || islogical(value)
                                javavalue = mat2str(value);
                            else
                                javavalue = char(value);
                            end
                    end
                case 'empty'
                    javavalue = [];
                otherwise
                    if isnumeric(value) || islogical(value)
                        javavalue = mat2str(value);
                    else
                        javavalue = char(value);
                    end
            end
        end
        
        function [value,stringconversion] = ConvertFromJava(self, javavalue)
        % Converts a Java object into the appropriate MatLab value.
            stringconversion = false;
            switch self.Shape
                case 'scalar'
                    switch self.PrimitiveType
                        case 'denserealdouble'
                            value = full(double(javavalue));
                        case 'sparserealdouble'
                            value = sparse(double(javavalue));
                        case 'denserealsingle'
                            value = full(single(javavalue));
                        case 'sparserealsingle'
                            value = sparse(single(javavalue));
                        case {'int8','uint8','int16','uint16','int32','uint32','int64','logical'}
                            value = cast(javavalue, self.PrimitiveType);
                        case {'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle'}
                            value = self.ConvertFromString(javavalue);
                            stringconversion = true;
                    end
                case {'row','column'}
                    switch self.PrimitiveType
                        case 'char'  % add a string property
                            if strcmp(self.Shape, 'row')
                                value = char(javavalue);
                            else
                                value = self.ConvertFromString(javavalue);
                                stringconversion = true;
                                return;
                            end
                        case 'cellstr'
                            value = strsplit(javavalue);
                        case 'logical'
                            if ~isempty(self.Domain)
                                value = strsetmatch(cell(javavalue), self.Domain);
                            else
                                value = self.ConvertFromString(javavalue);
                                stringconversion = true;
                                return;
                            end
                        case {...
                                'denserealdouble','sparserealdouble','denserealsingle','sparserealsingle',...
                                'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle',...
                                'int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
                            value = self.ConvertFromString(javavalue);
                            stringconversion = true;
                    end
                    switch self.Shape
                        case 'row'
                            value = reshape(value, 1, numel(value));
                        case 'column'
                            value = reshape(value, numel(value), 1);
                    end
                case 'matrix'
                    switch self.PrimitiveType
                        case {...
                                'denserealdouble','sparserealdouble','denserealsingle','sparserealsingle',...
                                'densecomplexdouble','sparsecomplexdouble','densecomplexsingle','sparsecomplexsingle',...
                                'char','int8','uint8','int16','uint16','int32','uint32','int64','uint64','logical'}
                            value = self.ConvertFromString(javavalue);
                            stringconversion = true;
                    end
            end
        end

        function value = ConvertFromString(self, text)
        % Converts matrix string representation to its numeric value.
        %
        % Input arguments:
        % text:
        %    a textual representation of a numeric scalar, vector or matrix
        %
        % Output arguments:
        % value:
        %    a numeric value satisfying all type constraints, or NaN if the
        %    string cannot be converted or the value cannot be coerced into
        %    the appropriate type
            text = char(text);
            switch self.Shape
                case 'scalar'
                    value = str2double(text);
                    status = ~isnan(value);
                otherwise
                    [value,status] = str2num(text); %#ok<ST2NM>
            end
            if all(int32(value(:)) == value(:))  % all numbers are integers
                value = int32(value);            % cast to integer
            end
            if ~status || ~self.CanAccept(value, true)
                value = NaN;
                return;
            end
            value = self.ConvertFromMatLab(value);
        end
        
        function value = ConvertFromMatLab(self, value)
        % Casts a MatLab value to the prescribed property type.
            if ~isreal(value) && isnumeric(self) && isreal(self)
                error('PropertyType:ArgumentTypeMismatch', 'Complex value cannot be coerced into a real type.');
            end
            if iscellstr(value) && strcmp(self.PrimitiveType, 'cellstr') || isobject(value) && strcmp(self.PrimitiveType, 'object')
                return
            end
            switch self.PrimitiveType
                case {'denserealdouble','densecomplexdouble'}
                    value = full(double(value));
                case {'sparserealdouble','sparsecomplexdouble'}
                    value = sparse(double(value));
                case {'denserealsingle','densecomplexsingle'}
                    value = full(single(value));
                case {'sparserealsingle','sparsecomplexsingle'}
                    value = sparse(single(value));
                case {'int8','uint8','int16','uint16','int32','uint32','int64','logical'}
                    value = cast(value, self.PrimitiveType);
                case {'char'}
                    value = char(value);
                otherwise
                    error('PropertyType:ArgumentTypeMismatch', 'Cannot coerce type %s into type %s.', class(value), self.PrimitiveType);
            end
        end
    end
    methods (Access = private)
        function tf = CanAcceptValue(self, value)
        % Whether the property constraints allow the value to be accepted.
            if islogical(self) && isvector(self)  % interpret logical vector as a set
                tf = isempty(self.Domain) || numel(value) == numel(self.Domain);  % size of logical vector corresponds to size of universe
            else  % test whether domain contains value
                tf = PropertyType.IsInDomain(self.Domain, value);
            end
        end
    end
    methods (Static)
        function obj = AutoDiscover(value)
        % Constructs a property type instance based on a value.
        % The property type is chosen to be the most specific possible that
        % fits the given value.
            if isobject(value)
                obj = PropertyType('object', PropertyType.AutoDiscoverShape(value));
                obj.ObjectType = class(value);
            else
                obj = PropertyType(PropertyType.AutoDiscoverType(value), PropertyType.AutoDiscoverShape(value));
            end
        end
        
        function type = AutoDiscoverType(value)
            clazz = class(value);
            switch clazz
                case {'logical','char','int8','uint8','int16','uint16','int32','uint32','int64','uint64'}
                    type = clazz;
                case {'single','double'}
                    if issparse(value)
                        sparsity = 'sparse';
                    else
                        sparsity = 'dense';
                    end
                    if isreal(value)
                        complexity = 'real';
                    else
                        complexity = 'complex';
                    end
                    type = [ sparsity complexity clazz ];
                case 'cell'
                    if iscellstr(value)
                        type = 'cellstr';
                    else
                        error('PropertyType:InvalidArgumentValue', ...
                            'Cell arrays other than cell array of strings are not supported.');
                    end
                otherwise
                    error('PropertyType:InvalidArgumentValue', ...
                        'Argument type "%s" is not supported.', class(value));
            end
        end
        
        function shape = AutoDiscoverShape(value)
            if ndims(value) > 2
                error('PropertyType:InvalidArgumentValue', ...
                    'Dimensions higher than 2 are not supported.');
            end
            if size(value,1) == 1 && size(value,2) == 1
                shape = 'scalar';
            elseif size(value,1) == 1
                shape = 'row';
            elseif size(value,2) == 1
                shape = 'column';
            elseif size(value,1) == 0 && size(value,2) == 0
                shape = 'empty';  % no dimensions
            else
                shape = 'matrix';
            end
        end
        
        function tf = CanAcceptType(generaltype, specifictype, precisionloss)
        % Type assignment check.
        %
        % Output arguments:
        % tf:
        %    true if a variable of the specific type could be assigned to a
        %    variable of the general type
            validateattributes(generaltype, {'char'}, {'nonempty','row'});
            validateattributes(specifictype, {'char'}, {'nonempty','row'});
            if nargin < 3
                precisionloss = false;
            end
            if strcmp(generaltype, specifictype)
                tf = true;
                return;
            end
            switch generaltype
                case {'densecomplexdouble','sparsecomplexdouble'}
                    type = {...
                        'int8','uint8','int16','uint16','int32','uint32','int64','uint64', ...
                        'denserealsingle','denserealdouble', ...
                        'densecomplexsingle','densecomplexdouble', ...
                        'sparserealsingle','sparserealdouble', ...
                        'sparsecomplexsingle','sparsecomplexdouble'};
                case {'denserealdouble','sparserealdouble'}
                    type = {...
                        'int8','uint8','int16','uint16','int32','uint32','int64','uint64', ...
                        'denserealsingle','denserealdouble', ...
                        'sparserealsingle','sparserealdouble'};
                case {'densecomplexsingle','sparsecomplexsingle'}
                    type = {...
                        'int8','uint8','int16','uint16','int32','uint32','int64','uint64', ...
                        'denserealsingle','densecomplexsingle','sparserealsingle','sparsecomplexsingle'};
                    if precisionloss  % permit loss of precision
                        type = [type {'denserealdouble','densecomplexdouble','sparserealdouble','sparsecomplexdouble'}];
                    end
                case {'denserealsingle','sparserealsingle'}
                    type = {...
                        'int8','uint8','int16','uint16','int32','uint32','int64','uint64', ...
                        'denserealsingle','sparserealsingle'};
                    if precisionloss  % permit loss of precision
                        type = [type {'denserealdouble','sparserealdouble'}];
                    end
                case {'int64','uint64'}
                    type = {'int8','uint8','int16','uint16','int32','uint32'};
                case {'int32','uint32'}
                    type = {'int8','uint8','int16','uint16'};
                case {'int16','uint16'}
                    type = {'int8','uint8'};
                case {'int8','uint8'}
                    type = {};
                case 'char'
                    type = {};
                case 'logical'
                    type = {'int8','uint8','int16','uint16','int32','uint32'};
                case 'cellstr'
                    type = {};
                otherwise
                    type = {};
            end
            tf = any(strmatch(specifictype, type, 'exact'));
        end
        
        function tf = CanAcceptShape(generalshape, specificshape)
        % Shape assignment check.
        %
        % Output arguments:
        % tf:
        %    true if a variable of the specific shape could be assigned to
        %    a variable of the general shape
            validateattributes(generalshape, {'char'}, {'nonempty','row'});
            validateattributes(specificshape, {'char'}, {'nonempty','row'});
            if strcmp(generalshape, specificshape)
                tf = true;
                return;
            end
            switch generalshape
                case 'matrix'
                    shape = {'scalar','row','column','empty'};
                case 'row'
                    shape = {'scalar','empty'};
                case 'column'
                    shape = {'scalar','empty'};
                case 'scalar'
                    shape = {};
                case 'empty'
                    shape = {};
                otherwise
                    shape = {};
            end
            tf = any(strmatch(specificshape, shape, 'exact'));
        end
        
        function tf = IsInDomain(domain, value)
            if isempty(domain)
                tf = true;
            elseif iscellstr(domain)
                tf = any(strcmp(value, domain)); 
            elseif iscell(domain)
                tf = any(cellfun(@(v) v==value, domain));
            elseif isnumeric(domain) && length(domain) == 2
                tf = value >= min(domain) && value <= max(domain);
            end
        end
    end
end
