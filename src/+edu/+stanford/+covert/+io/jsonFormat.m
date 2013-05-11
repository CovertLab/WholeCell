function result = jsonFormat(value)
%edu.stanford.covert.io.jsonFormat Serializes a value to JSON format.
%
%   See edu.stanford.covert.io.jsonParse
    result = jsonJava(value);
    result = char(result.toString());
end

function result = jsonJava(value)
    import com.twolattes.json.Json;

    sz = size(value);    
    cl = class(value);    
    switch cl
        case 'char'
            if numel(sz) == 2
                if ~any(sz)
                    result = Json.string('');
                    return;
                elseif sz(1) == 1 && sz(2)
                    result = Json.string(value);
                    return;
                end
            end
            result = jsonArray(...
                jsonSize(value),...
                Json.string(reshape(value, 1, numel(value))));
        case 'logical'
            if numel(sz) == 2 && all(sz == [1 1])
                if value
                    result = Json.TRUE;
                else
                    result = Json.FALSE;
                end
            else
                result = jsonNumberMatrix(Json.TRUE, sz, value);
            end
        case 'double'
            if any(isnan(value(:)))
                throw(MException('jsonFormat:NaN', '%d NaN(s) in [%s]-sized double',...
                    nnz(isnan(value)), num2str(size(value))));
            end
            value = min(max(value, -realmax), realmax);
            if numel(sz) == 2 
                if ~any(sz)
                    result = Json.fromString('[]');
                    return;
                elseif sz(1) == 1
                    switch sz(2)
                        case 0;
                        case 1, result = Json.number(value); return;
                        otherwise, result = jsonNumberArray(value); return;
                    end
                end
            end
            result = jsonDoubleMatrix(sz, value);
        case 'single'
            if any(isnan(value(:)))
                throw(MException('jsonFormat:NaN', '%d NaN(s) in [%s]-sized single',...
                    nnz(isnan(value)), num2str(size(value))));
            end
            result = jsonNumberMatrix(Json.number(int32(1)), sz, min(max(value, -realmax('single')), realmax('single')));
        case 'int8'
            result = jsonNumberMatrix(Json.number(int32(7)), sz, value);
        case 'uint8'
            result = jsonNumberMatrix(Json.number(int32(8)), sz, value);
        case 'int16'
            result = jsonNumberMatrix(Json.number(int32(15)), sz, value);
        case 'uint16'
            result = jsonNumberMatrix(Json.number(int32(16)), sz, int32(value));
        case 'int32'
            result = jsonNumberMatrix(Json.number(int32(31)), sz, value);
        case 'uint32'
            result = jsonNumberMatrix(Json.number(int32(32)), sz, int64(value));
        case 'int64'
            result = jsonNumberMatrix(Json.number(int32(63)), sz, value);
        case 'uint64'
            result = jsonUint64Matrix(Json.number(int32(64)), sz, value);
        case 'cell'
            result = jsonCellArray(sz, value);
        case 'struct'
            names = fieldnames(value);            
            nVals = numel(value);
            nNames = numel(names);
            if nVals == 1
                result = Json.object();
                for j = 1:nNames
                    name = names{j};
                    result.put(Json.string(name), jsonJava(value.(name)));
                end
            else
                result = javaArray('com.twolattes.json.Json$Value', nVals * nNames + 3);
                result(1) = Json.object();
                result(2) = jsonIntArray(sz);
                result(3) = jsonStringArray(names);
                offset = 3;
                for i = 1:nVals
                    v = value(i);
                    for j = 1:nNames
                        result(offset + j) = jsonJava(v.(names{j}));
                    end
                    offset = offset + nNames;                
                end            
                result = Json.array(result);
            end
        otherwise
            st = struct(value);
            mc = metaclass(value);
            for i = 1:length(mc.Properties)
                p = mc.Properties{i};
                if p.Dependent || p.Constant
                    st = rmfield(st, p.Name);
                end
            end
            object = Json.object();
            names = fieldnames(st);
            nNames = length(names);
            for j = 1:nNames
                name = names{j};
                object.put(Json.string(name), jsonJava(value.(name)));
            end            
            result = javaArray('com.twolattes.json.Json$Value', 2);
            result(1) = Json.string(cl);
            result(2) = object;
            result = Json.array(result);
    end
end

function result = jsonArray(varargin)
    import com.twolattes.json.Json;

    result = javaArray('com.twolattes.json.Json$Value', nargin);
    for i = 1:nargin
        result(i) = varargin{i};
    end
    result = Json.array(result);
end

function result = jsonDoubleMatrix(sz, numbers)
    import com.twolattes.json.Json;

    n = numel(numbers);
    result = javaArray('com.twolattes.json.Json$Value', n + 1);
    result(1) = jsonIntArray(sz);
    for i = 1:n
        result(i+1) = Json.number(numbers(i));
    end
    result = Json.array(result);
end

function result = jsonNumberMatrix(v1, sz, numbers)
    import com.twolattes.json.Json;

    n = numel(numbers);
    result = javaArray('com.twolattes.json.Json$Value', n + 2);
    result(1) = v1;
    result(2) = jsonIntArray(sz);
    for i = 1:n
        result(i+2) = Json.number(numbers(i));
    end
    result = Json.array(result);
end

function result = jsonUint64Matrix(v1, sz, numbers)
    import com.twolattes.json.Json;

    n = numel(numbers);
    result = javaArray('com.twolattes.json.Json$Value', n + 2);
    result(1) = v1;
    result(2) = jsonIntArray(sz);
    for i = 1:n
        result(i+2) = Json.number(java.math.BigDecimal(num2str(numbers(i))));
    end
    result = Json.array(result);
end

function result = jsonSize(value)
    result = jsonIntArray(size(value));
end

function result = jsonIntArray(value)
    result = jsonNumberArray(int32(value));
end

function result = jsonNumberArray(value)
    import com.twolattes.json.Json;

    n = numel(value);
    result = javaArray('com.twolattes.json.Json$Value', n);
    for i = 1:n
        result(i) = Json.number(value(i));
    end
    result = Json.array(result);
end

function result = jsonStringArray(value)
    import com.twolattes.json.Json;

    n = numel(value);
    if n > 0
        result = javaArray('com.twolattes.json.Json$Value', n);
        for i = 1:n
            result(i) = Json.string(value{i});
        end
        result = Json.array(result);
    else
        result = Json.fromString('[]');
    end
end

function result = jsonCellArray(sz, value)
    import com.twolattes.json.Json;

    n = numel(value);
    if n > 0
        result = javaArray('com.twolattes.json.Json$Value', 2);
        result(1) = jsonIntArray(sz);
        vals = javaArray('com.twolattes.json.Json$Value', n);
        for i = 1:n
            vals(i) = jsonJava(value{i});
        end
        result(2) = Json.array(vals);
    else
        result = javaArray('com.twolattes.json.Json$Value', 1);
        result(1) = Json.fromString('[]');
    end
    result = Json.array(result);
end
