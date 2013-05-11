function result = jsonParse(jsonString)
%jsonParse Deserializes a value from JSON format.
%
%   See edu.stanford.covert.io.jsonFormat

    try
        result = toMatlabValue(com.twolattes.json.Json.fromString(jsonString));
    catch ex
        if strcmp(ex.identifier, 'MATLAB:Java:GenericException')
            throw(MException('jsonParse:InvalidJSON', ex.message));
        else
            throw(ex);
        end
    end
end

function result = toMatlabValue(value)
    import com.twolattes.json.Json;

    switch class(value)
    case 'com.twolattes.json.Json$NumberImplBigDecimal'
        result = value.getNumber().doubleValue();
    case 'com.twolattes.json.Json$StringImpl'
        result = char(value.getString());
    case 'com.twolattes.json.Json$BooleanImpl'
        result = value.getBoolean();
    case 'com.twolattes.json.Json$ArrayImpl'
        if value.isEmpty()
            result = [];
            return;
        end
        iter = value.iterator();
        v1 = iter.next();
        switch class(v1)
        case 'com.twolattes.json.Json$BooleanImpl'
            sz = parseDoubleArray(iter.next());
            result = logical(parseNumberMatrix(iter, sz, 'int32', 'intValue'));
        case 'com.twolattes.json.Json$NumberImplBigDecimal'
            v2 = iter.next();
            if strcmp(class(v2), 'com.twolattes.json.Json$NumberImplBigDecimal')
                result = parseDoubleArray(value);
                return;
            end
            sz = parseDoubleArray(v2);
            switch v1.getNumber().intValue()
            case 1, result = parseNumberMatrix(iter, sz, 'single', 'floatValue');
            case 7, result = parseNumberMatrix(iter, sz, 'int8', 'shortValue');
            case 8, result = parseNumberMatrix(iter, sz, 'uint8', 'shortValue');
            case 15, result = parseNumberMatrix(iter, sz, 'int16', 'shortValue');
            case 16, result = parseNumberMatrix(iter, sz, 'uint16', 'intValue');
            case 31, result = parseNumberMatrix(iter, sz, 'int32', 'intValue');
            case 32, result = parseNumberMatrix(iter, sz, 'uint32', 'longValue');
            case 63, result = parseNumberMatrix(iter, sz, 'int64', 'longValue');
            case 64
                result = zeros(sz, 'uint64');
                for i = 1:numel(result)
                    result(i) = str2double(char(iter.next().getNumber().toString()));
                end
            end
        case 'com.twolattes.json.Json$ArrayImpl'
            if v1.isEmpty()
                result = cell.empty;
                return;
            end
            sz = parseDoubleArray(v1);
            if value.size() == 2
                v2 = iter.next();
                switch class(v2)
                case 'com.twolattes.json.Json$StringImpl'
                    if v2.isEmpty()
                        result = reshape('', sz);
                    else
                        result = reshape(char(v2.getString()), sz);
                    end
                case 'com.twolattes.json.Json$ArrayImpl'
                    result = parseCellArrayIterator(v2.iterator(), sz);
                end
            else
                result = parseNumberMatrix(iter, sz, 'double', 'doubleValue');
            end
        case 'com.twolattes.json.Json$ObjectImpl'
            sz = parseDoubleArray(iter.next());
            names = parseStringArray(iter.next());
            values = parseCellArrayIterator(iter, [length(names) prod(sz)]);
            result = reshape(cell2struct(values, names, 1), sz);
        case 'com.twolattes.json.Json$StringImpl'
            constructor = str2func(char(v1.getString()));
            properties = struct;
            entries = iter.next().entrySet().iterator();
            while entries.hasNext()
                e = entries.next();
                properties.(char(e.getKey().getString())) = ...
                    toMatlabValue(e.getValue());
            end
            result = constructor(properties);
        end
    case 'com.twolattes.json.Json$ObjectImpl'
        result = struct;
        entries = value.entrySet().iterator();
        while entries.hasNext()
            e = entries.next();
            result.(char(e.getKey().getString())) = ...
                toMatlabValue(e.getValue());
        end
    end
end

function result = parseDoubleArray(jsonArray)
    iter = jsonArray.iterator();
    result = zeros([1 jsonArray.size()]);
    for i = 1:numel(result)
        result(i) = iter.next().getNumber().doubleValue();
    end
end

function result = parseNumberMatrix(iter, sz, type, methodName)
    result = zeros(sz, type);
    for i = 1:numel(result)
        result(i) = iter.next().getNumber().(methodName)();
    end
end

function result = parseCellArrayIterator(iter, sz)
    result = cell(sz);
    for i = 1:numel(result)
        result{i} = toMatlabValue(iter.next());
    end
end

function result = parseStringArray(array)
    result = cell([array.size() 1]);
    iter = array.iterator();
    for i = 1:numel(result)
        result{i} = char(iter.next().getString());
    end
end
