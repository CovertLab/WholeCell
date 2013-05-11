function values = parseDoubles(delimiter, s)
%PARSEDOUBLES Splits a list of doubles at each occurrence of a delimiter.
%   Returns a row vector of double. Ignores a trailing delimiter, if present.

    if isempty(s)
        values = [];
    else
        slen = length(s);
        dlen = length(delimiter);
        k = [1-dlen strfind(s, delimiter)];
        if k(end) + dlen < slen
            k = [k slen+1];
        end
        n = length(k) - 1;
        values = zeros([1 n]);        
        for i = 1:n    
            values(i) = str2double(s(k(i)+dlen:k(i+1)-1));
        end
    end
end
