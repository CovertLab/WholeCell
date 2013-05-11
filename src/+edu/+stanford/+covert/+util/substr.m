function result = substr(s, i, j)
%SUBSTR Extracts a substring using s(i:j) when last > 0, or else s(i:end+j).

    if j > 0
        result = s(i:j);
    else
        result = s(i:end+j);
    end
end
