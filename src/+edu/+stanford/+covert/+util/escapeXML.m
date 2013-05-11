%XML excape a string
function str = escapeXML(str, fullyEscape)
if nargin < 2 || fullyEscape
    str = strrep(str, '&', '&amp;');
    str = strrep(str, '"', '&quot;');
    str = strrep(str, '''', '&#039;');
    str = strrep(str, '<', '&lt;');
    str = strrep(str, '>', '&gt;');
end

str = regexprep(str, ['[^a-zA-Z0-9&;:\.\- !@#\$%\^\*\(\)_\+=\{\}\[\]"''<>,?\\/|`~' sprintf('\r\n\t') ']'], '');