function [revision, differences] = revision()
%REVISION Gets the current Subversion revision and any uncommitted edits.

info = runSvnCommand('info');
[match, tokens] = regexp(info, 'Revision: (\d+)', 'match', 'tokens');
if isempty(match)
    revision = NaN;
    differences = '';
    return
else
    revision = str2double(tokens{1}{1});
    
    if nargout > 1
        differences = runSvnCommand('diff');
    end
end

function output = runSvnCommand(command)
[status, output] = system(['svn ' command]);
if status == 1
    [status, output] = system(['"lib/svn/svn" ' command]);
end
if status == 1
    throw(MException('revision:svnError', 'No subversion client available'));
end
