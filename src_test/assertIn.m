function assertIn(a, interval, interval_kind, custom_message)
%assertIn Assert that a scalar value is within an interval
%   assertIn(A, INTERVAL) throws an exception if scalar A is not in INTERVAL,
%   a two-element vector treated as a closed interval.
%
%   assertIn(A, INTERVAL, INTERVAL_KIND) behaves the same, but allows the caller
%   to specify whether the interval is closed ('[]'), open ('()'), or
%   half-closed ('[)' or '(]').
%
%   assertIn(a, INTERVAL, MESSAGE) and
%   assertIn(A, INTERVAL, INTERVAL_KIND, MESSAGE) prepend the string MESSAGE to
%   the assertion message when the assertion does not hold.
%
%   Examples
%   --------
%   % These calls return silently.
%   assertIn(3, [0 4]);
%   assertIn(pi, [0 4], 'example message');
%   assertIn(true, [false true], '(]', 'example message');
%
%   % These calls throw an error.
%   assertIn(3, [4 0]);
%   assertIn(-3, [0 4]);
%   assertIn(3, [0 3], '[)');
%
%   See also assertEqual, assertElementsAlmostEqual, assertVectorsAlmostEqual

% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/5/2010

switch nargin
    case 3
        if regexp(interval_kind, '[\[(][\])]')
            custom_message = '';
        else
            custom_message = interval_kind;
            interval_kind = '[]';
        end
    case 4
        if ~ischar(interval_kind) || ~regexp(interval_kind, '[\[(][\])]')
            throwAsCaller(MException('assertIn:invalidIntervalKind',...
                formatMessage(custom_message,...
                'The third argument must be ''[]'', ''[)'', ''(]'', or ''()''.')));
        end        
    otherwise
        interval_kind = '[]';
        custom_message = '';
end

validateattributes(a, {'numeric'}, {'vector', 'column'});
validateattributes(interval, {'numeric'}, {'size', [numel(a) 2]});

invalidValues = ...
    a < interval(:, 1) | interval_kind(1) == '(' & a == interval(:, 1) | ...
    a > interval(:, 2) | interval_kind(2) == ')' & a == interval(:, 2);

if any(invalidValues)
    msg = [];
    for i = 1:numel(invalidValues)
        if ~invalidValues(i)
            continue;
        end
        if ~isempty(msg)
            msg = [msg sprintf('\n')]; %#ok<AGROW>
        end
        msg = [msg sprintf('Element %d: %s outside interval %s', i, num2str(a(i)), formatInterval(interval(i, :), interval_kind))]; %#ok<AGROW>
    end
    throwAsCaller(MException('assertIn:outsideInterval',...
        formatMessage(custom_message, msg)));
end

function result = formatInterval(interval, kind)
result = [kind(1) num2str(interval(1)) ' ' num2str(interval(2)) kind(2)];

function msg = formatMessage(user_msg, msg)
if ~isempty(user_msg)
    msg = sprintf('%s\n%s', user_msg, msg);
end
