function assertAllEqual(a, B, custom_message)
%assertAllEqual Assert that all elements of a matrix equal a scalar value
%   assertAllEqual(a, B) throws an exception if any elements in B are not equal
%   to a, or if B is empty. Equality requires values to be of the same class.
%
%   assertAllEqual(a, B, MESSAGE) prepends the string MESSAGE to the assertion
%   message when the assertion does not hold.
%
%   Examples
%   --------
%   % This call returns silently.
%   assertAllEqual(1, ones([3 4]));
%
%   % This call throws an error.
%   assertAllEqual(false, [true false false]);
%
%   See also assertEqual, assertElementsAlmostEqual, assertVectorsAlmostEqual

% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/12/2010

if nargin < 3
    custom_message = '';
end

if ~isscalar(a)
    throwAsCaller(MException('assertAllEqual:nonScalarValue',...
        xunit.utils.comparisonMessage(custom_message, ...
        'The first argument must be a scalar.', a, B)));
end

if isempty(B)
    throwAsCaller(MException('assertAllEqual:emptyMatrix',...
        xunit.utils.comparisonMessage(custom_message, ...
        'The second argument must not be empty.', a, B)));
end

assertEqual(repmat(a, size(B)), B, custom_message);
