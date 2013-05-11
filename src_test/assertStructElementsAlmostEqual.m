function assertStructElementsAlmostEqual(varargin)
%assertElementsAlmostEqual Assert floating-point array elements almost equal.
%   assertElementsAlmostEqual(A, B, tol_type, tol, floor_tol) asserts that all
%   elements of floating-point arrays A and B are equal within some tolerance.
%   tol_type can be 'relative' or 'absolute'.  tol and floor_tol are scalar
%   tolerance values.
%
%   If the tolerance type is 'relative', then the tolerance test used is:
%
%       all( abs(A(:) - B(:)) <= tol * max(abs(A(:)), abs(B(:))) + floor_tol )
%
%   If the tolerance type is 'absolute', then the tolerance test used is:
%
%       all( abs(A(:) - B(:)) <= tol )
%
%   tol_type, tol, and floor_tol are all optional.  The default value for
%   tol_type is 'relative'.  If both A and B are double, then the default value
%   for tol and floor_tol is sqrt(eps).  If either A or B is single, then the
%   default value for tol and floor_tol is sqrt(eps('single')).
%
%   If A or B is complex, then the tolerance test is applied independently to
%   the real and imaginary parts.
%
%   Corresponding elements in A and B that are both NaN, or are both infinite
%   with the same sign, are considered to pass the tolerance test.
%
%   assertElementsAlmostEqual(A, B, ..., msg) prepends the string msg to the
%   output message if A and B fail the tolerance test.

%   Steven L. Eddins
%   Copyright 2008-2010 The MathWorks, Inc.

params = xunit.utils.parseFloatAssertInputs(varargin{:});

if ~isequal(size(params.A), size(params.B))
    message = xunit.utils.comparisonMessage(params.Message, ...
        'Inputs are not the same size.', ...
        params.A, params.B);
    throwAsCaller(MException('assertElementsAlmostEqual:sizeMismatch', ...
        '%s', message));
end

if ~(isstruct(params.A) && isstruct(params.B))
    message = xunit.utils.comparisonMessage(params.Message, ...
        'Inputs are not both structs.', ...
        params.A, params.B);
    throwAsCaller(MException('assertElementsAlmostEqual:notFloat', ...
        '%s', message));
end

%call recursively on each element of inputs
if numel(params.A) > 1
    for i = numel(params.A)
        assertStructElementsAlmostEqual(params.A(i), params.B(i), params.ToleranceType, params.Tolerance, params.FloorTolerance);
    end
    return;
end

%check that fields are almost equal
assertTrue(all(ismember(fieldnames(params.A), fieldnames(params.B))));
assertTrue(all(ismember(fieldnames(params.B), fieldnames(params.A))));
fieldNames = fieldnames(params.A);
for i = 1:numel(fieldNames)
    if isstruct(params.A.(fieldNames{i}))
        assertStructElementsAlmostEqual(params.A.(fieldNames{i}), params.B.(fieldNames{i}), params.ToleranceType, params.Tolerance, params.FloorTolerance);
    elseif isfloat(params.A.(fieldNames{i}))
        assertElementsAlmostEqual(params.A.(fieldNames{i}), params.B.(fieldNames{i}), params.ToleranceType, params.Tolerance, params.FloorTolerance);
    else
        assertEqual(params.A.(fieldNames{i}), params.B.(fieldNames{i}));
    end
end