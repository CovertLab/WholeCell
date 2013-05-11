%bilinearFit
% Fits data to bilinear function
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/11/2012
function [a, b, cutX, gof] = bilinearFit(x, y)

tf = ~isnan(x) & ~isnan(y);
x = x(tf);
y = y(tf);

[~, order] = sort(x);
x = x(order);
y = y(order);

a = NaN(numel(x), 2);
b = NaN(numel(x), 2);
sse = NaN(numel(x), 2);
rsquare = NaN(numel(x), 2);
adjrsquare = NaN(numel(x), 2);
for i = 2:numel(x)-2
    [fitobject1, gof1] = fit(x(1:i), y(1:i), 'poly1');
    [fitobject2, gof2] = fit(x(i+1:end), y(i+1:end), 'poly1');
    a(i, :) = [fitobject1.p2 fitobject2.p2];
    b(i, :) = [fitobject1.p1 fitobject2.p1];
    sse(i, :) = [gof1.sse gof2.sse];
    rsquare(i, :) = [gof1.rsquare gof2.rsquare];
    adjrsquare(i, :) = [gof1.adjrsquare gof2.adjrsquare];
end

[~, cutIdx] = nanmin(sum(sse, 2));
cutX = mean(x(cutIdx:cutIdx+1));

a = a(cutIdx, :);
b = b(cutIdx, :);
sse = sse(cutIdx, :);
rsquare = rsquare(cutIdx, :);
adjrsquare = adjrsquare(cutIdx, :);

sse0 = sum(sse);
rsquare0 = 1 - sse0 / sum((y-mean(y)).^2);
n = numel(x);
p = 3;
adjrsquare0 = 1 - (1 - rsquare0) * (n-1) / (n-p-1);

gof = struct(...
    'sse', [sse0 sse], ...
    'rsquare', [rsquare0 rsquare], ...
    'adjrsquare', [adjrsquare0 adjrsquare]);