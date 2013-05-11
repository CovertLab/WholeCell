function [x,nroot]=cubicfcn(a, b, c, d)
% ----------------------------------------------------------------------
% Usage:       [x,nroot]=cubicfcn(a, b, c, d)
%
% Solve a cubic equation where a, b, c, and d are real.
%   a*x^3 + b*x^2 + c*x + d = 0
%
% Public Variables
%   a, b, c, d  ... coefficients (input)
%   x           ... three (generally) complex solutions (output)
%   nroot       ... number of roots (output)
%
% Instructor: Nam Sun Wang
% ----------------------------------------------------------------------

% Local Variables:
%   y1, y2, y3  ... three transformed solutions
%
% Formula used are given in Tuma, "Engineering Mathematics Handbook", p7
%   (McGraw Hill, 1978).
%   Step 0: If a is 0. use the quadratic formula to avoid dividing by 0.
%   Step 1: Calculate p and q
%           p = ( 3*c/a - (b/a)**2 ) / 3
%           q = ( 2*(b/a)**3 - 9*b*c/a/a + 27*d/a ) / 27
%   Step 2: Calculate discriminant D
%           D = (p/3)**3 + (q/2)**2
%   Step 3: Depending on the sign of D, we follow different strategy.
%           If D<0, thre distinct real roots.
%           If D=0, three real roots of which at least two are equal.
%           If D>0, one real and two complex roots.
%   Step 3a: For D>0 and D=0,
%           Calculate u and v
%           u = cubic_root(-q/2 + sqrt(D))
%           v = cubic_root(-q/2 - sqrt(D))
%           Find the three transformed roots
%           y1 = u + v
%           y2 = -(u+v)/2 + i (u-v)*sqrt(3)/2
%           y3 = -(u+v)/2 - i (u-v)*sqrt(3)/2
%   Step 3b Alternately, for D<0, a trigonometric formulation is more convenient
%           y1 =  2 * sqrt(|p|/3) * cos(phi/3)
%           y2 = -2 * sqrt(|p|/3) * cos((phi+pi)/3)
%           y3 = -2 * sqrt(|p|/3) * cos((phi-pi)/3)
%           where phi = acos(-q/2/sqrt(|p|**3/27))
%                 pi  = 3.141592654...
%   Step 4  Finally, find the three roots
%           x = y - b/a/3
% ----------------------------------------------------------------------

%     pi = 3.141592654

% Step 0: If a is 0 use the quadratic formula. -------------------------
if (a == 0.)
    [x, nroot] = quadfcn(b, c, d);
    return
end

% Cubic equation with 3 roots
nroot = 3;

% Step 1: Calculate p and q --------------------------------------------
p  = c/a - b*b/a/a/3. ;
q  = (2.*b*b*b/a/a/a - 9.*b*c/a/a + 27.*d/a) / 27. ;

% Step 2: Calculate DD (discriminant) ----------------------------------
DD = p*p*p/27. + q*q/4. ;

% Step 3: Branch to different algorithms based on DD -------------------
if (DD < 0.)
    %       Step 3b:
    %       3 real unequal roots -- use the trigonometric formulation
    phi = acos(-q/2./sqrt(abs(p*p*p)/27.));
    temp1=2.*sqrt(abs(p)/3.);
    y1 =  temp1*cos(phi/3.);
    y2 = -temp1*cos((phi+pi)/3.);
    y3 = -temp1*cos((phi-pi)/3.);
else
    %       Step 3a:
    %       1 real root & 2 conjugate complex roots OR 3 real roots (some are equal)
    temp1 = -q/2. + sqrt(DD);
    temp2 = -q/2. - sqrt(DD);
    u = abs(temp1)^(1./3.);
    v = abs(temp2)^(1./3.);
    if (temp1 < 0.); u=-u; end
    if (temp2 < 0.); v=-v; end
    y1  = u + v;
    y2r = -(u+v)/2.;
    y2i =  (u-v)*sqrt(3.)/2.;
end

% Step 4: Final transformation -----------------------------------------
temp1 = b/a/3.;
y1 = y1-temp1;
if (DD < 0.)
    y2 = y2-temp1;
    y3 = y3-temp1;
else
    y2r=y2r-temp1;
end

% Assign answers -------------------------------------------------------
if (DD < 0.)
    x(1) = y1;
    x(2) = y2;
    x(3) = y3;
elseif (DD == 0.)
    x(1) =  y1;
    x(2) = y2r;
    x(3) = y2r;
else
    x(1) = y1;
    x(2) = y2r + y2i*1i;
    x(3) = y2r - y2i*1i;
end