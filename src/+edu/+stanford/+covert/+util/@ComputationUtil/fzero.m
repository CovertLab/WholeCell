function [b,fval,exitflag,output] = fzero(FunFcnIn,x,options,varargin)
%FZERO  Single-variable nonlinear zero finding. 
%   X = FZERO(FUN,X0) tries to find a zero of the function FUN near X0, 
%   if X0 is a scalar.  It first finds an interval containing X0 where the 
%   function values of the interval endpoints differ in sign, then searches 
%   that interval for a zero.  FUN is a function handle.  FUN accepts real 
%   scalar input X and returns a real scalar function value F, evaluated 
%   at X. The value X returned by FZERO is near a point where FUN changes 
%   sign (if FUN is continuous), or NaN if the search fails.  
%
%   X = FZERO(FUN,X0), where X0 is a vector of length 2, assumes X0 is a 
%   finite interval where the sign of FUN(X0(1)) differs from the sign of 
%   FUN(X0(2)). An error occurs if this is not true.  Calling FZERO with a
%   finite interval guarantees FZERO will return a value near a point where
%   FUN changes sign.
%
%   X = FZERO(FUN,X0), where X0 is a scalar value, uses X0 as a starting 
%   guess. FZERO looks for an interval containing a sign change for FUN and 
%   containing X0.  If no such interval is found, NaN is returned.  
%   In this case, the search terminates when the search interval 
%   is expanded until an Inf, NaN, or complex value is found. Note: if
%   the option FunValCheck is 'on', then an error will occur if an NaN or 
%   complex value is found.
%
%   X = FZERO(FUN,X0,OPTIONS) solves the equation with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, FunValCheck, OutputFcn, and PlotFcns. 
%
%   X = FZERO(PROBLEM) finds the zero of a function defined in PROBLEM. 
%   PROBLEM is a structure with the function FUN in PROBLEM.objective, 
%   the start point in PROBLEM.x0, the options structure in PROBLEM.options,
%   and solver name 'fzero' in PROBLEM.solver. The structure PROBLEM must have 
%   all the fields.
%
%   [X,FVAL]= FZERO(FUN,...) returns the value of the function described 
%   in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FZERO(...) returns an EXITFLAG that describes the 
%   exit condition of FZERO. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     1  FZERO found a zero X.
%    -1  Algorithm terminated by output function.
%    -3  NaN or Inf function value encountered during search for an interval
%         containing a sign change.
%    -4  Complex function value encountered during search for an interval 
%         containing a sign change.
%    -5  FZERO may have converged to a singular point.
%    -6  FZERO can not detect a change in sign of the function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FZERO(...) returns a structure OUTPUT
%   with the number of function evaluations in OUTPUT.funcCount, the
%   algorithm name in OUTPUT.algorithm, the number of iterations to
%   find an interval (if needed) in OUTPUT.intervaliterations, the
%   number of zero-finding iterations in OUTPUT.iterations, and the
%   exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fzero(@sin,3)
%     returns pi.
%        X = fzero(@sin,3,optimset('Display','iter')) 
%     returns pi, uses the default tolerance and displays iteration information.
%
%     FUN can be an anonymous function:
%        X = fzero(@(x) sin(3*x),2)
%
%     FUN can be a parameterized function.  Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) cos(c.*x);  % The parameterized function.
%        c = 2;                 % The parameter.
%        X = fzero(@(x) myfun(x,c),0.1)
%   
%   Limitations
%        X = fzero(@(x) abs(x)+1, 1) 
%     returns NaN since this function does not change sign anywhere on the 
%     real axis (and does not have a zero as well).
%        X = fzero(@tan,2)
%     returns X near 1.5708 because the discontinuity of this function near the 
%     point X gives the appearance (numerically) that the function changes sign at X.
%
%   See also ROOTS, FMINBND, FUNCTION_HANDLE.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 5.33.4.24 $  $Date: 2010/05/13 17:39:01 $

%  This algorithm was originated by T. J. Dekker.  An Algol 60 version,
%  with some improvements, is given by R. P. Brent in "Algorithms for
%  Minimization Without Derivatives", Prentice-Hall, 1973.  A Fortran
%  version is in Forsythe, Malcolm and Moler, "Computer Methods
%  for Mathematical Computations", Prentice-Hall, 1976.

%  Modified by Jonathan Karr, jkarr@stanford.edu 6/10/2011 for better error
%  reporting.

% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
procedure = ' ';

defaultopt = struct('Display','notify','TolX',eps,'FunValCheck','off','OutputFcn',[],'PlotFcns',[]);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(FunFcnIn,'defaults')
    b = defaultopt;
    return
end

% initialization
if nargin < 3, 
   options = []; 
end
% Detect problem structure input
if nargin == 1
    if isa(FunFcnIn,'struct')
        [FunFcnIn,x,options] = separateOptimStruct(FunFcnIn);
    else % Single input and non-structure.
        error('MATLAB:fzero:InputArg','The input should be either a structure with valid fields or at least two arguments to FZERO.');
    end
end

if nargin == 0 
    error('MATLAB:fzero:NotEnoughInputs',...
        'FZERO requires at least two input arguments.'); 
end

% Check for non-double inputs
if ~isa(x,'double')
  error('MATLAB:fzero:NonDoubleInput', ...
        'FZERO only accepts inputs of data type double.')
end

tol = optimget(options,'TolX',defaultopt,'fast');
funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');
printtype = optimget(options,'Display',defaultopt,'fast');
switch printtype
    case {'notify','notify-detailed'}
        trace = 1;
    case {'none', 'off'}
        trace = 0;
    case {'iter','iter-detailed'}
        trace = 3;
    case {'final','final-detailed'}
        trace = 2;
    otherwise
        trace = 1;
end
% Handle the output functions
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot functions
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

% Convert to function handle as needed.
[FunFcn,errStruct] = fcnchk(FunFcnIn,length(varargin));
if ~isempty(errStruct)
    error('MATLAB:fzero:InvalidFUN',errStruct.message)
end
% We know fcnchk succeeded if we got to here
if isa(FunFcn,'inline')      
    if isa(FunFcnIn,'inline')
        Ffcnstr = inputname(1);  % name of inline object such as f where f=inline('x*2');
        if isempty(Ffcnstr)  % inline('sin(x)')  
            Ffcnstr = formula(FunFcn);  % Grab formula, no argument name 
        end
        Ftype = 'inline object';
    else  % not an inline originally (string expression).
        Ffcnstr = FunFcnIn;  % get the string expression
        Ftype = 'expression';
    end
elseif isa(FunFcn,'function_handle') % function handle
    Ffcnstr = func2str(FunFcn);  % get the name passed in
    Ftype = 'function_handle';
else  % Not converted, must be m-file or builtin
    Ffcnstr = FunFcnIn;  % get the name passed in
    Ftype = 'function';
end

% Add a wrapper function to check for Inf/NaN/complex values
if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
    varargin = {FunFcn, varargin{:}};
    FunFcn = @checkfun;
end

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,[],'init',fcount,iter,intervaliter, ...
        [],procedure,[],[],[],[],varargin{:});
    if stop
        [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  trace > 0
            disp(output.message)
        end
        return;
    end
end

if  ~all(isfinite(x))
    error('MATLAB:fzero:Arg2NotFinite', 'Second argument must be finite.')
end

% Interval input
if (length(x) == 2) 
    if trace > 2
        disp(' ') %Initial blank line
    end
    a = x(1); savea=a;
    b = x(2); saveb=b;
    % Put first feval in try catch
    try
        fa = FunFcn(a,varargin{:});
    catch ME
        if ~isempty(Ffcnstr)
            ME.addCause(MException('MATLAB:fzero:InvalidFunctionSupplied', ...
                'FZERO cannot continue because user supplied %s ==> %s\nfailed with the error below.', ...
                Ftype, Ffcnstr)).rethrow();
        else
            ME.addCause(MException('MATLAB:fzero:InvalidFunctionSupplied', ...
                'FZERO cannot continue because user supplied %s\nfailed with the error below.', ...
                Ftype)).rethrow();
        end
        
    end
    
    fb = FunFcn(b,varargin{:});
    if any(~isfinite([fa fb])) || any(~isreal([fa fb]))
        error('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite',...
            'Function values at interval endpoints must be finite and real.')
    end
    fcount = fcount + 2;
    savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        msg = sprintf('Zero find terminated.');
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        msg = sprintf('Zero find terminated.');
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
        error('MATLAB:fzero:ValuesAtEndPtsSameSign',...
            'The function values at the interval endpoints must differ in sign.')
    end
    
    % Starting guess scalar input
elseif (length(x) == 1)
    if trace > 2 
        disp(' ')
        fprintf('Search for an interval around %g containing a sign change:\n',x);
        header = ' Func-count    a          f(a)             b          f(b)        Procedure';
    end
    % Put first feval in try catch
    try
        fx = FunFcn(x,varargin{:});
    catch ME
        if ~isempty(Ffcnstr)
            es = sprintf(['FZERO cannot continue because user supplied' ...
                    ' %s ==> %s\nfailed with the error below.'],  ...
                Ftype, Ffcnstr);
        else
            es = sprintf(['FZERO cannot continue because user supplied' ...
                    ' %s\nfailed with the error below.'],  ...
                Ftype);
        end
        ME.addCause(MException('MATLAB:fzero:InvalidFunctionSupplied', es)).rethrow();            
    end
    fcount = fcount + 1;  
    if fx == 0
        b = x;
        msg = sprintf('Zero find terminated.');
        if trace > 1
            disp(msg)
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fx;
        return
    elseif ~isfinite(fx) || ~isreal(fx)
        error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
            'Function value at starting guess must be finite and real.');
    end
    
    if x ~= 0, 
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;
    
    if trace > 2
        disp(header)
        procedure='initial interval';
        disp(sprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s',fcount,a,fa,b,fb, procedure));
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
            fx,procedure,a,fa,b,fb,varargin{:}); % a and b are x to start
        if stop
            [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  trace > 0
                disp(output.message)
            end
            return;
        end
    end

    while (fa > 0) == (fb > 0)
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  fa = FunFcn(a,varargin{:});
        fcount = fcount + 1;
        if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
            [exitflag,msg] = disperr(a,fa,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end

        if (fa > 0) ~= (fb > 0) % check for different sign
            % Before we exit the while loop, print out the latest interval
            if trace > 2
                procedure='search';
                disp(sprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s',fcount,a,fa,b,fb, procedure));
            end
            % OutputFcn and PlotFcns call
            if haveoutputfcn || haveplotfcn
                [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
                    fx,procedure,a,fa,b,fb,varargin{:});
                if stop
                    [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
                    if  trace > 0
                        disp(output.message)
                    end
                    return;
                end
            end
            break
        end
        
        b = x + dx;  fb = FunFcn(b,varargin{:});
        if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
            [exitflag,msg] = disperr(b,fb,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end
        fcount = fcount + 1;        
        if trace > 2
            procedure='search';
            disp(sprintf('%5.0f   %13.6g %13.6g %13.6g %13.6g   %s',fcount,a,fa,b,fb, procedure));
        end
        % OutputFcn and PlotFcns call
        if haveoutputfcn || haveplotfcn
            [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,'iter',fcount,iter,intervaliter, ...
                fx,procedure,a,fa,b,fb,varargin{:});
            if stop
                [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
                if  trace > 0
                    disp(output.message)
                end
                return;
            end
        end
    end % while
    
    if trace > 2
        disp(' ')
        fprintf('Search for a zero in the interval [%g, %g]:\n',a,b);
    end
    savea = a; savefa = fa; saveb = b; savefb = fb;
else
    error('MATLAB:fzero:LengthArg2', 'Second argument must be of length 1 or 2.');
end % if (length(x) == 2

fc = fb;
procedure = 'initial';
header2 = ' Func-count    x          f(x)             Procedure';
if trace > 2
    disp(header2)
end
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
    if trace > 2
        disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
    end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,b,'iter',fcount,iter,intervaliter, ...
            fb,procedure,savea,savefa,saveb,savefb,varargin{:});
        if stop
            [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  trace > 0
                disp(output.message)
            end
            return;
        end
    end
    
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value

% Output last chosen b
if trace > 2
    disp(sprintf('%5.0f   %13.6g %13.6g        %s',fcount, b, fb, procedure));
end

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,b,'iter',fcount,iter,intervaliter, ...
        fb,procedure,savea,savefa,saveb,savefb,varargin{:});
    if stop
        [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  trace > 0
            disp(output.message)
        end
        return;
    end
end

output.intervaliterations = intervaliter;
output.iterations = iter;
output.funcCount = fcount;
output.algorithm = 'bisection, interpolation';

if abs(fval) <= max(abs(savefa),abs(savefb))
    msg = sprintf('Zero found in the interval [%g, %g]',savea,saveb);
else
    exitflag = -5; 
    msg = sprintf([...
        'Current point x may be near a singular point. The interval [%g, %g] \n', ... 
        'reduced to the requested tolerance and the function changes sign in the interval,\n', ...
        'but f(x) increased in magnitude as the interval reduced.'],savea,saveb);
end
if trace > 1
    disp(' ')
    disp(msg)
end
output.message = msg;
% Outputfcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,b,'done',fcount,iter,intervaliter,fval,procedure,savea,savefa,saveb,savefb,varargin{:});
end


%------------------------------------------------------------------

function [exitflag,msg] = disperr(y, fy, trace)
%DISPERR Display an appropriate error message when FY is Inf, 
%   NaN, or complex.  Assumes Y is the value and FY is the function 
%   value at Y. If FY is neither Inf, NaN, or complex, it generates 
%   an error message.

if ~isfinite(fy)  % NaN or Inf detected
    exitflag = -3;
    msg = ...
        sprintf(['Exiting fzero: aborting search for an interval containing a sign change\n' ...
                 '    because NaN or Inf function value encountered during search.\n' ...
                 '(Function value at %g is %g.)\n' ...
                 'Check function or try again with a different starting value.'],y,fy);
    if trace > 0
        disp(msg)
    end
elseif ~isreal(fy) % Complex value detected
    exitflag = -4;
    msg = ...
        sprintf(['Exiting fzero: aborting search for an interval containing a sign change\n' ...
                 '    because complex function value encountered during search.\n' ...
                 '(Function value at %g is %s.)\n' ...
                 'Check function or try again with a different starting value.'],y,num2str(fy));
    if trace > 0
        disp(msg)        
    end
elseif ~isfinite(y) % Inf detected in bracketting stage
    exitflag = -6;
    msg = ...
        sprintf(['Exiting fzero: aborting search for an interval containing a sign change\n' ...
                 '    because no sign change is detected during search.\n' ...
                 'Function may not have a root.']);
    if trace > 0
        disp(msg)        
    end

else
    error('MATLAB:fzero:disperr:InvalidArg',...
        'DISPERR (in FZERO) called with invalid argument.')
end

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,state,fcount,iter,intervaliter,  ...
    f,procedure,a,fvala,b,fvalb,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter', or 'done'. 
% We do not handle the case 'interrupt' because we do not want to update
% xOutputfcn or optimValues (since the values could be inconsistent) before calling
% the outputfcn; in that case the outputfcn/plotfcns are called directly rather than
% calling it inside callOutputAndPlotFcns.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.funccount = fcount;
optimValues.iteration = iter;
optimValues.intervaliteration = intervaliter;
optimValues.fval = f;
optimValues.procedure = procedure;
optimValues.intervala = a;
optimValues.fvala = fvala;
optimValues.intervalb = b;
optimValues.fvalb = fvalb;

xOutputfcn = x;  % set xOutputfcn to be x
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fzero:InvalidState','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fzero:InvalidState','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end

%--------------------------------------------------------------------------
function [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization 
% is interrupted.

b = xOutputfcn;
fval = optimValues.fval;
exitflag = -1; 
output.intervaliterations = optimValues.intervaliteration;
output.iterations = optimValues.iteration;
output.funcCount = optimValues.funccount;
output.algorithm = 'bisection, interpolation';
output.message = 'Optimization terminated prematurely by user.';

%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FZERO handles it naturally.  ???
if isnan(f)
    error('MATLAB:fzero:checkfun:NaNFval', ...
        'User function ''%s'' returned NaN when evaluated at %g;\n FZERO cannot continue.', ...
        localChar(userfcn), x);  
elseif ~isreal(f)
    error('MATLAB:fzero:checkfun:ComplexFval', ...
        'User function ''%s'' returned a complex value when evaluated at %g;\n FZERO cannot continue.', ...
        localChar(userfcn),x);  
end

%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a string for printing

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch ME
        strfcn = '(name not printable)';
    end
end
