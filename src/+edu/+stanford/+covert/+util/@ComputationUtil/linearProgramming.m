function varargout = linearProgramming(maximizeFlag, f, A, b, lb, ub, constraintTypes, variableTypes, options)
% Provides a common interface to several linear programming solvers.
% Checks each solver for errors.
%
% That is it solves the problem:
%   max f'*x
%   subject to A*x~b, lb<=x<=ub
% where "~" can be equality, or double or single-sided inequality
% constraints if the chosen linear programming solver supports those
% types of constraints. Additionally the variables, x, can, if the chosen
% linear programming solver supports it, either continuous, integer, and
% boolean valued.
%
% Constraint Types: etiher a single value, or a column vector or values:
%   'F' Free (unbounded) variable (the constraint is ignored).
%   'U' Variable with upper bound ( A(i,:)*x <= b(i)).
%   'S' Fixed Variable (A(i,:)*x = b(i)).
%   'L' Variable with lower bound (A(i,:)*x >= b(i)).
%   'D' Double-bounded variable (A(i,:)*x >= -b(i) and A(i,:)*x <= b(i)).
%
% Variable Types: either a single value, or a column vector or values:
%   'C' Continuous variable.
%   'I' Integer variable
%   'B' Binary variable
%
% Output
%  x
%  lambda, 
%  fopt
%  errorFlag
%  errorMsg
%  extra
%
% Requirements: at least one of the following linear programing solvers
% - glpk -- GNU Linear programming kit, to solve linear programming problems
%   http://sourceforge.net/projects/glpkmex/
% - lp_solve -- open source solve linear programming solver
%   http://sourceforge.net/projects/lpsolve
% - clp -- open source solve linear programming solver
%   http://control.ee.ethz.ch/~joloef/mexclp.zip
% - bpmpd -- open source solve linear programming solver
%   http://www.pserc.cornell.edu/bpmpd/
% - qsopt -- open source solve linear programming solver
%   http://control.ee.ethz.ch/~joloef/mexqsopt.msql
% - Optimization Toolbox - has linprog, MATLAB's linear programming
%   solver. Note: linprog is not vevy good. We suggest you use another
%   solver instead.
%   http://www.mathworks.com/products/optimization/
%
% Author: Jonathan Karr
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/31/2008

%- A is non-empty
%- b is real column vector of length ncols(A)
%- lb, ub are real column vectors of length nrows(A)
validateattributes(A, {'numeric'}, {'nonempty'});
validateattributes(b, {'numeric'}, {'real', 'size', [size(A, 1) 1]});
if isempty(lb)
    lb = -Inf(size(A, 2), 1);
else
    validateattributes(lb, {'numeric'}, {'real', 'size', [size(A, 2) 1]});
end
if isempty(ub)
    ub = -Inf(size(A, 2), 1);
else
    validateattributes(ub, {'numeric'}, {'real', 'size', [size(A, 2) 1]});
end

%error check length sizes of contraint and variable types, expand if
%necessary
if isscalar(constraintTypes)
    if ~ismembc(constraintTypes, 'DFLSUdflsu')
        throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid constraints types'));
    end
    constraintTypes = constraintTypes(ones(length(b), 1), 1);
elseif length(constraintTypes) ~= length(b)
    throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid constraints types'));
elseif ~all(ismembc(constraintTypes, 'DFLSUdflsu'))
    throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid constraints types'));
end

if isscalar(variableTypes)
    if ~ismembc(variableTypes, 'BCIbci')
        throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid variable types'));
    end
    variableTypes = variableTypes(ones(size(A, 2), 1), 1);
elseif length(variableTypes) ~= length(lb)
    throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid variable types'));
elseif ~all(ismembc(variableTypes, 'BCIbci'))
    throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid variable types'));
end

%compute sense from maximizeFlag
%if maximize='maximize'->sense=-1
%otherwise->sense=1
sense = 1 - 2 *strcmp(maximizeFlag, 'maximize');

%solver linear programming problem using chosen solver
varargout = cell(max(nargout, 1), 1);
switch options.solver(1:2)
    case 'gl'; [varargout{:}] = runglpk(   f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'li'; [varargout{:}] = runlinprog(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);    
    case 'bp'; [varargout{:}] = runbpmpd(  f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'cl'; [varargout{:}] = runclp(    f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'lp'; [varargout{:}] = runlpsolve(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'qs'; [varargout{:}] = runqsopt(  f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'to'; [varargout{:}] = runtomlab( f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'cp'; [varargout{:}] = runcplex(  f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'gu'; [varargout{:}] = rungurobi( f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    case 'mo'; [varargout{:}] = runmosek(  f, A, b, lb, ub, sense, constraintTypes, variableTypes, options);
    otherwise; throw(MException('ComputationUtil:LinearProgramming', 'Support for %s not implemented', options.solver));
end

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runlinprog(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)

%setup equality, inequality constraints based on constraint types
Aeq = [];
beq = [];
if sum(constraintTypes == 'S') == length(constraintTypes)
    Aeq = A;
    A = [];
    beq = b;
    b = [];
elseif sum(constraintTypes == 'U') == length(constraintTypes)
elseif sum(constraintTypes == 'L') == length(constraintTypes)
    A = -A;
    b = -b;
else
    throw(MException('ComputationUtil:LinearProgramming', 'linprog: constraint types invalid'));
end

%throw error if integer or boolean variable types requested, linprog
%doesn't have this capability
if sum(variableTypes == 'C') ~= length(variableTypes)
    throw(MException('ComputationUtil:LinearProgramming', 'linprog: variable types invalid'));
end

%linprog options
if isfield(options.solverOptions, 'linprog')
    linprogoptions = options.solverOptions.linprog;
else
    linprogoptions = struct;
end

try
    [x, fopt, exitflag, ~, lambda] = linprog(sense * f, A, b, Aeq, beq, lb, ub, [], linprogoptions);
 
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', lambda);    
    fopt = fopt * sense;
    errorFlag = exitflag ~= 1;
    switch exitflag
        case  1; errorMsg = 'Function converged to a solution x.';
        case  0; errorMsg = 'Number of iterations exceeded options.MaxIter.';
        case -2; errorMsg = 'No feasible point was found.';
        case -3; errorMsg = 'Problem is unbounded.';
        case -4; errorMsg = 'NaN value was encountered during execution of the algorithm.';
        case -5; errorMsg = 'Both primal and dual problems are infeasible.';
        case -7; errorMsg = 'Search direction became too small. No further progress could be made.';
    end
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;    
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runglpk(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)

%make structure of glpk options
if isfield(options.solverOptions, 'glpk')
    glpkoptions = options.solverOptions.glpk;
else
    glpkoptions = struct;
end

%call glpk
try
    [x, fopt, status, extra] = glpkcc(f, A, b, lb, ub, constraintTypes, variableTypes, sense, glpkoptions);
    lambda = struct('reducedCosts', extra.redcosts, 'shadowPrices', extra.lambda);
    errorFlag = ~(status == 5 || status == 2);
    switch status
        %General Errors
        case 5; errorMsg = 'Solution is optimal';
        case 2; errorMsg = 'Solution is feasible';
        case 1; errorMsg = 'Solution is undefined';
        case 3; errorMsg = 'Solution is infeasible';
        case 4; errorMsg = 'No feasible solution exists';
        case 6; errorMsg = 'Solution is unbounded';
            
        %Simplex method Errors:
        case 101; errorMsg = 'Invalid basis';
        case 102; errorMsg = 'Singular matrix';
        case 103; errorMsg = 'Ill-conditioned matrix';
        case 104; errorMsg = 'Invalid bounds';
        case 105; errorMsg = 'Solver failed';
        case 106; errorMsg = 'Objective lower limit reached';
        case 107; errorMsg = 'Objective upper limit reached';
        case 108; errorMsg = 'Iteration limit exceeded';
        case 109; errorMsg = 'Time limit exceeded';
        case 110; errorMsg = 'No primal feasible solution';

        %Interior point method, mixed integer problem Errors:
        case 204; errorMsg = 'Unable to start the search.';
        case 205; errorMsg = 'Objective function lower limit reached.';
        case 206; errorMsg = 'Objective function upper limit reached.';
        case 207; errorMsg = 'Iterations limit exhausted.';
        case 208; errorMsg = 'Time limit exhausted.';
        case 209; errorMsg = 'No feasible solution.';
        case 210; errorMsg = 'Numerical instability.';
        case 211; errorMsg = 'Problems with basis matrix.';
        case 212; errorMsg = 'No convergence (interior).';
        case 213; errorMsg = 'No primal feasible solution (LP presolver).';
        case 214; errorMsg = 'No dual feasible solution (LP presolver).';
            
        otherwise, errorMsg = sprintf('Invalid error code %d', status);
    end
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runlpsolve(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)
%options
verbose = 3;
presolve = 0; %none
scaling = 4 + 64 + 128; %geometric + equilibrate + integers
if isfield(options.solverOptions, 'lp_solve')
    if isfield(options.solverOptions.lp_solve, 'verbose'), verbose = options.solverOptions.lp_solve.verbose; end
    if isfield(options.solverOptions.lp_solve, 'presolve'), presolve = options.solverOptions.lp_solve.presolve; end
    if isfield(options.solverOptions.lp_solve, 'scaling'), scaling = options.solverOptions.lp_solve.scaling; end    
end

%convert variable types to lp_solve format'
if any(variableTypes == 'B'); throw(MException('ComputationUtil:LinearProgramming', 'lp_solve: boolean variables are not allowed')); end;
xint = NaN(size(variableTypes));
xint(variableTypes == 'I') = 1;
xint(variableTypes == 'C') = 0;

%convert constraints types to lp_solve format
if any(constraintTypes == 'F'); throw(MException('ComputationUtil:LinearProgramming', 'lp_solve: free constraints are not allowed')); end;
if any(constraintTypes == 'D'); throw(MException('ComputationUtil:LinearProgramming', 'lp_solve: double-sided constraints are not allowed')); end;
con_types = NaN(size(constraintTypes));
con_types(constraintTypes == 'L') = 2;
con_types(constraintTypes == 'S') = 3;
con_types(constraintTypes == 'U') = 1;

try
    lp = mxlpsolve('make_lp', size(A, 1), size(A, 2));
    mxlpsolve('set_verbose', lp, verbose);
    mxlpsolve('set_mat', lp, A);
    mxlpsolve('set_rh_vec', lp, b);
    mxlpsolve('set_obj_fn', lp, f);
    mxlpsolve('set_sense', lp, (sense == -1) + 0);
    mxlpsolve('set_constr_type', lp, con_types);
    mxlpsolve('set_bounds', lp, lb, ub);
    mxlpsolve('set_int', lp, xint);    
    mxlpsolve('set_presolve', lp, presolve, mxlpsolve('get_presolveloops', lp));
    mxlpsolve('set_scaling', lp, scaling);
        
    status = mxlpsolve('solve', lp);
    errorFlag =~ status == 0 || status == 1 || status == 11 || status == 12;
    if ~errorFlag
        [fopt, x, duals] = mxlpsolve('get_solution', lp);
        reducedCosts = mxlpsolve('get_reduced_costs', lp);
    else
        fopt = NaN;
        x = NaN(size(lb));
        duals = NaN(size(A, 1), 1);
        reducedCosts = NaN(size(A, 2), 1);        
    end
    
    extra = struct;
    if nargin >= 6
        tmpFileName = ['lpsolve.' datestr(now, 30) '.txt'];
        
        mxlpsolve('set_outputfile', lp, tmpFileName);
        mxlpsolve('print_scales', lp);
        
        fid = fopen(tmpFileName, 'r');
        fgetl(fid);
        fgetl(fid);
        line = fgetl(fid);
        extra.objScaling = str2double(line(32:end));
        extra.rowScaling = zeros(size(A, 1), 1);
        extra.colScaling = zeros(size(A, 2), 1);
        for i = 1:size(A, 1)
            line = fgetl(fid);
            extra.rowScaling(i) = str2double(line(32:end));
        end
        for i = 1:size(A, 2)
            line = fgetl(fid);
            extra.colScaling(i) = str2double(line(32:end));
        end
        fclose(fid);                
    end
    
    mxlpsolve('delete_lp', lp);
    
    if nargin >= 6
        delete(tmpFileName);
    end
           
    lambda = struct('reducedCosts', reducedCosts, 'shadowPrices', duals);
    switch status
        case -2; errorMsg = 'Out of memory';
        case 0; errorMsg = 'An optimal solution was obtained';
        case 1; errorMsg = ['The model is sub-optimal. Only happens if there are integer variables and there is already an integer solution found. The solution is not guaranteed the most optimal one.\n'...
                '* A timeout occured (set via set_timeout or with the -timeout option in lp_solve)\n'...
                '* set_break_at_first was called so that the first found integer solution is found (-f option in lp_solve)\n'...
                '* set_break_at_value was called so that when integer solution is found that is better than the specified value that it stops (-o option in lp_solve)\n'...
                '* set_mip_gap was called (-g/-ga/-gr options in lp_solve) to specify a MIP gap\n'...
                '* An abort function is installed (put_abortfunc) and this function returned TRUE\n'...
                '* At some point not enough memory could not be allocated'];
        case 2; errorMsg = 'The model is infeasible';
        case 3; errorMsg = 'The model is unbounded';
        case 4; errorMsg = 'The model is degenerative';
        case 5; errorMsg = 'Numerical failure encountered';
        case 6; errorMsg = 'The abort routine returned TRUE. See put_abortfunc';
        case 7; errorMsg = 'A timeout occurred. A timeout was set via set_timeout';
        case 9; errorMsg = 'The model could be solved by presolve. This can only happen if presolve is active via set_presolve';
        case 10; errorMsg = 'The B&B routine failed';
        case 11; errorMsg = 'The B&B was stopped because of a break-at-first (see set_break_at_first) or a break-at-value (see set_break_at_value)';
        case 12; errorMsg = 'A feasible B&B solution was found';
        case 13; errorMsg = 'No feasible B&B solution found';
    end
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
    extra = struct;
end

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runqsopt(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)

%setup equality, inequality constraints based on constraint types
Aeq = [];
beq = [];
if sum(constraintTypes == 'S') == length(constraintTypes)
    Aeq = A;
    A = [];
    beq = b;
    b = [];
elseif sum(constraintTypes == 'U') == length(constraintTypes)
elseif sum(constraintTypes == 'L') == length(constraintTypes)
    A = -A;
    b = -b;
else
    throw(MException('ComputationUtil:LinearProgramming', 'qsopt: constraint types invalid'));
end

%throw error if integer or boolean variable types requested, qsopt
%doesn't have this capability
if sum(variableTypes == 'C') ~= length(variableTypes)
    throw(MException('ComputationUtil:LinearProgramming', 'qsopt: variable types invalid'));
end

%make structure of qsopt options
qsoptoptions = struct;
if isfield(options.solverOptions, 'qsopt')
    qsoptoptions = options.solverOptions.qsopt;
end

try
    [x, dual, status] = qsopt(sense * f, A, b, Aeq, beq, lb, ub, qsoptoptions);
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', dual);
    fopt = f' * x;
    errorFlag = status ~= 1;
    switch status
        case 1; errorMsg = 'optimal';
        case 2; errorMsg = 'infeasible';
        case 3; errorMsg = 'unbounded';
        case 4; errorMsg = 'iteration limit';
        case 5; errorMsg = 'time limit';
        case 6; errorMsg = 'other problem';
    end
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runclp(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)

%setup equality, inequality constraints based on constraint types
Aeq = [];
beq = [];
if sum(constraintTypes == 'S') == length(constraintTypes)
    Aeq = A;
    A = [];
    beq = b;
    b = [];
elseif sum(constraintTypes == 'U') == length(constraintTypes)
elseif sum(constraintTypes == 'L') == length(constraintTypes)
    A = -A;
    b = -b;
else
    throw(MException('ComputationUtil:LinearProgramming', 'clp: constraint types invalid'));
end

%throw error if integer or boolean variable types requested, clp
%doesn't have this capability
if sum(variableTypes == 'C') ~= length(variableTypes)
    throw(MException('ComputationUtil:LinearProgramming', 'clp: variable types invalid'));
end

%make structure of clp options
clpoptions = struct;
if isfield(options.solverOptions, 'clp')
    clpoptions = options.solverOptions.clp;
end

try
    [x, dual, status] = clp(zeros(length(f)), sense * f, A, b, Aeq, beq, lb, ub, clpoptions);
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', dual);
    fopt = f' * x;
    errorFlag = status ~= 0;
    switch status
        case 0; errorMsg = 'optimal';
        case 1; errorMsg = 'infeasible';
        case 2; errorMsg = 'unbounded';
    end
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runbpmpd(f, A, b, lb, ub, sense, constraintTypes, variableTypes, ~)

%throw error if integer or boolean variable types requested, bpmpd
%doesn't have this capability
if sum(variableTypes == 'C') ~= length(variableTypes)
    throw(MException('ComputationUtil:LinearProgramming', 'bpmpd: variable types invalid'));
end

%convert constraints types to bpmpd format
if any(constraintTypes == 'F'); throw(MException('ComputationUtil:LinearProgramming', 'bpmpd: free constraints are not allowed')); end;
if any(constraintTypes == 'D'); throw(MException('ComputationUtil:LinearProgramming', 'bpmpd: double-sided constraints are not allowed')); end;
e = NaN(size(constraintTypes));
e(constraintTypes == 'L') = 1;
e(constraintTypes == 'S') = 0;
e(constraintTypes == 'U') = -1;

%throw error if integer or boolean variable types requested, bpmpd
%doesn't have this capability
if sum(variableTypes == 'C') ~= length(variableTypes)
    throw(MException('ComputationUtil:LinearProgramming', 'bpmpd: variable types invalid'));
end

try
    llist = (1:length(lb))';
    ulist = (1:length(ub))';
    [x, dual, ~, w, errorMsg] = bp(zeros(length(f)), A, b, sense * f, e, llist, lb, ulist, ub, bpopt,0);
    lambda = struct('reducedCosts', w, 'shadowPrices', dual);
    fopt = f' * x;
    errorFlag =~ strcmp(errorMsg, 'optimal solution');
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runtomlab( f, A, b, lb, ub, sense, constraintTypes, ~, options)

b_L = -Inf(size(b));
b_U =  Inf(size(b));
b_L(constraintTypes == 'L') = b(constraintTypes == 'L');
b_U(constraintTypes == 'U') = b(constraintTypes == 'U');
b_L(constraintTypes == 'S') = b(constraintTypes == 'S');
b_U(constraintTypes == 'S') = b(constraintTypes == 'S');
b_L(constraintTypes == 'D') = -b(constraintTypes == 'D');
b_U(constraintTypes == 'D') =  b(constraintTypes == 'D');

prob = lpAssign(sense * f, A, b_L, b_U, lb, ub, [], [], [], [], [], [], [], [], []);

solver = 'minos';
prob.optParam = struct();
PriLev = 0;
if isfield(options.solverOptions, 'tomlab')
    if isfield(options.solverOptions.tomlab, 'solver'), solver = options.solverOptions.tomlab.solver; end    
    if isfield(options.solverOptions.tomlab, 'optParam'), prob.optParam = options.solverOptions.tomlab.optParam; end
    if isfield(options.solverOptions.tomlab, 'PriLev'), PriLev = options.solverOptions.tomlab.PriLev; end
end

try    
    result = tomRun(solver, prob, PriLev);
    errorFlag = result.ExitFlag ~= 0;
    errorMsg = result.ExitText;
    fopt = sense * result.f_k;
    x = result.x_k;
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runcplex(f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)
constraintTypes = upper(constraintTypes);
if ~ismembc(constraintTypes, 'LSU')
    throw(MException('ComputationUtil:LinearProgramming', 'linear programming: invalid constraints types'));
end

Aineq = [
    -A(constraintTypes == 'L', :)
    A(constraintTypes == 'U', :)
    ];
bineq = [
    -b(constraintTypes == 'L', :)
    b(constraintTypes == 'U', :)
    ];
    
Aeq = A(constraintTypes == 'S', :);
beq = b(constraintTypes == 'S', :);

try
    [x, fopt, exitFlag, output] =  cplexmilp(...
        sense * f, ...
        Aineq, bineq, ...
        Aeq, beq, ...
        [], [], [], lb, ub, variableTypes');
    errorFlag = exitFlag ~= 1;
    errorMsg = output.message;
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    fopt = NaN;
end
lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));

extra = struct('status', output.cplexstatus);

function [x, lambda, fopt, errorFlag, errorMsg, extra] = rungurobi( f, A, b, lb, ub, sense, constraintTypes, variableTypes, options)
b_L = -Inf(size(b));
b_U =  Inf(size(b));
b_L(constraintTypes == 'L') = b(constraintTypes == 'L');
b_U(constraintTypes == 'U') = b(constraintTypes == 'U');
b_L(constraintTypes == 'S') = b(constraintTypes == 'S');
b_U(constraintTypes == 'S') = b(constraintTypes == 'S');
b_L(constraintTypes == 'D') = -b(constraintTypes == 'D');
b_U(constraintTypes == 'D') =  b(constraintTypes == 'D');

xint = NaN(size(variableTypes));
xint(variableTypes == 'I') = 1;
xint(variableTypes == 'C') = 0;

PriLev = 0;
grbControl = struct();
if isfield(options.solverOptions, 'gurobi')
    if isfield(options.solverOptions.gurobi, 'PriLev'), PriLev = options.solverOptions.gurobi.PriLev; end
    if isfield(options.solverOptions.gurobi, 'grbControl'), grbControl = options.solverOptions.gurobi.grbControl; end
end

try    
    [x, ~, v, rc, fopt, ~, ~, status] = ...
        gurobi(sense * f, A, lb, ub, b_L, b_U, [], ...
        grbControl, PriLev, xint);
    fopt = fopt * sense;
    [errorMsg, exitFlag] = grbStatus(status);
    errorFlag = exitFlag ~= 0;
    lambda = struct('reducedCosts', rc, 'shadowPrices', v);
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;

function [x, lambda, fopt, errorFlag, errorMsg, extra] = runmosek( f, A, b, lb, ub, sense, constraintTypes, ~, options)
if sense == 1
    cmd = 'minimize';
else
    cmd = 'maximize';
end

b_L = -Inf(size(b));
b_U =  Inf(size(b));
b_L(constraintTypes == 'L') = b(constraintTypes == 'L');
b_U(constraintTypes == 'U') = b(constraintTypes == 'U');
b_L(constraintTypes == 'S') = b(constraintTypes == 'S');
b_U(constraintTypes == 'S') = b(constraintTypes == 'S');
b_L(constraintTypes == 'D') = -b(constraintTypes == 'D');
b_U(constraintTypes == 'D') =  b(constraintTypes == 'D');

lb = max(lb, -1e7);
ub = min(ub,  1e7);

opts = struct(...
    'MSK_IPAR_LOG', 0, 'MSK_IPAR_MAX_NUM_WARNINGS', 1e5, ...
    'MSK_IPAR_INTPNT_SCALING', 0, 'MSK_IPAR_SIM_SCALING', 0, 'MSK_IPAR_SIM_SCALING_METHOD', 1);
if isfield(options.solverOptions, 'mosek')
    opts = options.solverOptions.mosek;
end

try   
    result = msklpopt(f, A, b_L, b_U, lb, ub, opts, cmd);
    fopt = result.sol.itr.pobjval;
    x = result.sol.itr.xx;
    errorFlag = result.rcode ~= 0;
    errorMsg = result.rmsg;
    lambda = struct(...
        'reducedCosts', result.sol.itr.slx + result.sol.itr.sux, ...
        'shadowPrices', result.sol.itr.slc + result.sol.itr.suc);
catch exception
    errorFlag = 1;
    errorMsg = exception.message;
    x = NaN(size(lb));
    lambda = struct('reducedCosts', NaN(size(A, 2), 1), 'shadowPrices', NaN(size(A, 1), 1));
    fopt = NaN;
end

extra = struct;