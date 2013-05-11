function [x, objective, errorFlag, errorMsg] = quadraticProgramming(H, f, A, b, Aeq, beq, lb, ub, x0, options)
% Provides a common interface to several quadratic programming solvers.
% Checks each solver for errors.
%
% That is it solves the problem:
%   min 0.5 * x'*H*x + f'*x
%   subject to A*x<=b, Aeq*x=beq, lb<=x<=ub
%
% Requirements: at least one of the following quadratic programing solvers
% - clp -- open source solve quadratic programming solver
%   http://control.ee.ethz.ch/~joloef/mexclp.zip
% - minq
%   http://www.mat.univie.ac.at/~neum/software/minq/
% - qpip -- open source solve quadratic programming solver
%   http://sigpromu.org/quadprog/index.html
% - Optimization Toolbox - has quadprog, MATLAB's quadratic programming
%   solver. Note: quadprog is not very good. We suggest you use another
%   solver instead.
%   http://www.mathworks.com/products/optimization/
% - Wolf
%   http://www.mathworks.com/matlabcentral/fileexchange/27397-quadratic-programming-by-wolfs-method
%
% Author: Jonathan Karr
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/17/2011

solver = 'clp';
if exist('options', 'var') && isstruct(options)
    solver = options.solver;
end

solverOptions = struct;
if exist('options', 'var') && isstruct(options) && isfield(options, 'solverOptions') && isfield(options.solverOptions, options.solver)
    solverOptions = options.solverOptions.(options.solver);
end

switch solver
    case 'clp'
        [x, objective, exitFlag] = clp(H, f, A, b, Aeq, beq, lb, ub, solverOptions);
        switch exitFlag
            case 0, errorMsg = 'optimal';
            case 1, errorMsg = 'infeasible';
            case 2, errorMsg = 'unbounded';
        end
        errorFlag = exitFlag ~= 0;
    case 'minqsep'
        [m, n] = size(Aeq);
        printLevel = 1;
        if isfield(solverOptions, 'printLevel'), printLevel = solverOptions.printLevel; end;
        [x, objective, errorFlag] = minqsep(f, diag(H), [Aeq; eye(n); -eye(n)], ...
            [beq; lb; -ub], [true(m, 1); false(n, 1); false(n, 1)], printLevel, x0);
        switch errorFlag
            case 0, errorMsg = 'global minimizer found';
            case 1, errorMsg = 'approximate solution; feasible set probably empty';
            case 99, errorMsg = 'approximate solution; maxit exceeded';
        end
    case 'qpip'
        display = 0;
        mu = 0.0;
        method = 1;
        if isfield(solverOptions, 'display'), display = solverOptions.display; end;
        if isfield(solverOptions, 'mu'),      mu      = solverOptions.mu;      end;
        if isfield(solverOptions, 'method'),  method  = solverOptions.method;  end;
        [x, exitFlag] = qpip(H, f, A, b, Aeq, beq, lb, ub, display, mu, method);
        objective = x' *H * x + f' * x;
        errorFlag = exitFlag ~= 0;
        switch exitFlag
            case 0,    errorMsg = 'optimal';
            otherwise, errorMsg = 'other error';
        end
    case 'quadprog'
        [x, objective, exitFlag, output] = quadprog(H, f, A, b, Aeq, beq, lb, ub, x0, solverOptions);
        errorFlag = exitFlag <= 0;
        errorMsg = output.message;
    case 'wolf'
        [m, n] = size(Aeq);
        [x, objective] = wolf(H, f, [beq; lb; ub],[Aeq; eye(n); eye(n)], ...
            [zeros(m, 1); ones(n, 1); -ones(n, 1)], 1);
        errorFlag = 0;
        errorMsg = 'optimal';
    otherwise
        throw(MException('ComputationUtil:error', 'Invalid solver'));
end