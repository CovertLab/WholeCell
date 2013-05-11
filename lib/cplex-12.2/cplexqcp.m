function [x,fval,exitflag,output,lambda]=cplexqcp(H, f, Aineq, bineq, Aeq, beq, l, Q, r, lb, ub, x0, options)
%%
% Purpose
% Solve quadratically constrained linear/quadratic programming problems.
%
% Syntax
%    x = cplexqcp(H,f,Aineq,bineq)
%    x = cplexqcp(H,f,Aineq,bineq,Aeq,beq)
%    x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r)
%    x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub)
%    x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,x0)
%    x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,x0,options)
%    x = cplexqcp(problem)
%    [x,fval] = cplexqcp(...)
%    [x,fval,exitflag] = cplexqcp(...)
%    [x,fval,exitflag,output] = cplexqcp(...)
%    [x,fval,exitflag,output,lambda] = cplexqcp(...)
%
% Description
% Finds the minimum of a problem specified by
%    min      0.5*x'*H*x+f*x or f*x
%    st.      Aineq*x      <= bineq
%             Aeq*x         = beq
%             l*x + x'*Q*x <= r
%             lb <= x <= ub
%
% f, bineq, beq, l, lb and ub are column vectors.
% H, Aineq, Aeq, and Q are matrices.
%
% x = cplexqcp(H,f,Aineq,bineq) solves the quadratically constrained
% linear/quadratic programming problem min 1/2*x'*H*x + f*x subject to
% Aineq*x <= bineq. If no quadratic objective term exists, set H=[].
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq) solves the preceding problem while
% additionally satisfying the equality constraints Aeq*x = beq. If no
% inequalities exist, set Aineq=[] and bineq=[].
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r) solves the preceding problem
% while additionally satisfying the quadratic inequality constraints
% l*x + x'*Q*x <= r. If no equalities exist, set Aeq=[] and beq=[].
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub) defines a set of lower
% and upper bounds on the design variables, x, so that the solution is in
% the range lb <= x <= ub. If no quadratic inequalities exist, set l=[],
% Q=[] and r=[].
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,x0) sets the starting
% point to x0. If no bounds exist, set lb=[] and ub=[].
%
% x = cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,x0,options) minimizes
% with the default optimization options replaced by values in the structure
% options, which can be created using the function cplexoptimset. If you do
% not wish to give an initial point, set x0=[].
%
% x = cplexqcp(problem) where problem is a structure.
%
% [x,fval] = cplexqcp(...) returns the value of the objective function at
% the solution x: fval = 0.5*x'*H*x + f*x.
%
% [x,fval,exitflag] = cplexqcp(...) returns a value exitflag that describes
% the exit condition of cplexqcp.
%
% [x,fval,exitflag,output] = cplexqcp(...) returns a structure output that
% contains information about the optimization.
%
% [x,fval,exitflag,output,lambda] = cplexqcp(...) returns a structure
% lambda whose fields contain the Lagrange multipliers at the solution x.
%
% Input Arguments
% H         Double matrix for objective function
% f         Double column vector for objective function
% Aineq     Double matrix for linear inequality constraints
% bineq     Double column vector for linear inequality constraints
% Aeq       Double matrix for linear equality constraints
% beq       Double column vector for linear equality constraints
% l         Double column vector or matrix
%           Linear part of quadratic constraints
% Q         Double matrix or double matrix cell for quadratic
%           constraints
% r         Double or double row vector
%           Rhs of quadratic inequality constraints
% lb        Double column vector of lower bounds
% ub        Double column vector of upper bounds
% x0        Double column vector for initial point of x
% options   Options structure created with cplexoptimset
%
% problem   Structure containing the following fields:
%           H         Double matrix for objective function
%           f         Double column vector for objective function
%           Aineq     Double matrix for linear inequality constraints
%           bineq     Double column vector for linear inequality 
%                     constraints
%           Aeq       Double matrix for linear equality constraints
%           beq       Double column vector for linear equality constraints
%      	    qc 	      Struct vector
%  	        qc(i).a 	 Double column vector for linear part of the 
%                        quadratic constraint
%  	        qc(i).rhs    Double for righthand side for quadratic constraint
%  	        qc(i).Q 	 Double matrix for quadratic part of the
%                        quadratic constraint
%           lb        Double column vector of lower bounds
%           ub        Double column vector of upper bounds
%           x0        Double column vector for initial point of x
%           options   Options structure created with cplexoptimset
%
% Output Arguments
% x         Solution found by the optimization function. If exitflag > 0,
%           then x is a solution; otherwise, x is the value of the
%           optimization routine when it terminated prematurely.
% fval      Value of the objective function at the solution x
% exitflag  Integer identifying the reason the optimization algorithm
%           terminated
% output    Structure containing information about the optimization. The
%           fields of the structure are:
%           iterations         Number of iterations
%           algorithm          Optimization algorithm used
%           message            Exit message
%           time               Execution time of the algorithm
%           cplexstatus        Status code of the solution
%           cplexstatusstring  Status string of the solution
% lambda    Structure containing the Lagrange multipliers at the solution x
%           (separated by constraint type). It is only populated in the
%           case that the model has an empty set of quadratic constraints.
%           The fields of the structure are:
%           lower              Lower bounds lb
%           upper              Upper bounds ub
%           ineqlin            Linear inequalities
%           eqlin              Linear equalities
%
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexqcp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
