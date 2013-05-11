function [x,resnorm,residual,exitflag,output]=cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,x0,options)
%%
% Purpose
% Solve nonnegative, quadratically constrained least squares problems.
%
% Syntax
%    x = cplexlsqnonnegqcp(C,d)
%    x = cplexlsqnonnegqcp(C,d,Aineq,bineq)
%    x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq)
%    x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r)
%    x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,x0)
%    x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,x0,options)
%    x = cplexlsqnonnegqcp(problem)
%    [x,resnorm] = cplexlsqnonnegqcp(...)
%    [x,resnorm,residual] = cplexlsqnonnegqcp(...)
%    [x,resnorm,residual,exitflag] = cplexlsqnonnegqcp(...)
%    [x,resnorm,residual,exitflag,output] = cplexlsqnonnegqcp(...)
%    [x,resnorm,residual,exitflag,output] = cplexlsqnonnegqcp(...)
%
% Description
% Finds the minimum of a problem specified by
%
%    min      norm(C*x-d)^2
%    st.      Aineq*x      <= bineq
%             Aeq*x         = beq
%             l*x + x'*Q*x <= r
%             x >= 0
%
% d, bineq, beq, l are column vectors.
% C, Aineq, Aeq and Q are matrices.
%
% x = cplexlsqnonnegqcp(C,d) solves the least squares problem
% min norm(C*x-d)^2 such that x >= 0.
%
% x = cplexlsqnonnegqcp(C,d,Aineq,bineq) solves the preceding problem with
% the additional inequality constraints Aineq*x <= bineq.
%
% x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq) solves the preceding
% problem with the additional equality constraints Aeq*x = beq. If no
% inequalities exist, set Aineq=[] and bineq=[].
%
% x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r) solves the
% preceding problem while additionally satisfying the quadratic inequality
% constraints l*x + x'*Q*x <= r. If no equalities exist, set Aeq=[] and 
% and beq=[].
%
% x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,x0) sets the
% starting point for the algorithm to x0. If no quadratic inequalities
% exist, set l=[], Q=[] and r=[].
%
% x = cplexlsqnonnegqcp(C,d,Aineq,bineq,Aeq,beq,l,Q,r,x0,options)
% minimizes with the default optimization options replaced by values in the
% structure options, which can be created using the function cplexoptimset.
% If you do not wish to give an initial point, set x0=[].
%
% x = cplexlsqnonnegqcp(problem) where problem is a structure.
%
% [x,resnorm] = cplexlsqnonnegqcp(...) returns the value of the objective
% function at the solution x: resnorm = norm(C*x-d)^2.
%
% [x,resnorm,residual] = cplexlsqnonnegqcp(...) returns the residual at
% the solution: C*x-d.
%
% [x,resnorm,residual,exitflag] = cplexlsqnonnegqcp(...) returns a value
% exitflag that describes the exit condition of cplexlsqnonnegqcp.
%
% [x,resnorm,residual,exitflag,output] = cplexlsqnonnegqcp(...) returns a
% structure output that contains information about the optimization.
%
% Input Arguments
% C 	    Double matrix for objective function
% d 	    Double column vector for objective function
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
% ctype     String with possible char vales 'B','I','C','S','N'
%           ctype(j) to 'B', 'I','C', 'S', or 'N' to indicate
%           that x(j) should be binary, general integer,
%           continuous, semi-continuous or semi-integer
%           (respectively).
% x0        Double column vector for initial point of x
% options   Options structure created with cplexoptimset
%
% problem   Structure containing the following fields:
%           C 	    Double matrix for objective function
%           d 	    Double column vector for objective function
%           Aineq   Double matrix for linear inequality constraints
%           bineq   Double column vector for linear inequality
%                   constraints
%           Aeq     Double matrix for linear equality constraints
%           beq     Double column vector for linear equality constraints
%      	    qc 	    Struct vector
%  	        qc(i).a 	 Double column vector for linear part of the 
%                        quadratic constraint
%  	        qc(i).rhs    Double for righthand side for quadratic constraint
%  	        qc(i).Q 	 Double matrix for	the quadratic part of the
%                        quadratic constraint
%           x0      Double column vector for initial point of x
%           options Options structure created with cplexoptimset
%
% Output Arguments
% x         Solution found by the optimization function. If exitflag > 0,
%           then x is a solution; otherwise, x is the value of the
%           optimization routine when it terminated prematurely.
% resnorm 	Value of the objective function fun at the solution x
% residual 	Residual at the solution
% exitflag 	Integer identifying the reason the optimization algorithm
%           terminated
% output    Structure containing information about the optimization. The
%           fields of the structure are:
%           iterations         Number of iterations
%           algorithm          Optimization algorithm used
%           message            Exit message
%           time               Execution time of the algorithm
%           cplexstatus        Status code of the solution
%           cplexstatusstring  Status string of the solution
%
%
%  See also cplexoptimset
%

% ---------------------------------------------------------------------------
% File: cplexlsqnonnegqcp.m
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55
% Copyright IBM Corporation 2008, 2010. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with
% IBM Corp.
% ---------------------------------------------------------------------------
