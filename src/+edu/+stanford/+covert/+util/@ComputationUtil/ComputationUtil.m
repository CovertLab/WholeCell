% Utility functions for various computations:
% - linear programming
% - quadratic programming
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 5/12/2009
classdef ComputationUtil
    methods (Static = true)
        varargout = linearProgramming(...
            maximizeFlag, f, A, b, lb, ub, constraintTypes, variableTypes, options)
        
        [x, objective, errorFlag, errorMsg] = quadraticProgramming(...
            H, f, A, b, Aeq, beq, lb, ub, x0, options)
        
        [x, nroot] = cubicfcn(a, b, c, d)
        
        [b,fval,exitflag,output] = fzero(FunFcnIn,x,options,varargin)
        
        [a, b, cutX, gof] = bilinearFit(x, y)
        
        function [A, b] = calcIndependentLinearConstraints(A, b)
            %throw error if constraints are inconsistent
            if rank(A) - rank([A b]) > 0
                throw(MException('Simulation:error',...
                    'transcription unit constraints cannot be made consistent'));
            end
            
            %if constraints are dependent (eg. A is not full row rank) find submatrix
            %which is full row rank, and has rank equal to that of A
            [~, pivotRows] = rref(A');
            A = A(pivotRows, :);
            b = b(pivotRows);
        end
        
        function matrix = invertCompositionMatrix(matrix)
            matrix = diag(1 ./ sum(matrix)) * matrix';
        end
        
        function data = imputeMissingData(data, imputedValue)
            if ~exist('imputedValue', 'var')
                imputedValue = 'mean';
            end
            
            indicesToImpute = isnan(data) | isinf(data);
            switch imputedValue
                case 'mean', data(indicesToImpute) = mean(data(~indicesToImpute));
                case 0,      data(indicesToImpute) = 0;
            end
        end
        
        function val = roundHalfUp(val)
            val = round(val);
        end
        
        function val = roundHalfDown(val)
            val = floor(val) + (rem(val,1) > 0.5);
        end
        
        %http://www.mathworks.com/matlabcentral/newsreader/view_thread/299383
        function result = isdiag(A)
            
            if nnz(A) > size(A, 1) %<----Matt J added this pre-check
                result = false;
                return;
            end
            
            [I, J] = find(A);
            
            if ~isempty(I)
                result = all(I == J);
            else
                % make the simple choice that an all zero matrix
                % or an empty array is by definition diagonal, since
                % it has no non-zero off diagonals.
                result = true;
            end
        end
    end
end
