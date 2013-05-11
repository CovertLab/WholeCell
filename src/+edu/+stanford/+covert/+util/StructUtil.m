% Struct Utility functions
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/14/2010
classdef StructUtil
    methods (Static=true)
        function A = catstruct(varargin)
            % CATSTRUCT - concatenate structures
            %
            %   X = CATSTRUCT(S1,S2,S3,...) concates the structures S1, S2, ... into one
            %   structure X.
            %
            %   Example:
            %     A.name = 'Me' ;
            %     B.income = 99999 ;
            %     X = catstruct(A,B)
            %     % -> X.name = 'Me' ;
            %     %    X.income = 99999 ;
            %
            %   CATSTRUCT(S1,S2,'sorted') will sort the fieldnames alphabetically.
            %
            %   If a fieldname occurs more than once in the argument list, only the last
            %   occurence is used, and the fields are alphabetically sorted.
            %
            %   To sort the fieldnames of a structure A use:
            %     A = CATSTRUCT(A,'sorted') ;
            %
            %   To concatenate two similar array of structs use simple concatenation:
            %     A = dir('*.mat') ; B = dir('*.m') ; C = [A ; B] ;
            %
            %   When there is nothing to concatenate, the result will be an empty
            %   struct (0x0 struct array with no fields).
            %
            %   See also CAT, STRUCT, FIELDNAMES, STRUCT2CELL

            % for Matlab R13 and up
            % version 2.2 (oct 2008)
            % (c) Jos van der Geest
            % email: jos@jasen.nl

            % History
            % Created:  2005
            % Revisions
            %   2.0 (sep 2007) removed bug when dealing with fields containing cell
            %                  arrays (Thanks to Rene Willemink)
            %   2.1 (sep 2008) added warning and error identifiers
            %   2.2 (oct 2008) fixed error when dealing with empty structs (Thanks to
            %                  Lars Barring)

            N = nargin ;

            error(nargchk(1,Inf,N)) ;

            if ~isstruct(varargin{end}),
                if isequal(varargin{end},'sorted'),
                    sorted = 1 ;
                    N = N-1 ;
                    if N < 1,
                        A = struct([]) ;
                        return
                    end
                else
                    error('catstruct:InvalidArgument','Last argument should be a structure, or the string "sorted".') ;
                end
            else
                sorted = 0 ;
            end

            FN = cell(N,1) ;
            VAL = cell(N,1) ;

            for ii=1:N,
                X = varargin{ii} ;
                if ~isstruct(X),
                    error('catstruct:InvalidArgument',['Argument #' num2str(ii) ' is not a structure.']) ;
                end
                if ~isempty(X),
                    % empty structs are ignored
                    FN{ii} = fieldnames(X) ;
                    VAL{ii} = struct2cell(X) ;
                end
            end

            FN = cat(1,FN{:}) ;
            VAL = cat(1,VAL{:}) ;
            [UFN,ind] = unique(FN) ;

            %BEGIN added by jkarr
            UVAL = {VAL{ind}}';
            for i=1:length(UFN)
                idxs = find(strcmp(FN,UFN{i}));
                if(length(idxs)==1)
                    continue;
                end

                recurse=true;
                for j=1:length(idxs)
                    if(~isstruct(VAL{idxs(j)}))
                        recurse=false;
                        break;
                    end
                end

                if(recurse)
                    UVAL{i} = edu.stanford.covert.util.StructUtil.catstruct(VAL{idxs});
                end
            end
            VAL=UVAL;
            FN=UFN;
            %END added by jkarr

            %BEGIN removed by jkarr
            %if numel(UFN) ~= numel(FN),
            %    warning('catstruct:DuplicatesFound','Duplicate fieldnames found. Last value is used and fields are sorted') ;
            %    sorted = 1 ;
            %end
            %
            %if sorted,
            %    VAL = VAL(ind) ;
            %    FN = FN(ind) ;
            %end
            %END removed by jkarr

            if ~isempty(FN),
                % This deals correctly with cell arrays
                A = cell2struct(VAL, FN);
            else
                A = struct([]) ;
            end
        end
        
        function diff(struct1, struct2, indent)
            if ~exist('indent', 'var')
                indent = 0;
            end
            
            fieldNames1 = fieldnames(struct1);
            fieldNames2 = fieldnames(struct2);
            
            fieldsNotIn2 = setdiff(fieldNames1, fieldNames2);
            fieldsNotIn1 = setdiff(fieldNames2, fieldNames1);
            
            if ~isempty(fieldsNotIn2)
                fprintf([repmat('\t', 1, indent) 'Fields in struct1 not in struct2:\n' repmat('- %s\n', 1, numel(fieldsNotIn2))], fieldsNotIn2{:});
            end
            if ~isempty(fieldsNotIn1)
                fprintf([repmat('\t', 1, indent) 'Fields in struct2 not in struct1:\n' repmat('- %s\n', 1, numel(fieldsNotIn1))], fieldsNotIn2{:});
            end
            
            commonFieldNames = intersect(fieldNames1, fieldNames2);
            for i = 1:numel(commonFieldNames)
                val1 = struct1.(commonFieldNames{i});
                val2 = struct2.(commonFieldNames{i});
                class1 = class(val1);
                class2 = class(val2);
                if ~isequal(class1, class2)
                    fprintf([repmat('\t', 1, indent) 'Types of %s''s fields differ (%s v. %s)\n'], commonFieldNames{i}, class1, class2);
                    continue;
                end
                
                size1 = size(val1);
                size2 = size(val2);
                if ~isequal(size1, size2)
                    fprintf([repmat('\t', 1, indent) 'Sizes of %s''s fields differ (%s v. %s)\n'], commonFieldNames{i}, ...
                        sprintf(['%d' repmat(' x %d', 1, numel(size1)-1)], size1), ...
                        sprintf(['%d' repmat(' x %d', 1, numel(size2)-1)], size2));
                    continue;
                end
                
                switch class1                    
                    case 'cell'
                        idxs = find(cellfun(@(x, y) ~isequal(x, y), val1, val2));
                        if isempty(idxs); continue; end;
                        
                        fprintf([repmat('\t', 1, indent) 'Comparing %s fields\n'], commonFieldNames{i});
                        subs = cell(1, ndims(val1));
                        [subs{:}] = ind2sub(size1, idxs);
                        subs = cell2mat(subs);
                        val1 = val1(idxs);
                        for j = 1:numel(idxs)
                            fprintf([repmat('\t', 1, indent+1) '%s\t%s\n'], ...
                                sprintf(['(%d' repmat(', %d', 1, ndims(val1)-1)  ')'], subs(j, :)), ...
                                'differs');
                        end
                    case {'logical', 'int8', 'int16' ,'int32 ','int64', 'uint8', 'uint16', 'uint32', 'uint64', 'double', 'single', 'char'}
                        idxs = find(val1 ~= val2);
                        if isempty(idxs); continue; end;
                        
                        fprintf([repmat('\t', 1, indent) 'Comparing %s fields\n'], commonFieldNames{i});
                        subs = cell(1, ndims(val1));
                        [subs{:}] = ind2sub(size1, idxs);
                        subs = cell2mat(subs);
                        val1 = val1(idxs);
                        val2 = val2(idxs);
                        for j = 1:numel(idxs)
                            fprintf([repmat('\t', 1, indent+1) '%s\t%e\t%e\n'], ...
                                sprintf(['(%d' repmat(', %d', 1, ndims(val1)-1)  ')'], subs(j, :)), ...
                                val1(j), val2(j));
                        end
                    case {'edu.stanford.covert.util.SparseMat', 'edu.stanford.covert.util.CircularSparseMat'}
                        subs = find(val1 ~= val2);
                        if isempty(subs); continue; end;
                        
                        fprintf([repmat('\t', 1, indent) 'Comparing %s fields\n'], commonFieldNames{i});
                        val1 = val1(subs);
                        val2 = val2(subs);                        
                        for j = 1:numel(idxs)
                            fprintf([repmat('\t', 1, indent+1) '%s\t%e\t%e\n'], ...
                                sprintf(['(%d' repmat(', %d', 1, ndims(val1)-1)  ')'], subs(j, :)), ...
                                val1(j), val2(j)); 
                        end
                    otherwise
                        fprintf([repmat('\t', 1, indent) 'Comparing %s fields\n'], commonFieldNames{i});
                        edu.stanford.covert.util.StructUtil.diff(val1, val2, indent+1);
                end
            end
        end
    end
end