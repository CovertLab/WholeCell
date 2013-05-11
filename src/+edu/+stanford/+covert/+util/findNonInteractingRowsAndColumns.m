function [rowAssignments,colAssignments,blocks] = findNonInteractingRowsAndColumns(A)
% Factors a matrix into non-interacting blocks of rows and columns. Each
% block consists of the smallest set of rows and the smallest set of
% columns such that entries of the matrix outside the columns in the rows
% are zero and vice-versa.
%
% Author: Jonathan Karr
% Affiliation: Covert Lab, Depatment of Bioengineering, Stanford University
% Last updated: 11/5/2008

%Test Cases
% A = [1 0 1;
%      0 0 0;
%      1 0 0];
%
% A = [0 0 1;
%      0 1 0;
%      1 0 0];
%
% A = [1 0 1;
%      0 1 0;
%      1 0 0];

%assign rows, columns to blocks
nBlocks=0;
rowAssignments=zeros(size(A,1),1);
colAssignments=zeros(size(A,2),1);

%assign empty rows and columns to single block
for i=1:size(A,1)
    if(isempty(find(A(i,:),1)))
        nBlocks=1;
        rowAssignments(i)=nBlocks;
    end
end
for i=1:size(A,2)
    if(isempty(find(A(:,i),1)))
        nBlocks=1;
        colAssignments(i)=nBlocks;
    end
end

%assign remaining rows and columns
for i=1:size(A,1)
    if(rowAssignments(i)~=0); continue; end;
    nBlocks=nBlocks+1;
    rowAssignments(i)=nBlocks;
    [rowAssignments,colAssignments]=assignRecursively_row(A,i,nBlocks,rowAssignments,colAssignments);
end
for i=1:size(A,2)
    if(colAssignments(i)~=0); continue; end;
    nBlocks=nBlocks+1;
    colAssignments(i)=nBlocks;
    [rowAssignments,colAssignments]=assignRecursively_col(A,i,nBlocks,rowAssignments,colAssignments);
end

%assemble blocks
blocks=cell(nBlocks,1);
for i=1:nBlocks
    blocks{i}=A(rowAssignments==i,colAssignments==i);
end

% helper function
function [rowAssignments,colAssignments]=assignRecursively_row(A,row,nBlocks,rowAssignments,colAssignments)
for i=1:size(A,2)
    if(colAssignments(i)==0 && A(row,i)~=0)
        colAssignments(i)=nBlocks;
        [rowAssignments,colAssignments]=assignRecursively_col(A,i,nBlocks,rowAssignments,colAssignments);
    end
end

%helper function
function [rowAssignments,colAssignments]=assignRecursively_col(A,col,nBlocks,rowAssignments,colAssignments)
for i=1:size(A,1)
    if(rowAssignments(i)==0 && A(i,col)~=0)
        rowAssignments(i)=nBlocks;
        [rowAssignments,colAssignments]=assignRecursively_row(A,i,nBlocks,rowAssignments,colAssignments);
    end
end