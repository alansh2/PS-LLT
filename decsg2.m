function decsg2(gd,ns)
%DECSG2  Decompose constructive solid 2-D geometry into minimal regions
% Inputs:
%   gd: geometry description vector
%   ns: name-space matrix
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath,'/poly-test']);

% parse input
numreg = size(ns,2);
rgns = cell(1,numreg);
n1 = 1;
for i = 1:numreg
    n2 = n1 + 2*gd(n1);
    rgns{i} = reshape(gd(n1+1:n2),[gd(n1) 2]);
    n1 = n2 + 1;
end

out = rgns{1};
for i = 2:numreg
    n1 = size(out,1);
    n2 = size(rgns{i},1);
    e1 = [1:n1;2:n1 1].';
    e2 = [1:n2;2:n2 1].';
    [p,A] = lineline([out(e1(:,1),:) out(e1(:,2),:)],[rgns{i}(e2(:,1),:) rgns{i}(e2(:,2),:)]);
    if isempty(p), continue; end

    in1 = inpoly2(out,rgns{i},e2);
    in2 = inpoly2(rgns{i},out,e1);

    idx = find(A) - 1;
    l1 = [mod(idx,size(A,1))+1;0]; % row index with pad
    l2 = [floor(idx/size(A,1))+1;0]; % column index with pad

    % modify connectivity
    m1 = 0;
    m2 = 0;
    for j = 1:numel(idx)
        if in1(l1(j))
            e1(end+1,:) = [e1(l1(j)+m1,1) n1+j];
            e1(l1(j)+m1,1) = n1+j;
        else
            e1(end+1,:) = [n1+j e1(l1(j)+m1,2)];
            e1(l1(j)+m1,2) = n1+j;
        end
        if in2(l2(j))
            e2(end+1,:) = [e2(l2(j)+m2,1) n2+j];
            e2(l2(j)+m2,1) = n2+j;
        else
            e2(end+1,:) = [n2+j e2(l2(j)+m2,2)];
            e2(l2(j)+m2,2) = n2+j;
        end
        % multiplicity
        if l1(j) == l1(j+1)
            m1 = n1 + j - l1(j);
        else
            m1 = 0;
        end
        if l2(j) == l2(j+1)
            m2 = n2 + j - l2(j);
        else
            m2 = 0;
        end
    end

    % delete interior simplices
    e1(ismember(e1,find(in1))) = [];
    e2(ismember(e2,find(in2))) = [];
    e1 = reshape(e1,[],2);
    e2 = reshape(e2,[],2);

    disp(e1)
    disp('')
    disp(e2)
end