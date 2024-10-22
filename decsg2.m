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
    pop1 = [in1(e1(:,2));false(numel(idx),1)];
    pop2 = [in2(e2(:,2));false(numel(idx),1)];
    for j = 1:numel(idx)
        e1(end+1,:) = [n1+j e1(l1(j),2)];
        e1(l1(j),2) = n1+j;
        e2(end+1,:) = [n2+j e2(l2(j),2)];
        e2(l2(j),2) = n2+j;
        % update list of segments to delete
        pop1(n1+j) = pop1(l1(j));
        pop1(l1(j)) = false;
        pop2(n2+j) = pop2(l2(j));
        pop2(l2(j)) = false;
    end
    pop1(in1) = true;
    pop2(in2) = true;

    out = [out;p];
    rgns{i} = [rgns{i};p];

    % delete interior simplices
    e1(pop1,:) = [];
    e2(pop2,:) = [];
    e1 = reshape(e1,[],2);
    e2 = reshape(e2,[],2);




    figure(1); clf;
    plot(out(e1).',out(e1+size(out,1)).')
    hold on
    plot(rgns{i}(e2).',rgns{i}(e2+size(rgns{i},1)).')

    %EXAMPLE
    % t = 2*pi*(0:99).'/100;
    % gd = [5;0;0;0;2;2;-2;0;2;2;-2;100;cos(t);sin(t)];
    % ns = char('rect1','C1');
    % ns = ns.';
    % decsg2(gd,ns)
end