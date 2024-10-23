function pgon = decsg2(gd,sf,ns)
%DECSG2  Decompose constructive solid 2-D geometry into minimal regions
% Inputs:
%   gd: geometry description vector
%   sf: set formula
%   ns: name-space matrix
filename = mfilename('fullpath');
filepath = fileparts( filename );

addpath([filepath,'/poly-test']);

%-------------------------------------------------------------- parse input
% validate name-space and count regions
if ischar(ns)
    numreg = size(ns,2);
    ns = regexp(reshape(char(ns,''),[],1).','\w*','match');
elseif iscellstr(ns)
    numreg = numel(ns);
else
    error('decsg2:incorrectInputClass','Incorrect input class.');
end

% decompose set formula
sf = regexprep(sf,' ',''); % remove spaces
v  = regexp(sf,'[\+\-\*]','split');
op = regexp(sf,'[\+\-\*]','match');
[lchk,locns] = ismember(v,ns);
assert(all(lchk),'decsg2:undefinedString','Undefined string in set formula.');

% extract regions
rgns = cell(1,numreg);
n1 = 1;
for i = 1:numreg
    n2 = n1 + 2*gd(n1);
    rgns{i} = reshape(gd(n1+1:n2),[gd(n1) 2]);
    n1 = n2 + 1;
end

%--------------------------------------------- iteratively form the polygon
pgon = rgns{locns(1)};
counter = 0;
for i = locns(2:end)
    counter = counter + 1;

    n1 = size(pgon,1);
    n2 = size(rgns{i},1);
    e1 = [1:n1;2:n1 1].';
    e2 = [1:n2;2:n2 1].';
    % find intersections of boundaries
    [p,A] = lineline([pgon(e1(:,1),:) pgon(e1(:,2),:)],[rgns{i}(e2(:,1),:) rgns{i}(e2(:,2),:)]);
    if isempty(p), continue; end

    in1 = inpoly2(pgon,rgns{i},e2);
    in2 = inpoly2(rgns{i},pgon,e1);
    if op{counter}=='*', in1 = ~in1; end
    if op{counter}~='+', in2 = ~in2; end

    idx = find(A) - 1; % linear indices
    l1 = [mod(idx,size(A,1))+1;0]; % row index with pad
    l2 = [floor(idx/size(A,1))+1;0]; % column index with pad

    % modify connectivity
    pop1 = [in1(e1(:,2));false(numel(idx),1)];
    pop2 = [in2(e2(:,2));false(numel(idx),1)];
    e2 = e2 + n1 + numel(idx); % shift indices
    for j = 1:numel(idx)
        e1(end+1,:) = [n1+j e1(l1(j),2)];
        e1(l1(j),2) = n1+j;
        e2(end+1,:) = [n1+j e2(l2(j),2)]; % reference the same intersection
        e2(l2(j),2) = n1+j;
        % update list of segments to delete
        pop1(n1+j) = pop1(l1(j));
        pop1(l1(j)) = false;
        pop2(n2+j) = pop2(l2(j));
        pop2(l2(j)) = false;
    end
    pop1(in1) = true;
    pop2(in2) = true;

    % delete interior simplices
    e1(pop1,:) = [];
    e2(pop2,:) = [];
    e1 = reshape(e1,[],2);
    e2 = reshape(e2,[],2);

    % total region
    pgon = [pgon;p;rgns{i}];
    edge = [e1;e2];

    % sort vertices to form a simple polygon
    % create matrix M where M(i,j) = 1 if vertices i and j are connected; 0 otherwise
    M = sparse(edge,edge(:,[2 1]),true(numel(edge),1),size(pgon,1),size(pgon,1));
    ordering = zeros(size(edge,1),1);
    ordering(1) = edge(1,1); % start from any edge vertex
    skip = edge(1,2); % prevent revisiting of previous nodes
    for j = 2:size(edge,1)
        candidates = find(M(:,ordering(j-1)));
        ordering(j) = candidates(candidates~=skip);
        skip = ordering(j-1);
    end

    pgon = pgon(ordering,:); % sorted
end