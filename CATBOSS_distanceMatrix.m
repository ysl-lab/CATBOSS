chgcount = 0;
for i = 1:length(changes)
    if nnz(changes(i,2:size(changes,2))) > 0
        chgcount = chgcount + 1;
        chgindices(chgcount) = i;
    end
end

firstpt = chgindices(1);
lastpt = chgindices(2)-1;
lengths = zeros([chgcount-1,1]);
for i = 2:chgcount
    lengths(i-1) = 1+lastpt-firstpt;
    firstpt = chgindices(i);
    if i<chgcount
       lastpt = chgindices(i+1)-1;
    else
       lastpt = length(dihedrals);
    end
end

disp('Computing distance matrix');

dmtxcount = 1;
firstpt1=0;
firstpt2=0;
lastpt1=0;
lastpt2=0;
l1=0;
l2=0;
distmtx = zeros((chgcount-1)*(chgcount-2)/2,5);
for i=1:chgcount-2
   disp(i)
   for j=i+1:chgcount-1
       firstpt1 = chgindices(i);
       lastpt1 = chgindices(i+1)-1;
       firstpt2 = chgindices(j);
       lastpt2 = chgindices(j+1)-1;
       distmtx(dmtxcount,1)=i;
       distmtx(dmtxcount,2)=j;
       l1 = length(firstpt1:lastpt1);
       l2 = length(firstpt2:lastpt2);
       D = zeros(l1,l2);
       for k = 1:l1
           for l = 1:l2
               diffs = zeros(size(unshifted_dihed,2),1);
               for m = 1:size(unshifted_dihed,2)
                   diffs(m) = abs(unshifted_dihed(firstpt1+k-1,m)-unshifted_dihed(firstpt2+l-1,m));
                   diffs(m) = min(diffs(m),360-diffs(m));
               end               
               D(k,l) = sqrt(sum(diffs.^2));
           end
       end
       distmtx(dmtxcount,3) = emd(unshifted_dihed(firstpt1:lastpt1,1:2),unshifted_dihed(firstpt2:lastpt2,1:2)...
           ,[],[],'precomputed',D);
       distmtx(dmtxcount,4) = lengths(i);
       distmtx(dmtxcount,5) = lengths(j);
       dmtxcount=dmtxcount+1;
   end
end
fid = fopen('DISTANCE_MATRIX','w');
for i=1:size(distmtx,1)
fprintf(fid,'%d %d %15.6f %d %d\n',distmtx(i,1),distmtx(i,2),distmtx(i,3),distmtx(i,4),distmtx(i,5));
end
fclose all;

function [dist, flow] = emd(X, Y, X_weights, Y_weights, distance, D)
%Computes the Earth Mover's Distance between two weighted samples
% emd(X, Y[, X_weights, Y_weights, distance, D])
% X : First sample
% Y : Second sample
% X_weights : weights of elements in X (must sum to 1);
%    default: 1/|X| is used as the weight of each element
% Y_weights : weights of elements in Y (must sum to 1);
%    default: 1/|Y| is used as the weight of each element
% distance : valid distance type used by pdist2,
%    or "precomputed" if the pairwise distances D is supplied;
%    default: euclidean
% D : precomputed distance matrix; ignored unless distance="precomputed";
%    must be array of size |X|-by-|Y|
if nargin < 2
    error('Must at least provide two samples.');
end
if nargin < 3
    X_weights = [];
end
if nargin < 4
    Y_weights = [];
end
if nargin < 5
    distance = 'euclidean';
end
if nargin < 6
    D = [];
end

if ~strcmp(distance, 'precomputed')
    n = length(X);
    m = length(Y);
    D = pdist2(X, Y, distance);
    if isempty(X_weights)
        X_weights = ones(1, n)/n;
    elseif n ~= length(X_weights)
        error('Size mismatch of X and X_weights');
    else
        X_weights = reshape(X_weights, 1, []);
    end
    if isempty(Y_weights)
        Y_weights = ones(1, m)/m;
    elseif m ~= length(Y_weights)
        error('Size mismatch of Y and Y_weights');
    else
        Y_weights = reshape(Y_weights, 1, []);
    end
else
    if isempty(D)
        error('D must be supplied when distance=''precomputed''');
    end
    [n, m] = size(D);
    if isempty(X_weights)
        X_weights = ones(1, n)/n;
    elseif n ~= length(X_weights)
        error('Size mismatch of D and X_weights');
    else
        X_weights = reshape(X_weights, 1, []);
    end
    if isempty(Y_weights)
        Y_weights = ones(1, m)/m;
    elseif m ~= length(Y_weights)
        error('Size mismatch of D and Y_weights');
    else
        Y_weights = reshape(Y_weights, 1, []);
    end
end

if nargout > 1
    [dist, flow] = c_emd(X_weights, Y_weights, D);
else
    dist = c_emd(X_weights, Y_weights, D);
end

end

