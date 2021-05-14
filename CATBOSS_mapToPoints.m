%PART I: Pick out the changes

chgcount = 0;
chgindices = [];
for i = 1:length(changes)
    if nnz(changes(i,2:size(changes,2))) > 0
        chgcount = chgcount + 1;
        chgindices(chgcount) = i;
    end
end

%PART II: Find the means and STD's of all the traj segments

firstpt = chgindices(1);
lastpt = chgindices(2)-1;
grouping = zeros(length(dihedrals),1);
for i = 2:chgcount
    for j=firstpt:lastpt
        %change this 2 into 3 to enable halo control
       grouping(j) = CLUSTERASSIGNATION(i-1,2);
    end
    
    firstpt = chgindices(i);
  
    if i<chgcount
        lastpt = chgindices(i+1)-1;
    else
        lastpt = length(dihedrals);
    end
end