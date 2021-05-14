%PART I: Pick out the changes

chgcount = 0;
chgindices = [];
for i = 1:length(changes)
    if nnz(changes(i,2:size(changes,2))) > 0
        chgcount = chgcount + 1;
        chgindices(chgcount) = i;
    end
end

%PART II: Find the properties of trajectory segments

firstpt = chgindices(1);
lastpt = chgindices(2)-1;
means = zeros([chgcount-1,size(dihedrals,2)]);
stdevs = zeros([chgcount-1,size(dihedrals,2)]);
lengths = zeros([chgcount-1,1]);

sloped = zeros(chgcount-1,size(dihedrals,2));

for i = 2:chgcount
    means(i-1,:) = mean(dihedrals(firstpt:lastpt,1:size(dihedrals,2)),1);
    lengths(i-1) = 1+lastpt-firstpt;
    if firstpt ~= lastpt
        stdevs(i-1,:) = std(dihedrals(firstpt:lastpt,1:size(dihedrals,2)));
    else                                                                      
        stdevs(i-1,:) = ones(1,size(dihedrals,2));
    end

    for j = 1:size(dihedrals,2)
        mdl = fitlm(1:lengths(i-1),dihedrals(firstpt:lastpt,j));
        slope = mdl.Coefficients.Estimate(2);
        slopestd = std(mdl.Residuals.Raw)/sqrt(sum((dihedrals(firstpt:lastpt,j)-means(i-1,j)).^2));
        if abs(slope)>1.96*slopestd
            sloped(i-1,j)=1;
        end
    end
    
    firstpt = chgindices(i); 
    if i<chgcount
       lastpt = chgindices(i+1)-1;
    else
       lastpt = length(dihedrals);
    end
end
