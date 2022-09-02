function allRatios = AllLogRatio(W)
[K,D] = size(W);

no = D*(D-1)/2;
allRatios = zeros(K,no);

d = 1;

for i=1:D-1
    for j=i+1:D
        allRatios(:,d) = log(W(:,i) ./ W(:,j)); 
        d = d+1;
    end
end

end