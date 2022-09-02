function w = reverseLogRatio(r,n)

[K,N] = size(r);
w = [ones(K,1), exp(-r(:,1:n-1))];

w = w ./ sum(w,2);

end