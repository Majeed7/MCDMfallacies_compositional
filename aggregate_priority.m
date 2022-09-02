function [agg_weight,average_array,deviation_array,avgDev_array,info] = aggregate_priority(W,method)

if nargin < 2
    method = "mestimation";
end

ratios = AllLogRatio(W);
[dm_no,ratio_no] = size(ratios);
tol = 10e-10;
iter = 1000;
i=0;
ratio_agg = rand(ratio_no,1);

alpha = ones(1,dm_no);

if lower(method) == "mean"
    ratio_agg = mean(ratios);
    deviation = std(ratios);
    info.method = "mean";
elseif lower(method) == "median"
    ratio_agg = median(ratios);
    deviation = mad(ratios,1);
    info.method = "median";
elseif lower(method) == "mestimation"
    while i < iter
        
        ratio_agg_old = ratio_agg;
        lambda = alpha / sum(alpha);
        ratio_agg = ratios'*lambda';
        %lambda = sum(alpha) / size(alpha,1);
        %R_star = sum(alpha .* R,2) ./ sum(alpha,2);
        
        Error = ratios - repmat(ratio_agg,1,dm_no)';
        Error_norm = sqrt(sum(Error.^2,2));
        if(sum(Error_norm) == 0)
            break;
        end
        coef = length(Error_norm) / 5;
        %sigma = min((coef*norm(Error_norm,2)^2 / (length(Error_norm)^2)),100);
        sigma = (coef*norm(Error_norm,2)^2 / (length(lambda)^2));

        alpha_old = alpha;
        alpha = delta(Error_norm,sigma);
        
        
        if norm(alpha(:)-alpha_old(:),2) < tol && norm(ratio_agg-ratio_agg_old) < tol
            break;
        end
        
        i = i + 1;
    end
    
    temp = ((ratios - ratio_agg').^2) .* lambda';
    deviation = sum(temp).^.5;
    info.lambda = lambda;
    info.sigma = sigma;
    info.method = "mestimation";
end

%% computing average array and aggregated weight
average_array = tril(ones(size(W,2)),-1);
average_array(average_array==1) = ratio_agg;

average_array = exp(average_array - average_array')';

if lower(method) == "median"
    temp = average_array ./ sum(average_array);
    agg_weight = geomean(temp,2)';
else
    agg_weight = average_array(:,1)' ./ sum(average_array(:,1));
end

%% computing deviation array
deviation_array = tril(ones(size(W,2)),-1);
deviation_array(deviation_array==1) = deviation;

deviation_array = deviation_array + deviation_array';

%%
avgDev_array = tril(deviation_array) + triu((average_array),1);

end

function E = delta(v,sigma)
    K = length(v);
    E = zeros(1,K);
    for i=1:K
        E(i) = exp(-v(i)^2 / (2*sigma^2));
    end
end
