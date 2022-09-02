clc
clear 

load tourismData.mat


%% Aggregating priorities by different methods
[agg_mean,avgAry_mean,devAry_mean,avgDevAry_mean,info_mean] = aggregate_priority(Weights, 'mean');
[agg_median,avgAry_median,devAry_median,avgDevAry_median,info_median] = aggregate_priority(Weights, 'median');
[agg_hq,avgAry_hq,devAry_hq,avgDevAry_hq,info_hq] = aggregate_priority(Weights, 'mestimation');


deviant  = info_hq.lambda < 0.001;
deviants = Weights(deviant,:);
non_deviant = Weights(~deviant,:);

%% Clustering
[clus,cen] = kmeans(Weights,3,'Distance','cityblock');
sum(cen,2);
[clus_comps,cen_comps] = kmeans_compositional(Weights,3,'Distance','cityblock');
sum(cen_comps,2);
