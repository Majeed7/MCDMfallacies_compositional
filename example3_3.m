 clear
 clc

data = [0.2200    0.4347    0.2952    0.0501;
        0.2106    0.4336    0.3115    0.0442;
        0.3630    0.3122    0.1071    0.2178;
        0.2432    0.3861    0.3320    0.0386;
        0.2269    0.3811    0.3394    0.0526;
        ];
        
    
[agg_weight,average_array,deviation_array,avgDev_array,info] = aggregate_priority(data, "mean");
[agg_weight_hq,average_array_hq,deviation_array_hq,avgDev_array_hq,info_hq] = aggregate_priority(data, 'mestimation');

