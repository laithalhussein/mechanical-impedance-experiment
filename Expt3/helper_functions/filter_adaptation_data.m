function y = filter_adaptation_data(x, num_subjects, iqr_num)

%iqr_num = 3;
bad_idx = cell(num_subjects,1);
for ksub=1:num_subjects
    %%%% LIE data
    [y.LIE.V.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.LIE.V.P3.P(:,ksub), iqr_num);
    [y.LIE.V.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.LIE.V.P3.N(:,ksub), iqr_num);
    [y.LIE.V.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.LIE.V.P7.P(:,ksub), iqr_num);
    [y.LIE.V.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.LIE.V.P7.N(:,ksub), iqr_num);
    [y.LIE.V.P0.all(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.LIE.V.P0.all(:,ksub), iqr_num);
    
    %%%% MIX data
    [y.MIX.V.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.V.P3.P(:,ksub), iqr_num);
    [y.MIX.V.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.V.P3.N(:,ksub), iqr_num);
    [y.MIX.V.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.V.P7.P(:,ksub), iqr_num);
    [y.MIX.V.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.V.P7.N(:,ksub), iqr_num);
    [y.MIX.V.P0.all(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.V.P0.all(:,ksub), iqr_num);
    
    [y.MIX.C.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.C.P3.P(:,ksub), iqr_num);
    [y.MIX.C.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.C.P3.N(:,ksub), iqr_num);
    [y.MIX.C.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.C.P7.P(:,ksub), iqr_num);
    [y.MIX.C.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.MIX.C.P7.N(:,ksub), iqr_num);    
    
    %%%% HIE data
    [y.HIE.C.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HIE.C.P3.P(:,ksub), iqr_num);
    [y.HIE.C.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HIE.C.P3.N(:,ksub), iqr_num);
    [y.HIE.C.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HIE.C.P7.P(:,ksub), iqr_num);
    [y.HIE.C.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HIE.C.P7.N(:,ksub), iqr_num);
    [y.HIE.V.P0.all(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HIE.V.P0.all(:,ksub), iqr_num);
    
end

return


