function y = filter_adaptation_data(x, num_subjects, iqr_num)

%iqr_num = 3;
bad_idx = cell(num_subjects,1);
for ksub=1:num_subjects
    %%%% NH data
    [y.NH.V.P0.all(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.V.P0.all(:,ksub), iqr_num);
    [y.NH.V.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.V.P3.P(:,ksub), iqr_num);
    [y.NH.V.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.V.P3.N(:,ksub), iqr_num);
    [y.NH.V.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.V.P7.P(:,ksub), iqr_num);
    [y.NH.V.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.V.P7.N(:,ksub), iqr_num);
    
    [y.NH.C.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.C.P3.P(:,ksub), iqr_num);
    [y.NH.C.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.C.P3.N(:,ksub), iqr_num);
    [y.NH.C.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.C.P7.P(:,ksub), iqr_num);
    [y.NH.C.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.NH.C.P7.N(:,ksub), iqr_num);  
    
    %%%% HI data
    [y.HI.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HI.P3.P(:,ksub), iqr_num);
    [y.HI.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HI.P3.N(:,ksub), iqr_num);
    [y.HI.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HI.P7.P(:,ksub), iqr_num);
    [y.HI.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HI.P7.N(:,ksub), iqr_num);
    
    %%%% HD data
    [y.HD.P3.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HD.P3.P(:,ksub), iqr_num);
    [y.HD.P3.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HD.P3.N(:,ksub), iqr_num);
    [y.HD.P7.P(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HD.P7.P(:,ksub), iqr_num);
    [y.HD.P7.N(:,ksub), ~] = filter_iqr_vec_1118_2021a(x.HD.P7.N(:,ksub), iqr_num);
    
end

return


