home;



%take the Force issued on a FF trial and divide by velocity


f_sub.err_FF = get_sub_data(fr, repmat(idx.err_FF.all',num_subjects,1));
f_sub.err_EC = get_sub_data(fr, repmat(idx.err_EC.all',num_subjects,1));




vel_sub.err_FF = get_sub_data(vel, repmat(idx.err_FF.all',num_subjects,1));


% figure; hold on; 
% plot( squeeze(vel_sub.err_FF(1,1,:))); 
% plot( squeeze(f_sub.err_FF(1,1,:)));


visc_all = nan(num_subjects,1);
for i=1:num_subjects
   B_check_tmp = squeeze(f_sub.err_FF(i,1,:)) ./ squeeze(vel_sub.err_FF(i,1,:));
   B_check_tmp2 = round(B_check_tmp,2);
   if all(B_check_tmp2==B_check_tmp2(1)), visc_all(i) = B_check_tmp(1); end
   
end

visc_all