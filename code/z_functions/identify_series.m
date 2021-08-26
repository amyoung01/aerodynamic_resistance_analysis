function bool_gt_maxlength = identify_series(x,maxlength)
x = [x;0];
val_index = 1:length(x);
dx = [x(1);diff(x)];
length_of_series = val_index(dx==-1)-val_index(dx==1);
n_gt_maxlength = sum(length_of_series>maxlength);
bool_gt_maxlength = n_gt_maxlength>0;
end