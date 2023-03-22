s = 5; p = 3; q = 3; 
MthdName = sprintf('ExMthd_s%dp%dq%d.mat',s,p,q);
load(MthdName);

[OC_err,WSO_err] = WSO_and_order_cond_check(s,p,q,'non-stiff',A,b)
plot_stability_region(5,'non-stiff',A,b)