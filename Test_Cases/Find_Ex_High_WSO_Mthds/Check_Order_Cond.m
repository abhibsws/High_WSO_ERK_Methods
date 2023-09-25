s = size(A,1); p = 5; q = 5; 
%b = b';
[OC_err,WSO_err] = WSO_and_order_cond_check(s,p+1,q,'non-stiff',A,b)
plot_stability_region(5,'non-stiff',A,b)

norm(OC_err)
