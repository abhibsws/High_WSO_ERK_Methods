s = 9; p = 5; q = 4; 
fname_scheme = sprintf('Scheme_s%dp%dq%d_1',s,p,q);
load(fname_scheme)

Schemes = {};

for i = 1:length(Mthd)
    if isempty(Mthd{i}.A) == 0
        Schemes{end+1}.A = Mthd{i}.A;
        Schemes{end}.b = Mthd{i}.b;
        Schemes{end}.c = Mthd{i}.c;
    end
end

Err_norm = zeros(1,length(Schemes));
for i = 1:length(Schemes)
    A = Schemes{i}.A; b = Schemes{i}.b;
    [OC_err,WSO_err] = WSO_and_order_cond_check(s,p+1,q,'non-stiff',A,b);
    Err_norm(1,i) = norm(OC_err);
end



i = 2;
[OC_err,WSO_err] = WSO_and_order_cond_check(s,p,q,'non-stiff',Schemes{i}.A,Schemes{i}.b)
plot_stability_region(5,'non-stiff',Schemes{i}.A,Schemes{i}.b)
norm(OC_err)