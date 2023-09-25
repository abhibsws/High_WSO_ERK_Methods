s =5; p = 3; q = 3; no_it = 10;
Mthd= {};
for i = 1:5
    [A,b,c,exitflag] = FindOptSchemes(s,p,q,no_it);
    if exitflag == 1
        Mthd{i}.A = A; Mthd{i}.b = b; Mthd{i}.c = c;
    else
        Mthd{i}.A = []; Mthd{i}.b = []; Mthd{i}.c = [];
    end
end

save_data = 1;
if save_data 
    foldername_scheme = sprintf('SchemeData/');
    if exist(foldername_scheme,'dir') == 0
        mkdir(foldername_scheme);
    end
end
if save_data
    fname_scheme = sprintf('%sScheme_s%dp%dq%d',foldername_scheme,s,p,q);
    save(fname_scheme,'Mthd');
end

