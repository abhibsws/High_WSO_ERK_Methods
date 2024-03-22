function [A,b,c] = RK_Butcher_Tableau(s,p,q,scheme_no)
%=========================================================================%
% This code generates different Runge-Kutta (ERK) methods based on the
% parameters (s, p, q) and scheme_no. 
% ---Inputs---
    % s = number of stages of RK scheme.
    % p = order of RK scheme.
    % q = weak stage order of RK scheme.
    % scheme_no = unique serial number.
% ---Outputs---
    % Butcher tanleau (A,b,c)
%=========================================================================%
    if s == 3 && p == 3 && q == 1 && scheme_no == 1
        % SSPRK(3,3,1): principal error norm = 0.11
        A = [0,  0, 0;
             1,  0, 0;
           1/4,1/4, 0;];
        b = [1/6, 1/6, 2/3];
        c = sum(A,2);
    elseif s == 4 && p == 3 && q == 2 && scheme_no == 2
        %ExHighWSO-(4,3,2): principal error norm = 0.05
        A = [0,0,0,0;
            (3/10),0,0,0;
            (2/3),0,0,0;
            (-21/320),(45/44),(-729/3520),0];
        b = [(7/108),(500/891),(-27/44),(80/81)];
        c = [0,(3/10),(2/3),(3/4)]';
    elseif s == 5 && p == 3 && q == 3 && scheme_no == 3
        %ExHighWSO-(5,3,3): principal error norm = 0.07
        A = [0,0,0,0,0;
            (3/11),0,0,0,0;
            (285645/493487),(103950/493487),0,0,0;
            (3075805/5314896),(1353275/5314896),0,0,0;
            (196687/177710),(-129383023/426077496),(48013/42120),(-2268/2405),0];
        b = [(5626/4725),(-25289/13608),(569297/340200),(324/175),(-13/7)];
        c = [0,(3/11),(15/19),(5/6),1]';
    elseif s == 4 && p == 4 && q == 1 && scheme_no == 4
        % RK(4,4,1): principal error norm = 0.02
        A = [0  0  0  0;
               0.5 0  0  0;
               0  0.5 0  0;
               0  0  1  0];
        b = [1/6, 2/6, 2/6, 1/6];
        c = sum(A,2);
    elseif s == 6 && p == 4 && q == 3 && scheme_no == 5
        % ExHighWSO-(6,4,3): principal error norm = 0.01
        A = [0,0,0,0,0,0;
             1,0,0,0,0,0;
             (461/3920),(99/3920),0,0,0,0;
             (314/605),(126/605),0,0,0,0;
             (13193/197316),(39332/443961),(86632/190269),(-294151/5327532),0,0;
             (884721/773750),(52291/696375),(-155381744/135793125),(-53297233/355151250),(74881422/85499375),0];
        b = [(113/2880),(7/1296),(91238/363285),(-1478741/1321920),(147987/194480),(77375/72864)];
        c = [0,1,(1/7),(8/11),(5/9),(4/5)]';
    elseif s == 7 && p == 4 && q == 4 && scheme_no == 6
        %ExHighWSO-(7,4,4): principal error norm = 0.016
        A = [0,0,0,0,0,0,0;
            (13/15),0,0,0,0,0,0;
            (12603078031712033723970154667732315345979717833/24160502835995108267237194998715463575206403200),(1048907966089364624562691286403757878851144981/72481508507985324801711584996146390725619209600),0,0,0,0,0;
            (599677/612720),(1/185),(1/69),0,0,0,0;
            (424559865415888618629601372734144061495213187/3221400378132681102298292666495395143360853760),(2833374238559988231687781191318699943620201/644280075626536220459658533299079028672170752),(10939005/8358742409),0,0,0,0;
            (-57910884850960710216615685584299594389734738212701891/143200771200355662082106164157818074407670810332226816),(1058112267371143700865117923757427517682611842306739/12613539432155680079771009278149778574768931480040704),(-16209173194776291095/101264438213449707492),(-32817015/650336938),(19/34),0,0;
            (265884392436244759286808856348403413175819281494905197368751/270149059788838487928804328430790305367090530076692853841920),(-79382446353850270662474337268761271598625606923168650489267/55522863756255233270341476827744466906189590810235664261120),(60065272366553932534499312161/31525704659895735897804153776),(181308146225/261301556177),(-88/41),(51/64),0];
        b = [(-124664924851382288077728137/37822531451473304797511250),(60663184710162550730781989/9804158397668360944968750),(-276006970775403888708156064/46485365179396605438524625),(-265040017375280744373155776/74664519359574416777071875),(100421079686366171194956352/16321974107859175324078125),(57748716265592512/317670660091351875),(3694532608/2928669975)];
        c = [0,(13/15),(193/360),(719/720),(11/80),(1/36),(193/240)]';
    elseif s == 7 && p == 5 && q == 1 && scheme_no == 7
        % DP(7,5,1): principal error norm = 0.0005
        A = [0 0 0 0 0 0 0;
             0.2 0 0 0 0 0 0;
             0.075 0.225 0 0 0 0 0;
             0.9777777777777777  -3.7333333333333334   3.5555555555555554 0 0 0 0;
             2.9525986892242035 -11.595793324188385    9.822892851699436 -0.2908093278463649  0 0 0;
             2.8462752525252526 -10.757575757575758    8.906422717743473  0.2784090909090909  -0.2735313036020583 0 0;
             0.0911458333333333   0                    0.4492362982929021 0.6510416666666666  -0.322376179245283    0.130952380952381 0];
        b = [ 0.0911458333333333 0 0.4492362982929021 0.6510416666666666 -0.322376179245283  0.130952380952381 0];           
        c = sum(A,2);
    elseif s == 8 && p == 5 && q == 4 && scheme_no == 8
        %ExHighWSO-(8,5,4): principal error norm = 0.0077
        A = [0,0,0,0,0,0,0,0;
            0.64516129032258065E-1,0,0,0,0,0,0,0;
            0.11473736235193498E0,0.90390842776270145E-1,0,0,0,0,0,0;
            0.31436403860552034E0,(-0.86293863166923845E-1),0.16666666666666667E0,0,0,0,0,0;
            0.33125823507228792E1,(-0.40614484079007382E1),0.13541292150725958E1,0,0,0,0,0;
            0.34214855910511757E1,(-0.51131692112509222E1),0.20115703957091698E1,0.54165168602903831E0,(-0.66666666666666667E-1),0,0,0;
            0.83313507353076802E0,0.32560191616282976E0,(-0.19567245683807231E1),0.17191375642293256E1,0.11189486103529771E0,(-0.97560975609756098E-1),0,0;
            (-0.41537931248421012E1),0.95760165664257128E1,(-0.76571990866955658E1),0.30827696689952624E1,0.24877849404541579E0,0.10570339214133248E-1,(-0.10714285714285714E0),0];
        b = [(-0.52389485469257312E1),0.10367455564352238E2,(-0.71359557543211860E1),0.25587099588950165E1,0.39381932475820760E-1,0.19064003972189708E1,(-0.37628391975039718E1),0.22657956458088433E1];
        c = [0,0.64516129032258065E-1,0.20512820512820513E0,0.39473684210526316E0,0.60526315789473684E0,0.79487179487179487E0,0.93548387096774194E0,0.10000000000000000E1]';
    elseif s == 9 && p == 5 && q == 5 && scheme_no == 9
        %ExHighWSO-(9,5,5): principal error norm = 0.39
        A = [0,0,0,0,0,0,0,0,0;
            (1/19),0,0,0,0,0,0,0,0;
            (1/6),0,0,0,0,0,0,0,0;
            (5/16),0,0,0,0,0,0,0,0;
            (1/2),0,0,0,0,0,0,0,0;
            (11/16),0,0,0,0,0,0,0,0;
            (11448031/2850816),(-67411795275/16590798848),(51073011/43237376),(-23353/64148),(583825/8077312),(-1/116),0,0,0;
            (30521441823091/1986340257792),(-745932230071621375/35792226257928192),(42324456085/5966757888),(775674925/6453417096),(-38065236125/28020473856),(18388001255/24775053336),(-25/138),0,0;
            (544015925591990906117739018863/21097279127167116142731264000),(-51819957177912933732533469147783191/1292529408768612025127952939417600),(15141148893501140337719772533/769541606770966638202880000),(-22062343808701233885761491/5740046662014404900523000),(-180818957612953115541011736739/146721986657116762265358336000),(18393837528018836258241002593/22366927394951953576613895000),(-14372715851/701966192290),(-3316780581/34682124125),0];
        b = [(201919428075343316424206867/7205146638186855485778750),(-979811820279525173317561445351/23232888464237446713644747250),(-659616477161155066954978/262813990730721440278125),(10343523856053877739219144704/232857239079584284108576875),(-2224588357354685208355760476/50108519801935858605643125),(704220346724742597999572733952/31288349276326419946994221875),(-13778944/1751475),(92889088/11941875),(-714103988224/149255126145)];
        c = [0,(1/19),(1/6),(5/16),(1/2),(11/16),(5/6),(16/17),1]';
  %-------------------------------------------------------------------------%
    % elseif s == 10 && p == 5 && q == 5 && scheme_no == 9
    %     %ExHighWSO-(10,5,5): principal error norm = 0.39
    %     A = [                   0                   0                   0                   0                   0                   0                   0                   0                   0                   0;
    %             0.443729657631064                   0                   0                   0                   0                   0                   0                   0                   0                   0;
    %             0.557804963962067   0.143457925697194                   0                   0                   0                   0                   0                   0                   0                   0;
    %             0.641455431308224   0.934355258808585   0.802641520123623                   0                   0                   0                   0                   0                   0                   0;
    %             0.394029542407188   0.670426033704943  -0.210427444644478   0.004791755300044                   0                   0                   0                   0                   0                   0;
    %             0.776319730236621   0.628801841859860   0.980450309985540  -0.003894488709893  -0.000000000000000                   0                   0                   0                   0                   0;
    %             0.611092336957659  -0.065394116535949   1.164783118587564   1.068973527035473  -0.851884458348837  -1.048153473381065                   0                   0                   0                   0;
    %             1.022681461826575   0.020050660557230   0.705779439356984  -0.541885559124688   0.874702714370551   0.506461312765028  -0.185523190734235                   0                   0                   0;
    %             3.014659289782928  -0.211007163335989  -3.122649343194409  -4.842940867320532  -0.870988572773360   4.057116593028406   3.542497079835235   0.724131084055787                   0                   0;
    %             2.663589792919939  -1.647369842838934   0.183978958103717   0.321415453161979   1.292560440357945  -1.067440359347609  -0.061405984983368   0.439433955154278   0.275155588023367                   0];
    %     b = [   0.288638179781999  -0.547523242710040   2.093908348291901  -1.554217052826653   1.938968314440648   2.154807575503483  -2.821393479362506  -0.425027321382880   0.149322449885373   -0.277483771621322];
    %     c = sum(A,2);
    elseif s == 1 && p == 1 && q == 1 && scheme_no == 10
        % Explicit Euler 
        A = [0];
        b = [1];
        c = sum(A,2);
    end
end
%=========================================================================%






