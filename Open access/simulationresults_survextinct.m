load('present.mat')
%load("newwebsUT4Fa4Fm4H12Hb2data.mat")
%load("cutoffvalues.mat")
seasonlength = 5*365+1;
kl = length(present);
load("fishfreeavgbiomass")
%parfor i4 =1:length(cutoffvalues)
  %  Data = cell(1,400); %all data regarding variables for each network
   % cutofft = cutoffvalues(i4,1);
  %  cutoffb = cutoffvalues(i4,2);
    %cutofft = 0.6;
    %cutoffb = 0.45;
    observables = zeros(length(present),16); %target bycatch secextinct tyield(tot) byield(tot) tyield(5) byield(5) profits(tot) profits(5) avgfishseason
     medianobservables = zeros(length(present),16); %target bycatch secextinct tyield(tot) byield(tot) tyield(5) byield(5) profits(tot) profits(5) avgfishseason
    modeobservables = zeros(length(present),16);
    for i=1:length(present)
        avgintarbybiomass = Fishfreeavgbiomass{i};
        j = sprintf('openaccess_k%d_Emax1_cat2.000000e-03_cab2.000000e-03_costratio2.500000e-01_hill1.200000e+00_hsb2.000000e-01_fa4.000000e-01_fm4.000000e-01.mat',present(1,i));
        a = load(j);
        a = a.var;
        targetextinctions = zeros(length(a),1);  %capturing target extinctions
        bycatchextinctions = zeros(length(a),1); %capturing bycatch extinctions
        secondaryextinctions = zeros(length(a),1); %capturing secondary extinctions
        targetsurv = zeros(length(a),3);         %capturing target biomass loss
        bycatchsurv = zeros(length(a),3);        %capturing bycatch biomass loss
        biomasschange = zeros(length(a),3);             %capturing secondary extinctions
        avtaryield = zeros(length(a),5);                  %rolling averge (5 years, 4 years, 3 years) target yield of all iterations
        avbyyield = zeros(length(a),5);                   %avergae over 5 years bycatch yield of all iterations
        Avgfishseason = zeros(length(a),5);             %active fishing days in 1st, 2nd, 3rd, 4th and 5th years.
        tprofits = zeros(length(a),5);                   %profits
        fprofits = zeros(length(a),5); 
        
       % [i4 i]%final year profits
        for l=1:length(a)  %for every bycatch index

            b=a{1,l};
            TS = b{1,1};    %timeseries
            STS = b{1,2};   %stats for each iteration
            target = STS(1,1);
            bycatch = STS(1,2);
            spe = TS(:,1:30);
            efforts = TS(:,31:35);
            tarharvest = TS(:,36:40);
            byharvest = TS(:,41:45);
            %%
            % primary and secondary extinctions
            initialspeciesbiomass = spe(1,:);
            totalnonharvestspecies = nnz(initialspeciesbiomass) - 2;  %total species beside target and bycatch
            finalspeciesbiomass = spe(seasonlength,:);
            remainingspecies = nnz(finalspeciesbiomass);
            
           % secondaryextinctionindex = find(newlyextinct - alreadyextinct);     %indexes of those non-harvest species which got extinct after fishing
            if finalspeciesbiomass(target) ==0
                targetextinctions(l,1) = 1;        %primary extinctions target
            end
            if finalspeciesbiomass(bycatch) ==0
                bycatchextinctions(l,1) = 1;       %primary extinctions bycatch
            end
            secondaryextinctions(l,1) =1 - ((remainingspecies -2 + bycatchextinctions(l,1) + targetextinctions(l,1) )/totalnonharvestspecies) ;  %secondary extinctions
    
           %             %final/intial
            targetsurv(l,1) = finalspeciesbiomass(target)/initialspeciesbiomass(target);
            bycatchsurv(l,1) = finalspeciesbiomass(bycatch)/initialspeciesbiomass(bycatch);
            %finalyr 5 avg /last year avg before fishing happened
            targetsurv(l,3) = mean(spe(seasonlength-364:seasonlength,target))/avgintarbybiomass(l,1);
            bycatchsurv(l,3) = mean(spe(seasonlength-364:seasonlength,bycatch))/avgintarbybiomass(l,2);
            biomasschange(l,1) = sum(finalspeciesbiomass)/sum(initialspeciesbiomass);  %secondary extinctions
            %final/final-1year
            biomasschange(l,2) = sum(finalspeciesbiomass)/sum(spe(seasonlength-364,:));  %secondary extinctions
            if finalspeciesbiomass(target) ~=0
                targetsurv(l,2) = finalspeciesbiomass(target)/spe(seasonlength-364,target);
               % targetsurv(l,3) = mean(spe(seasonlength-364:seasonlength,target))/mean(spe(seasonlength - 2*365+1:seasonlength - 365,target));
            end
            if finalspeciesbiomass(bycatch) ~=0
                bycatchsurv(l,2) = finalspeciesbiomass(bycatch)/spe(seasonlength-364,bycatch);
               % bycatchsurv(l,3) = mean(spe(seasonlength-364:seasonlength,bycatch))/mean(spe(seasonlength - 2*365 +1:seasonlength - 365,bycatch));
            end
            
            %avg(finalyear)/avg(penultimateyear)
            
            biomasschange(l,3) = sum(mean(spe(seasonlength-364:seasonlength,:)))/sum(mean(spe(seasonlength - 2*365 +1:seasonlength - 365,:)));  %secondary extinctions
                      
                        
            %%
            %active fishing season per year
            totalef = sum(efforts,2);
            for t=1:5
                year = totalef(365*(t-1)+1:365*t);
                Avgfishseason(l,t) = nnz(year)/365;
            end
    
            %%
            %avg target and bycatch harvest
           for i1=1:5
                    avtaryield(l,i1) = mean((tarharvest(seasonlength,:) - tarharvest(seasonlength-(i1)*365+1,:))/i1);
                    avbyyield(l,i1) = mean((byharvest(seasonlength,:) - byharvest(seasonlength - i1*365+1,:))/i1);
           end
         
           %% profits

           %profits = (total yield - cost)(sum over all fishing days)
                targetts = TS(:,target);
                ftargets = TS(seasonlength - 364:seasonlength,target);
                N=5;
                de = 1:N;
                de = 0.65 + sort(de/(4*N)); %increasing sequence of detectibility
               tprofits(l,:)= sum(0.002*(de.*targetts - 0.25*spe(1,target)).*efforts);
               fprofits(l,:)= sum(0.002*(de.*ftargets - 0.25*spe(1,target)).*efforts(seasonlength - 364:seasonlength,:)); 
        end
         %mean
        qw = mean(Avgfishseason,1);
        ty = mean(avtaryield,1);
        ty2 = mean(avbyyield,1);
       % ty3 = mean(mean(tprofits));
        ty4 = mean(mean(fprofits));
         %median
        A5 = [mean(targetsurv,1) mean(bycatchsurv,1) mean(biomasschange,1) mean(targetextinctions,1) mean(bycatchextinctions,1) mean(secondaryextinctions,1) ty(1) ty2(1) ty4 qw(5)];
        observables(i,:) = A5; 
        qw = median(Avgfishseason,1);
        ty = median(avtaryield,1);
        ty2 = median(avbyyield,1);
       % ty3 = median(median(tprofits));
        ty4 = median(median(fprofits));
        A6 = [median(targetsurv,1) median(bycatchsurv,1) median(biomasschange,1) median(targetextinctions,1) median(bycatchextinctions,1) median(secondaryextinctions,1) ty(1) ty2(1) ty4 qw(5)];
        medianobservables(i,:) = A6;

        %mode
        qw = mode(Avgfishseason,1);
        ty = mode(avtaryield,1);
        ty2 = mode(avbyyield,1);
       % ty3 = median(median(tprofits));
        ty4 = mode(mode(fprofits));
        A6 = [mode(targetsurv,1) mode(bycatchsurv,1) mode(biomasschange,1) mode(targetextinctions,1) mode(bycatchextinctions,1) mode(secondaryextinctions,1) ty(1) ty2(1) ty4 qw(5)];
        modeobservables(i,:) = A6; 

    end
    j = sprintf('newmedian_openaccess_data.mat');
   var = cell(1,2);
    var{1}= mean(observables);
   % var{2} = median(medianobservables);
   % var{3} = mode(modeobservables);
   %median of the mean, instead of median of the median
   var{2} = median(observables);
    parsave(j,var);
%end