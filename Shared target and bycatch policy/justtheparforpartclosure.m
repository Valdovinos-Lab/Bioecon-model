function [] = justtheparforpartclosure(ty,ty2,SimConsLin,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,dataforrandomtarget)
%function [] = justtheparforpartclosure(ty,ty2,SimConsLin,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,Indexes)
  
parfor k=ty:ty2 
   % ui = Indexes{k}{2}; %total bycatch possible
   % targetindex = Indexes{k}{3};
    for i1 =1:length(cutofft)
        for j1=1:length(cutoffb)
           ui = dataforrandomtarget{k}; %total target/bycatch combinatirons possible
            Qw = cell(1,length(ui)); 
            targesc = cutofft(i1);
            byesc = cutoffb(j1);
            for i=1:length(ui)
               Qw{i} = discreteeffortssharedtargetbycatchquota(SimConsLin,k,N,costratio,cat,cab,targesc,byesc,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,ui(i,1),ui(i,2));
              % Qw{i} = discreteeffortssharedtargetbycatchquota(SimConsLin,k,N,costratio,cat,cab,targesc,byesc,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,targetindex,ui(i));
            end
            j = sprintf('sharedtargetbycatchquota_k%d_Emax%d_cutofft%d_cutoffb%d_cat%d_cab%d_costratio%d_hill%d_hsb%d_fa%d_fm%d.mat',k,inE,targesc,byesc,cat,cab,costratio,hill,hsb,fa,fm)
            parsave(j,Qw)
        end
    end
end
end