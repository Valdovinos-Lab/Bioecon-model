function [] = justtheparforpartsharedtarget(ty,ty2,SimConsLin,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,dataforrandomtarget)
   
parfor k=ty:ty2
    ui = dataforrandomtarget{k}; %total target/bycatch combinatirons possible
    Qw = cell(1,length(ui)); 
    for i1=1:4
        targesc = cutofft(i1);
            for i=1:length(ui)
               Qw{i} = discreteeffortssharedtarget(SimConsLin,k,N,costratio,cat,cab,targesc,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,ui(i,1),ui(i,2));
            end
         j = sprintf('sharedtargetquota_k%d_Emax%d_cutofft%d_cat%d_cab%d_costratio%d_hill%d_hsb%d_fa%d_fm%d.mat',k,inE,targesc,cat,cab,costratio,hill,hsb,fa,fm)
        parsave(j,Qw)
    end
    
end
end