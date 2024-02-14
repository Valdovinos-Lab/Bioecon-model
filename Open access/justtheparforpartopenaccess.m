function [] = justtheparforpartopenaccess(ty,ty2,SimConsLin,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,dataforrandomtarget)
   
parfor k=ty:ty2
    ui = dataforrandomtarget{k}; %total target/bycatch combinatirons possible
    Qw = cell(1,length(ui)); 
            for i=1:length(ui)
               Qw{i} = discreteeffortsopenaccess(SimConsLin,k,N,costratio,cat,cab,cutofft,cutoffb,nharvest,nbycatch,inE,nseasons,hill,hsb,fa,fm,ui(i,1),ui(i,2));
            end
    j = sprintf('openaccess_k%d_Emax%d_cat%d_cab%d_costratio%d_hill%d_hsb%d_fa%d_fm%d.mat',k,inE,cat,cab,costratio,hill,hsb,fa,fm)
    parsave(j,Qw)
end
end