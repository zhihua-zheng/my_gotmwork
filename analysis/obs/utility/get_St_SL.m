function USt_sl = get_St_SL(uSt,vSt,zSt,sld)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% find index for surface layer depth in zSt

nz  = size(uSt,1);
ntm = size(uSt,2);

[SL,Z]   = meshgrid(sld,zSt);
where_sl = Z + SL;
sl_inx   = nan(1,ntm)*nz;

for j = 1:ntm
    
    if ~isnan(where_sl(end,j))
        
        % first positive value
        sl_inx(j) = find(where_sl(:,j) >= 0,1,'first'); 
    end
end

%% average Stokes drift within surface layer

uSt_sl = nan(1,ntm);
vSt_sl = nan(1,ntm);

for j = 1:ntm
    
    slj = sl_inx(j);
    
    if ~isnan(slj)
        
        uSt_sl(j) = trapz(zSt(slj:end),uSt(slj:end,j),1);
        vSt_sl(j) = trapz(zSt(slj:end),vSt(slj:end,j),1);
    end
end

USt_sl = sqrt(uSt_sl.^2 + vSt_sl.^2);

end

