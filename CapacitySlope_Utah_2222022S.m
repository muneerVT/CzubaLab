%% BMS_CapacitySlope
% This function determines the capacity of an individual link and updates
% slope accordingly.

% Jon Czuba
% February 17, 2015

%%
%cycle through each link
for i=1:LinkNum
    %only do this check capacity if parcels are in link
    if ~isempty(P_vol{t,i})

        %porosity of all sediment in the link: active and inactive
        Dg_Lp(t,i)=(2.^(sum(P_vol{t,i}./sum(P_vol{t,i}).*...
            (log10(P_d{t,i})./log10(2)))));
        Lp(t,i)=0.13+0.21./(((Dg_Lp(t,i).*1000)+0.002).^0.21);%Wu and Wang 2006

        %First In Last Out
        %compute cumulative volume in link beginning with last in
            locidx1=find(P_loc{t,i}<=0.5);
            cvol1=fliplr(cumsum(fliplr(P_vol{t,i}(locidx1))));
            exc1=find(cvol1>0.7*capacity(i,1),1,'last');
            x1=locidx1(1:exc1);
            
            locidx2=find(P_loc{t,i}>0.5 & P_loc{t,i}<=1);
            cvol2=fliplr(cumsum(fliplr(P_vol{t,i}(locidx2))));
            exc2=find(cvol2>0.3*capacity(i,1),1,'last');
            x2=locidx2(1:exc2);
           %determine which parcels that were the first to enter are above capacity         
            exc=sort(cat(2,x1,x2));

        %if parcels have been identified above capacity then
        if ~isempty(exc1)
            %set their status to inactive, 1
            P_storage{t,i}(1,exc)=1;
            %set their travel time to 0
            P_tt{t,i}(1,exc)=0;
            
            %UPDATE Elevations beginning here
            %compute volume inactive
            vstor=sum(P_vol{t,i}(exc))./(1-Lp(t,i));%m3 vol in storage, porosity Lp=0.4
            
            %determine upstream links
            usid=find(Connect(:,2)==i);
            %include current link and upstream links that are not lakes
            elevid=cat(1,i,usid(Lake(usid)==0));
            %update elevation at upstream end of link
            %volume is placed at upstream end of link and along current link
            %length and upstream link lengths to compute new elev
           mxRelev=mnelev+Slope_raw.*Length;
%            if t==1
%             mxRelev=mnelev+Slope(i,t).*Length;
%             
%            else
%             mxRelev=mnelev+Slope(i,t-1).*Length;
%            end
         Elev(t,i)=mxRelev(i,1)+(2.*vstor)./sum(B(elevid,1).*Length(elevid,1));
            %UPDATE Elevations complete
            
            clear usid elevid exc exc1 cvol1 x1 exc2 cvol2 x2 locidx2 locidx1
        end 
    end
end

%Elev(t,202)= mxelev(202)+2;
USlink = setdiff(GridID,ToNode);
for i=1:LinkNum
    if ismember(i,USlink,'legacy')
        Elev(t,i)= mxelev(i)+1;
    end
end

%compute/update slope
for i=1:LinkNum
    if i==OutletLinkID
        Slope(i,t)=(Elev(t,i)-mnelev(i,1))./Length(i,1);
    else
        Slope(i,t)=(Elev(t,i)-Elev(t,Connect(i,2)))./Length(i,1);
    end
end
Slope(Slope<1e-3)=1e-3;
%H(1:timesteps,1:LinkNum)=NaN;

%%

for i=1:LinkNum
    if isempty(P_loc{t,i})
        continue
    end
grvidx=P_d{t,i}>0.002;
actidx=P_storage{t,i}==0;
actvol=sum(P_vol{t,i}(actidx));

grvactidx=and(actidx,grvidx);
grvactvol=sum(P_vol{t,i}(grvactidx));

% Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
% (-log10(P_d{t,i}(actidx).*1000)./log10(2)))));
% if Dg(t,i)==0.5
%     Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
%         (log10(P_d{t,i}(actidx).*1000)./log10(2)))))./1000;
% end

% Dg(t,i)=((sum(P_vol{t,i}(actidx)./actvol.*...
% ((P_d{t,i}(actidx))))));

Dg(t,i)=(2.^(sum(P_vol{t,i}(actidx)./actvol.*...
    (log10(P_d{t,i}(actidx))./log10(2)))));
DG(t,i)=(2.^(sum(P_vol{t,i}(grvactidx)./grvactvol.*...
    (log10(P_d{t,i}(grvactidx))./log10(2)))));

%Lp(t+1,i)=0.13+0.21./(((Dg(t,i).*1000)+0.002).^0.21);%Wu and Wang 2006

% Dg=(2.^(sum(Fdfpsd.*...
%     (log10(Dpsd)./log10(2)))));

actsandidx=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(1));
actsandvol=sum(P_vol{t,i}(actsandidx));
Fs=actsandvol./actvol;

taursg=rho.*R.*g.*Dg(t,i).*(0.021+0.015.*exp(-20.*Fs));

% cumulative PSD all sizes
clear cpsd
cpsd(1,1:length(Dcsd),1)=NaN;
for j=1:length(Dcsd)
    actjidx=and(P_storage{t,i}==0,P_d{t,i}<=Dcsd(j));
    actftjvol=sum(P_vol{t,i}(actjidx));
    cpsd(1,j)=actftjvol./actvol;
end

for j=1:length(Dpsd)
    actjidxD=and(P_storage{t,i}==0,P_d{t,i}<=Dpsd(j));
    actftjvolD=sum(P_vol{t,i}(actjidxD));
    Dcpsd(i,j)=actftjvolD./actvol;
end

% compute D84
Dxx=0.84;
for j=1:length(Dcsd)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D84(t,i)=exp(interp1(cpsd(1,j:j+1),log(Dcsd(j:j+1,1)),Dxx));
    end
end
% compute D50
Dxx=0.50;
for j=1:length(Dcsd)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D50(t,i)=exp(interp1(cpsd(1,j:j+1),log(Dcsd(j:j+1,1)),Dxx));
    end
end
% compute D16
Dxx=0.16;
for j=1:length(Dcsd)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D16(t,i)=exp(interp1(cpsd(1,j:j+1),log(Dcsd(j:j+1,1)),Dxx));
    end
end

% cumulative PSD gravel only, no sand
clear cpsd
cpsd(1,1:length(DcsdG),1)=NaN;
for j=1:length(DcsdG)
    actjidx=and(grvactidx,P_d{t,i}<=DcsdG(j));
    actftjvol=sum(P_vol{t,i}(actjidx));
    cpsd(1,j)=actftjvol./grvactvol;
end

% compute D84
Dxx=0.84;
for j=1:length(DcsdG)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D84G(t,i)=exp(interp1(cpsd(1,j:j+1),log(DcsdG(j:j+1,1)),Dxx));
    end
end
% compute D50
Dxx=0.50;
for j=1:length(DcsdG)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D50G(t,i)=exp(interp1(cpsd(1,j:j+1),log(DcsdG(j:j+1,1)),Dxx));
    end
end
% compute D16
Dxx=0.16;
for j=1:length(DcsdG)-1
    if cpsd(1,j)<Dxx && cpsd(1,j+1)>=Dxx
        D16G(t,i)=exp(interp1(cpsd(1,j:j+1),log(DcsdG(j:j+1,1)),Dxx));
    end
end

%DXX(isnan(DXX))=Dpsd(1);
% 
% % Schneider 2016 WRR
% a1=6.5;
% a2=3.97;
% e=1.5;
% 
% num=a2.*(H(t,i)./D84).^(5/6);
% den=sqrt(a1.^2+a2.^2.*(H(t,i)./D84).^(5/3));
% 
% Sred=Slope(i,1).*(num./den).^e;
% tau=rho.*g.*H(t,i).*Sred;


n=.045;
H2(i)=((Q2(i).*n)./(B(i,1).*sqrt(Slope(i,t)))).^(3/5); %Bankfull H
tau2=rho.*g.*H2(i).*Slope(i,1);
% if Slope(i,1)>0.04
%     n=n/2;
% end
%H(t,i)=((Q(t,i).*(2.*Dg(t,i)).^(1/6))./(B(i,1).*sqrt(g.*Slope(i,1)).*8.1)).^(3/5);
H(t,i)=((Q(t,i).*n)./(B(i,1).*sqrt(Slope(i,t)))).^(3/5);


if H(t,i)>H2(i)
    H(t,i)=H2(i); % highest depth set as bankfull
end


tauSS(t,i)=rho.*g.*H(t,i).*Slope(i,1);
  
%H(H> bH(i))= bH(i);
  
% Rickemann and Recking SP
%a1=6.5;
%a2=2.5;
%d84=5/1000;

%CHECK parenthesis
clear SP ftotal fo
SP=H(t,i)./D84(t,i);
ftotal=8.*(a1.^2+a2.^2.*SP.^(5/3))./(a1.*a2.*SP).^2;
fo=8./(6.5.*SP.^(1/6)).^2;

So(i,t)=Slope(i,t).*(sqrt(fo./ftotal)).^1.5; %t/1

tauSP=rho.*g.*H(t,i).*So(i,t); %t/1
tau=rho.*g.*H(t,i).*Slope(i,t);


 tau=tauSP; % stress partitioning ON



for k=1:gsclass
    % because ==Dpsd this formulation will not work when there is abrasion!
    actidxj=and(P_storage{t,i}==0,P_d{t,i}==Dpsd(k));
    actvolj=sum(P_vol{t,i}(actidxj));
    Fj=actvolj./actvol;
    
    b=0.67./(1+exp(1.5-Dpsd(k)./Dg(t,i)));
    
    %v moved out of k loop
    %tau=rho.*g.*H(t,i).*Slope(i,1);
    taurj=taursg.*(Dpsd(k)./Dg(t,i)).^b;
    tautaurj=tau./taurj; %tau or tauSP
    
    if tautaurj<1.35
        Wj=0.002.*tautaurj.^7.5;
    else
        Wj=14.*(1-0.894./sqrt(tautaurj)).^4.5;
    end
    
%     if k==1
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta(i)./Wj./tau.^(3/2)./Fs;
%     else
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta(i)./Wj./tau.^(3/2)./(1-Fs)./Fj;
%     end
    
%     if k==1
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./Fs;
%     else
%         P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./(1-Fs)./Fj;
%     end
    
    P_tt{t,i}(actidxj)=rho.^(3/2).*R.*g.*Length(i).*theta./Wj./tau.^(3/2)./Fj;
    
    clear actidxj actvolj Fj b taurj tautaurj Wj
    
end

clear actidx actvol actsandidx actsandvol Fs taursg tau e num den Sred
end

%%
% limit=ceil(capacity./pvol);
% 
% % determine which inputs are active during current timestep
% n = cellfun(@length,P_loc(t,:));
% 
% excess=find((n'-limit)>0);
% 
% if ~isempty(excess)
%     for k = 1:numel(excess)
%         
%         % First In First Out
%         [srt, sind]=sort(P_loc{t,excess(k)},'descend');
%         %include storage set to 0; no storage set to 1
%         %P_active{t,excess(k)}(1,sind(1:limit(excess(k))))=0;
%         P_active{t,excess(k)}(1,sind(limit(excess(k))+1:end))=0;
%         
%         % First In Last Out
% %        P_active{t,excess(k)}(1,limit(1:excess(k)))=0;
%         
%         
%         clear srt sind
%     end
% end

