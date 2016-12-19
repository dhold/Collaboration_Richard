function [coup_com_save,coup_acom_save,const_factor,nn]=...
        multi_site_drude_and_brownian_HOM_4(cc_com,cc_acom,vv,Kap1,Kap2,flg)
%This version only generating the matricies which can be used to generate
%the HEOM
 % cc_com is coefficient of the commutator terms
 % cc_acom is the coefficient of the anti commutator terms
 % N is the number of sites / Chromophores
 N = length(cc_com); 

 % Kap1 e.g. = 4*omega_0 (frequency of significant transition), this
 % truncation parameter is used instead of / as well as Kap2, this
 % truncates the heirarchy based on cc dot vv via the approximation
 % exp(-cc.vv*t) ~ delta(t)/(cc.vv)
% Kap2 e.g.= 4; %second truncation parameter, there are Kappa + 1 explicitly 
 %treated exponentials, this number is the maximum order the sum of the
 %number of powers can take, e.g. n_1 + n_2 +... + n_Kappa <= Kap2

totpoles = sum(cellfun(@length,cc_com));
%tic

%% Construct labels for auxilliary density matrix

if isempty(Kap1)
    Kap1 = inf; %no truncation based on Kap1, I frequency cut off truncation
end
    
tot_poles = sum(cellfun(@length,cc_com));
if Kap2>0 && tot_poles>0 %else trivial and no hierarchy

numwithn = zeros(1,Kap2+1); %size(cc2,2) is number of poles in J(omega)
    for kk = 0:Kap2
        numwithn(kk+1) = nchoosek(tot_poles-1+kk,kk);
%number of density operators with each order, equal to size of 
%irreducible nth power rep of sym group U(N*(Kappa+size(cc2,2))) 
    end

if length(Kap1)>1 %this is an overloading to allow you to just pass nn 
    %straight to the function if it has already been calculated
   nn = Kap1;
   totfreq = horzcat(vv{:}).';
else
    
nn = zeros(sum(numwithn),tot_poles); %if Kap1 <inf then less elements may be present
totfreq = horzcat(vv{:}).';
count1=0; tmp3 = 0;
for k =2:Kap2+1

        count2 = count1 + size(tmp3,1); 
        
        tmp3 = zeros(length(count1+1:count2)*size(nn,2),size(nn,2));
       for j = count1+1:count2
        tmp = repmat(nn(j,:),size(nn,2),1); %select element from tier below
        tmp2 = tmp+eye(size(nn,2)); %add one to every position
        tmp3(1+(j-(count1+1))*size(nn,2):(j-count1)*size(nn,2) ,:) = tmp2;       
       end
       %remove any elements of tmp3 which fail to satisfy the high
       %frequency cut off condition
       lg = true(size(tmp3,1),1);
       if isfinite(Kap1)
            for j =1:length(count1+1:count2)*size(nn,2)
                lg(j) = tmp3(j,:)*real(totfreq) < Kap1;
            end
       end
       tmp3 = tmp3(lg,:); %remove states that don't satisfy the cut off
       
       %next remove duplicate entries
       jj=1;
       while jj < size(tmp3,1)-1
           lg = any(tmp3(jj+1:end,:)  ~= repmat(tmp3(jj,:),size(tmp3,1)-jj,1),2);
           tmp3 = tmp3([true(jj,1);lg],:);
           jj=jj+1;
       end
       count1 = count2;
       nn(count1+1:count1 + size(tmp3,1),:) = tmp3;
end
end

%would be good to include a complicated truncation condition such that
%weakly coupled modes could be taken only to lower orders

%clean any nn values which are all zero besides the first
nn = nn([true;any(nn(2:end,:)~=0,2)],:);
tierofnn = sum(nn,2);
for k  =2:Kap2+1    
    numwithn(k) = sum(tierofnn == (k-1)); %set actual number included
end
else
   nn=0; numwithn=1; %trivial case, DON'T ENTER Kap1 as nn and also Kap2 =0;
end

%toc
%tic

%% Calculate the coupling from the heirarchy

%general up and down coupling terms
coup_com = zeros(size(nn,1));  %coefficient will be the sqrt of number it couples to (upward)
%coup_p1_save= zeros(size(nn,1),size(nn,1),N);  %takes up too much memory

coup_com_save{N}= sparse(coup_com);

%must also have the term that accounts for the fact that down coupling
%terms are different except for the matsubara terms
%coup_m1 = coup_com ;
coup_m1_acom = coup_com ;  %anti commutator part
coup_acom_save= coup_com_save;

rng_j =0; %full_rng = 1:totpoles; 
if nargin ~=6 %calculate for forward acting direction
 
for j = 1:N
    if ~isempty(cc_com{j}) %no bath at this mode for some reason
    rng_j = rng_j(end)+1:length(cc_com{j})+rng_j(end); 
    rng_rest = [1:rng_j(1)-1 , rng_j(end)+1:totpoles];
    cc = cc_com{j};  ccI = cc_acom{j};
for k =1:sum(numwithn(1:end-1))
    
    currenttier = tierofnn(k);
    tierabove = currenttier+1;

    tierpickr = abs(tierofnn-tierabove)<eps(10); %picks correct tier

    nnsec = nn(k,rng_j); %section of interest
    nnsec2 = nn(k,rng_rest); %rest
    %which obviously must still match, but must be element of j section
    %changed

    temp0 = repmat(nnsec,numwithn(tierabove+1),1); 
    %make same size as nn in tier above in dimension 1
    temp02 = repmat(nnsec2,numwithn(tierabove+1),1); 
    
    temp = temp0 - nn(tierpickr,rng_j); 
    %take away elements in the tier above

    temp4 = temp02 - nn(tierpickr,rng_rest);
    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    
    temp3 = sum(temp,2) < -1/2 & sum(temp,2) > -3/2; %ones with a -1 net 
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    tierpickr(tierpickr) = comb_lg ; 
    %use to put elements in the right row position in coupling matrix,
    %comb_lg is displaced by the number of tiers below
    
    if any(comb_lg)
        
    [~,temp44] = find(temp(comb_lg,:)); %find row position, i.e. which coefficient
    temp7 = sqrt(sqrt(abs(cc(temp44)).^2+abs(ccI(temp44)).^2));
    coup_com(k,tierpickr)=temp7.*sqrt(1+nnsec(temp44));
    
    %elements with minus are in the positions which are the transpose of this 
    %due to the rescaling the coefficients should be the same up and down
    %for the matsubara terms and the commutator part of the drude terms
    %now second down coupling operator
    
    coup_com(tierpickr,k) = cc(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
       % ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));
     %  if any(cc_acom{j}(temp44))
    coup_m1_acom(tierpickr,k) = ccI(temp44).*sqrt(1+nnsec(temp44))./temp7;%...
      %  ./ sqrt(sqrt(abs(cc(j,temp44)).^2+abs(ccI(j,temp44)).^2));    
     %  end
    end

end
    %Note I am scaling things w.r.t. (abs(a)^2+abs(b)^2)^(1/4)
    %for a term with operator a V^O + b V^X
    
    end

coup_com_save{j} = sparse(coup_com); coup_com = 0*coup_com;
coup_acom_save{j}= sparse(coup_m1_acom); coup_m1_acom=coup_com;
  %down coupling from real frequency and residue modes is just  coup_p1.'
end
else %if sixth arg given do the same for left acting (operator t_dep)
    
    
    
for j = 1:N
    if ~isempty(cc_com{j}) %no bath at this mode for some reason
    rng_j = rng_j(end)+1:length(cc_com{j})+rng_j(end); 
    rng_rest = [1:rng_j(1)-1 , rng_j(end)+1:totpoles];
    cc = cc_com{j};  ccI = cc_acom{j};
for k =2:sum(numwithn(1:end))
    
    currenttier = tierofnn(k);
    tierbelow = currenttier-1;

    tierpickr = abs(tierofnn-tierbelow)<eps(10); %picks tier below

    nnsec = nn(k,rng_j); %section of interest
    nnsec2 = nn(k,rng_rest); %rest

    temp0 = repmat(nnsec,numwithn(tierbelow+1),1); 
    temp02 = repmat(nnsec2,numwithn(tierbelow+1),1); 
    
    temp = temp0 - nn(tierpickr,rng_j);    
    temp2 = sum(temp~=0,2) ==1; %only one element +/- 1 diff
    temp3 = sum(temp,2) > +1/2 & sum(temp,2) < +3/2; %ones with a +1 net 
    temp4 = temp02 - nn(tierpickr,rng_rest);
    temp4 = all(temp4==0,2); %match rest of the terms
    
    comb_lg = temp2 & temp3 & temp4; %full condition required
    tierpickr(tierpickr) = comb_lg ; 
    %use to put elements in the right row position in coupling matrix,
    %comb_lg is displaced by the number of tiers below

    if any(comb_lg)
        
    [~,temp44] = find(temp(comb_lg,:)); %find row position, i.e. which coefficient
    temp7 = sqrt(sqrt(abs(cc(temp44)).^2+abs(ccI(temp44)).^2));
    %down coupling terms
    coup_com(tierpickr,k) = cc(temp44).*sqrt(1+nnsec(temp44))./temp7;
    coup_m1_acom(tierpickr,k) = ccI(temp44).*sqrt(1+nnsec(temp44))./temp7;
 
    %up coupling term
    coup_com(k,tierpickr)=temp7.*sqrt(1+nnsec(temp44));

    end

end
    end

coup_com_save{j} = sparse(coup_com); coup_com = 0*coup_com;
coup_acom_save{j}= sparse(coup_m1_acom); coup_m1_acom=coup_com;
end    
    
    
end



% calculate the sum_{j=1}^N sum_{k=0}^Kappa n_{jk} v_{ik}
if size(nn,1)==1
  const_factor = 0; 
else
    %reshape(vv.',numel(vv),1)
    %([vv(1,:).';vv(2,:).';vv(3,:).';vv(4,:).'])
    
const_factor = nn*totfreq;%DECAY FACTOR
end
%toc
