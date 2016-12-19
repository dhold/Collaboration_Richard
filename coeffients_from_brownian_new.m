function [cc1,cc2R,cc2I,vv1,vv2,truncation_correction,cc] = ...
    coefficients_from_brownian(lambda,gamma,om_0,Temp,Kappa,lam_dru,gam_dru)
%This function generates the coefficients for every
%  \tilde{V}_{j,m} = \int dtau exp(-gamma_{v,m} + i \lambda)(t-tau) 
                            % c_R V(tau)^x - i c_I V(tau)^o
% This assumes the spectrum is made up of under and overdamped modes, the
% set of frequencies is given as an input
% Also allows the inclusion of Drude modes
% cc1 is the commutator coefficient from the matsubara modes
% cc2R is the commutator coefficient from the J(w) modes
% cc2I is the anticommutator coefficient from the J(w) modes
% vv1 is the exponential coefficient from matsubara terms
% vv2 is the (potentially complex) exponential coefficients
% cc 
% Tested to give the correct results for Drude and underdamped 
if isempty(lam_dru) && isempty(lambda)
    cc1=[]; cc2R=[]; cc2I=[]; vv1=[]; vv2=[];
    truncation_correction=0; cc =[]; return
end


light_speed = 299792458; %units m s^{-1}, NOT cm 
    length_unit = 100; % units m, this is 1/[1(cm)^{-1}]
    hbar = 1.05457173*10^(-34); %units Js
    boltz_const = 1.3806488 * 10^(-23); 
Beta = (2 * pi * hbar * light_speed * length_unit)/ ( Temp * boltz_const);

  %first calculate the poles from (1-exp(-beta *omega))^(-1)
 vv1 = 2*pi*(1:Kappa)/Beta; %matsubara frequencies
  om_rng = -1i*vv1; %pole locations
 cc1 = om_rng*0;   truncation_correction = 0;
    cnt =1;
if ~isempty(om_0)    
lg = om_0./gamma/2 > 1; %under damped
lg2 = om_0./gamma/2==1; %critically damped
lg3 = om_0./gamma/2 < 1; %over damped


cc2R = zeros(1,2*sum(lg)+2*sum(lg3)+sum(lg2)+length(lam_dru)); 
cc2I = cc2R; vv2 = cc2R;
cc = [vv1,vv2]*0;
xi = zeros(size(lambda));

rnj1 = find(lg); rnj2 =find(lg2); rnj3 = find(lg3);
 for lp = 1:length(rnj1) %first underdamped most complicated due to the complex frequency
     k = rnj1(lp);
     xi(k) = sqrt(om_0(k)^2-gamma(k)^2/4); 
     %seen this also denoted gam/2, this is just a diff definition
     prefc = lambda(k)*om_0(k)^2/(2*xi(k));
     
     %polepos1 = xi(k) - 1i*gamma(k)/2; %exp(-(gamma/2 + i xi)) 
     %polepos2 = -xi(k) - 1i*gamma(k)/2; %exp(-(gamma/2 - i xi))
     vv2([cnt,cnt+1]) = gamma(k)/2 +[1,-1]* 1i*xi(k); %frequency, pole loc at -i*vv2
     cc2I([cnt,cnt+1]) = -1i*prefc*[1,-1]; % V^o coeff   
  cc2R([cnt,cnt+1]) = prefc*(-sinh(Beta*xi(k)) + 1i*[1,-1].*sin(Beta*gamma(k)/2))...
                        ./(cos(Beta*gamma(k)/2)-cosh(Beta*xi(k))); %V^x coefficient   
  cc([cnt,cnt+1]) = [1i,-1i].*prefc.*(cot(Beta*(gamma(k)/2 + [1i,-1i].*xi(k))/2)-1i);
                    
  cnt=cnt+2;
 end

  for lp = 1:length(rnj2)  %crit damped %not tested yet
        k = rnj2(lp);      
    cc2R(cnt) = lambda(k)*om_0(k)*Beta/2/(sin(Beta*om_0(k)/2))^2;   %V^x coedfficient  
     vv2(cnt) =  om_0(k);            %frequency
        cnt=cnt+1;
                
  end
 
  for lp = 1:length(rnj3)  %overrrrdamped, poles imaginary %not tested yet
     k = rnj3(lp);
     xi(k) = sqrt(gamma(k)^2/4-om_0(k)^2);
     polepos1 = sqrt(om_0(k)^2 - gamma(k)^2/2 + gamma(k)*xi(k));  
     polepos2 = sqrt(om_0(k)^2 - gamma(k)^2/2 - gamma(k)*xi(k)); 
     
    vv2(cnt) =  -1i*polepos1;
     tmp = lambda(k)*gamma(k)*om_0(k)^2/(2*om_0(k)^2-gamma(k)^2);
     %tmp = tmp*(1+1i*cot(Beta*polepos1/2)); 
    cc2R(cnt) = tmp;   
    cc2I(cnt) = -tmp*cot(Beta*polepos1/2);  %V^o coedfficient  
                %frequency
      cnt=cnt+1;
    cc2R(cnt) = real(tmp);   
    cc2I(cnt) = -tmp*cot(Beta*polepos2/2);     
     cnt=cnt+1;       
                
  end
 for k = 1:length(lambda)
 cc1 = cc1 + (lambda(k).*gamma(k).*om_0(k).^2 ).*(vv1./( (om_0(k).^2 ...
        - om_rng.^2).^2 + gamma(k).^2.*om_rng.^2));
  
    
 truncation_correction = truncation_correction + lambda(k)*(...
          (sin(Beta*gamma(k)/2) + gamma(k)*sinh(Beta*xi(k))/(2*xi(k)))/(...
          cos(Beta*gamma(k)/2)-cosh(Beta*xi(k))) + 2*gamma(k)/Beta/om_0(k)^2);

 end
else
    cc = zeros(1,length(lam_dru)+Kappa);
end
if ~isempty(gam_dru)
  for k = 1:length(gam_dru) %so overdamped treated with drude approx
    cc(cnt) = lam_dru(k)*gam_dru(k)*(cot(Beta*gam_dru(k)/2)-1i);       
    cc2R(cnt) = lam_dru(k)*gam_dru(k)*cot(Beta*gam_dru(k)/2);   
    cc2I(cnt) = -lam_dru(k)*gam_dru(k);  %V^o coefficient 
     vv2(cnt) =  gam_dru(k);            %frequency
        cnt=cnt+1;

 truncation_correction = truncation_correction + lam_dru(k)*(...
         2-Beta*gam_dru(k)*cot(Beta*gam_dru(k)/2))/Beta/gam_dru(k);     
     if ~isempty(cc1)
             cc1 = cc1 + lam_dru(k).*gam_dru(k).*vv1./( gam_dru(k).^2-vv1.^2)/Beta;
     end
  end
end

 cc1 = cc1*4; %pole residue, 
 %truncation_correction is  dot(cc1(j,Kappa+1:infinity),1./vv1(Kappa+1,2:infinity))
 
 cc(cnt:end) = cc1; %these are the actual pole resides, the same for all but the complex ones
 
  truncation_correction = truncation_correction - dot(cc1,1./(vv1));