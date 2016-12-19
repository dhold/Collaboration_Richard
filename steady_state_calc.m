
Temperature_in_Kelvin = 300;
[t_sc,B] = inv_cm_unit_sys_standard(Temperature_in_Kelvin); %B = 1/kB T

mode_inc = 3;
switch mode_inc
    case 0 %all
lg = zeros(48,1);    
    case 1 %based on if BOmodes_tmp(:,3)>500 & BOmodes_tmp(:,1)<15 then exclude
lg = [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,1,1,1,1,1 ...
        ,0,1,1,1,1,1,1,0,1,1,1,1,1,1]'; %value of 1 removes from HEOM
    case 2 %BOmodes_tmp(:,3)>600 | BOmodes_tmp(:,1)<20
lg = [1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 ...
      ,1,1,1,1,1,1,1,1,1,1,1,1,1,1]';   %takes only 2 modes     
    case 3 %to take all out   
lg = ones(48,1) ; 
end
lg2 = true(48,1); %local vector to select particular modes
max_tier = 3; Kappa = 2;
specden;  %call to get all the info including spec den and HEOM and Hamiltonian

[M_e,E1] = eig(H_1); %eigenstate basis
%%
 %calculate some coupling rates with my own code
int_calc_point = 100000;
%calculate "Marcus" rate
[R_12,R_21,g_broad1,g_broad_ex1,lam_ex1,MM1] = Forster_outcouple(B,H_2, ...
            [true;false],Drude_modes_ex(7:8),BO_modes_ex(7:8),int_calc_point);
%check detail balance
EI = MM1{2}{2}-lam_ex1(1); Ealpha = MM1{3}{2}-lam_ex1(2);
DE1 = EI - Ealpha; DBC = exp(-B*DE1)*R_12/R_21; discrep1 = DBC-1;
        
%[out_rates,M_mrf,g_broadtest] = mod_red_outcouple(B,H_2,[true;false],...
%    Drude_modes_ex(7:8),BO_modes_ex(7:8));
%calculate gen Forster rates
[M_12,M_21,g_broad,g_broad_ex,lam_ex,MM] = Forster_outcouple(B,H_12, ...
            [true(6,1);false],Drude_modes_ex(1:7),BO_modes_ex(1:7));
Eex = MM{2}{2}-lam_ex(1:end-1); EI = MM{3}{2}-lam_ex(end);
 DE =    Eex-EI; DBC = exp(-B*DE).*M_12./(M_21.');  discrep = DBC-1;

%{
 Richard params
forward exciton -> CT: [  4.46166470e-01   2.79877206e-03   1.22445590e+01   7.79503568e-01   
4.32373274e-01   1.94641875e-01] (lowest to highest exciton, i think including the reorg energies)

backward CT -> exciton: [  8.00480726e-02   5.08904514e-04   1.16400212e+00   1.14238066e-01   
4.40222557e-02   5.38981883e-03] (again lowest to highest exciton)

forward I -> alpha: 0.361221381544

backward alpha -> I: 0.0417627640778        
%}        
 %% parameters from Richard Stones
Rich1 = [  4.46166470e-01   2.79877206e-03   1.22445590e+01   7.79503568e-01   4.32373274e-01   1.94641875e-01];
Rich2 = [  8.00480726e-02   5.08904514e-04   1.16400212e+00   1.14238066e-01   4.40222557e-02   5.38981883e-03] ;
RI_a = 0.361221381544; Ra_I = 0.0417627640778  ; 
%predict from DBC on ALL states included here
Ra_IDBC = exp(-B*DE1)*RI_a;
Rich2DBC = exp(-B*DE.').*Rich1 ;

%% Derive Liouvillian Terms

L_so = @(aa) kron((aa').',aa) - kron(speye(size(aa)),(aa')*aa)/2 ...
               - kron(((aa')*aa).',speye(size(aa)))/2  ; %Lindblad form

 n_elec = 60000; %number of photons or some shit
 gamma_ex = 0.00125; %transfer rate to lowest exciton
bb = zeros(10); bb(1,5) = 1; bb=sparse(bb); %annihilation of lowest ex   
%derive the rate of transfer between ground and lower exciton state
L_pump =  gamma_ex*sparse(( (n_elec+1)*L_so(bb) + n_elec*L_so(bb'))); 

%next the operator that deals with the generalised Forster transfer
L_genF = L_pump*0;
for k = 1:6
bb = zeros(10); bb(4,4+k) = 1;
L_genF = L_genF + Rich1(k)*L_so(bb) + Rich2(k)*L_so(bb');
end
%and the Marcus transfer
bb = zeros(10); bb(3,4) = 1; 
L_Marcus = RI_a*L_so(bb) +  Ra_I*L_so(bb');
%next dealing with alpha -> beta decay
GammaR = 201.638; bb = zeros(10); bb(2,3) = 1; 
L_ab = GammaR*L_so(bb); %only one way
% finally beta -> g decay
GammaL = 201.638; bb = zeros(10); bb(1,2) = 1; 
L_bg = GammaL*L_so(bb); %only one way

L_total = L_pump + L_genF + L_Marcus + L_ab + L_bg; %only acts on first tier?
%This assumes that essentially the correlations with the bath don't effect
%these transitions rates, which they will for L_genF... but hey

H_tot = diag([zeros(4,1);diag(E1)]); %Hamiltonian part, diagonal 
%L_tot = L_total-1i*(kron(speye(10),H_tot)-kron(H_tot.',speye(10)));
L_tot = -1i*(kron(speye(10),H_tot)-kron(H_tot.',speye(10)));
%% Generate full prop op for the dynamics
V_coup_tot = zeros(10,10,N); V_coup_sys = zeros(N,N,N);
for j = 1:N %project each term into the exciton basis
    tmp = zeros(N); tmp(j,j) = 1;    V_1 = M_e'*tmp*M_e;
V_coup_tot(:,:,j) = blkdiag(zeros(4),V_1); %zeros are g, alpha, beta, I
V_coup_sys(:,:,j) = V_1; %zeros are g, alpha, beta, I
end                

%prop_op_full =H_prop_gen3(L_tot,[],V_coup_tot,V_coup_tot,...
%                QQ_topass,const_factor,coup_com_save,coup_acom_save);
prop_op_sys =H_prop_gen3(E1,E1,V_coup_sys,V_coup_sys,...
                QQ_topass,const_factor,coup_com_save,coup_acom_save);            
            
%prop_op_full(1:100,1:100) = prop_op_full(1:100,1:100)+L_total; %HEOM + rate
%this is one mother of a matrix!

%% Calculate system only dynamics

propeq_sys = @(t,v)  prop_op_sys*v;

rho_vec = zeros(N^2*size(nn,1),1);  rho_vec(1)=1;
t_range = linspace(0,0.05,100); 

rnj = 1:N^2; pop_lg = false(N^2,1); pop_lg(1:N+1:N^2)=true;

%options = odeset('OutputFcn',@ODE_output_trim);   
%ODE_output_trim(t_range,rnj,'input'); ode45(propeq_sys,t_range,rho_vec,options);
%[tout,rho_out] = ODE_output_trim([],[],'get_data');
[tout,rho_out] = ode45(propeq_sys,t_range,rho_vec);
%%
figure; plot1 = plot(tout,rho_out(:,pop_lg)); 
for k =1:sum(pop_lg)
set(plot1(k),'DisplayName',strcat('exciton ',num2str(k)));
end
xlabel('time (ps)'); ylabel('exciton population')
%%
tmp = reshape(rho_out(:,1:length(pop_lg)).' ,[N,N,length(tout)]);
coherence_set = zeros(sum(~pop_lg)/2,length(tout)); cnt = 0;
for k1 = 1:N; for k2=k1+1:N;
cnt = cnt+1;
coherence_set(cnt,:) = tmp(k1,k2,:)./sqrt(tmp(k1,k1,:).*tmp(k2,k2,:));
end; end     
        
figure; plot2 = plot(tout,abs(coherence_set)); cnt=0;
for k1 = 1:N; for k2=k1+1:N;
cnt = cnt+1;
set(plot2(cnt),'DisplayName',strcat('\rho_{',num2str(k1),num2str(k2),'}'));
end; end     

xlabel('time (ps)'); ylabel('normalised coherence $|\rho_{ab}|/\sqrt{\rho_{aa}\rho_{bb}}$')
 %%  Trim it down so it doesn't include coherences that can never form
 %at least I think they can't
 
 lg_mat = blkdiag(eye(4),ones(6));
 lg_mat = reshape(lg_mat,numel(lg_mat),1); 
 lg_matHEOM = logical(repmat(lg_mat,[size(nn,1),1]));
 pop_lg = logical([reshape(eye(10),100,1);zeros((size(nn,1)-1)*100,1)]);
 pop_lg = pop_lg(lg_matHEOM);
 
 prop_op_red = prop_op_full(lg_matHEOM,lg_matHEOM);
 clear prop_op_full %this actually takes up a fair bit of memory
 %% Calculate steady state
 tol  = max(2*max(size(prop_op_red)) * norm(prop_op_red,1) * eps,200*eps);
 %[~, rho_ss] = spspaces(prop_op_red,2);
[rho_ss,V]= eigs(prop_op_red,1,'SM');

rho_sys = zeros(10);  rho_sys(logical(lg_mat)) = rho_ss(1:sum(lg_mat(:)));
rho_sys = rho_sys/trace(rho_sys);
figure; bar(diag(abs(rho_sys)))

%% Can calculate dynamics as well I guess

%initial state -> take single quanta in the lowest exciton state
%bath I guess is still doing whatever it was in the ground state
rho_init = zeros(10); rho_init(5,5) = 1; rho_init = reshape(rho_init,100,1);
rho_vec = zeros(sum(lg_mat(:))*size(nn,1),1); 
rho_vec(1:sum(lg_mat(:)))= rho_init(logical(lg_mat));

propeq = @(t,v)  prop_op_red*v;

rnj1 = 1:sum(lg_mat(:)); %full primary system 
rnj2 = 5:7:sum(lg_mat(:)); %exciton populations

options = odeset('OutputFcn',@ODE_output_trim);   

t_range1 = linspace(0,0.1,500);   t_range2 = linspace(0.1,1,500); 
ODE_output_trim(t_range1,rnj1,'input'); ode45(propeq,t_range1,rho_vec,options);
[tout1,rho_out1] = ODE_output_trim([],[],'get_data');
%[tout1,rho_out1] = ode45(propeq,t_range,rho_vec);
ODE_output_trim(t_range2,rnj1,'input'); ode45(propeq,t_range2,rho_out1(:,end),options);
[tout2,rho_out2] = ODE_output_trim([],[],'get_data');
%[tout2,rho_out2] = ode45(propeq,t_range,rho_out1(end,:).');

%figure; plot([tout1;tout2(2:end)],[rho_out1(:,rnj);rho_out2(2:end,rnj)])
figure; plot([tout1;tout2(2:end)],[rho_out1,rho_out2(:,2:end)])
