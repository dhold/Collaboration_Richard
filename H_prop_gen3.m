function [prop_op]=H_prop_gen3(H_1,H_2,V1,V2,QQ,const_factor,coup_com,coup_acom,flg)
   
%Generates HEOM given coupling coefficients and decay factors and
%what type of state there is e.g. |e><g|, |e><e| or |f><g>
%Hamiltonian 1 is the left state Ham hamiltonian 2 is the right state one
% so for |e><g| , H_1 = H_e,  H_2 = H_g,  V_1 = delta_{ej} V_2 = 0;
% one calls    H_prop_gen2(H_e,H_g,fock_rep1,zeros(1,N),QQ,...)
%
% for rho^(uv) -> H_1 = H_uu , H_2 = H_vv,  V1 = V_uu, V2 = V_vv;
% for operator V^{vu} -> H_1 = H_uu, .. same

%% Construct Ltot, the total Louville for self coupling density matricies

% For NN_1 X NN_1 A and NN_2 X NN_2 B we have the following
% reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
% reshape(C * B, NN_1*NN_2,1) = kron(B.',eye(NN_1))*reshape(C,NN_1*NN_2,1)
% reshape(A * C * B, NN_1*NN_2,1) = kron(B.',A)*reshape(C, NN_1*NN_2,1)
% Hence H rho - rho H can be expressed with a flattened rho as L *rho

N = length(coup_com); HL = length(coup_acom{1});

if nargin <9 %normal -i* (H_uu O_uv - O_uv H_vv)
    if ~isempty(H_2)    
 L = -1i*(kron(eye(length(H_2)),H_1)-kron(H_2.',eye(length(H_1))));          
    else
 L = H_1; %pass directly, use if additional decoherence terms included        
    end
else %assume operator time dep -i* (O_vu H_uu  - H_vv O_vu )
    if ~isempty(H_2)    
 L = -1i*(kron(H_1.',eye(length(H_2)))-kron(eye(length(H_1)),H_2) );          
    else
 L = H_1; %pass directly, use if additional decoherence terms included        
    end
end
      %include decay factor first  
temp = kron(const_factor,ones(length(L),1));  
prop_op = -sparse(1:length(temp),1:length(temp),temp);

if size(QQ,1) == 1 %assume same for every site
    QQ = repmat(QQ,N*viblvls,1);
end
 Q = sparse(size(L,1),size(L,2));

for j = 1:N %site loop 
    %(usually away, unless you have some strange system bath interaction)
    
    if ndims(V1) == 3 %can give general form for exciton basis etc
    V_1 = V1(:,:,j);
    V_2 = V2(:,:,j);
    else %assume diagonal
        V_1 = diag(V1(:,j)); V_2 = diag(V2(:,j)); 
    end
% reshape(A * C, NN_1*NN_2,1) = kron(eye(NN_2),A)*reshape(C, NN_1*NN_2,1)
% reshape(C * B, NN_1*NN_2,1) = kron(B.',eye(NN_1))*reshape(C,NN_1*NN_2,1)  
 %reshape(A * B * A, N^2,1) = kron(A.',A)*reshape(B, N^2,1)
 
    %[Q_i,[Q_i,rho_n]] = (1 X Q_i*Q_i + Q_i^T*Q_i X 1 - 2 Q_i^T X Q_i ) rho_vec
    % with rho_vec the flattened (vector) version of rho_n (matrix)
    %[Q_uv ,rho_uv] = Quu rho_uv - rho_uv Q_vv 
    if nargin <9 %(c_k V_uu p_uv - c_k^* p_uv V_vv)
    %Qcom = kron(V_1.',eye(size(V_2))) - kron(eye(size(V_1)),V_2); %commtator
    Qcom =  kron(speye(size(V_2)),sparse(V_1))-kron(sparse(V_2.'),speye(size(V_1)));
    
    Qacom = kron(speye(size(V_2)),sparse(V_1))+kron(sparse(V_2.'),speye(size(V_1))) ; %anti commutator
    %keep factor of -1i from -1i tilde(V)^X term
    %factor of -1i cancels with +1i in front of imag part
    else %(c_k O_vu p_uu  - c_k^* p_vv O_vu) = 
        %re(c_k) (O_vu p_uu - p_vv O_vu)+ i*im(c_k) (O_vu p_uu + p_vv O_vu)  
     
        
    Qcom = kron(sparse(V_1.'),speye(size(V_2))) - kron(speye(size(V_1)),sparse(V_2)); %left acting commutator
    
    Qacom = kron(sparse(V_1.'),speye(size(V_2))) + kron(speye(size(V_1)),sparse(V_2)); %anti commutator    
   % Q1 =   Qcom*Qcom
   % Q2 = kron((V_1^2).',eye(size(V_2))) + kron(eye(size(V_1)),V_2^2)...
   %         - kron(V_1.',V_2) - kron(V_1,V_2.')
    end

    Q = Q + QQ(j,1).*Qcom*Qcom + QQ(j,2).*Qacom*Qcom;
  
    %kron(A,Qcom) places the matrix "Qcom" at every single position in
    %superoperator corresponding to a mixing between different tiers, the
    %coefficient is given by kron A, Qcom
    prop_op = prop_op  -1i * kron(sparse(coup_com{j}),Qcom) +...
                        kron(sparse(coup_acom{j}),Qacom);
    
end

     prop_op = prop_op + kron(sparse(eye(HL)),sparse(L-Q));
    
