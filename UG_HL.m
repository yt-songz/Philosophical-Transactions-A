% clear;
% clc;

N_P = 100;
N_R = 100;

% h=0.5;
% l=0.1;
% M_P = 5;
% M_R = 5;

X=[]; %HH (Proposor, receiver)
Y=[]; %HL
Z=[]; %LH 
W=[]; %LL

theta_PH = 1;
theta_PL = 0;
theta_RH = 1;
theta_RL = 0;
MP_vec = [];
MR_vec = [];
s_vec = [];
h_vec = [];
l_vec = [];
% s = 1;
startValue=0;
endValue=100;
step=5;
% M_R = 0;

for aaa=0:1:4
    s = 0.1 * 10^aaa;
    for h=0.4:0.01:0.6
        for l=0.1:0.01:0.3
            for M_R=startValue:step:endValue
                for M_P = startValue:step:endValue
                
                mp = M_P;
                mr = M_R;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ת�Ƹ��ʾ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TransitionProbabilityMatrix=zeros(4,4);
            
                %X-Y --- 1,2 HH and HL (H proposer population, and transition between H and L receivers)
                %k_P = N_P --- all proposers are HP
                ggg=0;
                for k=0:N_R-1
                    hhh=1;
                    P_xy=0;
                    P_yx=0;
                    for m=1:k
            
                        x=m;
                        P_xy=h*(N_P+M_P*theta_PH)/(N_P+M_P);
                        P_yx=(h*(N_P+M_P*theta_PH)+l*(N_P-N_P+M_P*(1-theta_PH)))/(N_P+M_P);
            
                        T = x/(N_R-x)*(N_R-x+theta_RL*M_R)/(x+(1-theta_RL)*M_R);
                        hhh=hhh*T*(1+exp(s*(P_yx-P_xy)))/(1+exp(s*(P_xy-P_yx)));
            
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(1,2)=1/ggg;
                
            
                %X-Z --- 1,3 HH and LH (H receiver population, and transition between H and L proposers)
                %k_R = N_R --- all recerivers are HR
                ggg=0;
                for k=0:N_P-1
                    hhh=1;
                    P_xz = 0;
                    P_zx = 0;
                    for m=1:k
                        x=m;
                        P_xz=1-h;
                        P_zx=(1-l)*(M_R*theta_RL)/(N_R+M_R); % H receiver is considered
            
                        T=x*(N_P-x + (1-theta_PH)*M_P)/(N_P-x)/(x+theta_PH*M_P);
                        hhh=hhh*T*(1+exp(s*(P_zx-P_xz)))/(1+exp(s*(P_xz-P_zx)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(1,3)=1/ggg;
                
                %X-W --- 1,4 HH and LL (no transition)
            
                TransitionProbabilityMatrix(1,4)=0;
               
               
                %Y-X --- 2,1 HL and HH (H prposer population, and transition between L and H receivers)
                %k_P = N_P --- all recerivers are HR
                %number of H is N_R-m
                ggg=0;
                for k=0:N_R-1
                    hhh=1;
                    P_yx=0;
                    P_xy=0;
                    for m=1:k
                        x=N_R-m;
            
                        P_yx=(h*(N_P+M_P*theta_PH)+l*(M_P*(1-theta_PH)))/(N_P+M_P);
                        P_xy=h*(N_P+M_P*theta_PH)/(N_P+M_P);
            
                        T =(N_R-x)/x*(x+(1-theta_RL)*M_R)/(N_R-x+theta_RL*M_R);      
            
                        hhh=hhh*T*(1+exp(s*(P_xy-P_yx)))/(1+exp(s*(P_yx-P_xy)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(2,1)=1/ggg;
                
            
                %Y-Z --- 2,3 HL and LH
                
                TransitionProbabilityMatrix(2,3)=0;
            
            
                %Y-W --- 2,4 HL and LL (L receiver population, and transition between H and L proposers)
                %k_R = 0 --- all recerivers are LR
                %number of H is m
                ggg=0;
                for k=0:N_P-1
                    hhh=1;
                    P_yw=0;
                    P_wy=0;
                    for m=1:k
                        x = m;
            
                        P_yw=1-h;
                        P_wy=(1-l)*(N_R+M_R*(1-theta_RH))/(N_R+M_R);
            %           
                        T=x*(N_P-x + (1-theta_PH)*M_P)/(N_P-x)/(x+theta_PH*M_P);
            
                        hhh=hhh*T*(1+exp(s*(P_wy-P_yw)))/(1+exp(s*(P_yw-P_wy)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(2,4)=1/ggg;
                
                
                %Z-X --- 3,1 LH and HH (H receiver population, and transition between L and H proposers)
                %k_R = N_R --- all recerivers are HR
                %number of H is N_P-m
                ggg=0;
                for k=0:N_P-1
                    hhh=1;
                    P_zx=0;
                    P_xz=0;
                    for m=1:k
                        x = N_P-m;
            
                        P_zx=(1-l)*(M_R*(1-theta_RH))/(N_R+M_R);
                        P_xz=1-h; % H receiver is considered
            
                        T=(N_P-x)*(x+theta_PH*M_P)/x/(N_P-x + (1-theta_PH)*M_P);
                        hhh=hhh*T*(1+exp(s*(P_xz-P_zx)))/(1+exp(s*(P_zx-P_xz)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(3,1)=1/ggg;
                
                %Z-Y --- 3,2 LH and HL
                
                TransitionProbabilityMatrix(3,2)=0;
                
                %Z-W --- 3,4 LH and LL (L proposer population, and transition between H and L receiver)
                %k_P = 0 --- all proposers are LP
                %number of H is m
                ggg=0;
                for k=0:N_R-1
                    hhh=1;
                    P_zw=0;
                    P_wz=0;
                    for m=1:k
                        x=m; 
                        P_zw=h*(0+M_P*theta_PH)/(N_P+M_P);
                        P_wz=(h*(0+M_P*theta_PH)+l*(N_P-0+M_P*(1-theta_PH)))/(N_P+M_P);
                        T = x/(N_R-x)*(N_R-x+theta_RL*M_R)/(x+(1-theta_RL)*M_R);
                        hhh=hhh*T*(1+exp(s*(P_wz-P_zw)))/(1+exp(s*(P_zw-P_wz)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(3,4)=1/ggg;
                
                
                %W-X --- 4,1 LL and HH
                
                TransitionProbabilityMatrix(4,1)=0;
                
                %W-Y --- 4,2 LL and HL (L receiver population, and transition between L and H proposer)
                %k_R = 0 --- all receivers are LR
                %number of H is N_P-m
                ggg=0;
                for k=0:N_P-1
                    hhh=1;
                    P_wy=0;
                    P_yw=0;
                    for m=1:k
                        x=N_P-m;
                        P_wy=(1-l)*(N_R+M_R*(1-theta_RH))/(N_R+M_R);
                        P_yw=1-h;
                        T=(N_P-x)*(x+theta_PH*M_P)/x/(N_P-x + (1-theta_PH)*M_P);
                        hhh=hhh*T*(1+exp(s*(P_yw-P_wy)))/(1+exp(s*(P_wy-P_yw)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(4,2)=1/ggg;
                
                %W-Z --- 4,3 LL and LH (L proposer population, and transition between L and H receiver)
                %k_P = 0 --- all proposers are LP
                %number of H is x=N_R-m
                ggg=0;
                for k=0:N_R-1
                    hhh=1;
                    P_wz=0;
                    P_zw=0;
                    for m=1:k
                        x=N_R-m; 
                        P_wz=(h*(0+M_P*theta_PH)+l*(N_P+M_P*(1-theta_PH)))/(N_P+M_P);
                        P_zw=h*(0+M_P*theta_PH)/(N_P+M_P);
            
                        % T1=(N_R-m)/N_R * (m+(1-theta_PH)*M_R)/(N_R+M_R)/(1+exp(s*(P_zw-P_wz))); %L+
                        % T2=m/N_R * (N_R-m+theta_PH*M_R)/(N_R+M_R)/(1+exp(s*(P_wz-P_zw))); % L-
                        % hhh=hhh*T2/T1;
            
            
                        T =(N_R-x)/x*(x+(1-theta_RH)*M_R)/(N_R-x+theta_RL*M_R);
                        hhh=hhh*T*(1+exp(s*(P_zw-P_wz)))/(1+exp(s*(P_wz-P_zw)));
                    end
                    ggg=ggg+hhh;
                end
                TransitionProbabilityMatrix(4,3)=1/ggg;
              
                
                for mm=1:4
                    TransitionProbabilityMatrix(mm,mm)=1-sum(TransitionProbabilityMatrix(:,mm))/2;
                end
            
            
                %%%ֻ��C�˳�
            %     TransitionProbabilityMatrix=TransitionProbabilityMatrix(1:3,1:3);
                %%%ֻ��D�˳�
            %    TransitionProbabilityMatrix=TransitionProbabilityMatrix([1,2,4],[1,2,4]);
                %%%%%both �����ı�
                [xx,yy]=eig(TransitionProbabilityMatrix);%x��ÿһ��ֵ��ʾ����A��һ������������y�ĶԽ�Ԫ��ֵ��������A������ֵ
                [row,list]=find(yy==max(max(yy)));
                
                eigenvector_normalized=xx(:,list)./sum(xx(:,list));%����ֵΪ1������������һ��
            %     MMM=[MMM;M];
            %     sss=[sss;s];
            %     sigma=[sigma;sigma2];
                X=[X;double(eigenvector_normalized(1))];
                Y=[Y;double(eigenvector_normalized(2))];
                Z=[Z;double(eigenvector_normalized(3))];
                W=[W;double(eigenvector_normalized(4))];
                MP_vec = [MP_vec; mp];
                MR_vec = [MR_vec; mr];
                s_vec = [s_vec; s];
                l_vec = [l_vec; l];
                h_vec = [h_vec; h];
                Results=[s_vec h_vec l_vec MP_vec MR_vec X Y Z W];
            
                end 
            end
        end
    end
end
       
