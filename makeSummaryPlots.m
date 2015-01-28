
% makeSummaryPlots - This script runs each of the different criterions 
%   and creates plots.  Chane the switch value below to run different 
%   criterion.  Change the q values (pre-scaler) and p (exponents) 
%   as needed.


clpear, close all, clc

%Muscle Moment Arms (m)
D=[0.050  -0.0620 -0.0720     0.0340      0       0       0       0;...
    0       0       -0.0340     0.0500      0.0420  -0.0200 0       0;...
    0       0       0           0           0       -0.0530 -0.0530 0.037];
%Muscle Max Isometric Force (N)
fMax=[1917  1967    3878       1718        8531    2596    3734    1233];

%Target Joint Moments (N-m)
mTarget=[-50;50;-50];

%Initial Guess for Muscle Forces
fo=zeros(1,8);



switch 0   %Change this value to run different static optmizations
    
    case 0   % Perform SO with Stress Directly
        q=1;  %The pre-scaling coeffecient
        for i=1:9
            [f(:,i),m(:,i),s,info]=z2_polyOptStressScaling(D,fMax,...
                mTarget,fo,i+1,q);
        end
        legText={'p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10'};
        
    case 1  % Perform SO with polynomial criterion
        q=10;  %The pre-scaling coeffecient
        for i=1:9
            [f(:,i),m(:,i),s,info]=z2_polyOpt(D,fMax,mTarget,fo,i+1,q);
        end
        legText={'p=2','p=3','p=4','p=5','p=6','p=7','p=8','p=9','p=10'};
        
    case 2  %Perform SO with higher p and with minmax
        q=10; %The pre-scaling coeffecient
        
        p=[2,3,5,9,16,28,50];
        for i=1:length(p)
            [f(:,i),m(:,i),s,info]=z2_polyOpt(D,fMax,mTarget,fo,p(i),q);
        end 
        [fMinmax,mMinmax,s,info]=minMaxOpt(D,fMax,mTarget,fo);
        f=[f,fMinmax];
        m=[m,mMinmax];
        legText={'p=2','p=3','p=5','p=9','p=16','p=28','p=50','minmax'};
        
        
    case 3  %Perform SO with activation directly (OpenSim like)
        q=10;
        p=[2,3,5,9,16,28,50];
        for i=1:length(p)
            [f(:,i),m(:,i),s,info]=z2_polyOptActivation(D,fMax,mTarget,...
                fo,p(i),q)
        end
        [fMinmax,mMinmax,s,info]=minMaxOpt(D,fMax,mTarget,fo);
        f=[f,fMinmax];
        m=[m,mMinmax];
        legText={'p=2','p=3','p=5','p=9','p=16','p=28','p=50','minmax'};
        
    case 4  %Perfom SO with maxnorm poly
        q=1;
        p=[2,3,5,9,16,28,50];
        for i=1:length(p)
            [f(:,i),m(:,i),s,info]=z2_polyOptWithOuterExponent(D,fMax,...
                mTarget,fo,p(i),q)
        end
        legText={'p=2','p=3','p=5','p=9','p=16','p=28','p=50'};
end

% Plot Muscle Force
figure
bar(f,'grouped')
set(gca,'fontsize',22)
a={'Ilipos', 'Glut',  'Hams',    'Rect',   'Vast',  'Gast', 'Soleus','Tib'};
set(gca,'xTickLabel',a)
ylabel('Muscle Force (N)')
set(gcf,'Position',[671 300 929 664])
l=legend(legText)
set(l,'fontsize',20)

%Plot Muscle Activation
for i=1:size(f,2)
    act(:,i)=f(:,i)./fMax';
end
figure
bar(act,'grouped')
set(gca,'ylim',[0 0.3]);
set(gca,'fontsize',22)
set(gca,'xTickLabel',a)
ylabel('Activation (N/N)')
set(gcf,'Position',[671 300 929 664])
l=legend(legText)
set(l,'fontsize',20)

%Plot Resulting Moments
figure
bar(m,'grouped')
set(gca,'fontsize',22)
a={'Hip','Knee','Ankle'};
set(gca,'xTickLabel',a)
ylabel('Joint Moment (N-m)')
l=legend(legText)
set(l,'fontsize',20)
set(gcf,'Position',[671 300 929 664])






