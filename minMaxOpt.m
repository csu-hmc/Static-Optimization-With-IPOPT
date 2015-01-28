function [fOpt,mResult,s,infoOut]=minMaxOpt(D,fMax,mTarget,fo)

%minMaxOpt - Perform Static Optimization with Minmax
%
%
%[fOpt,mResult,s,infoOut]=minMaxOpt(D,fMax,mTarget,fo)
%
%       Inputs:
%           D - Matrix (Joints by Muscle) of moments arms
%           fMax - Max Isometric Force for Each  (row vector)
%           mTarget - Moments at Each Joint (column vector)
%           fo - Initial guess for each force (row vector)
%
%       Outputs:
%           fOpt - Result Forces
%           mResults - result moments (should match mTarget)
%           s - value of objective function
%           infoOut - info returned from IPOPT

global c

c.D=D;
c.fMax=fMax;

fo=[fo 0];  %Append the slack variable

% Set the bounds and constraints
options.lb = zeros(1,size(fo,2));  % Lower bound on the variables.
options.ub = [c.fMax inf];  % Upper bound on the variables.

nMuscles=length(fo)-1;
infV(1:nMuscles)=-Inf;
z(1:nMuscles)=0;
options.cl = [mTarget;infV'];  % Lower bounds on the constraint functions.
options.cu = [mTarget;z'];   % Upper bounds on the constraint functions.

% Set the IPOPT options
%options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
%options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-7;

% The callback functions.
funcs.objective         = @objFunc;
funcs.constraints       = @constrFunc;
funcs.gradient          = @gradObjFunc;
funcs.jacobian          = @constraintJac;
sparseMat=constraintJac(fo);
funcs.jacobianstructure = @() sparseMat;


[fOpt infoOut] = ipopt(fo,funcs,options);

mResult=fOpt(1:end-1)*c.D';
mResult=mResult';
s=constrFunc(fOpt);

fOpt=fOpt(1:end-1)';
%---------------------------------------------------
function s=objFunc(f)
%The objective to minimize 
%   Input  f: muscle forces as provided by IPOPT

s=f(end);


%---------------------------------------------------
function ds_df=gradObjFunc(f)
%The gradient of the objective
%   Input  f: muscle forces as provided by IPOPT

ds_df=zeros(size(f,1),size(f,2));
ds_df(end)=1;

%---------------------------------------------------
function m=constrFunc(f)
%The constraints
%   Input  f: muscle forces as provided by IPOPT
%   Output m: moments at the joints and muscle stress (column vector)

global c

s=f(end);  %Max activation slack variable
f=f(1:end-1);  % Forces activation

mAct=f./c.fMax-s;

mMoments=c.D*f';
m=[mMoments;mAct'];

%---------------------------------------------------
function dm_df=constraintJac(f)
%The jacobian (sparse gradient) of the constraints
%   Input  f: muscle forces as provided by IPOPT
%   Output dm_df: change in moments with change in forces

global c
n=length(f);  % Number of variables 



% Stress equations

dm_dfAct=diag(1./c.fMax);
sVect=-ones(n-1,1);
dm_dfAct=[dm_dfAct sVect];


% Moment equations
nMoments=size(c.D,1);
dm_dfMoments=[c.D zeros(nMoments,1)];

dm_df=sparse([dm_dfMoments;dm_dfAct]);
