function  [fOpt,m,s,infoOut]=z2_polyOptActivation(D,fMax,mTarget,fo,p,q)

%z2_polyOptActivation - Perform Static Optimization with Activation Units
%   and Polynomial Criterion similar to OpenSim
%
%[fOpt,m,s,infoOut]=z2_polyOptActivation(D,fMax,mTarget,fo,p,q)
%
%       Inputs:
%           D - Matrix (Joints by Muscle) of moments arms
%           fMax - Max Isometric Force for Each  (row vector)
%           mTarget - Moments at Each Joint (column vector)
%           fo - Initial guess for each force (row vector)
%           p - polynomial power to use
%           q - prescaler value
%
%       Outputs:
%           fOpt - Result Forces
%           mResults - result moments (should match mTarget)
%           s - value of objective function
%           infoOut - info returned from IPOPT

if nargin<6
    q=1;
end

global c
c.D=D;
c.fMax=fMax;
c.pCsa=c.fMax/25e4;
c.q=q;
c.n=p;

% Set the bounds and constraints
options.lb = zeros(1,length(fo));  % Lower bound on the variables.
options.ub = ones(1,length(fo));   % Upper bound on the variables.
options.cl = mTarget;   % Lower bounds on the constraint functions.
options.cu = mTarget;   % Upper bounds on the constraint functions.

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
funcs.jacobianstructure = @() sparse(c.D);


ao=fo./c.fMax;  %Convert forces to activation

[a infoOut] = ipopt(ao,funcs,options);  %Run IPOPT

f=a.*c.fMax;  %Convert activation to force
fOpt=f'  ;
m=constrFunc(a);
s=objFunc(a);

%---------------------------------------------------
function s=objFunc(a)
%The objective to minimize
%   Input  f: muscle forces as provided by IPOPT

global c
s=sum((c.q.*a).^c.n);


%---------------------------------------------------
function ds_df=gradObjFunc(a)
%The gradient of the objective
%   Input  f: muscle forces as provided by IPOPT

global c
ds_df=  c.n.*  ( (c.q.*a).^(c.n-1) );


%---------------------------------------------------
function m=constrFunc(a)
%The constraints
%   Input  f: muscle forces as provided by IPOPT
%   Output m: moments at the joints

global c
f=a.*c.fMax;
m=c.D*f';  %Calculate the moments

%---------------------------------------------------
function dm_df=constraintJac(a)
%The jacobian (sparse gradient) of the constraints
%   Input  f: muscle forces as provided by IPOPT
%   Output dm_df: change in moments with change in forces

global c

for i=1:size(c.D,1)
    fD(i,:)=c.D(i,:).*c.fMax;
end
dm_df=sparse(fD);  %dm_df is just the moment arm