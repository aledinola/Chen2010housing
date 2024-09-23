% Chen (2010) - A life-cycle analysis of social security with housing
%
% Two endogenous states: a (financial assets), h (housing)
% One exogenous state: z (stochastic earnings)

% I renamed two things: 
%   Chen calls the stochastic process on earnings eta, I call it z.
%   Chen calls the deterministic earnings as function of age epsilonj, I call it kappaj

clear;clc;close all
addpath(genpath('C:\Users\aledi\Documents\GitHub\VFIToolkit-matlab\VFIToolkit-matlab'))
%addpath(genpath('C:\Users\aledi\OneDrive\Documents\GitHub\VFIToolkit-matlab'))

%% 
n_d=0; % share of time to invest in new human capital
%n_a=[201,31]; % Assets, Housing
n_a=[301,11]; % Assets, Housing
n_z=7; % Chen (2010) uses 7 points and Tauchen. 

N_j=66; % total number of periods
Params.agejshifter=19; % starting age 20
% So model is annual for ages 20 to 85
Params.J=N_j;

%% Parameters

% Discount factor
Params.beta=0.9622;

% Preferences
Params.sigma=2; % curvature of utility
Params.theta=0.8954; % share of nondurable consumption in the utility function

% Demographics
Params.Jr=46; % retire at age 65
Params.agej=1:1:Params.J; % agej is the period
Params.n=0.01; % Population growth rate

% Deterministic productivity growth
Params.g=0.015;

% Production function
Params.alpha=0.2743; % share of capital in Cobb-Douglas production fn
% Depreciation rates
Params.delta_r=0.0254; % annual depreciation rate for rental housing
Params.delta_o=0.013; % annual depreciation rate for owner-occupied housing
Params.delta_k=0.0951; % annual depreciation rate for physical capital (financial assets)

% Housing
Params.phi=0.05; % transaction cost for changing amount of housing
Params.gamma=0.2; % downpayment ratio (when buying a house)

% Earnings
Params.rho_z=0.96; % autocorrelation coeff of AR(1) on earnings
Params.sigma_z_epsilon=0.045; % standard deviation of innovations to the AR(1) on earnings
Params.sigma_z1=0.38; % standard deviation of initial earnings (at age 20, period 1)
% Note: asymptotic std dev of z is sqrt((0.045^2)/(1-0.96^2))=0.1607
% Which means that the std dev of earnings is going to decrease over age.
% Seems weird.
% Actually, looking more closely, sigma_z1 (Chen calls it sigma_eta1) is
% described as the "variance of labor-income shocks at initial age" with
% no mention of 'log'.

% Pension
Params.vartheta=0.483; % Not used for anything (Chen2010 reports this a the 'replacement rate', but is not used here)

% Taxes
Params.tau_p=0.107;

%First the raw data of Hansen 1993:
%Table II. Weights assigned to age-sex
%         Males       Females
%Age      Weight      Weight
%16-19    0.56        0.52
%20-24    0.78        0.69
%25-34    1.14        0.89
%35-44    1.37        0.90
%45-54    1.39        0.87
%55-64    1.33        0.84
%65 +     0.89        0.66
%14-17    0.56        0.52
%14-19    0.56        0.52
%18-24    0.78        0.69
%25-44    1.24        0.89
%45-64    1.37        0.86
% kappaj is set based on the first entries of the data for males
% Rather than use this, Chen (2010) codes "./benchmark/demographics.f90"
% contains "Age-efficiency Units from Hansen (1991)" which I copy here
Params.kappaj=[1.0000, 1.0719, 1.1438, 1.2158, 1.2842, 1.3527, 1.4212, 1.4897, 1.5582, 1.6267, 1.6952, 1.7217, 1.7438, 1.7748, 1.8014, 1.8279, 1.8545, 1.8810, 1.9075, 1.9341, 1.9606, 1.9623, 1.9640, 1.9658, 1.9675, 1.9692, 1.9709, 1.9726, 1.9743, 1.9760, 1.9777, 1.9700, 1.9623, 1.9546, 1.9469, 1.9392, 1.9315, 1.9238, 1.9161, 1.9084, 1.9007, 1.8354, 1.7701, 1.7048, 1.6396];
%[0.918005352232869, 0.928276694799422,0.937322109949312, 0.945265376137042, 0.952226660882508, 0.958318429643330, 0.963643605752677, 0.968295020073741, 0.972355583361627, 0.975898844621107, 0.978989735486561, 0.981685382871588, 0.984035921985104, 0.986085272134275, 0.987871856136990, 0.989429255229458, 0.990786797877776, 0.991970084616345, 0.993001453018585, 0.993900387831011, 0.994683881593424, 0.995366750990553, 0.995961913899121, 0.996480631710744, 0.996932721087412, 0.997326738879082, 0.997670143523069, 0.997969435863394, 0.998230281979946, 0.998457620303540, 0.998655755012837, 0.998828437460621, 0.998978937157669, 0.999110103649428, 0.999224420451403, 0.999324052060782, 0.999410884932041, 0.999486563190821, 0.999552519761282, 0.999610003495669, 0.999660102819339, 0.999703766338693, 0.999741820802016, 0.999774986753218, 0.999803892174763];
Params.kappaj=[Params.kappaj,zeros(1,Params.J-Params.Jr+1)]; % zeros for retirement

% Survival probabilities
% F. C. Bell and M. L. Miller (2005), Life Tables for the United States Social Security Area 1900-2100, Actuarial Study No. 120, Office of the Chief Actuary
% http://www.socialsecurity.gov/oact/NOTES/s2000s.html
% Table 8 — Period Probabilities of Death Within One Year (qx) at Selected Exact Ages, by Sex and Calendar Year (Cont.)
% The raw data from there is
%          Sex and Exact Age
%     |  Male                                                Female
% Year| [0 30 60 65 70 100]                                  [0 30 60 65 70 100]
% 2010| [0.00587,0.00116,0.01086,0.01753,0.02785,0.39134]    [0.00495,0.00060,0.00734,0.01201,0.01912,0.34031]
% I just take the numbers for Males, and then set my actual values based on a linear interpolation of the data.
dj_temp=interp1([0,30,60,65,70,100],[0.00587,0.00116,0.01086,0.01753,0.02785,0.39134],0:1:100,'linear');
Params.sj=1-dj_temp(20:85);
Params.sj(1)=1;
Params.sj(end)=0;
% I use linear interpolation to fill in the numbers inbetween those reported by Bell & Miller (2005).
% I have additionally imposed that the prob of death at age 20 be zero and that prob of death at age 85 is one.
% Chen (2010) replication codes contains a "./benchmark/surv.txt" contains
% the survival probabilities he used, they are quite different, so I now provide those
Params.sj=[0.9985995, 0.9985575, 0.998498, 0.9985075, 0.9986075, 0.998641, 0.998743, 0.9987235, 0.9987385, 0.998705, 0.9986825, 0.9986935, 0.9986795, 0.9986615, 0.998573, 0.9984555, 0.998305, 0.9981125, 0.9979675, 0.997858, 0.9977755, 0.997652, 0.997384, 0.99714, 0.9969315, 0.9967195, 0.99666, 0.996468, 0.9962415, 0.9959035, 0.995637, 0.9953435, 0.9950225, 0.9946615, 0.9942525, 0.9937795, 0.993245, 0.9926615, 0.992031, 0.9913415, 0.9906035, 0.9897745, 0.9887855, 0.987602, 0.9862635, 0.9847635, 0.9832185, 0.981763, 0.9804655, 0.979231, 0.977836, 0.976205, 0.974443, 0.9725355, 0.9704165, 0.967929, 0.9650565, 0.9619165, 0.9585285, 0.954784, 0.9506335, 0.945818, 0.940017, 0.933017, 0.9249375, 0.9160355];
Params.sj(end)=0; % I need this for the accidental bequests calculations. Not sure how Chen (2010) treated it.

%% Grids and exogenous shocks

%% Grid for housing
% Chen used 20 as maximum, 0 as minimum. Here 0 means that you are a renter
minh   = 0;
maxh   = 20; % maximum value of housing: 
h_curv = 3.0;
h_grid = minh + (maxh-minh)*linspace(0,1,n_a(2)).^h_curv';

%% Grid for financial assets
maxassets=100; % maximum value of assets. Chen used 100
% we want possibility of negative assets
a_curv = 3.0;

% First, the negative part of the grid: equally spaced from 
% -(1-Params.gamma)*maxh up to 0. Note that the max value of borrowing ever 
% possible is (1-Params.gamma)*maxh
n_a_neg = round(0.1*n_a(1));
a_grid_neg = linspace(-(1-Params.gamma)*maxh,0,n_a_neg)';

% Second, the positive part of the grid, NOT equally spaced but denser the
% closer it is to 0.
n_a_pos = n_a(1)-n_a_neg+1;
a_grid_pos = 0+(maxassets-0)*linspace(0,1,n_a_pos).^a_curv';
asset_grid = [a_grid_neg(1:n_a_neg-1);a_grid_pos];
% Note: -(1-Params.gamma)*maxh is the minimum possible value of assets

if ~isequal(length(asset_grid),n_a(1))
    error('asset_grid must have exactly n_a(1) points')
end

%% Combined grids for endogenous state variables
a_grid = [asset_grid; h_grid]; % stacked column vector
% size(a_grid) = [sum(n_a),1]
% kfttoptions.initialj1sigmaz=Params.sigma_z1;
% kfttoptions.nSigmas=2; % plus/minus two std deviations as the max/min grid points
% [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs] = discretizeLifeCycleAR1_KFTT(0,Params.rho_z,Params.sigma_z_epsilon,n_z,N_j,kfttoptions);
% z_grid_J=exp(z_grid_J);

% Just for comparison, use Tauchen to see what the grid Chen (2010) would have had looks like
tauchenoptions=struct();
Tauchen_q=sqrt(7/2); % plus/minus two std deviations as the max/min grid points % DONT KNOW WHAT Chen2010 actually used
[z_grid_log,pi_z]=discretizeAR1_Tauchen(0,Params.rho_z,sqrt(Params.sigma_z_epsilon),n_z,Tauchen_q,tauchenoptions);
z_grid=exp(z_grid_log); % Ranges from 0.7 to 1.4
% Stationary distribution of pi_z
aux = pi_z^1000;
prob_stat = aux(1,:)';
z_mean = dot(z_grid,prob_stat);
z_grid = z_grid/z_mean;

%% Initial distribution of earnings "z" for age 1
jequaloneDistz = discretize_normal(z_grid_log,sqrt(Params.sigma_z1));

% -- As a check, compute the variance of log earnings over age

% Step 1: Compute the distribution of z for each age, starting from the
% distribution of log(z1). This is NOT the stationary distribution, of
% course
pilab = zeros(n_z,Params.Jr-1);
pilab(:,1) = jequaloneDistz;
for j=1:Params.Jr-1
    pilab(:,j+1)=pi_z'*pilab(:,j);
end

% Step2: Compute mean and standard deviation of log(z) by age
var_log_z = zeros(Params.Jr-1,1);
for j=1:Params.Jr-1
    log_z_mean = dot(z_grid_log,pilab(:,j));
    var_log_z(j) = dot((z_grid_log-log_z_mean).^2,pilab(:,j));
end

figure
plot(1:Params.Jr-1,var_log_z)
xlabel('Age, j')
title('Variance of log(z_t)')
print('var_logs','-dpng')

%% Age distribution
AgeWeightParamNames={'mewj'};
Params.mewj=ones(1,N_j);
for jj=1:N_j-1
    Params.mewj(jj+1)=Params.mewj(jj)*Params.sj(jj)/(1+Params.n); % iteratively apply formula that mewj is determined by conditional survival probabilities (sj) and population growth rate (n)
end
Params.mewj=Params.mewj/sum(Params.mewj); % Normalize total mass to one
fracret = sum(Params.mewj(Params.Jr:end));

%% Aggregate labor in efficiency units, N_agg
N_agg = 0;
for jj=1:Params.Jr-1
    N_agg_jj = sum(z_grid.*pilab(:,jj));
    N_agg   = N_agg + Params.kappaj(jj)*N_agg_jj*Params.mewj(jj);
end

d_grid=[]; % no d variable
jequaloneDistz = gpuArray(jequaloneDistz);

%% Return fn
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr)...
    Chen2010_ReturnFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,sigma,gamma,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);

%% General eqm parameters
Params.r=0.06354260; % This is an initial guess, but it is also the value that Chen (2010) reports as being the initial stationary general eqm
Params.Tr=0.04988135; % accidental bequests of assets and housing

% Three other parameters are being determined in general equilibrium, but
% we can shortcut for these. They are b, p and w.

% p and w are both determined in general eqm. But since both are just
% simple functions of other general eqm parameters, we can just hard-code
% them (inside return function and elsewhere).
Params.p=(Params.r+Params.delta_r)/(1+Params.r); % price of housing
Params.w=(1-Params.alpha)*((Params.r+Params.delta_k)/Params.alpha)^(Params.alpha/(Params.alpha-1)); % wage
% Because labor supply is exogenous you can actually figure out of b is without having to
% solve general eqm (point 4 in Definition 1 of his Appendix A.1).
%Params.b=0.1;

% price of housing
p=(Params.r+Params.delta_r)/(1+Params.r);
% wage
w=(1-Params.alpha)*((Params.r+Params.delta_k)/Params.alpha)^(Params.alpha/(Params.alpha-1));

% Note that the pension balance condition is
% b*fracret = tau_p*w*N
% where fracret, N are exogenous
Params.b = Params.tau_p*Params.w*N_agg/fracret;

%% Solve the value function
vfoptions.divideandconquer=1;
vfoptions.level1n=[9,n_a(2)];
vfoptions.verbose=1;
tic
[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
time_vfi=toc;
% [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);

% vfoptions.divideandconquer=0;
% [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
% 
% max(abs(V(:)-V2(:)))
% max(abs(Policy(:)-Policy2(:)))

%% Initial age 20 (period 1) distribution
jequaloneDist=zeros([n_a,n_z],'gpuArray');
% Everyone is born with zero assets and zero housing
[~,zeroassetindex]=min(abs(asset_grid));
jequaloneDist(zeroassetindex,1,:)=shiftdim(jequaloneDistz,-2); % initial dist of z, with zero assets and zero housing

%% Agent distribution
simoptions=struct(); % defaults
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% AggVars
FnsToEvaluate.A=@(aprime,hprime,a,h,z) a;
FnsToEvaluate.N=@(aprime,hprime,a,h,z,kappaj) kappaj*z;
FnsToEvaluate.H=@(aprime,hprime,a,h,z) h; % housing
FnsToEvaluate.Hr=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b) (hprime==0)*Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b); % rental housing=housing services used by renters
FnsToEvaluate.pensiontaxrevenue=@(aprime,hprime,a,h,z,tau_p,w,kappaj) tau_p*w*kappaj*z;
FnsToEvaluate.pensionspend=@(aprime,hprime,a,h,z,agej,Jr,b) (agej>=Jr)*b;
FnsToEvaluate.accidentalbeqleft=@(aprime,hprime,a,h,z,r,sj,delta_o) (1+r)*aprime*(1-sj)+(1-delta_o)*hprime*(1-sj);

FnsToEvaluate.Homeownership=@(aprime,hprime,a,h,z) (hprime>0);
FnsToEvaluate.TotalWealth=@(aprime,hprime,a,h,z) a+h; % NOT SURE IF THIS IS CORRECT DEFINITION OF TOTAL WEALTH
FnsToEvaluate.LoanToValueRatio=@(aprime,hprime,a,h,z) (hprime>0)*abs(aprime/hprime);
FnsToEvaluate.earnings=@(aprime,hprime,a,h,z,w,kappaj) w*kappaj*z;
FnsToEvaluate.consumption=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr) Chen2010_ConsumptionFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b,Tr);
FnsToEvaluate.housingservices=@(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b) Chen2010_HousingServicesFn(aprime,hprime,a,h,z,kappaj,r,tau_p,theta,phi,alpha,delta_k,delta_o,delta_r,agej,Jr,b);


% [V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid_J, pi_z_J, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
% StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptions);
% AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,[],simoptions);
% AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid_J,simoptions);
%[V,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid, a_grid, z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
%StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

tic
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,[],simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy, FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
time_stats=toc;

K=AggVars.A.Mean-AggVars.Hr.Mean*(1-Params.p); % physical capital
Y=(K^Params.alpha)*(AggVars.N.Mean^(1-Params.alpha)); % output

% Chen2010, bottom of pg 603 says that the Gini coeff for labor income is 0.40
tic
simoptions.agegroupings=[1,Params.Jr]; % working age and retirees
AgeConditionalStats=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AgeConditionalStats.earnings.Gini(1) % Substantially below 0.4
AllStats.earnings.Gini % 0.3, so still not 0.4

% Quick look at the replacement rate
Params.vartheta % 0.48
Params.b/AgeConditionalStats.earnings.Mean(1) % 0.46, so is correct

% Plot some life-cycle profiles to see more about what is going on
simoptions=struct(); % back to defaults
simoptions.whichstats(1) = 1; % only the mean
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);

time_life=toc;

figure(1)
subplot(3,1,1); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.earnings.Mean)
title('Earnings (age conditional mean of)')
subplot(3,1,2); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.A.Mean)
title('Financial Wealth (age conditional mean of)')
subplot(3,1,3); plot(Params.agejshifter+Params.agej,AgeConditionalStats2.H.Mean)
title('Housing (age conditional mean of)')

figure(2)
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.consumption.Mean)
hold on
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.housingservices.Mean)
plot(Params.agejshifter+Params.agej,AgeConditionalStats2.earnings.Mean)
hold off
legend('consumption','housing services', 'earnings')
title('Consumption, Housing Services and Earnings')

Homeownership_5y = mean_k_width(AgeConditionalStats2.Homeownership.Mean,5);
age_5y = Params.agejshifter+(1:5:N_j)';
figure
plot(age_5y(1:end-1),Homeownership_5y,'LineWidth',2)
xlabel('Age')
title('Homeownership rate (conditional on age)')
print('own','-dpng')

%% Make table

% Gini for financial wealth
% StationaryDist(a,h,z,j). Compute the marginal distribution of assets "a"
% and call it p
StationaryDist_a = sum(gather(StationaryDist),[2,3,4]);
[~, ~, ~, giniw] = lrzcurve(StationaryDist_a, asset_grid);
GiniA2 = giniw;

% Chen (2010) Table 2 on page 605 reports the following
fprintf('Quantitative properties of the benchmark economy \n')
% Comments at end of each line are values in the paper for comparison
fprintf('Targeted Variables \n')
fprintf('Payroll tax rate is %1.3f \n',Params.tau_p) % 0.107
fprintf('r is %1.2f%% \n',100*Params.r) % 6.35%
fprintf('K/Y is %1.3f \n',K/Y) % 1.729
fprintf('H/Y is %1.3f \n',AggVars.H.Mean/Y) % 1.075
fprintf('Homeownership rate is %2.1f%% \n', 100*AggVars.Homeownership.Mean) % 65.0%
fprintf('Nontargeted Variables \n')
% fprintf('pH/(C+pH) is %2.1f \% \n',100*Params.p*AggVars.H.Mean/()) % 10.7% I SKIPPED AS CANT BE BOTHERED CREATING CONSUMPTION FN
fprintf('Ho/(A+Ho) is %2.1f%% \n',100*AggVars.H.Mean/(AggVars.A.Mean+AggVars.H.Mean)) % 32.2%  % I think this is the correct calculation, not sure
fprintf('Gini for total wealth is %1.3f \n',AllStats.TotalWealth.Gini) % 0.73
fprintf('Gini for financial wealth is %1.3f \n',GiniA2) % 0.93
fprintf('Gini for housing is %1.3f \n',AllStats.H.Gini) % 0.52 % NOT SURE HOW RENTAL HOUSING SERVICES ARE TREATED HERE, GUESSING JUST AS ZEROS?
fprintf('Mean loan-to-value ratio (for borrowers) is %2.1f \n',100*AggVars.LoanToValueRatio.Mean/AggVars.Homeownership.Mean) % 0.93

disp('RUNNING TIMES')
fprintf('Time VFI: %f \n',time_vfi)
fprintf('Time aggvars and stats: %f \n',time_stats)
fprintf('Time life profiles: %f \n',time_life)





