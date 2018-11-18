function RunWarwickVLmodelFittingAndPredictions()
% RUNWARWICKVLMODELFITTINGANDPREDICTIONS MATLAB code to set input 
% parameters and perform the fitting and predictions for the Warwick VL 
% model (model W). See Supplementary File 1 for a full description of the
% model.
% 
% [NVSTAR,DSTRCTNGTVELL,TOTALNGTVELL] = FITVLMODEL(...) fits the model to 
% the district monthly numbers of VL cases in the CARE data using maximum 
% likelihood estimation to estimate the district baseline sandfly-to-human 
% ratios (SHRs) NVSTAR. DSTRCTNGTVELL is a vector of the district negative
% log-likelihoods for the fitted model and TOTALNGTVELL is the total 
% negative log-likelihood of the model. The inputs and outputs of the 
% function are saved in the MAT file 'FittingRslts.mat'. The sub-function 
% PLOTMODELFIT(...) plots the predicted monthly numbers of VL cases for 
% each district for model W for the fitted SHR (as shown in Figure S2 in 
% Supplementary File 3) and saves them in the file 'EstdMnthlyNumCases.mat'.
%
% [NVSTARCENS,CENSDSTRCTNGTVELL,TOTALCENSNGTVELL] = GGRPHCLCROSSVLDTN(...)
% performs the geographical cross-validation of the model (i.e. censoring 
% of the data for one district and estimation of the SHR for that district 
% from the data for the other districts). NVSTARCENS is a vector of the 
% estimated baseline SHRs for the censored districts, CENSDSTRCTNGTVELL is 
% a vector of the censored district negative log-likelihoods for the 
% estimated SHRs and TOTALCENSNGTVELL the overall negative log-likelihood
% of the model for the estimated SHRs. The inputs and outputs of the
% function are saved in the MAT file 'GgrphclCrossVldtnRslts.mat'. The
% predicted monthly number of cases for each censored district with the
% estimated SHRs (as shown in Figure 3 in the main text) are plotted and 
% saved (in 'EstdMnthlyNumCasesCensDstrctsGgrphclCrossVldtn.mat') by the 
% sub-function PLOTMODELFIT(...).
%
% [DSTRCTDEVS,TOTALDEV,CENSDSTRCTDEVS,TOTALCENSDEV] = CALCDEVIANCE(DATA,
% DSTRCTNGTVELL,CENSDSTRCTNGTVELL) calculates the district and overall 
% deviances DSTRCTDEVS and TOTALDEV of the fitted model (without censoring)
% from a saturated model that fits the data exactly, and the individual and 
% total deviances CENSDSTRCTDEVS and TOTALCENSDEV for the censored 
% districts in the geographical cross-validation, and saves them in
% 'ngtveLLsAndDevs.mat'.
%
% PREDICTVLINCDNCE(...) predicts and plots the district VL incidences up to
% 2020 with the continuation of current average intervention levels, and 
% incidences for sub-districts with different pre-control endemicities for
% different intervention strategies. Predicted incidences are saved in 
% 'PrdctdDstrctVLIncdncesCrrntIntvtns.mat' and 
% 'PrdctdSubdstrctVLIncdncesAltntveIntvtns.mat' respectively.

clear

%% SET PARAMETER VALUES
% Estimated district populations on 1st January 2012 from average growth 
% rate between 2001 and 2011 Indian censuses (district order: Saharsa, 
% E. Champaran, Samastipur, Gopalganj, Begusarai, Khagaria, Patna, 
% W. Champaran)
N=[1941389; 5221406; 4351328; 2603227; 3003668; 1707689; 5886089; 2842528];
% Number of districts
Nd=numel(N);

%%% Transmission parameters (all rates per day):
% Sandly biting rate
b=1/4;
% Probability of transmission from infected sandfly to susceptible human 
pH=1;  
% Rate determining mean duration of asymptomatic infection (1/gamma) 
gamma=1/150;

% Shape and scale parameters for Erlang distributions fitted to district 
% onset-to-treatment (OT) times
shape=[2;2;2;1;1;2;2;1];
scale=[16.652078781485827  17.079326237413952;
       26.714961555964166  22.027950735360275;
       19.104695873444257  17.059782774328347;
       47.491316303932365  42.751110852896730;
       21.520690102540208  22.624998488174622;
       18.300518310127735  16.620370938536908;
       22.446629145069924  21.490566045050578;
       63.640623181263237  60.600000368982123];
% Rate determining mean OT duration, calculated from shape and scale 
% parameters
delta=1./bsxfun(@times,scale,shape);
% Number of infected states (human + sandfly) 
% = # clinical VL sub-compartments + asymptomatics + treated 1 + treated 2  
%   + latently infected flies + infectious flies                                                    
n=shape+5;

% Rate determining duration of 1st and 2nd treatment for VL
tau=1/28;
% Rate determining mean duration of immunity
kappa=1/(5*365);
% District-specific death rates
death_rate=[7.6;7.8;6.7;6.3;6.2;9.3;5.0;8.7]/1000;
% Human death rate
mu=-log(1-death_rate)/365;
% Excess mortality rate due to VL
muK=1/150;
% Excess mortality rate due to 1st or 2nd treatment
muT=1/600;
% Overall mortality rates for untreated VL and VL patients undergoing 
% treatment
mu1=mu+muK;
mu2=mu+muK+muT;
% Proportion of asymptomatic individuals who develop clinical VL
f1=0.03;
% Proportion of VL patients who have 2nd treatment
f2=[102/1595; 99/1361; 127/909; 165/752; 52/483; 56/387; 27/286; 10/157];
% Relative amplitude of seasonal variation in sandfly birth rate
a1=0.3;
% Phase shift of seasonal variation in sandfly birth rate
a2=2*pi/3;
% Frequency of oscillations in sandfly birth rate (annual)
omega=2*pi/365;
% Relative infectivity of asymptomatic individuals
p1=0.025;
% Infectivity of clinical VL cases (reference value)
p2=1;
% Relative infectivity of patients undergoing 1st or 2nd treatment
p3=0.5;
% Rate determining mean duration of latent infection in sandfly
sigma=1/8;
% Rate determining mean sandfly life expectancy
muV=1/14;

%%% IRS parameters:
% District 2012 IRS coverage levels (as proportions)
cvrg=[72.46;63.16;56.68;44.42;62.09;58.14;58.70;53.66]/100;
% IRS efficacy factor
eff=0.006;

% Number of years before 2012 to run the model for to reach equilibrium
eqlbrtn_time=200;
% Time of reduction in OT
tOT=(eqlbrtn_time+1)*365;
% Time of start of IRS
tIRS=(eqlbrtn_time-1)*365;

%% SET INITIAL CONDITIONS
% Population proportions of human stages (susceptible (S), asymptomatic (A), 
% clinical VL 1 (K1), treated 1 (T1), treated 2 (T2))
PropH=[0.9,0.01,4e-5,5e-5,5e-6];
% Proportions of sandfly stages (susceptible (SV), infectious (IV))
PropV=[0.99,0.005];
% Cumulative number of VL cases (C) and proportion of clinical VL 2 (K2)
C0andPropK2=[0.001,4e-5];
% Combined proportions when there are 2 sub-compartments for clinical VL
initial=[PropH,PropV,C0andPropK2]; % initial=[S0,A0,K10,T10,T20,SV0,IV0,C0,K20] 
                                   % (see transmssnODEs(...) below)
% Combined proportions when there is 1 compartment for clinical VL
initial2=[initial(1:2),initial(3)+initial(end),initial(4:8)];

%% INPUT DATA
% Numbers of VL cases per month (columns) for each district (rows) (with 
% approximated onset dates for 37 cases with missing onset dates)
data = [128	121	102	95	81	82	72	66	73	64	44	60	73	82	93	96	74	44;
        138	107	107	86	90	57	67	58	39	43	33	58	86	91	74	42	30	12;
        57	55	87	63	50	58	60	38	20	28	22	24	49	68	53	58	50	29;
        54	51	47	56	46	38	26	26	27	24	15	17	37	44	45	44	51	21;
        31	45	29	33	32	29	24	20	15	18	20	15	25	33	33	21	25	11;
        44	29	25	22	7	18	10	14	16	19	5	14	22	22	32	17	21	20;
        18	24	21	15	17	26	11	13	9	16	10	12	10	13	16	13	12	6;
        7	8	8	5	5	10	10	8	3	3	10	6	9	5	13	7	9	5];
% Number of days in each month
monthdays=[31,28,31,30,31,30,31,31,30,31,30,31];
% District names
dstrctNames={'Saharsa','EChamparan','Samastipur','Gopalganj','Begusarai','Khagaria','Patna','WChamparan'};

%% SET SIMULATION PARAMETERS
% Times at which to evaluate number of VL cases
tint=[0,cumsum([repmat(monthdays,1,eqlbrtn_time+1),monthdays(1:6)])];
% Time span for simulation
tspan=[0,max(tint)];
% Relative tolerance for ODE solver
optionsDE=odeset('RelTol',1e-5);

%% FIT DISTRICT SANDFLY-TO-HUMAN RATIOS (SHRs) 
% Guess baseline SHRs (nVstar=nV(0)*exp(-muV*a1/omega*sin(a2)))
nVstar0=[0.4;0.4;0.35;0.275;0.4;0.4;0.4;0.35];
% Set lower and upper bounds for SHRs
lb=0.01;
ub=1;
% Set function tolerance for likelihood optimisation
options=optimoptions('fmincon','TolFun',1e-2);
% Fit SHRs
[nVstar,dstrctNgtveLL,totalNgtveLL]=FitVLModel(nVstar0,lb,ub,options,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,data,tspan,initial,initial2,optionsDE,tint,cvrg,eff,tOT,tIRS,Nd,dstrctNames);

%% PERFORM GEOGRAPHICAL CROSS VALIDATION OF MODEL
% Estimated district populations on 1st July 2012
N2012=[1963972; 5289197; 4401084; 2625987; 3039108; 1730364; 5949098; 2879280];
% District numbers of cases in 2012
Ncases2012=sum(data(:,1:12),2);
% Calculate average district incidences in 2012
dstrctIncdnces=10000*Ncases2012./N2012;

% Perform cross-validation: censor data for one district, fit SHRs for 
% other 7 districts, and use to estimate SHR for censored district. N.B. 
% since SHR is fitted separately for each district, this is equivalent to 
% fitting all district SHRs as above, censoring SHR for one district and 
% using SHRs for other 7 districts to predict censored district SHR. Here
% we do the latter for speed and convenience.
[nVstarCens,censDstrctNgtveLL,totalCensNgtveLL]=GgrphclCrossVldtn(nVstar,dstrctIncdnces,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,data,tspan,initial,initial2,optionsDE,tint,cvrg,eff,tOT,tIRS,Nd,dstrctNames);

%% CONVERT ESTIMATED BASELINE SHRs TO ANNUAL MEAN SHR
nVmean=convertBaselineSHRtoMeanSHR(nVstar,muV,a1,a2,omega);
save('meanSHR.mat','nVmean')

%% CALCULATE DEVIANCES OF MODEL WITH AND WITHOUT CENSORING
[dstrctDevs,totalDev,censDstrctDevs,totalCensDev]=calcDeviance(data,dstrctNgtveLL,censDstrctNgtveLL);

%% PREDICT FUTURE VL INCIDENCE UNDER CURRENT AND ALTERNATIVE INTERVENTIONS
% Set number of years for which to predict incidence
fut_yrs=8;
% Predict district incidences up to 2020 under current intervention levels
% (average IRS coverage = 60%, average OT time = 40 days) and incidence for 
% sub-districts with different pre-control endemicities under current and 
% alternative interventions (80% IRS coverage, 20-day OT time)
PredictVLIncdnce(fut_yrs,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,initial,initial2,cvrg,eff,eqlbrtn_time,tOT,tIRS,Nd,dstrctNames)

function [nVstar,dstrctNgtveLL,totalNgtveLL]=FitVLModel(nVstar0,lb,ub,options,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,data,tspan,initial,initial2,optionsDE,tint,cvrg,eff,tOT,tIRS,Nd,dstrctNames)

nVstar=zeros(Nd,1);
dstrctNgtveLL=zeros(Nd,1);

% Fit district SHRs
parfor i=1:Nd
% for i=1:Nd % uncomment this line and comment above line to switch to
% normal for-loop if running MATLAB version without Parallel Computing 
% Toolbox
    % Select initial conditions depending on number of clinical VL 
    % compartments
    if n(i)==7 % 2 VL sub-compartments
        initiali=initial;
    else % 1 VL compartment
        initiali=initial2;
    end
    
    % Maximise likelihood to estimate district SHRs
    try
        [nVstar(i),dstrctNgtveLL(i)]=fmincon(@(nVstar)calcNgtveLogLikelihood(nVstar,n(i),b,pH,gamma,delta(i,:),tau,kappa,mu(i),mu1(i),mu2(i),f1,f2(i),a1,a2,omega,p1,p2,p3,sigma,muV,N(i),data(i,:),tspan,initiali,optionsDE,tint,cvrg(i),eff,tOT,tIRS),nVstar0(i,:),[],[],[],[],lb,ub,[],options);
    catch % return warning if R0<=1 for initial parameter values
        warning('R0<=1 for initial values of parameters. Reset initial parameter values such that R0>=1.')
        nVstar(i)=NaN;
        dstrctNgtveLL(i)=NaN;
    end    
end

% Calculate overall negative log-likelihood for model
totalNgtveLL=sum(dstrctNgtveLL);
% Save fitting results
save FittingRslts
% Plot model fit
PlotModelFit(n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,tint,initial,initial2,optionsDE,data,cvrg,eff,tOT,tIRS,Nd,dstrctNames,'fitting')

function [nVstarCens,censDstrctNgtveLL,totalCensNgtveLL]=GgrphclCrossVldtn(nVstar,dstrctIncdnces,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,data,tspan,initial,initial2,optionsDE,tint,cvrg,eff,tOT,tIRS,Nd,dstrctNames)

nVstarCens=zeros(Nd,1);
censDstrctNgtveLL=zeros(Nd,1);

% Perform cross-validation
for i=1:Nd
    % Censor baseline SHR for district i
    nVstari=[nVstar(1:i-1);nVstar(i+1:end)];
    
    % Estimate SHR for censored district from SHRs for uncensored districts
    nVstarCens(i)=estCensDstrctSHR(nVstari,dstrctIncdnces,i,dstrctNames{i},muV,a1,a2,omega);
    if n(i)==7 % 2 VL sub-compartments
        initiali=initial;
    else % 1 VL compartment
        initiali=initial2;
    end
    % Calculate negative log-likelihood for censored district
    censDstrctNgtveLL(i)=calcNgtveLogLikelihood(nVstarCens(i),n(i),b,pH,gamma,delta(i,:),tau,kappa,mu(i),mu1(i),mu2(i),f1,f2(i),a1,a2,omega,p1,p2,p3,sigma,muV,N(i),data(i,:),tspan,initiali,optionsDE,tint,cvrg(i),eff,tOT,tIRS);  
end

% Calculate total negative log-likelihood for censored districts
totalCensNgtveLL=sum(censDstrctNgtveLL(~isnan(censDstrctNgtveLL)));
% Save cross-validation results
save GgrphclCrossVldtnRslts
% Plot model fit
PlotModelFit(n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstarCens,tint,initial,initial2,optionsDE,data,cvrg,eff,tOT,tIRS,Nd,dstrctNames,'cross-validation')

function PredictVLIncdnce(fut_yrs,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,initial,initial2,cvrg,eff,eqlbrtn_time,tOT,tIRS,Nd,dstrctNames)

%% Set simulation parameters
% Times at which to evaluate incidence
tint=1:365*(eqlbrtn_time+fut_yrs);
% Time span for simulation
tspan=[0,tint(end)];
% Relative tolerance for ODE solver
optionsDE=odeset('RelTol',1e-6);

%% Run simulations forward for districts
incdnce=zeros(Nd,numel(tint));
for i=1:Nd
    if n(i)==7 % 2 VL sub-compartments
        initiali=initial;
    else % 1 VL compartment
        initiali=initial2;
    end
    % Solve transmission ODEs
    sol=ode15s(@(t,Y)transmssnODEs(t,Y,b,pH,gamma,delta(i,:),tau,kappa,mu(i),mu1(i),mu2(i),f1,f2(i),a1,a2,omega,p1,p2,p3,sigma,muV,N(i),nVstar(i),cvrg(i),eff,tOT,tIRS),tspan,initiali,optionsDE);
    % Evaluate solution at times in tint
    Y=deval(sol,tint)';
    % Plot solutions
    [incdnce(i,:),~,~]=calcAndPlotIncdnce(tint,Y,eqlbrtn_time,fut_yrs,N(i),dstrctNames{i},true,'crrnt');
    saveas(gcf,['PrdctdDstrctVLIncdnceCrrntIntvtns' dstrctNames{i}])
end
save('PrdctdDstrctVLIncdncesCrrntIntvtns','tint','tspan','tOT','tIRS','incdnce')

%% Predict sub-district-level incidence for alternative intervention strategies
%%% Set pre-control parameter values:
% Times at which to evaluate incidence
tint1=1:365*eqlbrtn_time;
% Time span for simulation for fitting pre-control incidence
tspan1=[0,tint1(end)];
% Times of reduction in OT and start of IRS (at end of simulation time)
tOT1=eqlbrtn_time*365;
tIRS1=eqlbrtn_time*365;
% Endemicity levels
incdnces=[2,5,10];
% IRS coverage level
cvrg=0.6;
% Mean onset-treatment time
OT=40;
% Death rates (use means of district values)
mu=mean(mu);
mu1=mean(mu1);
mu2=mean(mu2);
% Proportion of VL patients who have 2nd treatment (use mean of district
% values)
f2=mean(f2);
% Sub-district population (approx. average of Bihar sub-district
% populations)
N1=2e5;

% Initial guesses for pre-control SHRs in different endemicity sub-districts
nVstar1_0=[0.28,0.3,0.32];
% Make vectors for storing pre-control SHRs
Nincdnces=numel(incdnces);
nVstar1=NaN(Nincdnces,1);
incdnceDiff=NaN(Nincdnces,1);
% Set function tolerance for fitting pre-control SHRs 
options=optimset('TolFun',1e-3);
% Fit pre-control SHRs
for i=1:Nincdnces
    % Minimise difference between incidence predicted by model and chosen 
    % pre-control incidence
    [nVstar1(i),incdnceDiff(i)]=fminsearch(@(nVstar)calcIncdnceDiff(nVstar,n(1),b,pH,gamma,1/OT,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N1,incdnces(i),tspan1,initial,optionsDE,tint1,cvrg,eff,tOT1,tIRS1),nVstar1_0(i),options);
end

%%% Set control parameter values
% IRS coverage levels
cvrgs=[0.6,0.8];
% Mean OT times
OTs=[40,40;40,20];

% Predict sub-district incidence for different IRS coverage and OT 
% reduction combinations
Ncvrgs=numel(cvrgs);
N_OTs=size(OTs,1);
incdnce1=NaN(Nincdnces*(Ncvrgs+N_OTs),numel(tint));
lgtxt=cell(Ncvrgs*N_OTs,1);
for i=1:Nincdnces
    for j=1:Ncvrgs
        for k=1:N_OTs
            % Solve transmission ODEs
            sol=ode15s(@(t,Y)transmssnODEs(t,Y,b,pH,gamma,1./OTs(k,:),tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N1,nVstar1(i),cvrgs(j),eff,tOT1,tIRS1),tspan,initial,optionsDE);
            % Evaluate solution at times in tint
            Y=deval(sol,tint)';
            % Make new figure for each pre-control endemicity level
            if j==1 & k==1
                newfig=true;
            else
                newfig=false;
            end
            % Plot solutions
            [incdnce1((Ncvrgs+N_OTs)*(i-1)+Ncvrgs*(j-1)+k,:),ax,h1(Ncvrgs*(j-1)+k,:)]=calcAndPlotIncdnce(tint,Y,eqlbrtn_time,fut_yrs,N1,['Pre-control endemicity = ' num2str(incdnces(i)) '/10,000/yr'],newfig,'altntve');
            lgtxt{Ncvrgs*(j-1)+k}=['IRS cvrge = ' num2str(100*cvrgs(j)) '%, OT = ' num2str(OTs(k,2)) ' days'];
            hold(ax,'on')
        end
    end
    legend(h1(:,1),lgtxt)
    hold(ax,'off')
    saveas(gcf,['PrdctdSubdstrctVLIncdnceAltntveIntvtnsEndmcty=' num2str(incdnces(i)) 'per10000perYr'])
end
save('PrdctdSubdstrctVLIncdncesAltntveIntvtns','tint','tspan','incdnces','cvrg','OT','tIRS1','tOT1','nVstar1','cvrgs','OTs','Ncvrgs','N_OTs','incdnce1')

function dstrctNgtveLL=calcNgtveLogLikelihood(nVstar,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,data,tspan,initial,optionsDE,tint,cvrg,eff,tOT,tIRS)

% Calculate R0 with seasonality
R0=calcSeasonalR0(n,nVstar,b,pH,gamma,delta(1),tau,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV);
if R0<=1 % i.e. no stable endemic equilibrium
    dstrctNgtveLL=NaN; % don't run simulation & set ngtve log-likelihood to NaN
else
    % Solve transmission ODEs
    sol=ode15s(@(t,Y)transmssnODEs(t,Y,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,cvrg,eff,tOT,tIRS),tspan,initial,optionsDE);
    % Evaluate solution at times in tint
    Y=deval(sol,tint)';
    % Calculate negative log-likelihood
    dstrctNgtveLL=-LogLikelihood(Y,data);
end

function dYdt=transmssnODEs(t,Y,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,cvrg,eff,tOT,tIRS)

dYdt=zeros(numel(Y),1);

% Set step-change in OT time at t=tOT
if t<=tOT
    delta=delta(1);
else
    delta=delta(2);
end

% Human stages, as proportions of the total human population
S=Y(1); % susceptible
A=Y(2); % asymptomatic
K1=Y(3); % clinical VL 1
T1=Y(4); % treated 1
T2=Y(5); % treated 2

if numel(Y)==9 % 2 VL sub-compartments
    K2=Y(end); % clinical VL 2
    % Calculate force of infection (FOI) towards flies
    lambdaV=b*(p1*A+p2*(K1+K2)+p3*(T1+T2));
    % ODE for K1
    dYdt(3)=f1*gamma*A-(2*delta+mu1)*K1;
    % ODE for K2
    dYdt(end)=2*delta*K1-(2*delta+mu1)*K2;
    % ODE for T1
    dYdt(4)=2*delta*K2-(tau+mu2)*T1;
    % Set human birth rate to match population loss due to death and excess 
    % mortality from VL and treatment
    alpha=mu+(mu1-mu)*(K1+K2)+(mu2-mu)*(T1+T2);
    % Proportion recovered
    R=1-sum(Y(1:5))-K2;
elseif numel(Y)==8 % 1 VL compartment
    % Calculate FOI towards flies
    lambdaV=b*(p1*A+p2*K1+p3*(T1+T2));
    % ODE for K1
    dYdt(3)=f1*gamma*A-(delta+mu1)*K1;
    % ODE for T1
    dYdt(4)=delta*K1-(tau+mu2)*T1;
    % Set human birth rate to match population loss due to death and excess 
    % mortality from VL and treatment
    alpha=mu+(mu1-mu)*K1+(mu2-mu)*(T1+T2);
    % Proportion recovered
    R=1-sum(Y(1:5));    
end

% Sandfly stages, as proportions of the total sandfly population
SV=Y(6); % susceptible
IV=Y(7); % infectious
EV=1-SV-IV; % latently infected

% Proportional increase in sandfly death rate due to IRS
epsilon=cvrg*eff;
% Start IRS at t=tIRS
if t<=tIRS
    % SHR before start of IRS
    nV=nVstar*exp(muV*a1*sin(omega*t+a2)/omega);
else
    % SHR after start of IRS
    nV=nVstar*exp(muV*(a1*sin(omega*t+a2)/omega-epsilon*(t-tIRS)));
end

% Calculate FOI towards humans
lambdaH=b*pH*nV*IV;

% Calculate sandfly birth rate
alphaV=muV*(1+a1*cos(omega*t+a2));

% ODEs
% Humans
dYdt(1)=alpha-(lambdaH+mu)*S+kappa*R;
dYdt(2)=lambdaH*S-(gamma+mu)*A;
dYdt(5)=f2*tau*T1-(tau+mu2)*T2;

% Sandflies
dYdt(6)=alphaV-(lambdaV+alphaV)*SV;
dYdt(7)=sigma*EV-alphaV*IV;

% Cumulative incidence of VL cases
dYdt(8)=f1*gamma*N*A;

function LL=LogLikelihood(Y,data)
% Cumulative number of VL cases
C=Y(:,8);
% Calculate monthly numbers of cases
mnthlyNumCases=C((end-17):end)-C((end-18):(end-1));
% Calculate log-likelihood assuming cases are Poisson distributed with rate
% parameter mnthlyNumCases
LL=sum(data'.*log(mnthlyNumCases)-mnthlyNumCases-log(factorial(data')));

function R0=calcSeasonalR0(n,nVstar,b,pH,gamma,delta,tau,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV)
% Calculate R0 with seasonal variation in sandfly population using Floquet 
% theory approach of Bacaer 2007 (see Sec 3.1.2 of Supplementary File 1)

% Make transmissions and transitions matrices, T and Sigma, for infection 
% sub-system (equations for infected states)
T=zeros(n,n);
if n==7 % 2 VL sub-compartments
    T(6,1:5)=b*[p1,p2,p2,p3,p3];
    Sigma=diag([gamma+mu,2*delta+mu1,2*delta+mu1,tau+mu2,tau+mu2,0,0],0);
    Sigma(3,2)=-2*delta;
    Sigma(4,3)=-2*delta;
    Sigma(5,4)=-f2*tau;
elseif n==6 % 1 VL compartment
    T(5,1:4)=b*[p1,p2,p3,p3];
    Sigma=diag([gamma+mu,delta+mu1,tau+mu2,tau+mu2,0,0],0);
    Sigma(3,2)=-delta;
    Sigma(4,3)=-f2*tau;
end
Sigma(2,1)=-f1*gamma;
Sigma(end,end-1)=-sigma;

% Initial guess for R0
lambda0=1.5;
options=optimset('Display','off');
% Calculate R0
R0=fsolve(@(lambda)calcSpctrlRdsFndmntlMtrxSoln(lambda,n,T,Sigma,nVstar,b,pH,sigma,muV,a1,a2,omega),lambda0,options);

function F=calcSpctrlRdsFndmntlMtrxSoln(lambda,n,T,Sigma,nVstar,b,pH,sigma,muV,a1,a2,omega)

% Calculate fundamental matrix solution of linearised infection subsystem
% eigenvalue problem (Eqns (25)-(26) in Sec 3.1.2 of Supplementary File 1)
tspan=[0,2*pi/omega];
psi=zeros(n,n);
for i=1:n
   X0=zeros(n,1);
   X0(i)=1;
   [~,X]=ode45(@(t,X)LnrsdInfctnSubsystem(t,X,lambda,T,Sigma,nVstar,b,pH,sigma,muV,a1,a2,omega),tspan,X0);
   psi(:,i)=X(end,:)';
end

if sum(isnan(psi))==0 & sum(isinf(psi))==0
    F=max(eig(psi))-1;
else
    F=1e10;
end

function dXdt=LnrsdInfctnSubsystem(t,X,lambda,T,Sigma,nVstar,b,pH,sigma,muV,a1,a2,omega)
% Define linearised infection sub-system
nV=nVstar*exp(muV*a1*sin(omega*t+a2)/omega);
alphaV=muV*(1+a1*cos(omega*t+a2));
T(1,end)=b*pH*nV;
Sigma(end-1,end-1)=sigma+alphaV;
Sigma(end,end)=alphaV;
dXdt=(T/lambda-Sigma)*X;

function nVstarCens=estCensDstrctSHR(nVstar,dstrctIncdnces,i,dstrctName,muV,a1,a2,omega)
% Make vector of 2012 average district incidences excluding incidence for
% censored district
dstrctIncdncesi=[dstrctIncdnces(1:i-1);dstrctIncdnces(i+1:end)];
% Perform least-squares linear regression of district incidences on SHRs
X=[ones(numel(dstrctIncdncesi),1) dstrctIncdncesi];
b=X\nVstar;
nVstarCens=[1 dstrctIncdnces(i)]*b;
nVstar_est=X*b;
rsqrd=1-sum((nVstar-nVstar_est).^2)/sum((nVstar-mean(nVstar)).^2);
% Convert SHRs to mean annual SHRs
nVmean=convertBaselineSHRtoMeanSHR(nVstar,muV,a1,a2,omega);
nVmeanCens=convertBaselineSHRtoMeanSHR(nVstarCens,muV,a1,a2,omega);
nVmean_est=convertBaselineSHRtoMeanSHR(nVstar_est,muV,a1,a2,omega);
% Plot linear regression and estimated censored district SHR
figure; plot(nVmean,dstrctIncdncesi,'bx',nVmean_est,dstrctIncdncesi,'k-',nVmeanCens,dstrctIncdnces(i),'ro');
xlabel('Average SHR'); ylabel('2012 average VL incidence (per 10,000/yr)'); legend('Fitted SHRs',['r^2 = ' num2str(rsqrd)],'Estd censrd SHR','Location','NorthWest')
saveas(gcf,['DstrctIncdncesVsSHRsGgrphclCrossVldtnCensDstrct' dstrctName])

function nVmean=convertBaselineSHRtoMeanSHR(nVstar,muV,a1,a2,omega)
nVmean=nVstar*integral(@(t)exp(muV*a1/omega*sin(omega*t+a2)),0,365)/365;

function incdnceDiff=calcIncdnceDiff(nVstar,n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,incdnce,tspan,initial,optionsDE,tint,cvrg,eff,tOT,tIRS)

% Calculate R0 with seasonality
R0=calcSeasonalR0(n,nVstar,b,pH,gamma,delta,tau,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV);
if R0<=1 % no endemic equilibrium
    incdnceDiff=NaN; % don't run simulation & set error in incidence to NaN
else
    % Solve transmission ODEs
    sol=ode15s(@(t,Y)transmssnODEs(t,Y,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,cvrg,eff,tOT,tIRS),tspan,initial,optionsDE);
    % Evaluate solution at times in tint
    Y=deval(sol,tint)';
    % Calculate annual incidence from model and difference from required incidence
    C=Y(:,8);
    AnnlIncdnce=1e4*(C(end)-C(end-365))/N;
    incdnceDiff=abs(AnnlIncdnce-incdnce);
end

function [incdnce,ax,h]=calcAndPlotIncdnce(t,Y,eqlbrtn_time,fut_yrs,N,dstrctName,newfig,str)

if newfig
    figure;
    ax=gca;
else
    ax=gca;
end
C=Y(:,8);
% Calculate VL incidence (per 10,000/yr)
incdnce=1e4*365*(C-[0;C(1:end-1)])'/N;
% Make index vector for times to plot incidence (including 2 years prior to
% start of interventions)
Nt=numel(t);
indx=(Nt-(fut_yrs+2)*365+1):Nt;

if strcmp(str,'crrnt') % plot predicted district incidences
h=plot(2012+t(indx)/365-eqlbrtn_time,incdnce(indx),2012+t(indx)/365-eqlbrtn_time,ones(1,numel(t(indx))),'r--');
xlabel('Year'); ylabel('VL Incidence (per 10,000/yr)'); title(dstrctName)
elseif strcmp(str,'altntve') % plot predicted incidence at sub-district level for alternative interventions
h=plot(t(indx)/365-eqlbrtn_time,incdnce(indx),t(indx)/365-eqlbrtn_time,ones(1,numel(t(indx))),'r--');
xlabel('Years after start of intervention'); ylabel('VL Incidence (per 10,000/yr)'); title(dstrctName)
end

pause(0.001)

function [dstrctDevs,totalDev,censDstrctDevs,totalCensDev]=calcDeviance(data,dstrctNgtveLL,censDstrctNgtveLL)

% Calculate negative log-likelihood of saturated model
dstrctNgtveLLsat=SaturatedNgtveLogLikelihood(data);
% Calculate deviance of model for each district
dstrctDevs=Deviance(dstrctNgtveLL,dstrctNgtveLLsat);
% Calulate total deviance
totalDev=sum(dstrctDevs);

% Calculate deviance of model for each censored district in
% cross-validation
censDstrctDevs=Deviance(censDstrctNgtveLL(~isnan(censDstrctNgtveLL)),dstrctNgtveLLsat(~isnan(censDstrctNgtveLL)));
% Calculate total model deviance with censoring
totalCensDev=sum(censDstrctDevs);
save('ngtveLLsAndDevs.mat','dstrctNgtveLL','dstrctNgtveLLsat','dstrctDevs','totalDev','censDstrctNgtveLL','censDstrctDevs','totalCensDev')

function dev=Deviance(ngtveLL,ngtveLLsat)
dev=2*(ngtveLL-ngtveLLsat);

function ngtveLLsat=SaturatedNgtveLogLikelihood(data)
ngtveLLsat=-sum(data.*log(data)-data-log(factorial(data)),2);

function PlotModelFit(n,b,pH,gamma,delta,tau,kappa,mu,mu1,mu2,f1,f2,a1,a2,omega,p1,p2,p3,sigma,muV,N,nVstar,tint,initial,initial2,optionsDE,data,cvrg,eff,tOT,tIRS,Nd,dstrctNames,str)
mnthlyNumCases=zeros(Nd,18);
% Plot model fit for each district
for i=1:Nd
    if n(i)==7 % 2 VL sub-compartments
        initiali=initial;
    else % 1 VL compartment
        initiali=initial2;
    end
    % Solve transmission ODEs for fitted SHR for district i
    [~,Y]=ode15s(@(t,Y)transmssnODEs(t,Y,b,pH,gamma,delta(i,:),tau,kappa,mu(i),mu1(i),mu2(i),f1,f2(i),a1,a2,omega,p1,p2,p3,sigma,muV,N(i),nVstar(i),cvrg(i),eff,tOT,tIRS),tint,initiali,optionsDE);
    % Cumulative number of VL cases
    C=Y(:,8);
    % Monthly numbers of cases 
    mnthlyNumCases(i,:)=(C((end-17):end)-C((end-18):(end-1)))';
    % Plot monthly numbers of cases
    figure;
    mnths=1:18;
    plot(mnths,mnthlyNumCases(i,:),mnths,data(i,:))
    xlim([1 18])
    ax=gca;
    ax.XTick=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
    ax.XTickLabel={'Jan-12','Feb-12','Mar-12','Apr-12','May-12','Jun-12','Jul-12','Aug-12','Sep-12','Oct-12','Nov-12','Dec-12','Jan-13','Feb-13','Mar-13','Apr-13','May-13','Jun-13'};
    ax.XTickLabelRotation=45;
    xlabel('Month'); ylabel('No. of cases'); title(dstrctNames{i})
    legend('Model','Data')
    if strcmp(str,'fitting')
        saveas(gcf,['VLModelFit' dstrctNames{i}])
    elseif strcmp(str,'cross-validation')
        saveas(gcf,['VLModelFitGgrphclCrossVldtn' dstrctNames{i}])
    end
end
% Save estimated monthly numbers of cases
if strcmp(str,'fitting')
    save('EstdMnthlyNumCases','mnthlyNumCases')
elseif strcmp(str,'cross-validation')
    save('EstdMnthlyNumCasesCensDstrctsGgrphclCrossVldtn','mnthlyNumCases')
end
    