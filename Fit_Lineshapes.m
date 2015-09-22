# -*-Octave-*-

function [] = Fit_Lineshapes(MUT="WT", KD=76e-6, kiso=10.0, freq_bound=[7865.5,7450.3], Output_File="temp.out",\
              Cyp_Concentration=[0,5,10,20,5,10,20,50,100], FGP_Concentration=[1000,500,500,500,1000,1000,1000,1000,1000],\
              iterations=200, acquisition_time=1, acquisition_points=14045, zero_fill=2^16)

%An octave function, will fit lineshape data, along with dissociation constant, isomerization rate, and bound chemical shift values to determine microscopic rate constants kab, kba, kbc, kcb, kcd, and kdc
%Data are 2D extracted hsqcs comprising lineshapes that will be summed over nitrogen to generate a 1D proton lineshape to fit. Data must be stored in:
%  ./dat/MUT/FGPconcentration_CYPconcentration/hsqc_RES#_RESconformation.txt
%  e.g.: ./data/WT/1000_50/hsqc_6_trans.txt
%  Data folder must contain data for residues 6 and 7, cis and trans, and also contain an ‘hsqc_11_ref.txt’ file that contains reference data for GSW peptide reference peak
%  Output file is in Octave format, with parameters as listed below
%
%  function [] = Fit_Lineshapes(MUT, KD, kiso, freq_bound, F_Out, Cyp_Conc, FGP_Conc, iterations, acquisition_time, acquisition_points, zero_fill)
%
%  MUT = string protein/mutant name (default, “WT”)
%  KD = Dissociation constant to fit (M), if set to 0, will not use in fit
%  kiso = isomerization rate, at 1 mM peptide, 20 uM protein (s^-1), if set to 0, will not use in fit
%  freq_bound = vector of proton frequencies, in Hz, of bound residues, if set to 0, will not use in fit
%  F_out = output file (Default “temp.out”)
%  Cyp_Conc = Cyp concentrations to use, in uM, first value must be free peptide, with Cyp_Conc==0 (default = [0, 5, 10, 20, 5, 10, 20, 50, 100])
%  FGP_Conc = Peptide concentrations to use, in uM (default = [1000, 500, 500, 500, 1000, 1000, 1000, 1000, 1000])
%  iterations = number of fits to perform (default = 200)
%  acquisition_time = collected FID time (in s, default 1)
%  acquisition_points = number of points acquired in FID (default 14045)
%  zero_fill = if zero filling is applied, number of points (default 2^16)
%
%  Output File Parameters (N = number of iterations)
%  PARS_OUT = fit parameters, a 10xN matrix of parameters, with each column corresponding to [kab; kba; kbc; kcb; kcd; kdc; w_6_bound_trans, w_6_bound_cis, w_7_bound_trans, w_7_bound_cis]
%  r2_OUT = r^2 values for each fit
%  f_OUT = XxN matrix of fits, where X is the number of data points fit
%  full_data_vector = linearized 1D proton data set to which fit is performed
%  data_freq_linear = frequencies corresponding to each data point in full_data_vector
%  Populations_OUT = 5xN matrix of the populations occupied given 1 mM peptide and 20 uM Cyp, given the fit parameters, with each column corresponding to [free_trans; bound_trans; bound_cis; free_cis; free_Cyp]


global full_data_vector GLOBAL

more off %Prevent hang-ups on paging

% Load all data
GLOBAL.RES_to_USE=[6,6,7,7];
GLOBAL.CONFORMATION={"trans","cis","trans","cis"};

%Constraints to use flags
if(freq_bound(1) != 0) GLOBAL.constrain_w=1; end
GLOBAL.w_scale=1e-3; %scale bound w's so they're in the ~10 range
if(KD != 0) GLOBAL.constrain_Kd=1; end
GLOBAL.Kd_scale=1e5; %scale Kd's so they're in the ~10 range
if(kiso != 0) GLOBAL.constrain_ZZ=1; end %No scaling needed for ZZ - already in the ~10 range
GLOBAL.constrain_lineshape=1;
GLOBAL.R20bound=40.0*ones(size(GLOBAL.RES_to_USE));

%Normalize data flag
GLOBAL.norm_factor=1000; %Normalize area under each condition to this value to keep all parameters in the same range
GLOBAL.normalize=1; %Normalize data from each experiment to the same area under the curve

GLOBAL.MUT=MUT;
GLOBAL.Cyp_Conc=Cyp_Concentration;
GLOBAL.FGP_Conc=FGP_Concentration;

%Initialize Output file;
fout=strcat(Output_File);

%Parameters
Kd_match=KD;
ZZ_match=kiso;

data_freq_linear=[];

%Determine time pts
np=acquisition_points;
if(zero_fill>np) np=zero_fill; end
GLOBAL.t_app=acquisition_time*np./acquisition_points;
GLOBAL.t_step=t=(0:GLOBAL.t_app./np:GLOBAL.t_app-GLOBAL.t_app./np);
GLOBAL.full_freq=[0:1./t(2)./size(t,2):1./t(2)-1./t(2)./size(t,2)];

% Load/Process Data

%Load reference data
use_ref="GSW";
%use_ref="10";

if(strcmp(use_ref,"10"))
  %load and process data for reference residue #10 to deterine reference peak shifts:
  reference_shift=zeros(size(GLOBAL.Cyp_Conc,2),2);
  for n=1:size(GLOBAL.Cyp_Conc,2)
    dat_temp=load(strcat("dat/",GLOBAL.MUT,"/",int2str(GLOBAL.FGP_Conc(n)),"_",int2str(GLOBAL.Cyp_Conc(n)),"/hsqc_","10_both.txt")); %load data
    nptsH=(find(dat_temp(:,2)!=dat_temp(1,2))-1)(1); %find # H points
    nptsN=size(dat_temp,1)./nptsH; %find # N points
    dat_temp_lineshape=sum(reshape(dat_temp(:,3),nptsH,nptsN),2)'; %sum over N
    dat_temp_freq=dat_temp(1:nptsH); 
    reference_shift(n,:)=[find(dat_temp_lineshape==max(dat_temp_lineshape(1:floor(nptsH./2)))), \
    find(dat_temp_lineshape==max(dat_temp_lineshape(floor(nptsH./2):end)))]; %Position of peak 10, as the mean of the two N-split peaks
  end
  peak_shift=mean(reference_shift.-repmat(reference_shift(1,:),size(reference_shift,1),1),2)*   \
  (max(dat_temp_freq)-min(dat_temp_freq))/(size(dat_temp_freq,2)-1); %shift to apply, with the free peptide spectrum as the reference
end

if(strcmp(use_ref,"GSW"))
  "GSW";
  %Load GSW reference peak (labeled as residue 11) for referencing
  get_ref=1;
  if(get_ref==1)
    reference_shift=zeros(size(GLOBAL.Cyp_Conc,2),2);
    ic_freq=[6954,7322];
    ref_res=[11,12];
    for m=1
      for n=1:size(GLOBAL.Cyp_Conc,2)
        dat_temp=load(strcat("dat/",GLOBAL.MUT,"/",int2str(GLOBAL.FGP_Conc(n)),"_",int2str(GLOBAL.Cyp_Conc(n)),"/hsqc_",int2str(ref_res(m)),"_ref.txt")); %load data
        nptsH=(find(dat_temp(:,2)!=dat_temp(1,2))-1)(1); %find # H points
        nptsN=size(dat_temp,1)./nptsH; %find # N points
        dat_temp_lineshape=sum(reshape(dat_temp(:,3),nptsH,nptsN),2)'; %sum over N
        dat_temp_freq=dat_temp(1:nptsH); 

        %Define frequencies
        freq_range_ref=[find(GLOBAL.full_freq>=min(dat_temp_freq))(1),find(GLOBAL.full_freq<=max(dat_temp_freq))(end)]; 
        ref_freq=GLOBAL.full_freq(freq_range_ref(1):freq_range_ref(2));
        %interpolate data to resolve rounding differences
        data_interp=interp1(dat_temp_freq,dat_temp_lineshape,ref_freq);
        %define freq ranges

        %fit to lineshape to get freqency
        weights=ones(size(data_interp));

        ic=[1.7e4,0.14,93,8,ic_freq(m)];
        incr = [0.001,0.001,0.001,0.001,1e-6];

        Fit_Single_Lineshape_Handle=@fit_single_lineshape;
        [f,ref_pars,cvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(freq_range_ref,data_interp,ic,Fit_Single_Lineshape_Handle,1e-8,100,weights,incr);
        reference_shift(n,m)=ref_pars(5);

      end
    end
    peak_shift=reference_shift(1,1)-reference_shift(:,1);
  else
    %peak_shift=-[0.00000;-0.30148;-0.16430;-0.16436;0.31815;0.04772;0.23559;0.19515;0.17605];
  end
end


% Load remaining residues
data_free=full_data_vector=[];
for n=1:size(GLOBAL.Cyp_Conc,2)
  data_single_cond=[];
  for m=1:size(GLOBAL.RES_to_USE,2)
    dat_temp=load(strcat("dat/",GLOBAL.MUT,"/",int2str(GLOBAL.FGP_Conc(n)),"_",int2str(GLOBAL.Cyp_Conc(n)),"/hsqc_",int2str(GLOBAL.RES_to_USE(m)),"_",GLOBAL.CONFORMATION{m},".txt")); %load data
    %reshape and sum data over nitrogen
    nptsH=(find(dat_temp(:,2)!=dat_temp(1,2))-1)(1);
    nptsN=size(dat_temp,1)./nptsH;
    dat_temp_lineshape=sum(reshape(dat_temp(:,3),nptsH,nptsN),2)';
    dat_temp_lineshape-=min(dat_temp_lineshape);
    dat_temp_freq=dat_temp(1:nptsH);
    %Define frequencies
    GLOBAL.freq_range{n,m}=[find(GLOBAL.full_freq>=min(dat_temp_freq))(1),find(GLOBAL.full_freq<=max(dat_temp_freq))(end)];
    GLOBAL.data_freq{n,m}=GLOBAL.full_freq(GLOBAL.freq_range{n,m}(1):GLOBAL.freq_range{n,m}(2));
    if(n!=1) data_freq_linear=[data_freq_linear,GLOBAL.data_freq{n,m}]; end
    %interpolate data to resolve rounding differences
    data_interp=interp1(dat_temp_freq,dat_temp_lineshape,GLOBAL.data_freq{n,m});
    %define freq ranges of each data block 
    %Shift data to reference  
    data_interp=interp1(GLOBAL.data_freq{n,m}+peak_shift(n),data_interp,GLOBAL.data_freq{n,m},"extrap");
    %Zero data based on 50 pts in the middle of the spectrum
    data_interp-=mean(data_interp(floor(size(data_interp,2)./2)-25:(floor(size(data_interp,2)./2)+25)));
    data_single_cond=[data_single_cond,data_interp];

  end
  %Put data into vectors
  if(n==1)
    data_free=[data_free,data_single_cond];
  else
    if(GLOBAL.normalize) %normalize to GLOBAL.norm_factor
      data_temp{n,m}=data_single_cond./sum(data_single_cond)*GLOBAL.norm_factor;
      full_data_vector=[full_data_vector,data_single_cond./sum(data_single_cond)*GLOBAL.norm_factor];
    else
      full_data_vector=[full_data_vector,data_single_cond];
    end
  end
end

if(GLOBAL.normalize==0) full_data_vector=full_data_vector./sum(full_data_vector).*size(GLOBAL.Cyp_Conc,2).*GLOBAL.norm_factor; end

%Determine parameters in the absence of enzyme
%Once determined, will hard-code in to 10 decimal pts. to speed up troubleshooting 
calc_noenzyme_parameters=1;
if(calc_noenzyme_parameters==1)
  
  weights=ones(size(data_free));
  ic=[1.7e4,0.14,93,7,8,8,8,8,7637,7800,7442,7620];
  incr = [0.001,0.001,ones(1,6)*0.001,ones(1,4)*1e-6];
  
  Fit_all_noenzyme_Handle = @fit_all_noenzyme;
  [f,No_enzyme_pars,cvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(1,data_free,ic,Fit_all_noenzyme_Handle,1e-8,100,weights,incr);

  for n=1:size(No_enzyme_pars,1)
    printf("%0.10f,",No_enzyme_pars(n));
  end 
else
%  No_enzyme_pars=[35484.0883760160,0.1589869669,93.1888223180,6.9918387095,8.7497376647,10.6510708857,9.8458157527,9.4057383508,7637.9298789454,7802.5299809480,7444.3629700334,7620.2797614272];
end
SCALE=No_enzyme_pars(1);
GLOBAL.CT_RATIO=No_enzyme_pars(2);
GLOBAL.JNH=No_enzyme_pars(3);
GLOBAL.JHH=No_enzyme_pars(4);
GLOBAL.R_FREE=No_enzyme_pars(5:8);
GLOBAL.W_FREE=No_enzyme_pars(9:12);

%Set initial conditions/ranges
ic=[5e6,500.0,2000.0,2000.0,500.0,5.0e6,freq_bound(1),freq_bound(1),freq_bound(2),freq_bound(2)];
ic_range=[1e6,5e6;100,1000;100,2000;100,2000;10,500;5e6,2e7]
%set steps in fit
incr=[ones(1,6)*0.001,ones(1,4)*1e-3];
data_size=size(full_data_vector,2);

%Generate weighting vector:
data_weights=[];
for n=2:size(GLOBAL.Cyp_Conc,2)
  data_single_cond=[];
  for m=1:size(GLOBAL.RES_to_USE,2)
    %Weight cis heavier than trans, as the values are lower
    data_weights=[data_weights,ones(1,size(GLOBAL.data_freq{n,m},2))./(1+(GLOBAL.CT_RATIO-1)*strcmp(GLOBAL.CONFORMATION{m},"cis"))];
  end
end
data_weights /= mean(data_weights);% normalize average values of data weights to 1

weights=[];
wt_scale=100;
%Scale wt. by ~average value and secondary scaling factor
if(GLOBAL.constrain_ZZ) weights=[weights,data_size./7.0./wt_scale]; end
if(GLOBAL.constrain_Kd) weights=[weights,data_size./10.0./wt_scale]; end
if(GLOBAL.constrain_w) weights=[weights,10,10]; end
if(GLOBAL.constrain_lineshape) weights=[weights,data_weights]; end

options.bounds=[5e5,5e7;1,5e5;1,1e5;1,5e5;1,1e5;5e5,5e7; \
ic(7)-900*2,ic(7)+900*2;ic(8)-900*2,ic(8)+900*2;ic(9)-900*2,ic(9)+900*2;ic(10)-900*2,ic(10)+900*2];
ic_save=ic;
PARS_OUT=ic_OUT=r2_OUT=f_OUT=cvg_OUT=corp_OUT=covr_OUT=Z_OUT=Populations_OUT=scale_OUT=[];
PARS_OUT_1=r2_OUT_1=f_OUT_1=[];

Full_Data_Set_Handle=@full_data_set;
for n=1:iterations
  cvg=r2=0;
  while(cvg !=1 || r2<0.95)
    %randomize ICs
    ic_rand=ic_range(:,1)+rand(size(ic_range,1),1).*(ic_range(:,2)-ic_range(:,1));
    ic=[ic_rand',ic(7:end)];
    try
      GLOBAL.global_counter=0;  
      %Generate data vector to fit
      fit_data=[];
      if(GLOBAL.constrain_ZZ) fit_data=[fit_data,ZZ_match]; end
      if(GLOBAL.constrain_Kd) fit_data=[fit_data,Kd_match*GLOBAL.Kd_scale]; end
      if(GLOBAL.constrain_w) fit_data=[fit_data,freq_bound*GLOBAL.w_scale]; end
      if(GLOBAL.constrain_lineshape) fit_data=[fit_data,full_data_vector]; end

      [f,PARS,cvg,iter,corp,covp,covr,stdresid,Z,r2]=leasqr(GLOBAL,fit_data,ic,Full_Data_Set_Handle,1e-6,100,weights,incr,'dfdp',options);
      "----------------------------------------------------------------------------------------------------------------------"
      n
      r2

    catch
      "------------------------------ERROR--------------------------"
      lasterr
      "------------------------------ERROR--------------------------"
    end_try_catch
  end

  PARS_OUT=[PARS_OUT,PARS];
  r2_OUT=[r2_OUT,r2];
  f_OUT=[f_OUT,f];
  Populations_OUT=[Populations_OUT,return_EQ(PARS(1:6),1e-3,20e-6)];
  save(fout,"PARS_OUT","r2_OUT","f_OUT","full_data_vector","data_freq_linear","Populations_OUT");
end

end
%----------------------------------------------------------------------------------------------------%
function [kex,FGP_tc_ratio]=calc_ZZ(kab_in,kba_in,kbc_in,kcb_in,kcd_in,kdc_in,CypA_Conc_ZZ_in,FGP_Conc_ZZ_in=0.001);
%Input microcopic rate constants (in s-1 or Ms-1) and CypA conc (in M)
%return exchange rate for given CypA and FGP concentrations
%Initiate Global Variables
% states a:Free-trans
%        b:Enzyme-trans
%        c:Enzyme-cis
%        d:Free-cis

kab=kab_in; kba=kba_in; kbc=kbc_in; kcb=kcb_in; kcd=kcd_in; kdc=kdc_in;
CypA_Conc_ZZ=CypA_Conc_ZZ_in; FGP_Conc_ZZ=FGP_Conc_ZZ_in;

K_zz=[kab,kba,kbc,kcb,kcd,kdc];

%Determine Equilibrium Concentrations:
EQ_options.bounds=[0,1;0,1;0,1;0,1;0,1];
Solve_EQ_Handle=@solve_EQ;
Partial_EQ_Handle=@partial_EQ;
[f,EQ_conc_ZZ,cvg,iter]=leasqr([kab,kba,kbc,kcb,kcd,kdc],[0,0,0,CypA_Conc_ZZ,FGP_Conc_ZZ],[0.001*kba./(kba+kcd),1e-6,1e-6,0.001*kcd./(kba+kcd),CypA_Conc_ZZ./2.0],Solve_EQ_Handle,1e-6,1e6,ones(1,5),0.001*ones(1,5),Partial_EQ_Handle,EQ_options);
if(cvg != 1) error("Equillibrium Calculation did not converge"); end

FGP_tc_ratio=EQ_conc_ZZ(1)/EQ_conc_ZZ(4);

E=EQ_conc_ZZ(5);

Mtc=[-E*kab,   kba  ,    0   , 0     ;\
E*kab,-kba-kbc,   kcb  , 0     ;\
0    ,   kbc  ,-kcb-kcd, 0     ;\
0    ,    0   ,   kcd  , 0     ];

Mct=[ 0    ,   kba  ,    0   , 0     ;\
0    ,-kba-kbc,   kcb  , 0     ;\
0    ,   kbc  ,-kcb-kcd, E*kdc ;\
0    ,    0   ,   kcd  ,-E*kdc ];

%Calculate Eigenvalues/vectors
[V,lam_ct]=eig(Mct);
[V,lam_tc]=eig(Mtc);

ktc=min(abs(diag(lam_tc)(abs(diag(lam_tc))>1e-5)));
kct=min(abs(diag(lam_ct)(abs(diag(lam_ct))>1e-5)));

kex=(kct+ktc);

end
%----------------------------------------------------------------------------------------------------%
function output=solve_EQ(X,A)
%For Solving 4 state exhchange at equillibrium
%A(1:4) = concentration of states a-d
%A(5) = concentration of free CypA
%X = [kab,kba,kbc,kcb,kcd,kdc]
%x1-x3 are flux through each of the a-b, b-c, c-d pathways - Set to Zero
%x4 is total CypA concentration
%x5 is total FGP concentration

kab=X(1);
kba=X(2);
kbc=X(3);
kcb=X(4);
kcd=X(5);
kdc=X(6);

x1=A(1)*A(5)*kab-A(2)*kba;
x2=A(2)*kbc-A(3)*kcb;
x3=A(3)*kcd-A(4)*A(5)*kdc;
x4=A(5)+A(2)+A(3);
x5=sum(A(1:4));

output=[x1,x2,x3,x4,x5];
end
%----------------------------------------------------------------------------------------------------%
function output=partial_EQ(X,f,A,dp,func)
%Calculates Jacobian of equation in function solve_EQ 
kab=X(1);
kba=X(2);
kbc=X(3);
kcb=X(4);
kcd=X(5);
kdc=X(6);

output = [A(5)*kab,0,0,0,1 ;\
-kba,kbc,0,0,1   ;\
0,-kcb,kcd,1,1   ;\
0,0,-A(5)*kdc,0,1;\
A(1)*kab,0,-A(4)*kdc,1,0]';
end
%----------------------------------------------------------------------------------------------------%
function EQ_conc=return_EQ(K,FGP_Conc_in,CypA_Conc_in)
Solve_EQ_Handle=@solve_EQ;
Partial_EQ_Handle=@partial_EQ;
[f,EQ_conc,cvg,iter]=leasqr(K,[0,0,0,CypA_Conc_in,FGP_Conc_in],[FGP_Conc_in/4,FGP_Conc_in/4,FGP_Conc_in/4,FGP_Conc_in/4,CypA_Conc_in/2],Solve_EQ_Handle,1e-8,10000,ones(1,5),0.001*ones(1,5),Partial_EQ_Handle);
end
%----------------------------------------------------------------------------------------------------%
function [lineshape,freq]=Calc_NoEnzyme_lineshape(c_ratio,w,JNH,JHH,R,peak="trans")  
%Calculates linshapes for the peptide in the absense of exchange btwn cis and trans
global GLOBAL

t=GLOBAL.t_step;
M=0;
for nn=1:4
  J=-JNH./2+JNH*(nn<=2)-JHH./2+JHH*rem(nn,2);
  if(strcmp(peak,"trans")==1)  M+=exp((2*pi()*i*(w(1)+J)-R(1))*t);end
  if(strcmp(peak,"cis")==1) M+=c_ratio*exp((2*pi()*i*(w(1)+J)-R(1))*t);end
end
lineshape=real(fft(M));
freq=[0:1/t(2)/(size(lineshape,2)):1/t(2)-t(2)/(size(lineshape,2))];
end
%----------------------------------------------------------------------------------------------------%
function output=fit_all_noenzyme(X,A)
%Use with leqsqr to determine c/t ratio, J coupling constants, and free relaxation parameters
% X=[scaling_factor,c_ratio,JNH,JHH,Rs,ws]
global GLOBAL
output=[]; 
Rstart=5;
wstart=Rstart+size(unique(GLOBAL.RES_to_USE),2)*2;
for m=1:size(GLOBAL.RES_to_USE,2)
  Rall = A(Rstart);
  x=A(1).*Calc_NoEnzyme_lineshape(A(2),A(wstart),A(3),A(4),Rall,GLOBAL.CONFORMATION{m});
  Rstart+=1;
  wstart+=1;
  output=[output,x(GLOBAL.freq_range{1,m}(1):GLOBAL.freq_range{1,m}(2))];
end
end
%----------------------------------------------------------------------------------------------------%
function output=fit_single_lineshape(X,A)
%Use with leqsqr to determine c/t ratio, J coupling constants, and free relaxation parameters
% A=[scaling_factor,JNH,JHH,R,w]
% X=freq_range
global GLOBAL
x=A(1).*Calc_NoEnzyme_lineshape(1,A(5),A(2),A(3),A(4));
output=x(X(1):X(2));
end
%----------------------------------------------------------------------------------------------------%
function [lineshape,freq_out,shift_w]=Calc_lineshape(k,w,R,Cyp,FGP,EQ_conc,peak="both")
%Determine lineshapes in the presence of enzyme
%w = [w_trans,w_trans_bound,w_cis_bound,w_cis]
%R = [R_trans,R_trans_bound,R_cis_bound,R_cis]
%k=[kab,kba,kbc,kcb,kcd]
%peak = "cis" or "trans" - will return only this spectrum or "both" will return a sum of both
global GLOBAL

%Shift frequencies to reduce FT cost
freq_window=450; %window of freqencies to simulate about 0 after shift
if(strcmp(peak,"trans")) w_ref=[w(1),w(1)]; end
if(strcmp(peak,"cis")) w_ref=[w(4),w(4)]; end
if(strcmp(peak,"bound")) w_ref=[w(2),w(3)]; end
freq=GLOBAL.full_freq; %Frequencies

shift_w=[0,0];
shift_w(1)=freq(find(abs(freq-w_ref(1))==min(abs(freq-w_ref(1))))); %find closest measured freqency to shift by
shift_w(2)=freq(find(abs(freq-w_ref(2))==min(abs(freq-w_ref(2))))); %find closest measured freqency to shift by
shift_w=mean(shift_w);
w-=shift_w; %shift frequencies

%Set new time points
N_measure=size(freq,2);
N_tpts=size(find((freq<(max(w_ref)+freq_window)).*(freq>(min(w_ref)-freq_window))),2);

t=(0:GLOBAL.t_app./N_tpts:GLOBAL.t_app-GLOBAL.t_app/N_tpts);

w_trans_free=w(1);  w_trans_bound=w(2);  w_cis_bound=w(3);  w_cis_free=w(4);
R2_trans_free=R(1);  R2_trans_bound=R(2);  R2_cis_bound=R(3);  R2_cis_free=R(4);
kab=k(1);  kba=k(2);  kbc=k(3);  kcb=k(4);  kcd=k(5);  kdc=k(6);

zf=size(t,2);
lineshape=zeros(1,zf);
fid_4states=zeros(4,zf);
M=0;
for nn=1:4
  J=-GLOBAL.JNH./2+GLOBAL.JNH*(nn<=2)-GLOBAL.JHH./2+GLOBAL.JHH*rem(nn,2);
  if(strcmp(peak,"bound")) J=0; end
  %Evolution matrix
  M=[2*pi()*i*(w_trans_free+J)-R2_trans_free-kab*EQ_conc(5),  kba , 0 , 0;\
  kab*EQ_conc(5), 2*pi()*i*(w_trans_bound+J)-R2_trans_bound-kba-kbc, kcb, 0;\
  0, kbc, 2*pi()*i*(w_cis_bound+J)-R2_cis_bound-kcb-kcd, kdc*EQ_conc(5);\
  0, 0, kcd, 2*pi()*i*(w_cis_free+J)-R2_cis_free-kdc*EQ_conc(5)];
  %Calculate Eigenvalues/vectors
  [V,lam]=eig(M);
  %Set initial conditions/Calculate pre-exponential terms for specific solution
  M0=[EQ_conc(1:4)];
  C=inv(V)*M0;
  %generate fid
  fid_4states+=V*(exp(diag(lam)*t).*repmat(C,1,size(t,2)));
end
if(strcmp(peak,"trans"))
  lineshape=real(fft(fid_4states(1,:),zf));
end
if(strcmp(peak,"cis"))
  lineshape=real(fft(fid_4states(4,:),zf));
end
if(strcmp(peak,"bound"))
  lineshape=real(fft(fid_4states(2,:),zf));
  lineshape+=real(fft(fid_4states(3,:),zf));
end
lineshape*=N_measure./N_tpts;

freq_out=[0:1/t(2)/(size(lineshape,2)):1/t(2)-t(2)/(size(lineshape,2))];
freq_out-=freq_out(ceil(((size(freq_out,2)+1)./2)));
end
%----------------------------------------------------------------------------------------------------%
function output=full_data_set(GLOBAL,A)
%A=[kba,kcb,kcd,R_bound,w_bound,Scale_factor]
%R_bound and w_bound are [bound_trans,bound_cis] for each residue
%Constrain_w is a flag to constrain the max of the bound state, 1=constrain (output 3) 0=don't constrain
global full_data_vector GLOBAL

K=[A(1),A(2),A(3),A(4),A(5),A(6)];

R_bound=GLOBAL.R20bound;
w_bound=A(7:end);

ls_output=[];

EQ_conc=1e-6.*[GLOBAL.FGP_Conc(1)./(1+GLOBAL.CT_RATIO),1,1,GLOBAL.FGP_Conc(1)*GLOBAL.CT_RATIO/(1+GLOBAL.CT_RATIO),GLOBAL.Cyp_Conc(1)];
EQ_options.bounds=[0,1;0,1;0,1;0,1;0,1;];
for n=2:size(GLOBAL.Cyp_Conc,2)
  start_count=1;
  single_cond=[];

  Solve_EQ_Handle=@solve_EQ;
  Partial_EQ_Handle=@partial_EQ;
  [f,EQ_conc,cvg,iter]=leasqr(K,[0,0,0,1e-6*GLOBAL.Cyp_Conc(n),1e-6*GLOBAL.FGP_Conc(n)],EQ_conc,Solve_EQ_Handle,1e-6,1e6,ones(1,5),0.001*ones(1,5),Partial_EQ_Handle,EQ_options);

  for m=1:size(GLOBAL.RES_to_USE,2)
    R_all = [GLOBAL.R_FREE(start_count),R_bound(start_count),R_bound(start_count+1),(GLOBAL.R_FREE(start_count+1))];
    w_all = [GLOBAL.W_FREE(start_count),w_bound(start_count),w_bound(start_count+1),(GLOBAL.W_FREE(start_count+1))];
    [x,freq_temp,freq_shift]=Calc_lineshape(K,w_all,R_all,GLOBAL.Cyp_Conc(n)*1e-6,GLOBAL.FGP_Conc(n)*1e-6,EQ_conc,GLOBAL.CONFORMATION{m});
    if(strcmp(GLOBAL.CONFORMATION{m},"cis")) start_count+=2; end

    x=fftshift(x); 
    x-=min(x);
    freq_temp+=freq_shift; %shift back subtracted frequency
    ls_temp=interp1(freq_temp,x,GLOBAL.data_freq{n,m}); %interpolate back onto data frequencies (corrects for rounding)

    %subtract off baseine, defined by middle 50 pts.
    ls_temp-=mean(ls_temp(floor(size(ls_temp,2)./2)-25:(floor(size(ls_temp,2)./2)+25)));
    single_cond=[single_cond,ls_temp];
  end
  if(GLOBAL.normalize) single_cond/=sum(single_cond)/GLOBAL.norm_factor; end  %normalize area to GLOBAL.norm_factor
  ls_output=[ls_output,single_cond];
end
if(GLOBAL.normalize==0) ls_output=ls_output./sum(ls_output).*size(GLOBAL.Cyp_Conc,2).*GLOBAL.norm_factor; end  %normalize area to GLOBAL.norm_factor

%if(GLOBAL.constrain_w) %Constrain bound chemical shift values
freq_output=[0,0];
[f,EQ_conc,cvg,iter]=leasqr(K,[0,0,0,4e-3,1e-3],[1e-6,0.5e-3,0.5e-3,1e-6,3e-3],Solve_EQ_Handle,1e-6,1e6,ones(1,5),0.001*ones(1,5),Partial_EQ_Handle);
for m=1:2
  R_all = [GLOBAL.R_FREE((m-1)*2+1),R_bound((m-1)*2+1),R_bound((m-1)*2+2),GLOBAL.R_FREE((m-1)*2+2)];
  w_all = [GLOBAL.W_FREE((m-1)*2+1),w_bound((m-1)*2+1),w_bound((m-1)*2+2),(GLOBAL.W_FREE((m-1)*2+2))];
  [x,freq_temp,freq_shift]=Calc_lineshape(K,w_all,R_all,4e-3,1e-3,EQ_conc,"bound");
  freq_temp=round(freq_temp.*10)./10;
  x=fftshift(x);
  x-=min(x);
  freq_temp+=freq_shift;
  freq_temp=round(freq_temp.*10)./10;
  freq_output(m)=freq_temp(x==max(x));
end  
%end

%Calcualate kexeff
[zz_output,zz_ratio]=calc_ZZ(K(1),K(2),K(3),K(4),K(5),K(6),20e-6);

%Calc Kdapp
Keq=zz_ratio;
Kdc=K(5)/K(6);
Kdt=K(2)/K(1);
kd_output=Kda=(Keq+1)/(Keq/Kdt+1/Kdc);

GLOBAL.global_counter++;

%Cut out if ratio is too low for multiple rounds
if(GLOBAL.global_counter==1) GLOBAL.low_frac=zeros(50,1); end
GLOBAL.low_frac(2:end)=GLOBAL.low_frac(1:49);
GLOBAL.low_frac(1)=(zz_ratio<1);
if(sum(GLOBAL.low_frac)==50) error("ZZ Ratio too low"); end

if(rem(GLOBAL.global_counter,10)==0) %plot out during fit
  hold off
  plot(full_data_vector,"2")
  hold on
  plot(ls_output,"3")
  print(strcat(GLOBAL.MUT,".eps"),"-deps","-color","-S900,300")
end

%Print out parameters during run
printf("%s %6.2e %6.0f %6.0f %6.0f %6.0f %6.2e / %5.0f %5.0f %5.0f %5.0f \n   %6.3e / %4.2f / %4.2f / %4.0f %4.0f / %4.3f /%4.3f \n",GLOBAL.MUT,K(1),K(2),K(3),K(4),K(5),K(6),w_bound(1),w_bound(2),w_bound(3),w_bound(4),kd_output,zz_output,zz_ratio,freq_output(1),freq_output(2),R_bound(1),(corr(ls_output,full_data_vector)).^2);
fflush(stdout);

%Determine Output:
output=[];
if(GLOBAL.constrain_ZZ) output=[output,zz_output]; end
if(GLOBAL.constrain_Kd) output=[output,kd_output.*GLOBAL.Kd_scale]; end
if(GLOBAL.constrain_w)  output=[output,freq_output.*GLOBAL.w_scale]; end
if(GLOBAL.constrain_lineshape)  output=[output,ls_output]; end 

end
%----------------------------------------------------------------------------------------------------%
