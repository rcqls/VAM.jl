struct FamilyModel
    dHR::Vector{Float64}
    dHL::Vector{Float64}
    dhR::Vector{Float64}
    dhL::Vector{Float64}
    d2HR::Matrix{Float64}
    d2HL::Matrix{Float64}
    d2h::Matrix{Float64}
    dhd::Matrix{Float64}
    nb_params::Int
end

#   virtual NumericVector get_params() = 0;

# 	virtual void set_params(NumericVector par) = 0;

# 	virtual double hazardRate(double x) = 0;

# 	virtual double cumulative_hazardRate(double x) = 0;

# 	virtual double inverse_cumulative_hazardRate(double x) = 0;

# 	virtual double hazardRate_derivative(double x) = 0;

#   virtual double inverse_hazardRate(double x) = 0;

# 	virtual double* hazardRate_param_derivative(double x,bool R) = 0;//Size:nb_params_-1

# 	virtual double* cumulative_hazardRate_param_derivative(double x,bool R) = 0;//Size:nb_params_-1

#   virtual double* hazardRate_derivative_param_derivative(double x) = 0;//Size:nb_params_-1

#   virtual double hazardRate_2derivative(double x) = 0;

#   virtual double* hazardRate_param_2derivative(double x) = 0;//Size:nb_params_*(:nb_params_-1)/2

#   virtual double* cumulative_hazardRate_param_2derivative(double x,bool R) = 0;//Size:nb_params_*(:nb_params_-1)/2


#   void init_Familiy(){
#     dHR = new double[nb_params_-1];
#     dHL = new double[nb_params_-1];
#     dhR = new double[nb_params_-1];
#     dhL = new double[nb_params_-1];
#     dhd = new double[nb_params_-1];

#     d2HR = new double[(nb_params_-1)*nb_params_/2];
#     d2HL = new double[(nb_params_-1)*nb_params_/2];
#     d2h = new double[(nb_params_-1)*nb_params_/2];
#   }

#   int nb_params() {
#     return nb_params_;
#   }
# };

struct WeibullFamilyModel <: FamilyModel
    α::Float64
    β::Float64
end
# public:
#   WeibullFamilyModel(NumericVector par) {
#     nb_params_=2;
#   	alpha=par[0];beta=par[1];
#     init_Familiy();
#   }

#   double alpha, beta;

#   NumericVector get_params() {
#   	NumericVector out(2);
#   	out[0]=alpha;out[1]=beta;
#   	return out;
#   }

#   void set_params(NumericVector par) {
#   		alpha=par[0];beta=par[1];
#   }

  

#   double hazardRate(double x) {
#   	return (x<=0 ? 0 : alpha*beta*pow(x,beta-1));
#   }

#   double inverse_hazardRate(double x) {
#     return (x<=0 ? 0 : pow(x/alpha/beta,1/(beta-1)));
#   }

#   double cumulative_hazardRate(double x) {
#   	return (x<=0 ? 0 : alpha*pow(x,beta));
#   }

#   double inverse_cumulative_hazardRate(double x) {
#   	 return (x<=0 ? 0 : pow(x/alpha,1/beta));

#   }

#   double hazardRate_derivative(double x) {
#   	return (x<=0 ? 0 : alpha*beta*(beta-1)*pow(x,beta-2));
#   }

#   double* hazardRate_param_derivative(double x,bool R) {
#     double *dh;
#     if(R) dh=dhR; else dh=dhL;
#   	dh[0]=(x<=0 ? 0 : alpha*(1+beta*log(x))*pow(x,beta-1));
#     return dh;
#   }

#   double* cumulative_hazardRate_param_derivative(double x,bool R) {
#     double *dH;
#     if(R) dH=dHR; else dH=dHL;
#     dH[0]= (x<=0 ? 0 : alpha*log(x)*pow(x,beta));
#     return dH;
#   }

#   double* hazardRate_derivative_param_derivative(double x) {
#     dhd[0]= (x<=0 ? 0 : alpha*(2*beta-1+beta*(beta-1)*log(x))*pow(x,beta-2));
#     return dhd;
#   }

#   double hazardRate_2derivative(double x) {
#     return (x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x,beta-3));
#   }

#   double* hazardRate_param_2derivative(double x) {
#     d2h[0]= (x<=0 ? 0 : alpha*(2+beta*log(x))*log(x)*pow(x,beta-1));
#     return d2h;
#   }

#   double* cumulative_hazardRate_param_2derivative(double x,bool R) {
#     double *d2H;
#     if(R) d2H=d2HR; else d2H=d2HL;
#     d2H[0]= (x<=0 ? 0 : alpha*pow(log(x),2)*pow(x,beta));
#     return d2H;
#   }
 
# };


struct LogLinearFamilyModel <:  FamilyModel
    α::Float64
    β::Float64
end
#   public:
#     LogLinearFamilyModel(NumericVector par) {
#       nb_params_=2;
#       alpha=par[0];beta=par[1];
#       init_Familiy();
      
#       //see just above the explanation
#       //Why not this value (and also some visual check), may be that can be changed by user in further version ?
#       LDorder=5;
#       bxLim=0.000001;
#     }
  
#   double alpha, beta;
#   //With the cumulative_hazardRate and its param derivaatives, 
#   //  theres is a numerical problem (due to undertermined limits) when beta*x tends to 0.
#   //  The problem is solved by replacing the functions values by their limited development in beta*x=0.
#   // LDorder corresponds to the order of the corresponding develomment
#   // bxLim corresponds to the limiting value above which the limited development is used.
#   int LDorder;
#   double bxLim;

#   NumericVector get_params() {
#     NumericVector out(2);
#     out[0]=alpha;out[1]=beta;
#     return out;
#   }
  
#   void set_params(NumericVector par) {
#     alpha=par[0];beta=par[1];
#   }
  
#   double hazardRate(double x) {
#     return (x<0 ? 0 : alpha*exp(beta*x));
#   }

#   double inverse_hazardRate(double x) {
#     return log(x/alpha)/beta;
#   }
  
#   double cumulative_hazardRate(double x) {
#     double res;
#     if(abs(beta*x)<bxLim){
#       double prec=beta*x/2;
#       res=1+prec;
#       for (int i=1; i< LDorder ;i++){
#         prec=prec*beta*x/(i+1);
#         res=res+prec;
#       }
#       res=alpha*res*x;
#     } else {
#       res=alpha*(exp(beta*x)-1)/beta;
#     }
#     return res;
#   }
  
#   double inverse_cumulative_hazardRate(double x) {
#     return log(1+x*beta/alpha)/beta;
    
#   }
  
#   double hazardRate_derivative(double x) {
#     return (x<0 ? 0 : alpha*beta*exp(beta*x));
#   }
  
#   double* hazardRate_param_derivative(double x,bool R) {
#     double *dh;
#     if(R) dh=dhR; else dh=dhL;
#     dh[0]= alpha*x*exp(beta*x) ;
#     return dh;
#   }
  
#   double* cumulative_hazardRate_param_derivative(double x,bool R) {
#     double *dH;
#     double res;
#     double prec;
#     if(R) dH=dHR; else dH=dHL;

#     if(abs(beta*x)<bxLim){
#       prec=beta*x/6;
#       res=0.5+2*prec;
#       for (int i=1; i< LDorder ;i++){
#         prec=prec*beta*x/(i+3);
#         res=res+(i+2)*prec;
#       }
#       res=alpha*res*pow(x,2);
#     } else {
#       res=alpha*(x*exp(x*beta)/beta-(exp(beta*x)-1)/pow(beta,2));
#     }

#     dH[0]=  res;
#     return dH;
#   }
  
#   double* hazardRate_derivative_param_derivative(double x) {
#     dhd[0]= alpha*exp(beta*x)*(1+beta*x) ;
#     return dhd;
#   }

#   double hazardRate_2derivative(double x) {
#     return (x<=0 ? 0 : alpha*pow(beta,2)*exp(beta*x));
#   }

#   double* hazardRate_param_2derivative(double x) {
#     d2h[0]=  alpha*pow(x,2)*exp(beta*x) ;
#     return d2h;
#   }
  
#   double* cumulative_hazardRate_param_2derivative(double x,bool R) {
#     double *d2H;
#     double res;
#     double prec;
#     if(R) d2H=d2HR; else d2H=d2HL;
    
#     if(abs(beta*x)<bxLim){
#       prec=beta*x/24;
#       res=(2/3)+6*prec;
#       for (int i=1; i< LDorder ;i++){
#         prec=prec*beta*x/(i+4);
#         res=res+(i+2)*(i+3)*prec;
#       }
#       res=alpha*res*pow(x,3);
#     } else {
#       res=alpha*(pow(x,2)*exp(x*beta)/beta-2*x*exp(x*beta)/pow(beta,2)+2*(exp(beta*x)-1)/pow(beta,3));
#     }

#     d2H[0]= res;
#     return d2H;
#   }

# };

struct Weibull3FamilyModel <: FamilyModel
    α::Float64
    β::Float64
    c::Float64
end
# public:
#   Weibull3FamilyModel(NumericVector par) {
#     nb_params_=3;
#     alpha=par[0];beta=par[1];c=par[2];
#     init_Familiy();
#   }

#   double alpha, beta, c;

#   NumericVector get_params() {
#     NumericVector out(3);
#     out[0]=alpha;out[1]=beta;out[2]=c;
#     return out;
#   }

#   void set_params(NumericVector par) {
#       alpha=par[0];beta=par[1];c=par[2];
#   }

#   double hazardRate(double x) {
#     return (x<0 ? 0 : alpha*beta*pow(x+c,beta-1));
#   }

#   double inverse_hazardRate(double x) {
#     return (x<=0 ? 0 : pow(x/alpha/beta,1/(beta-1))-c);
#   }

#   double cumulative_hazardRate(double x) {
#     return (x<=0 ? 0 : alpha*(pow(x+c,beta)-pow(c,beta)));
#   }

#   double inverse_cumulative_hazardRate(double x) {
#      return (x<=0 ? 0 : pow(pow(c,beta)+x/alpha,1/beta)-c);

#   }

#   double hazardRate_derivative(double x) {
#     return (x<=0 ? 0 : alpha*beta*(beta-1)*pow(x+c,beta-2));
#   }

#   double* hazardRate_param_derivative(double x,bool R) {
#     double *dh;
#     if(R) dh=dhR; else dh=dhL;
#     dh[0]=(x==0 ? 0 : alpha*(1+beta*log(x+c))*pow(x+c,beta-1));
#     dh[1]=(x<=0 ? 0 : alpha*beta*(beta-1)*pow(x+c,beta-2));
#     return dh;
#   }

#   double* cumulative_hazardRate_param_derivative(double x,bool R) {
#     double *dH;
#     if(R) dH=dHR; else dH=dHL;
#     dH[0]= (x==0 ? 0 : alpha*(log(x+c)*pow(x+c,beta)-log(c)*pow(c,beta)));
#     dH[1]= (x<=0 ? 0 : alpha*beta*(pow(x+c,beta-1)-pow(c,beta-1)));
#     return dH;
#   }

#   double* hazardRate_derivative_param_derivative(double x) {
#     dhd[0]= (x==0 ? 0 : alpha*(2*beta-1+beta*(beta-1)*log(x+c))*pow(x+c,beta-2));
#     dhd[1]=(x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x+c,beta-3));
#     return dhd;
#   }

#   double hazardRate_2derivative(double x) {
#     return (x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x+c,beta-3));
#   }

#   double* hazardRate_param_2derivative(double x) {
#     d2h[0]= (x==0 ? 0 : alpha*(2+beta*log(x+c))*log(x+c)*pow(x+c,beta-1));
#     d2h[1]=(x==0 ? 0 : alpha*pow(x+c,beta-2)*(2*beta-1+beta*(beta-1)*log(x+c)));
#     d2h[2]=(x<=0 ? 0 : alpha*beta*(beta-1)*(beta-2)*pow(x+c,beta-3));
#     return d2h;
#   }

#   double* cumulative_hazardRate_param_2derivative(double x,bool R) {
#     double *d2H;
#     if(R) d2H=d2HR; else d2H=d2HL;
#     d2H[0]= (x==0 ? 0 : alpha*(pow(log(x+c),2)*pow(x+c,beta)-pow(log(c),2)*pow(c,beta)));
#     d2H[1]=(x==0 ? 0 : alpha*(pow(x+c,beta-1)*(beta*log(x+c)+1)-pow(c,beta-1)*(beta*log(c)+1)));
#     d2H[2]=(x<=0 ? 0 : alpha*beta*(beta-1)*(pow(x+c,beta-2)-pow(c,beta-2)));
#     return d2H;
#   }
 
# };


# FamilyModel* newFamilyModel(List family);
