mutable struct Model <: AbstractModel
	
    k::Int # current system
    nb_system::Int #number of system

	nbPM::Int
    idMod::Int
	nb_paramsMaintenance::Int
    nb_paramsFamily::Int
    nb_paramsCov::Int

	data::Vector{DataFrame} #::List # TODO

	#Additional Covariates stuff
	data_cov #::DataFrame
	params_cov::Vector{Float64}
	sum_cov::Float64 #to save the computation

	time::Vector{Float64}
	type::Vector{Int}

    indType::Float64 #TODO: remove S4 if unused!

    Vleft::Float64
    Vright::Float64
    hVleft::Float64


	dVleft::Vector{Float64}
    dVright::Vector{Float64}
    
	d2Vleft::Matrix{Float64}
    d2Vright::Matrix{Float64}

	A::Float64
	dA::Vector{Float64}
	d2A::Matrix{Float64}

	mu::Int
	VR_prec::Float64
	dVR_prec::Float64
    d2VR_prec::Matrix{Float64}

    comp::Compute

	models::Vector{AbstractMaintenanceModel}

	family::FamilyModel

	maintenance_policy::AbstractMaintenancePolicy
end

# 	FamilyModel* get_family() {
# 	 	return family;
# 	}

# 	//Ununsed and can create bugs
# 	//NumericVector get_dB(int k) {
# 	// 	NumericVector dBkR(nb_paramsMaintenance);
# 	// 	for (int i=0;i<nb_paramsMaintenance;i++){
# 	// 		dBkR[i]=dB[k*nb_paramsMaintenance+i];
# 	// 	}
# 	// 	return dBkR;
# 	// }

# 	// NumericMatrix get_d2B(int k){
# 	// 	int j;
# 	// 	int di=k*nb_paramsMaintenance*(nb_paramsMaintenance+1)/2;
# 	// 	NumericMatrix d2BkR(nb_paramsMaintenance,nb_paramsMaintenance);
# 	// 	for (int i=0;i<nb_paramsMaintenance;i++) {
# 	// 		d2BkR(i,i)=d2B[di+i*(i+1)/2+i];
# 	// 		for (j=0;j<i;j++) {
# 	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
# 	// 			d2BkR(i,j)=d2B[di+i*(i+1)/2+j];
# 	// 			d2BkR(j,i)=d2BkR(i,j);
# 	// 		}
# 	// 	}
# 	// 	return d2BkR;
# 	// }

# 	// List get() {
# 	// 	int j;
# 	// 	int n_params=nb_paramsMaintenance+nb_paramsFamily-1;

# 	// 	List ret;
# 	// 	ret["S1"]=NumericVector::create(S1);ret["S2"]=NumericVector::create(S2);ret["S0"]=NumericVector::create(S0);ret["S3"]=NumericVector::create(S3);
# 	// 	ret["Vright"]=NumericVector::create(Vright);ret["Vleft"]=NumericVector::create(Vleft);
# 	// 	//printf("S0=%f, S1=%f, S2=%f, S3=%f, Vright=%f, Vleft=%f\n",S0,S1,S2,S3,Vright,Vleft);
# 	// 	NumericVector dS1R(n_params),dS2R(n_params), dS3R(nb_paramsMaintenance);
# 	// 	NumericMatrix d2S1R(n_params,n_params),d2S2R(n_params,n_params), d2S3R(nb_paramsMaintenance,nb_paramsMaintenance);
# 	// 	ret["dS1"]=dS1R;ret["dS2"]=dS2R;ret["dS3"]=dS3R;
# 	// 	ret["d2S1"]=d2S1R;ret["d2S2"]=d2S2R;ret["d2S3"]=d2S3R;

# 	// 	for (int i=0;i<n_params;i++) {
# 	// 		dS1R[i]=dS1[i];
# 	// 		dS2R[i]=dS2[i];
# 	// 		//printf("dS1[%d]=%f, d2S2=%f\n",i,dS1[i],dS2[i]);
# 	// 		d2S1R(i,i)=d2S1[i*(i+1)/2+i];
# 	// 		d2S2R(i,i)=d2S2[i*(i+1)/2+i];
# 	// 		for (j=0;j<i;j++) {
# 	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
# 	// 			d2S1R(i,j)=d2S1[i*(i+1)/2+j];
# 	// 			d2S2R(i,j)=d2S2[i*(i+1)/2+j];
# 	// 			//printf("d2S1[%d,%d][%d]=%f, d2S2=%f\n",i,j,i*(i+1)/2+j,d2S1[i*(i+1)/2+j],d2S2[i*(i+1)/2+j]);
# 	// 			d2S1R(j,i)=d2S1R(i,j);
# 	// 			d2S2R(j,i)=d2S2R(i,j);
# 	// 			//printf("d2S1[%d,%d][%d]=%f, d2S2=%f\n",j,i,i*(i+1)/2+j,d2S1R(i,j),d2S2R(i,j));
# 	// 		}
# 	// 	}
# 	// 	NumericVector dVrightR(nb_paramsMaintenance),dVleftR(nb_paramsMaintenance);
# 	// 	NumericMatrix d2VrightR(nb_paramsMaintenance,nb_paramsMaintenance),d2VleftR(nb_paramsMaintenance,nb_paramsMaintenance);
# 	// 	NumericVector dAR(nb_paramsMaintenance);
# 	// 	NumericMatrix d2AR(nb_paramsMaintenance,nb_paramsMaintenance);
# 	// 	NumericVector dCR(nb_paramsMaintenance);
# 	// 	NumericMatrix d2CR(nb_paramsMaintenance,nb_paramsMaintenance);

# 	// 	ret["dVright"]=dVrightR;ret["dVleft"]=dVleftR;
# 	// 	ret["d2Vright"]=d2VrightR;ret["d2Vleft"]=d2VleftR;
# 	// 	ret["A"]=NumericVector::create(A);
# 	// 	ret["dA"]=dAR;
# 	// 	ret["d2A"]=d2AR;
# 	// 	ret["C"]=NumericVector::create(C);
# 	// 	ret["dC"]=dCR;
# 	// 	ret["d2C"]=d2CR;
# 	// 	for (int i=0;i<nb_paramsMaintenance;i++) {
# 	// 		dVrightR[i]=dVright[i];
# 	// 		dVleftR[i]=dVleft[i];
# 	// 		dS3R[i]=dS3[i];
# 	// 		dAR[i]=dA[i];
# 	// 		dCR[i]=dC[i];
# 	// 		d2VrightR(i,i)=d2Vright[i*(i+1)/2+i];
# 	// 		d2VleftR(i,i)=d2Vleft[i*(i+1)/2+i];
# 	// 		d2S3R(i,i)=d2S3[i*(i+1)/2+i];
# 	// 		d2AR(i,i)=d2A[i*(i+1)/2+i];
# 	// 		d2CR(i,i)=d2A[i*(i+1)/2+i];
# 	// 		for (j=0;j<i;j++) {
# 	// 			//i and j(<=i) respectively correspond to the line and column indices of (inferior diagonal part of) the hessian matrice
# 	// 			d2VrightR(i,j)=d2Vright[i*(i+1)/2+j];
# 	// 			d2VleftR(i,j)=d2Vleft[i*(i+1)/2+j];
# 	// 			d2S3R(i,j)=d2S3[i*(i+1)/2+j];
# 	// 			d2AR(i,j)=d2A[i*(i+1)/2+j];
# 	// 			d2CR(i,j)=d2C[i*(i+1)/2+j];
# 	// 			d2VrightR(j,i)=d2VrightR(i,j);
# 	// 			d2VleftR(j,i)=d2VleftR(i,j);
# 	// 			d2S3R(j,i)=d2S3R(i,j);
# 	// 			d2AR(j,i)=d2AR(i,j);
# 	// 			d2CR(j,i)=d2CR(i,j);
# 	// 		}
# 	// 	}
# 	// 	List BR(max_mem);
# 	// 	List dBR(max_mem);
# 	// 	List d2BR(max_mem);
# 	// 	ret["B"]=BR;
# 	// 	ret["dB"]=dBR;
# 	// 	ret["d2B"]=d2BR;
# 	// 	for (j=0;j<max_mem;j++){
# 	// 		BR[j]=B[j];
# 	// 		dBR[j]=get_dB(j);
# 	// 		d2BR[j]=get_d2B(j);
# 	// 	}

# 	// 	return ret;
# 	// }

# 	void set_data(List data_);

# 	void select_data(int i);

# 	DataFrame get_selected_data(int i);

# 	NumericVector get_params();

#     void set_params(NumericVector pars);

#     double virtual_age(double x) ;

#     double virtual_age_inverse(double x);

#     void update_Vleft(bool with_gradient,bool with_hessian);

#     List get_virtual_age_infos(double by, double from, double to);

# 	void init_computation_values();

# 	//Covariates related

# 	int current_system;

#     void select_current_system(int i,bool compute) {
#         //Covariates related
#         current_system=i;
# 		//simulation: compute=false since only computation in c++ and set_current_system in R
# 		//mle: compute=true since both computation and select_current_system in c++
# 		if(compute) compute_covariates();
#     }

# 	double compute_covariates(); //output maybe useful inside R

# 	double get_covariate(int j);

# 	void set_covariates(List covariates_);


# private:
# 	void set_models(List models_);

#     void set_family(List family_);

#     void set_maintenance_policy(List maintenance_policy_);

# 	void init(List model_);

# 	void init_virtual_age_infos();

# 	DataFrame get_virtual_age_info(double from,double to,double by, double expCov);

# };