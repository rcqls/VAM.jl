mutable struct Model <: AbstractModel
	
    k::Int # current system
    nb_system::Int #number of system
	nbPM::Int
    id_mod::Int # current maintenance model
	id_params::Int # current index in maintenance model parameters
	nb_params_maintenance::Int
    nb_params_family::Int
    nb_params_cov::Int
	mu::Int

	data::Vector{DataFrame}
	models::Vector{AbstractMaintenanceModel}
	family::FamilyModel
	maintenance_policy::AbstractMaintenancePolicy

	#Additional Covariates stuff
	data_cov::DataFrame
	params_cov::Vector{Float64}
	sum_cov::Float64 #to save the computation

	time::Vector{Float64}
	type::Vector{Int}

    indType::Float64

    Vleft::Float64
    Vright::Float64
    hVleft::Float64


	dVleft::Vector{Float64}
    dVright::Vector{Float64}
    
	d2Vleft::Vector{Float64}
    d2Vright::Vector{Float64}

	A::Float64
	dA::Vector{Float64}
	d2A::Vector{Float64}

	VR_prec::Vector{Float64}
	dVR_prec::Vector{Float64}
    d2VR_prec::Vector{Float64}

	comp::Compute

	Model() = new()
end

function init!(m::Model)
		m.k = 0  # current system
		m.nb_system = 1 #number of system
		if isdefined(m, :models) && length(m.models) > 0
			m.nbPM = length(m.models) - 1
		else
			m.nbPM = 0
		end
		m.id_mod = 0
		m.nb_params_maintenance=0
		for (id, mm) in enumerate(m.models)
			m.nb_params_maintenance += nb_params(mm)
		end
		m.nb_params_family = nb_params(m.family)
		m.nb_params_cov = 0
		m.mu = max_memory(m)

		m.data=DataFrame[]

		#Additional Covariates stuff
		# data_cov::DataFrame
		# params_cov::Vector{Float64}
		# sum_cov::Float64 #to save the computation

		## internal
		m.time = Float64[]
		m.type = Int[]

		m.indType = 0

		m.Vleft = 0
		m.Vright = 0
		m.hVleft = 0


		m.dVleft = zeros(m.nb_params_maintenance)
		m.dVright = zeros(m.nb_params_maintenance)
		
		nb2d = m.nb_params_maintenance * (m.nb_params_maintenance + 1) ÷ 2

		m.d2Vleft = zeros(nb2d)
		m.d2Vright = zeros(nb2d)

		m.A = 0
		m.dA = zeros(m.nb_params_maintenance)
		m.d2A = zeros(nb2d)
		if m.mu > 0
			m.VR_prec = zeros(m.mu)
			m.dVR_prec = zeros(m.mu * m.nb_params_maintenance)
			m.d2VR_prec = zeros(m.mu * nb2d)
		end
		m.comp = Compute(m)
		return nothing
end

nb_params(m::Model)::Int = m.nb_params_family + m.nb_params_maintenance + n.nb_params_cov
params(m::Model)::Vector{Float64} = cat(params(m.family),(map(m.models) do mm;params(mm); end)...,dims=1)
function params!(m::Model, p::Vector{Float64})
	from, to = 1, nb_params(m.family)
	params!(m.family,p[from:to])
	for mm in m.models
		if nb_params(mm) > 0
			from = to + 1
			to = from + nb_params(mm) - 1
			params!(mm, p[from:to])
		end
	end
end

function init_compute!(m::Model)
	init!(m.comp)
	for mm in m.models
		init!(mm)
	end
end

virtual_age(m::Model, x::Float64)::Float64 = m.Vright + (x  - m.time[m.k]) * m.A
virtual_age_inverse(m::Model, x::Float64) = (x - m.Vright) / m.A + m.time[m.k]

function update_Vleft!(m::Model;with_gradient::Bool=false, with_hessian::Bool=false)
	# /*if(model->k < 10) printf("Vleft:%lf\n", model->Vleft);*/
	m.Vleft = virtual_age(m, m.time[m.k + 1])
	# //printf("Vleft:%lf\n", model->Vleft);
	if with_hessian 
		for i in 0:(m.nb_paramsMaintenance - 1) 
            m.dVleft[i] = m.dVright[i] + (m.time[m.k+1] - m.time[m.k]) * m.dA[i]
            for j in 0:i
                m.d2Vleft[i * (i + 1) / 2 + j]= m.d2Vright[i * (i + 1) / 2 + j] + (m.time[m.k+1]  - m.time[m.k]) * m.d2A[i * (i + 1) / 2 + j]
			end
		end
	elseif with_gradient
		for i in 0:(m.nb_paramsMaintenance - 1)
            m.dVleft[i]= m.dVright[i] + (m.time[m.k+1]  - m.time[m.k]) * m.dA[i]
		end
	end
end

# void VamModel::set_data(List data_) {
# 	data=data_;
# 	nb_system=data.size();
# 	//printf("Number of systems: %d\n",nb_system);
# 	select_data(0);//default when only one system no need to
# }

# void VamModel::select_data(int i) {
# 	//In particular, if no data the following is skipped!
# 	if(data.size() > i) {
# 		List data2=data[i];
# 		//OLD CODE before gcc6: time = data2[0]; type = data2[1];//0 stand for Time and 1 for Type
# 		time = as<std::vector<double> >(data2[0]); type = as<std::vector<int> >(data2[1]); //Thanks to Lea, seems to be related to introduction of gcc Version 6 (see also https://github.com/apache/incubator-mxnet/issues/2185)
# 	}
# }

# DataFrame VamModel::get_selected_data(int i) {
# 	select_data(i);//Skipped if data is unset (see above)
# 	return DataFrame::create(_["Time"]=time,_["Type"]=type);
# };

function init_virtual_age_infos(m::Model)
		# int i;
    	# k=0;
    	# idMod=0; //id of current model
    	# S1 = 0;
    	# Vright=0;
    	# A=1;
    	# for(i=0;i<nbPM + 1;i++) models->at(i)->init();
end

function get_virtual_age_info(m::Model, from::Float64, to::Float64, by::Float64, expCov::Float64)
# 	double s=ceil((to-from)/by);
# 	int n=static_cast<int>(s);
# //printf("ici=%d,%lf (%lf,%lf,%lf)\n",n,s,to,from,by);
# 	std::vector<double> t(n+1);
# 	std::vector<double> v(n+1);
# 	std::vector<double> h(n+1); //i as intensity
# 	std::vector<double> H(n+1); //I for cumulative intensity
# 	std::vector<double> F(n+1); //F for conditional cumulative distribution function
# 	std::vector<double> S(n+1); //S for conditional survival function
# 	std::vector<double> f(n+1); //S for conditional survival function

# 	t[0]=from;t[n]=to;
# 	v[0]=virtual_age(from);v[n]=virtual_age(to);
# 	h[0]=expCov*A*family->hazardRate(v[0]);h[n]=expCov*A*family->hazardRate(v[n]);
# 	H[0]=S1;H[n]=S1+expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0]));
# 	F[0]=0;F[n]=1-exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	S[0]=1;S[n]=exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	f[0]=expCov*A*family->hazardRate(v[0]);f[n]=expCov*A*family->hazardRate(v[n])*exp(-expCov*(family->cumulative_hazardRate(v[n])-family->cumulative_hazardRate(v[0])));
# 	double by_t=(t[n]-t[0])/s;
# 	double by_v=(v[n]-v[0])/s;

# 	for(int i=1;i<n;i++) {
# 		t[i]=t[i-1]+by_t;//printf("t[%d]=%lf\n",i,t[i]);
# 		v[i]=v[i-1]+by_v;
# 		h[i]=expCov*A*family->hazardRate(v[i]);
# 		H[i]=S1+expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0]));
# 		F[i]=1-exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 		S[i]=exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 		f[i]=expCov*A*family->hazardRate(v[i])*exp(-expCov*(family->cumulative_hazardRate(v[i])-family->cumulative_hazardRate(v[0])));
# 	}

# 	return DataFrame::create(
# 		_["t"]=NumericVector(t.begin(),t.end()),
# 		_["v"]=NumericVector(v.begin(),v.end()),
# 		_["i"]=NumericVector(h.begin(),h.end()),
# 		_["I"]=NumericVector(H.begin(),H.end()),
# 		_["F"]=NumericVector(F.begin(),F.end()),
# 		_["S"]=NumericVector(S.begin(),S.end()),
# 		_["f"]=NumericVector(f.begin(),f.end())
# 	);
end

function get_virtual_age_infos(m::Model, from::Float64, to::Float64, by::Float64)

	# // Only one system first!
	# init_virtual_age_infos();
	# double expCov=1;
	# if (nb_paramsCov>0) expCov=exp(compute_covariates());
	# int n=time.size() - 1;
	# List res(n);
	# while(k < n) {
	# 	//printf("k=%d/n=%d,(%lf,%lf)\n",k,n,time[k],time[k+1]);
	# 	update_Vleft(false,false);
	# 	if(from > time[k] || time[k+1] > to ) res[k] = R_NilValue;
	# 	else res[k]=get_virtual_age_info(time[k],time[k+1],by,expCov);
	# 	S1 += expCov*(family->cumulative_hazardRate(Vleft) - family->cumulative_hazardRate(Vright));
	# 	//gradient_update_for_current_system();
	# 	int type2=type[k + 1];
	# 	if(type2 < 0) type2=0;
	# 	models->at(type2)->update(false,false);
	# }
	# return res;
end

# //Covariates related
# void VamModel::set_covariates(List model) {
# 	sum_cov=0.0;
# 	nb_paramsCov=0;
# 	if(model.containsElementNamed("covariates")) {
# 		List covariates_=model["covariates"];
# 		data_cov=covariates_["data"];
# 		params_cov=covariates_["params"];
# 		nb_paramsCov=params_cov.size();
# 	}
# }

function compute_covariates(m::Model)::Float64
	sum_cov=0.0
	# for j in 0:(m.nb_params_cov - 1)
	# 	# NumericVector var=data_cov[j];
	# 	# sum_cov += params_cov[j] * var[current_system];
	# 	# //printf("syst=%d,j=%d,th=%lf,params_cov=%lf\n",current_system,j,params_cov[j],var[current_system]);
	# end
	return sum_cov;
end

has_maintenance_policy(m::Model)::Bool = isdefined(m,:maintenance_policy) || isnothing(m.maintenance_policy)

function max_memory(m::Model)::Int
	maxmem = 1
	for mm in m.models
		if isdefined(mm,:m)
			if mm.m > maxmem
				maxmem = mm.m
			end
		end
	end
	return maxmem
end

function save_id_mod(m::Model, id_mod::Int)
	m.id_mod = id_mod
	if id_mod == 0
		m.id_params = 0
	else
		m.id_params += nb_params(m.models[id_mod]) 
	end
end