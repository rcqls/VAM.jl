mutable struct MaintenanceModel
    model::Model
    id::Int
    id_params::Int
end


    # MaintenanceModel(VamModel* model_) {
    # 	model = model_;
    # }

    # virtual ~MaintenanceModel() {};

    # virtual NumericVector get_params() = 0;

    # virtual  void set_params(NumericVector par, int ind) = 0;//ind indicates the indice of vector par at which the parameters to set begin

    # virtual  void init() = 0;

    # virtual int nb_params() = 0;

    # virtual void update(bool with_gradient,bool with_hessian) = 0;



    # void set_id(int id_) {
    # 	id=id_;
    # }

    # void set_id_params(int id_params_) {
    #     id_params=id_params_;
    # }

mutable struct ARA1 <: MaintenanceModel
    ρ::Float64
end

    # ARA1(NumericVector rho_,VamModel* model_) : MaintenanceModel(model_) {
    # 	rho = rho_[0];
    # }

    # NumericVector get_params() {
    # 	NumericVector out(1);
    # 	out[0]=rho;
    # 	return out;
    # }

    # void set_params(NumericVector par, int ind) {
    # 	rho=par[ind];
    # }

    # void init(){
    # }

    # int nb_params(){
    #     return 1;
    # }

    # void update(bool with_gradient,bool with_hessian);

mutable struct ARAInf <: MaintenanceModel 
	ρ::Float64
end

    # ARAInf(NumericVector rho_,VamModel* model_) : MaintenanceModel(model_) {
    # 	rho = rho_[0];
    # }

    # NumericVector get_params() {
    # 	NumericVector out(1);
    # 	out[0]=rho;
    # 	return out;
    # }

    # void set_params(NumericVector par, int ind) {
    # 	rho=par[ind];
    # }

    # void init(){
    # }

    # int nb_params(){
    #     return 1;
    # }

    # void update(bool with_gradient,bool with_hessian);


struct AGAN <: MaintenanceModel
end
# public:

#     AGAN(VamModel* model_) : MaintenanceModel(model_) {
#     }

#     NumericVector get_params() {
#         NumericVector out(0);
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 0;
#     }

#     void update(bool with_gradient,bool with_hessian);

# };

struct ABAO <: MaintenanceModel
end

# public:

#     ABAO(VamModel* model_) : MaintenanceModel(model_) {
#     }

#     NumericVector get_params() {
#         NumericVector out(0);
#         return out;
#     }

#     void set_params(NumericVector par,int ind) {
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 0;
#     }

#     void update(bool with_gradient,bool with_hessian);

# };

struct AGAP <: MaintenanceModel
end

# public:

#     AGAP(VamModel* model_) : MaintenanceModel(model_) {
#     }

#     NumericVector get_params() {
#         NumericVector out(0);
#         return out;
#     }

#     void set_params(NumericVector par,int ind) {
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 0;
#     }

#     void update(bool with_gradient,bool with_hessian);

# };

mutable struct QAGAN <: MaintenanceModel
    ρ::Float
end

# public:

#     QAGAN(VamModel* model_) : MaintenanceModel(model_) {
#     }

#     NumericVector get_params() {
#         NumericVector out(0);
#         return out;
#     }

#     void set_params(NumericVector par,int ind) {
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 0;
#     }

#     void update(bool with_gradient,bool with_hessian);

# };

mutable struct QR <: MaintenanceModel
    ρ::Float64
end

# public:

#     QR(NumericVector rho_,VamModel* model_) : MaintenanceModel(model_) {
#         rho = rho_[0];
#     }

#     NumericVector get_params() {
#         NumericVector out(1);
#         out[0]=rho;
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#         rho=par[ind];
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 1;
#     }

#     void update(bool with_gradient,bool with_hessian);

struct F_GQR
end
# public:
#     virtual double eval(double x) = 0;
# };

struct ID_GQR <:  F_GQR
end
# public:
#     double eval(double x){
#         return x;
#     }

# };

struct LOG_GQR <: F_GQR
end
# public:
#     double eval(double x){
#         return log(x+1);
#     }

# };

struct SQRT_GQR <: F_GQR
end
# public:
#     double eval(double x){
#         return sqrt(x);
#     }

# };

mutable struct GQR <: MaintenanceModel
    ρ::Float64
    K::Float64
    f::F_GQR
end
# public:

#     GQR(NumericVector rho_, std::string extra, VamModel* model_) : MaintenanceModel(model_) {
#         rho = rho_[0];
#         K=0;
#         if(extra.compare("identity")==0){
#             f=new id_GQR();
#         } else if(extra.compare("log")==0){
#             f=new log_GQR();
#         }  else if(extra.compare("sqrt")==0){
#             f=new sqrt_GQR();
#         } else {
#             std::cout<<"Undefined argument"<< extra<< "for GQR model: replaced by identity function\n";
#             f=new id_GQR();
#         }
#     }

#     NumericVector get_params() {
#         NumericVector out(1);
#         out[0]=rho;
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#         rho=par[ind];
#         K=0;
#     }

#     void init(){
#         K=0;
#     }

#     int nb_params(){
#         return 1;
#     }

#     void update(bool with_gradient,bool with_hessian);

mutable struct GQR_ARA1 <:  MaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
    # GQR_ARA1(NumericVector rho_, std::string extra, VamModel* model_) : MaintenanceModel(model_) {
    #     rho_QR = rho_[0];
    #     rho_ARA = rho_[1];
    #     K=0;
    #     if(extra.compare("identity")==0){
    #         f=new id_GQR();
    #     } else if(extra.compare("log")==0){
    #         f=new log_GQR();
    #     }  else if(extra.compare("sqrt")==0){
    #         f=new sqrt_GQR();
    #     } else {
    #         std::cout<<"Undefined argument"<< extra<< "for GQR model: replaced by identity function\n";
    #         f=new id_GQR();
    #     }
    # }

    # NumericVector get_params() {
    #     NumericVector out(2);
    #     out[0]=rho_QR;
    #     out[1]=rho_ARA;
    #     return out;
    # }

    # void set_params(NumericVector par, int ind) {
    #     rho_QR=par[ind];
    #     rho_ARA=par[ind+1];
    #     K=0;
    # }

    # void init(){
    #     K=0;
    # }

    # int nb_params(){
    #     return 2;
    # }

    # void update(bool with_gradient,bool with_hessian);


mutable struct GQR_ARAInf <:  MaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
end
# public:

#     GQR_ARAInf(NumericVector rho_, std::string extra, VamModel* model_) : MaintenanceModel(model_) {
#         rho_QR = rho_[0];
#         rho_ARA = rho_[1];
#         K=0;
#         if(extra.compare("identity")==0){
#             f=new id_GQR();
#         } else if(extra.compare("log")==0){
#             f=new log_GQR();
#         }  else if(extra.compare("sqrt")==0){
#             f=new sqrt_GQR();
#         } else {
#             std::cout<<"Undefined argument"<< extra<< "for GQR model: replaced by identity function\n";
#             f=new id_GQR();
#         }
#     }

#     NumericVector get_params() {
#         NumericVector out(2);
#         out[0]=rho_QR;
#         out[1]=rho_ARA;
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#         rho_QR=par[ind];
#         rho_ARA=par[ind+1];
#         K=0;
#     }

#     void init(){
#         K=0;
#     }

#     int nb_params(){
#         return 2;
#     }

#     void update(bool with_gradient,bool with_hessian);

mutable struct ARAm <: MaintenanceModel
    ρ::Float64
    m::Int
end

# public:

#     ARAm(NumericVector rho_, int m_, VamModel* model_) : MaintenanceModel(model_) {
#         rho = rho_[0];
#         m = m_;
#     }

#     NumericVector get_params() {
#         NumericVector out(1);
#         out[0]=rho;
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#         rho=par[ind];
#     }

#     void init(){
#     }

#     int nb_params(){
#         return 1;
#     }

#     void update(bool with_gradient,bool with_hessian);


mutable struct GQR_ARAm <: MaintenanceModel
    ρQR::Float64
    ρARA::Float64
    K::Float64
    f::F_GQR
    m::Int
end

# public:

#     GQR_ARAm(NumericVector rho_, std::string extra, int m_, VamModel* model_) : MaintenanceModel(model_) {
#         rho_QR = rho_[0];
#         rho_ARA = rho_[1];
#         m = m_;
#         K=0;
#         if(extra.compare("identity")==0){
#             f=new id_GQR();
#         } else if(extra.compare("log")==0){
#             f=new log_GQR();
#         }  else if(extra.compare("sqrt")==0){
#             f=new sqrt_GQR();
#         } else {
#             std::cout<<"Undefined argument"<< extra<< "for GQR model: replaced by identity function\n";
#             f=new id_GQR();
#         }
#     }

#     NumericVector get_params() {
#         NumericVector out(2);
#         out[0]=rho_QR;
#         out[1]=rho_ARA;
#         return out;
#     }

#     void set_params(NumericVector par, int ind) {
#         rho_QR=par[ind];
#         rho_ARA=par[ind+1];
#         K=0;
#     }

#     void init(){
#         K=0;
#     }

#     int nb_params(){
#         return 2;
#     }

#     void update(bool with_gradient,bool with_hessian);

# private:
#     double rho_QR;
#     double rho_ARA;
#     int m;
#     double K;
#     f_GQR *f;
# };

# MaintenanceModel* newMaintenanceModel(List maintenance,VamModel* model);
