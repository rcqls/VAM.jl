abstract type StopPolicy end

# StopPolicy* newStopPolicy(SimVam* sim,List policy);

struct SizeGreaterThanStopPolicy <: StopPolicy
    size::Int
end
    # public:
#     SizeGreaterThanStopPolicy(SimVam* sim_,int size_): StopPolicy(sim_) {
#         size=size_;
#     }

#     ~SizeGreaterThanStopPolicy() {};

    # bool ok();

	# 	void first() {};

    # int size;


# };

# //M is for CM or PM, type=-1 or 1,2,...

struct SizeOfTypeGreaterThanStopPolicy <:StopPolicy
    size::Int
    type::Int
    count::Int
end
#     public:
#     SizeOfTypeGreaterThanStopPolicy(SimVam* sim_,int type_,int size_): StopPolicy(sim_) {
#         type=type_;
#         size=size_;
#         count=0;
#     }

#     ~SizeOfTypeGreaterThanStopPolicy() {};

#     bool ok();

# 		void first() {};

#     int size;

#     int type;

#     int count;

# };

struct TimeGreaterThanCensorshipStopPolicy <: StopPolicy
    to_init::Bool
    time::Float64
    #expr::Language
    #env::Environment
end
#     public:
#     TimeGreaterThanCensorshipStopPolicy(SimVam* sim_,double time_,Language expr_,Environment env_): StopPolicy(sim_) {
#         time=time_;expr=expr_;env=env_;
# 				to_init=(time<0);
#     }

#     ~TimeGreaterThanCensorshipStopPolicy() {};

#     bool ok();

# 		void first();

# 		bool to_init;
		
#     double time;

# 		Language expr;

# 		Environment env;


# };

struct TimeGreaterThanStopPolicy <: StopPolicy
    time::Float64
end
# public:
#     TimeGreaterThanStopPolicy(SimVam* sim_,double time_): StopPolicy(sim_) {
#         time=time_;
#     }

#     ~TimeGreaterThanStopPolicy() {};

#     bool ok();

# 		void first() {};

#     double time;


# };

struct AndStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end

# public:
#     AndStopPolicy(SimVam* sim_,List policies_): StopPolicy(sim_) {
#         for(
#             List::iterator it=policies_.begin();
#             it != policies_.end();
#             ++it
#         ) {
#             List policy=*it;
#             StopPolicy*  sp=newStopPolicy(sim_,policy);
#             if(!(sp == NULL)) policies.push_back(sp);
#         }
#     }

#     ~AndStopPolicy() {
#         for(
#             std::vector<StopPolicy*>::iterator it=policies.begin();
#             it != policies.end();
#             ++it
#         ) {
#             delete *it;
#         }
#     };

#     bool ok();

# 		void first();

#     std::vector<StopPolicy*> policies;


# };

struct OrStopPolicy <: StopPolicy
    policies::Vector{StopPolicy}
end
# public:
#     OrStopPolicy(SimVam* sim_,List policies_): StopPolicy(sim_) {
#         for(
#             List::iterator it=policies_.begin();
#             it != policies_.end();
#             ++it
#         ) {
#             List policy=*it;
#             StopPolicy*  sp=newStopPolicy(sim_,policy);
#             if(!(sp == NULL)) policies.push_back(sp);
#         }
#     }

#     ~OrStopPolicy() {
#         for(
#             std::vector<StopPolicy*>::iterator it=policies.begin();
#             it != policies.end();
#             ++it
#         ) {
#             delete *it;
#         }
#     };

#     bool ok();

# 		void first();

#     std::vector<StopPolicy*> policies;


# };