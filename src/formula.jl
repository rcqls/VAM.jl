# Time & Type ~ (ABAO()|Weibull(1,3)
# Time & Type ~ (ARA1(~Beta(1.08,0.108)) | Weibull(~NonInform(),~Gamma(32,0.097)))

# ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4))))

# Systeme & Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4)))

function parse_model(ex_f::Expr)
    m = Model()
    if Meta.isexpr(ex_f, :call)
        ex_m = ex_f
        def_names = ["Time", "Type"] # default names
        ## model detection
        if ex_f.args[1] == :&
            ## No names given on the left side of :~ => def_names used
            ex_m.args[2] = ex_m.args[2].args[2] # remove unused the tilde
        elseif ex_f.args[1] == :~
            ## Left part is names and right part is the model
            def_names = map(ex_f.args[2].args[2:3]) do n
                string(n)
            end
            ex_m = ex_f.args[3]
        end
        ## parsing model (ex_m)
        m.models = AbstractMaintenanceModel[]
        #print(ex_m.args)
        if ex_m.args[2].args[1] == :|
            ex_cm = ex_m.args[2]
            parse_cm(m, ex_cm)
            ## PMs (Preventive Maintenances) and MPs (Maintenance Policies)
            ex_pm = ex_m.args[3]
            if ex_pm.args[1] == :|
                # PMs
                ex_pms = ex_pm.args[2]
                if ex_pms.args[1] == :+
                    # several PMs
                    for pm in ex_pms.args[2:end]
                        push!(m.models,eval(pm))
                    end
                else
                    # only 1 PM
                    push!(m.models,eval(ex_pms))
                end
                # Maintenance policies
                ex_mps = ex_pm.args[3]
                if ex_mps.args[1] == :*
                    # several MPs
                    maintenance_policies = AbstractMaintenancePolicy[]
                    for mp in ex_mps.args[2:end]
                        push!(maintenance_policies,eval(complete_name!(mp,1,"MaintenancePolicy")))
                    end
                    m.maintenance_policy = MaintenancePolicyList(maintenance_policies)
                else
                    # only 1 MP
                    m.maintenance_policy = eval(complete_name!(ex_mps,1,"MaintenancePolicy"))
                end
            end
        else
            ##No PM
            parse_cm(m, ex_m)
        end
        return m
    end
end

function parse_cm(m::AbstractModel,ex_cm::Expr)
    if ex_cm.args[1] == :|
        push!(m.models,eval(ex_cm.args[2])) # CM (Corrective Maintenance)
        fm = ex_cm.args[3]
        #fm.args[1] = complete_name(fm.args[1], "FamilyModel")
        if isa(fm.args[2], Expr) && fm.args[2].head == :parameters
            # covariates
            fm2=Expr(:call, fm.args[1], fm.args[3:end]..., fm.args[2].args[1].args)
            println(fm2)
            m.family = eval(complete_name!(fm2, 1, "FamilyModel"))
        else
            m.family = eval(complete_name!(fm, 1, "FamilyModel"))
        end
    end
end

macro vam(ex_f)
    return parse_model(ex_f)
end

## TO REMOVE
# macro sim(ex_f, ex_s)
#     mod = parse_model(ex_f)
#     sim = Sim(mod, formula_translate(ex_s))
#     init!(sim.model)
#     return sim
# end

# macro sim(ex_f)
#     println(ex_f.args)
#     mod = parse_model(ex_f)
#     sim = Sim(mod, nothing)
#     init!(sim.model)
#     return sim
# end

macro stop(ex_s)
    return formula_translate(ex_s).args
end

function formula_translate(ex_f::Expr)
    Meta.parse(replace(string(ex_f), "size" => "s", "time" => "t"))
end

function complete_name!(ex::Expr, i::Int, append::String)::Expr
    ex.args[i] = Symbol(string(ex.args[i]) * append)
    return ex
end

