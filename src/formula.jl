# Time & Type ~ (ABAO()|Weibull(1,3)
# Time & Type ~ (ARA1(~Beta(1.08,0.108)) | Weibull(~NonInform(),~Gamma(32,0.097)))

# ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4))))

# Systeme & Temps & Type ~ (ARA1(.5) | Weibull(0.01,2.5)) & (ARAInf(.7)+ARAInf(.3)|Periodic(12,prob=c(0.6,0.4)))

macro model(ex_f)
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
            push!(m.models,eval(ex_m.args[2].args[2])) # corrective maintenance
            fm = ex_m.args[2].args[3]
            fm.args[1] = Symbol(string(fm.args[1]) * "FamilyModel")
            m.family = eval(fm)
        end 
        ## 
        return m
        # m.models.push(ex_f.args[])
    end
    # return m
end

# parse_formula()

# parse.vam.formula <- function(formula) {
# 	## Needs to have this envir to evaluate params (otherwise, beta was found in baseenv() first before globalenv() for example!)
# 	envir.eval <- parent.frame(4)
# 	eval.vam <- function(e) eval(e,envir.eval)
# 	if(formula[[1]] != as.name("~")) stop("Argument has to be a formula")
# 	if(length(formula) == 2) {
# 		response <- NULL
# 		cm <- formula[[2]]
# 	} else {
# 		tmp <- formula[[2]]
# 		## simplify parenthesis
# 		while(tmp[[1]] == as.name("(")) tmp <- tmp[[2]]
# 		if(tmp[[1]] != as.name("&") && length(tmp) != 3) stop("Left part of formula of the form 'Time & Type'!")
# 		if(length(tmp[[2]])==3 && tmp[[2]][[1]]==as.name("&")) {
# 			response <- c(as.character(tmp[[2]][[2]]),as.character(tmp[[2]][[3]]),as.character(tmp[[3]]))
# 		} else response <- c(as.character(tmp[[2]]),as.character(tmp[[3]]))
# 		cm <- formula[[3]]
# 	}
# 	## simplify parenthesis
# 	while(cm[[1]] == as.name("(")) cm <- cm[[2]]
# 	pms <- list()
# 	policy <- NULL
# 	if(there.is.pm <- (cm[[1]] == as.name("&"))) { # there is a PM part
# 		pm <- cm[[3]]
# 		cm <- cm[[2]]
# 		# deal with PM part
# 		if(pm[[1]] == as.name("(")) {
# 			pm <- pm[[2]]
# 			if(pm[[1]] != as.name("|")) {
# 				## Case: No maintenance policy
# 				#stop("Need a policy to manage Preventive Maintenance")
# 				policy <- NULL
# 			} else {
# 				policy <- pm[[3]]
# 				if(policy[[1]] == as.name("*")) {
# 					## Case: Composition of maintenance policies
# 					# recursive function to detect maintenance policies
# 					run.over.policies<-function(p) {
# 						if(p[[1]] == as.name("*")) {
# 							run.over.policies(p[[2]])
# 							run.over.policies(p[[3]])
# 						} else if(is.name(p[[1]])) {
# 							p[[1]] <- as.name(paste0(as.character(p[[1]]),".maintenance.policy"))
# 							policies <<- c(policies,list(p))
# 						}
# 					}
# 					## init policies and
# 					policies <- list()
# 					run.over.policies(policy)
# 					## print(policies)
# 					policy <- policies ##[[1]]
# 				} else if(is.name(policy[[1]])) {
# 					## Case: One maintenance policy
# 					policy[[1]] <- as.name(paste0(as.character(policy[[1]]),".maintenance.policy"))
# 				}

# 				# PMs
# 				pm <- pm[[2]]
# 			}
# 			# parser for pm
# 			parse.pm <- function(pm) {
# 				if(is.name(pm[[1]])) {
# 					pm[[1]] <- as.name(paste0(as.character(pm[[1]]),".va.model"))
# 				}
# 				pm
# 			}
# 			cpt.pms <- 0
# 			while(pm[[1]] == as.name("+") ) {
# 				if(length(pm) == 3) {
# 					pms[[cpt.pms <- cpt.pms + 1]] <- parse.pm(pm[[3]])
# 					pm <- pm[[2]]
# 				}
# 			}
# 			pms[[cpt.pms <- cpt.pms + 1]] <- parse.pm(pm)
# 		} else stop("Need parenthesis around the Preventive Maintenance terms")
# 	}
# 	# deal with CM PART
# 	cms <- list()

# 	# parser for cm
# 	parse.cm <- function(cm) {
# 		# print(there.is.pm)
# 		# print(cm)
# 		if(there.is.pm) {
# 			if(cm[[1]] == as.name("(")) cm <- cm[[2]]
# 			else stop("CM needs a family!")
# 		}
# 		if(cm[[1]] != as.name("|")) stop("CM needs a family!")
# 		family <- cm[[3]]
# 		if(is.name(family[[1]])) {
# 			family[[1]] <- as.name(paste0(as.character(family[[1]]),".family.cm"))
# 		}
# 		cm <- cm[[2]]
# 		if(is.name(cm[[1]])) {
# 			cm[[1]] <- as.name(paste0(as.character(cm[[1]]),".va.model"))
# 		}
# 		list(model=cm,family=family)
# 	}
# 	cpt.cms <- 0
# 	while( cm[[1]] == as.name("+") ) {
# 		if(length(cm) == 3) {
# 			cms[[cpt.cms <- cpt.cms + 1]] <- parse.cm(cm[[3]])
# 			cm <- cm[[2]]
# 		}
# 	}
# 	cms[[cpt.cms <- cpt.cms + 1]] <- parse.cm(cm)

# 	## Parse covariates
# 	parse.covariates <- function(expr) {
# 		form<-list()
# 		params <- c()
# 		add_term <- function(term,sign) {
# 			##print(term)
# 			if(term[[1]]==as.name("*")) {
# 				form <<- c(as.character(term[[3]]),form)
# 				if(sign == as.name("-")) {
# 				 	if(!(is.numeric(eval(term[[2]])))) stop("Only + is admitted between covariates definition in Bayes case.") else param_expr <- eval.vam(parse(text=paste0(sign,as.character(eval(term[[2]])))))
# 				} else param_expr <- as.vector(eval.vam(term[[2]]))
# 				params<<- c(param_expr,params)
# 			}
# 			##print(list(form=form,params=params))
# 		}
# 		while(expr[[1]]==as.name("+") || expr[[1]]==as.name("-")) {
# 			add_term(expr[[3]],as.character(expr[[1]]))
# 			expr <- expr[[2]]
# 		}
# 		add_term(expr,"")
# 		list(formula=eval(parse(text=paste0("~",paste(form,collapse="+")))),params=params)
# 	}

# 	convert.family <- function(fam) {
# 		# eval.vam is here to evaluate the value if it is a symbol!
# 		# We want to detect if there is covariates or not. Normally the operator | delimitating covariates has priority, expect possibly in Baysian case.
# 		if(has.covariates <- (length(fam[[length(fam)]])>1 && fam[[length(fam)]][[1]] == as.name("|"))) {
# 			covariates_expr <- fam[[length(fam)]][[3]]
# 			fam[[length(fam)]] <- fam[[length(fam)]][[2]] # first argument of last terms becomes last argument of family
# 		} else {
# 		#In Bayesian case the ~ operator of the description of the prior distribution of the last argument of the family can possibly has priority on the operator | delimitating covariates
# 			if(length(fam[[length(fam)]])>1 && fam[[length(fam)]][[1]] == as.name("~") && length(fam[[length(fam)]][[2]])>1 && fam[[length(fam)]][[2]][[1]] == as.name("|")){
# 				has.covariates <- TRUE
# 				covariates_expr <- fam[[length(fam)]][[2]][[3]]
# 				fam[[length(fam)]][[2]] <- fam[[length(fam)]][[2]][[2]]
# 			}
# 		}

# 		res<-list(
# 				name=as.character(fam[[1]]),
# 				params=sapply(fam[-1],function(e) as.vector(eval.vam(e)))
# 				## instead of : params=sapply(cm$family[-1],as.vector)
# 				## which does not work with negative real since element of tmp[-1] interpreted as call!
# 		)
# 		if(has.covariates) {
# 			res$covariates <- parse.covariates(covariates_expr)
# 		}
# 		return(res)
# 	}
# 	convert.pm <- function(pm) {
# 		n_pip<-c()
# 		if(length(pm)>1){
# 			for(i in 2:length(pm)){
# 				if((length(pm[[i]])==3)&&(pm[[i]][[1]]==as.name("|"))) {
# 					n_pip<-c(n_pip,i)
# 				}
# 			}
# 		}
# 		if(length(n_pip)==0) {
# 			list(
# 				name=as.character(pm[[1]]),
# 				params=as.vector(if(length(pm)==1) numeric(0) else sapply(pm[2:length(pm)],function(e) as.vector(eval.vam(e))))
# 			)
# 		} else if(length(n_pip)==1) {
# 			if(n_pip<(length(pm)-1)) {
# 				stop("Maximum two arguments after a | in a maintenance effect!")
# 			} else if(n_pip == length(pm)) {
# 				if( typeof(tryCatch( as.double(eval.vam(pm[[length(pm)]][[3]])) ,error=function(e){FALSE},finally=function(e){TRUE}))!="logical"){
# 	  				if((round(eval.vam(pm[[length(pm)]][[3]])) != eval.vam(pm[[length(pm)]][[3]]))||(round(eval.vam(pm[[length(pm)]][[3]]))<=0)) {
# 	  					stop("Memory argument of a maintenance model has to be a strictly positive integer!")
# 	  				} else {
# 	  	  				list(
# 									name=as.character(pm[[1]]),
# 									params=as.vector(if(length(pm)==2) as.vector(eval.vam(pm[[2]][[2]])) else c(sapply(pm[2:(length(pm)-1)],function(e) as.vector(eval.vam(e))),as.vector(eval.vam(pm[[length(pm)]][[2]])))),
# 									m=as.integer(eval.vam(pm[[length(pm)]][[3]]))
# 		  					)
# 						}
# 	  			} else {
# 	  				list(
# 							name=as.character(pm[[1]]),
# 							params=as.vector(if(length(pm)==2) as.vector(eval.vam(pm[[2]][[2]])) else c(sapply(pm[2:(length(pm)-1)],function(e) as.vector(eval.vam(e))),as.vector(eval.vam(pm[[length(pm)]][[2]])))),
# 							extra=as.character(pm[[length(pm)]][[3]])
# 						)
# 	  			}
# 			} else {
# 				if( typeof(tryCatch( as.double(eval.vam(pm[[length(pm)-1]][[3]])) ,error=function(e){FALSE},finally=function(e){TRUE}))!="logical"){
#   				if((round(eval.vam(pm[[length(pm)-1]][[3]]))!=eval.vam(pm[[length(pm)-1]][[3]]))||(round(eval.vam(pm[[length(pm)-1]][[3]]))<0)) {
#   					stop("Memory argument of a maintenance model has to be a positive integer!")
#   				} else {
#   	  				list(
# 								name=as.character(pm[[1]]),
# 								params=as.vector(if(length(pm)==3) as.vector(eval.vam(pm[[2]][[2]])) else c(sapply(pm[2:(length(pm)-2)],function(e) as.vector(eval.vam(e))),as.vector(eval.vam(pm[[length(pm)-1]][[2]])))),
# 								m=as.integer(eval.vam(pm[[length(pm)-1]][[3]])),
# 								extra=as.character(pm[[length(pm)]])
# 	  					)
# 					}
# 				} else {
# 					if( typeof(tryCatch( as.double(eval.vam(pm[[length(pm)]])) ,error=function(e){FALSE},finally=function(e){TRUE}))=="logical"){
# 						stop("At least one of the two argument of maintenance model after a | must be a memory that is to say a non negative positive integer!")
# 					} else {
# 						if((round(eval.vam(pm[[length(pm)]]))!=eval.vam(pm[[length(pm)]]))||(round(eval.vam(pm[[length(pm)]]))<0)) {
# 	  						stop("Memory argument of a maintenance model has to be a positive integer!")
# 	  					} else {
# 	  	  					list(
# 								name=as.character(pm[[1]]),
# 								params=as.vector(if(length(pm)==3) as.vector(eval.vam(pm[[2]][[2]])) else c(sapply(pm[2:(length(pm)-2)],function(e) as.vector(eval.vam(e))),as.vector(eval.vam(pm[[length(pm)-1]][[2]])))),
# 								m=as.integer(eval.vam(pm[[length(pm)]])),
# 								extra=as.character(pm[[length(pm)-1]][[3]])
# 		  					)
# 	  	  				}
# 					}

# 				}
# 			}
# 		} else {
# 			stop("Maximum one | in a maintenance effect!")
# 		}



# 	 #  if((length(pm)==1)||(pm[[length(pm)]][[1]]!=as.name("|"))) {
# 		# list(
# 		# 	name=as.character(pm[[1]]),
# 		# 	params=as.vector(if(length(pm)==1) numeric(0) else sapply(pm[2:length(pm)],function(e) as.vector(eval(e))))
# 		# )
# 	 #  } else if ( typeof(tryCatch( as.double(eval(pm[[length(pm)]][[3]])) ,error=function(e){FALSE},finally=function(e){TRUE}))!="logical"){
# 	 #  	if((round(eval(pm[[length(pm)]][[3]]))!=eval(pm[[length(pm)]][[3]]))||(round(eval(pm[[length(pm)]][[3]]))<0)) {
# 	 #  		stop("Memory argument of a maintenance model has to be a positive integer!")
# 	 #  	} else {
# 	 #  	  list(
# 		# 	name=as.character(pm[[1]]),
# 		# 	params=as.vector(if(length(pm)==2) pm[[2]][[2]] else c(sapply(pm[2:(length(pm)-1)],function(e) as.vector(eval(e))),as.vector(eval(pm[[length(pm)]][[2]])))),
# 		# 	m=as.integer(eval(pm[[length(pm)]][[3]]))
# 		#   )
# 		# }
# 	 #  }	else {
# 	 #  	list(
# 		# 	name=as.character(pm[[1]]),
# 		# 	params=as.vector(if(length(pm)==2) pm[[2]][[2]] else c(sapply(pm[2:(length(pm)-1)],function(e) as.vector(eval(e))),as.vector(eval(pm[[length(pm)]][[2]])))),
# 		# 	extra=as.character(pm[[length(pm)]][[3]])
# 		# )
# 	 #  }
# 	}
# 	convert.mp <- function(mp) {#maintenance policy
# 		if(is.null(mp)) list(name="None")
# 		else if(is.list(mp)) {
# 			list(name="MaintenancePolicyList",policies=lapply(mp,convert.mp))
# 		}
# 		else {

# 			## The function defining the maintenance policy
# 			## (registered in maintenance-policy-register.R or in any other R file)
# 			mp.fct <- eval.vam(mp[[1]])
# 			## params used in the call mp
# 			pars <- as.list(match.call(mp.fct,mp))[-1]

# 			## Default values are then automatically completed using declaration of maintenance policy
# 			pars.default <- (as.list(mp.fct)->tmp)[-length(tmp)]
# 			pars.default <- pars.default[sapply(pars.default,function(e) nchar(as.character(e)))!=0]
# 			for(e in names(pars.default)) if(is.null(pars[[e]])) pars[[e]] <- pars.default[[e]]

# 			##print(list(pars=pars))

# 			## deal with model parameter which has a specific treatment
# 			mod <- NULL
# 			if(!is.null(pars[["model"]])) {
# 				mod <- rcpp(eval.vam(pars[["model"]]))
# 				pars[["model"]] <- NULL
# 			}

# 			res <- list(
# 				name=as.character(mp[[1]]),
# 				params=lapply(pars,eval.vam)
# 			)
# 			res[["with.model"]] <- !is.null(mod)
# 			if(!is.null(mod)) res[["model"]] <- mod
# 			res
# 		}
# 	}


# 	res<-list(
# 		response=response,
# 		models=c(list(convert.pm(cms[[1]]$model)),lapply(pms[rev(seq(pms))],convert.pm)),
# 		family=convert.family(cms[[1]]$family),
# 		pm.policy=convert.mp(policy)
# 	)
# 	## covariates direct acces
# 	res$covariates <- res$family$covariates
# 	res$family$covariates <- NULL

# 	res$max_memory <- max(1,unlist(sapply(res$models,function(e) e$m)),na.rm=TRUE)
# 	res

# }