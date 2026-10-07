try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE)
suppressPackageStartupMessages(library(lme4))
fits <- list()
for (id in c("p5p95_lower","p5p95_upper","sire_subset_base_upper")) {
  path <- paste0("results/models/",id,".rds")
  original_path <- paste0("results/models/",id,"_initial_unstable.rds")
  m <- readRDS(if(file.exists(original_path)) original_path else path)
  dat <- model.frame(m)
  cat("Checking",id,"using adaptive quadrature and alternative starts\n")
  agq <- glmer(formula(m),data=dat,family=binomial(),nAGQ=9,
        control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=200000)))
  crosscheck <- glmer(formula(m),data=dat,family=binomial(),nAGQ=9,
         start=list(theta=getME(agq,"theta"),fixef=fixef(agq)),
         control=glmerControl(optimizer="nloptwrap",
                  optCtrl=list(maxeval=200000,ftol_abs=1e-10,xtol_abs=1e-8)))
  fit <- agq
  delta <- max(abs(fixef(fit)-fixef(crosscheck)))
  stopifnot(is.null(fit@optinfo$conv$lme4$messages),fit@optinfo$conv$opt==0,delta < .01)
  file.copy(path,paste0("results/models/",id,"_initial_unstable.rds"),overwrite=FALSE)
  saveRDS(agq,paste0("results/models/",id,"_AGQ9_checked.rds"))
  saveRDS(crosscheck,paste0("results/models/",id,"_AGQ9_secondary_optimizer.rds"))
  saveRDS(fit,path)
  fits[[id]] <- data.frame(model=id,initial_logLik=as.numeric(logLik(m)),
       initial_variance=as.data.frame(VarCorr(m))$vcov[1],
       selected_logLik=as.numeric(logLik(fit)),
       selected_variance=as.data.frame(VarCorr(fit))$vcov[1],
       AGQ9_variance=as.data.frame(VarCorr(agq))$vcov[1],
       selected_nAGQ=9,
       max_fixed_effect_change_between_optimizers=delta,
       secondary_optimizer_messages=paste(crosscheck@optinfo$conv$lme4$messages,collapse="; "),
       messages=paste(fit@optinfo$conv$lme4$messages,collapse="; "))
  print(fits[[id]])
}
write.csv(do.call(rbind,fits),"results/numerical_checks.csv",row.names=FALSE)
