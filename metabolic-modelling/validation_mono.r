library(data.table)
library(readxl)
library(cobrar)
library(ggplot2)
library(stringr)

setwd("src/OMM12-modelling_github/")

model.file1      <- "dat/OMM_munich-aam_20241205.RDS"
model.file2      <- "dat/OMM_munich-csm-aam_20241205_xml/"
model.file3      <- "~/uni/OMM_munich/dat/OMM_munich-rfastcormics2-aam_20261009.RDS"
validation.file <- "dat/OMM12-uptake-product-matrix.xlsx"
medium.file     <- "dat/AAM_nomuc_rev5.csv"


check_growth <- function(model, medium=NULL, gs.medium=NULL, obj=NULL, uptake=-100, printCosts=F, printExchanges=T, returnSol=F, useNames=T, verbose=TRUE, mtf=FALSE, cutoff=1e-6, solver=NA, fva=F, sampling.dat=NULL, sampling.cor=0.99, cffba=F, full=F, returnEx=F, returnMod=F, uncertainty=F, cutoff.print=0, normalize=F){
  if( full ){
  if( cobrar::COBRAR_SETTINGS("SOLVER") == "cplex" ) solver_ok=1 else if( cobrar::COBRAR_SETTINGS("SOLVER") == "glpk" ) solver_ok=5
    fva=T; mtf=T; printCosts=T; useNames=T
  }
  ex.idx <- grep("EX_", model@react_id)
  if(is.null(medium)){
    medium <- model@react_id[ex.idx][which(model@lowbnd[ex.idx]<0)]
    uptake <- model@lowbnd[ex.idx][which(model@lowbnd[ex.idx] < 0)]
  } 
  if(length(uptake) == 1) uptake <- rep(uptake, length(medium))
  if(!is.null(gs.medium)){
    medium <- paste0("EX_",gs.medium$compounds,"_e0")
    uptake <- -gs.medium$maxFlux
    uptake <- ifelse(uptake < -1000, -1000, uptake)
  }
  if(!is.null(obj)){
    model <- changeObjFunc(model, react=obj)
    stop("not missing cobrar")
  }else if( all(grepl("^cpd", medium) )){
    medium <- paste0("EX_", medium, "_e0")
  }
  ub <- model@uppbnd
  lb <- model@lowbnd
  lb[ex.idx] <- 0
  idx <- match(medium, model@react_id)
  if(sum(is.na(idx))>0 & verbose) cat("Not found in model:", medium[is.na(idx)], "\n")
  lb[na.omit(idx)] <- uptake[!is.na(idx)]
  model.constraint <- cobrar::changeBounds(model, react=1:length(lb), ub=ub, lb=lb)
  sol.fba <- cobrar::fba(model.constraint)
  if( verbose ) message("state:",sol.fba@stat==solver_ok, "\tr=", round(sol.fba@obj,-log10(cutoff)))
  if(mtf & sol.fba@stat==solver_ok){ # minimize total flux
    sol.mtf <- cobrar::pfba(model.constraint)
    oldfluxsum <- sum(abs(sol.fba@fluxes))
    newfluxsum <- sum(abs(sol.mtf@fluxes))
    if( verbose ) message("mtf state:", "\told flux sum=", oldfluxsum, "\tnew flux sum=", newfluxsum)
  }
  if( fva ){
    sol.fva <- cobrar::fva(model.constraint, react=ex.idx)
    fva_max <- sol.fva$max.flux
    fva_min <- sol.fva$min.flux
  }
  if(normalize==T) norm <- sol.fba@obj else norm <- 1
  sol.out <- round(sol.fba@fluxes[ex.idx]/norm,2)
  out.dt <- data.table(rxn=model@react_id[ex.idx])
  if(useNames) out.dt[,name:=model@react_name[ex.idx]]
  out.dt[,fba:=sol.out]
  if(mtf) out.dt[,mtf:=round(sol.mtf@fluxes[ex.idx]/norm,2)]
  if(fva){
      out.dt[,fva_min:=round(fva_min/norm,2)]
      out.dt[,fva_max:=round(fva_max/norm,2)]
      out.dt[,fva_var:=ifelse(round(fva_max-fva_min,2)>0, "*", "")]
  }
  if(printCosts){ 
      out.dt[,red_costs:=round(sol.fba@redCosts[ex.idx],2)]
  }
  if(printExchanges){
    print(out.dt[abs(rowMeans(out.dt[, .SD, .SDcols=is.numeric])) > cutoff.print][order(fba)])
  }
  return(out.dt)
}

exp.mono.dat <- read_excel(validation.file, skip=1, sheet=1)
exp.mono.dt <- data.table(exp.mono.dat[-which(apply(exp.mono.dat,1,function(x)sum(is.na(x)))>2),]) # remove empty rows

exp.mono.dat2 <- read_excel(validation.file, skip=1, sheet=3)
exp.mono.dt2 <- data.table(exp.mono.dat2[-which(apply(exp.mono.dat2,1,function(x)sum(is.na(x)))>2),]) # remove empty rows



models1 <- cobrar::readSybilmod(model.file1)
xml.files2 <- list.files(model.file2, pattern="xml$", full.names=T)
xml.files2 <- xml.files2[grep("M10|AAM|mouse|batch|strict",xml.files2,invert=T)]
models2 <- lapply(xml.files2, function(x) cobrar::readSBMLmod(x))
names(models2) <- gsub(".xml","",basename(xml.files2))
models3 <- cobrar::readSybilmod(model.file3)
models3 <- models3[grep("mono", names(models3))]
names(models3) <- str_remove(names(models3), "_mono")
for(i in 1:length(models3)) models3[[i]]@mod_id <- names(models3)[i]

add_validation <- function(models, model.file, exp.mono.dt){
  # simulate models
  sim.dt <- data.table()
  for(mod in models){
    if(is.null(medium.file)){
      mod.ex.dt <- check_growth(mod, full=T, returnEx=T, printExchanges=F, verbose=F, normalize=T)
    }else{
      medium <- fread(medium.file)
      mod.ex.dt <- check_growth(mod, gs.medium=medium, full=T, returnEx=T, printExchanges=F, verbose=F, normalize=T)
    }
    sim.dt <- rbind(sim.dt, data.table(org=mod@mod_id,mod.ex.dt))
  }
  sim.dt[,sim:=mtf]
  sim.dt[,sim:=round(sim,1)]
  
  sim.dt[,ID:=str_extract(rxn, "(?<=EX_)cpd[0-9]+")]
  
  suppressWarnings(exp.mono.melt <- melt(exp.mono.dt[!is.na(ID),-"Comment"], id.vars=c("Metabolite","ID", "Category"), variable.name="org", value.name="exp"))
  
  exp.mono.melt <- exp.mono.melt[!(Category=="AA" & exp!=-1)] # exclude amino acid production/inactivity because of uncertainty of protein degradation
  
  # consider unclear metabolites in experimental data (e.g. L/D lactate) by avaraging simulated fluxes 
  unclear.mets <- exp.mono.melt[,length(unlist(str_split(ID,"/"))),by=ID][V1>1,ID]
  sim2.dt <- sim.dt[!ID %in% unlist(str_split(unclear.mets,"/")), .(org,sim,rxn,ID)]
  for(met in unclear.mets) sim2.dt <- rbind(sim2.dt, sim.dt[grepl(gsub("/","|",met),ID),.(sim=mean(sim), rxn=NA,ID=met),by=org])
  
  # merge experimental and simulated data
  sel.merge.dt <- merge(exp.mono.melt, sim2.dt, by=c("ID","org"), all.x=T)
  sel.merge.dt[is.na(sim), sim:=0]
  suppressWarnings(sel.merge.dt[,exp2:=factor(as.numeric(exp))])
  sel.merge.dt <- sel.merge.dt[!is.na(exp2)] # do not consider entries with unclear experimental state
  sel.merge.dt[,sim2:=factor(ifelse(sim>0,1,ifelse(sim<0,-1,0)))] # convert simulated fluxes to three-state variable (1 production, -1 consumption, 0 nothing)
  
  # compare predicted with experimental data
  val.dt <- data.table()
  err.dt <- data.table()
  for(cat in c("all", unique(sel.merge.dt$Category))){
    if(cat=="all"){
      sel.merge.dt2 <- sel.merge.dt
    }else{
      sel.merge.dt2 <- sel.merge.dt[Category==cat]
    }
    sel.cm <- caret::confusionMatrix(sel.merge.dt2$exp2, sel.merge.dt2$sim2) # three classes (-1, +1, 0)
    val.tmp <- data.table(model=basename(model.file),medium=ifelse(!is.null(medium.file),basename(medium.file),""),cases=sum(sel.cm$table), t(sel.cm$overall))
    err.tmp <- sel.merge.dt2[exp2!=sim2, .(Metabolite,ID,org,exp2,sim2)]
    err.tmp[,error:=ifelse(exp2!=0&sim2==0, "FN", ifelse((exp2==-1&sim2==1)|(exp2==1&sim2==-1),"FD",ifelse(exp2==0&sim2!=0,"FP",NA)))]
    cor.dt <- sel.merge.dt2[exp2==sim2, .(Metabolite,ID,org,exp2,sim2)]
    cor.dt[,correct:=ifelse(exp2==0&sim2==0, "TN", ifelse((exp2==-1&sim2==-1)|(exp2==1&sim2==1),"TP","?"))]
    if(cat=="all"){
      print(transpose(head(err.tmp[error=="FP",list(error=unique(error),.N),by=Metabolite][order(N,decreasing=T)]),make.names="Metabolite"))
      print(transpose(head(err.tmp[error=="FN",list(error=unique(error),.N),by=Metabolite][order(N,decreasing=T)]),make.names="Metabolite"))
    }
    val.tmp2 <- data.table(val.tmp[,list(model,medium,category=cat,cases, Accuracy, Kappa)], as.data.frame.matrix(t(table(err.tmp$error))), as.data.frame.matrix(t(table(cor.dt$correct))))
    val.dt <- rbind(val.dt, val.tmp2, fill=T)
    err.dt <- rbind(err.dt, err.tmp, fill=T)
  }
  return(list(val.dt, err.dt))
}

validation1 <- add_validation(models1, model.file1, exp.mono.dt)
validation1a <- add_validation(models1, model.file1, exp.mono.dt2)
validation2 <- add_validation(models2, model.file2, exp.mono.dt)
validation2a <- add_validation(models2, model.file2, exp.mono.dt2)
validation3 <- add_validation(models3, model.file3, exp.mono.dt)
validation3a <- add_validation(models3, model.file3, exp.mono.dt2)
val.dt <- rbind(validation1a[[1]], validation2a[[1]], validation3a[[1]], fill=T)

val.dt[category=="all"][order(Accuracy)]

val.sel <- val.dt[model=="OMM_munich-aam_20241205.RDS", .(model, category, Accuracy, Kappa)]

val.sel[category=="all", category:="All"]
val.sel[, category:=factor(category, levels=c("All", "AA", "CH", "FA", "OA", "Other"), ordered=F)]

ggplot(val.sel, aes(x=Accuracy, y=Kappa, fill=category)) + geom_point(size=5, shape=21) + theme_bw(base_size=14)
ggsave("img/validation_mono.pdf", width=4, height=3)
ggsave("img/validation_mono.png", width=4, height=3)

val.dt[, id:=ifelse(grepl("csm", model), "csm", ifelse(grepl("rfastcormics2", model),"csm2", "org")) ]
ggplot(val.dt, aes(x=Accuracy, y=Kappa, fill=category)) + geom_point(size=5, shape=21) + theme_bw(base_size=14) + ggrepel::geom_text_repel(aes(label=id), size=3, box.padding = 0.5)
ggsave("img/validation_mono_model-cmp.pdf", width=4.5, height=3)
ggsave("img/validation_mono.model-cmp.png", width=4.5, height=3)
