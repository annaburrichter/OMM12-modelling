library(cobrar)
library(readxl)
library(Biostrings)

files <- list.files("../dat/", pattern="*.ffn.gz", full.names=T)
proteome1 <- data.table(readxl::read_excel("../dat/reanalysed_single_vs_community.xlsx", sheet=1, na=c("","NA"))) # comparison of batch community metaproteome only because batch culture was also used for proteome2 data set
proteome1.I48 <- data.table(readxl::read_excel("../dat/reanalysed_single_vs_community.xlsx", sheet=2, na=c("","NA")))
proteome2 <- data.table(readxl::read_excel("../dat/reanalysed_mouse_vs_invitro.xlsx", sheet=1, na=c("","NA")))
setnames(proteome1, old=c("community(batch)_1", "community(batch)_2", "community(batch)_3"),new=c("batch_1","batch_2","batch_3"))
setnames(proteome1.I48, old=c("community(batch)_1", "community(batch)_2", "community(batch)_3"),new=c("batch_1","batch_2","batch_3"))
proteome1[is.na(proteome1)] <- 0; proteome1.I48[is.na(proteome1.I48)] <- 0; proteome2[is.na(proteome2)] <- 0
proteome1[,strain:=str_extract(species, "[A-Z]+[0-9]+$")]
proteome1.I48[,strain:=str_extract(species, "[A-Z]+[0-9]+$")]
proteome2[,strain:=str_extract(species, "[A-Z]+[0-9]+$")]

models <- readRDS("../dat/OMM_munich-aam_20241205.RDS")

stat.dt <- data.table()
rxn.dt <- data.table()
gene_prot_rxn.dt <- data.table()
for(f in files){
  	id <- str_remove(basename(f), "\\.ffn.gz")
  	seq <- readDNAStringSet(f)
	match.dt <- data.table(protein_id=str_extract(names(seq),"(?<=_cds_).*?(?=_[0-9]+)"), locus_tag=str_extract(names(seq),"(?<=locus_tag\\=).*?(?=\\])"))
	
	mod <- models[[id]]
	mod.genes <- unique(unlist(mod@genes))
	anno.shared <- unique(match.dt[protein_id %in% mod.genes, locus_tag]) # proteins shared between model and protein fasta
	mod.gene_rxn <- sapply(mod.genes, function(g) paste(mod@react_id[grep(g, mod@genes)], collapse=","))
	gene_prot_rxn.dt <- rbind(gene_prot_rxn.dt, data.table(id, proteom=match.dt$locus_tag[match(mod.genes, match.dt$protein_id)], model=mod.genes, reaction=mod.gene_rxn))
	
	# in vitro
	for(col in c(id, "batch")){ # only batch comparison 
		# median of all replicates
		if(id=="I48"){
			proteins1 <- unique(na.omit(proteome1.I48[strain==id,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome1.I48)), by=`locus tag_old`][V1>0,`locus tag_old`])) # I48 specific proteome
		}else{
		  proteins1 <- unique(na.omit(proteome1[strain==id,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome1)), by=`locus tag_old`][V1>0,`locus tag_old`])) # proteins in proteome	
		}
		shared1 <- intersect(anno.shared, proteins1)
		mod.cov1 <- length(shared1)/length(mod.genes) # % of model genes with proteins in proteome
		prot.cov1 <- length(shared1)/length(proteins1) # % of proteome shared with model
		stat.dt <- rbind(stat.dt, data.table(id, col, mod.cov=mod.cov1, prot.cov=prot.cov1))
		if(id=="I48"){
			prot.active <- proteome1.I48[strain==id & `locus tag_old` %in% shared1,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome1.I48)), by=`locus tag_old`][V1>0,`locus tag_old`]			
		}else{
			prot.active <- proteome1[strain==id & `locus tag_old` %in% shared1,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome1)), by=`locus tag_old`][V1>0,`locus tag_old`]			
		}
		genes.inactive <- match.dt[!locus_tag %in% prot.active & protein_id %in% unlist(mod@genes), protein_id]
		idx.rxn.inactive <- match(geneDel(mod, genes.inactive), mod@react_id)
		rxn.active <- mod@react_id[-idx.rxn.inactive]
		rxn.has.gene <- ifelse(mod@gprRules[-idx.rxn.inactive]!="",T,ifelse(mod@react_id[-idx.rxn.inactive]=="bio1",T,F))
		fwrite(data.table(rxn.active=rxn.active[rxn.has.gene]), paste0("../dat/csm/",id,"-",col,"_strict.csv"))
		fwrite(data.table(rxn.active=rxn.active), paste0("../dat/csm/",id,"-",col,"_all.csv"))
		rxn.dt <- rbind(rxn.dt, data.table(id, col, rxn.active[rxn.has.gene]))
	}
	# in vivo
	for(col in c("AAM","M10", "mouse")){
		# median of all replicates
		proteins2 <- unique(na.omit(proteome2[strain==id,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome2)), by=`locus tag_old`][V1>0,`locus tag_old`])) # proteins in proteome
		shared2 <- intersect(anno.shared, proteins2)
		mod.cov2 <- length(shared2)/length(mod.genes) # % of model genes with proteins in proteome
		prot.cov2 <- length(shared2)/length(proteins2) # % of proteome shared with model
		stat.dt <- rbind(stat.dt, data.table(id, col, mod.cov=mod.cov2, prot.cov=prot.cov2))
		prot.active <- proteome2[strain==id & `locus tag_old` %in% shared2,apply(.SD,1,median),.SDcols=grep(col,colnames(proteome2)), by=`locus tag_old`][V1>0,`locus tag_old`]
		genes.inactive <- match.dt[!locus_tag %in% prot.active & protein_id %in% unlist(mod@genes), protein_id]
		idx.rxn.inactive <- match(geneDel(mod, genes.inactive), mod@react_id)
		rxn.active <- mod@react_id[-idx.rxn.inactive]
		rxn.has.gene <- ifelse(mod@gprRules[-idx.rxn.inactive]!="",T,ifelse(mod@react_id[-idx.rxn.inactive]=="bio1",T,F))
		fwrite(data.table(rxn.active=rxn.active[rxn.has.gene]), paste0("../dat/csm/",id,"-",col,"_strict.csv"))
		fwrite(data.table(rxn.active=rxn.active), paste0("../dat/csm/",id,"-",col,"_all.csv"))		
		rxn.dt <- rbind(rxn.dt, data.table(id, col, rxn.active[rxn.has.gene]))
	}
}
saveRDS(list(stat.dt, rxn.dt), "../dat/core-rxn.RDS")
fwrite(gene_prot_rxn.dt, "../dat/gene-prot-rxn.csv")
