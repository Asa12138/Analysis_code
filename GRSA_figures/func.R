
geo_list2=tibble::tribble(
         ~GEO_id,                    ~Disease,    ~PathID, ~Normal, ~Case,   ~Pubmed,                              ~Tissue,
  "GSE14924_CD4",    "Acute Myeloid Leukemia", "hsa05221",     10L,   10L, 19710498L,                         "CD4 T Cell",
  "GSE14924_CD8",    "Acute Myeloid Leukemia", "hsa05221",     11L,   10L, 19710498L,                         "CD8 T Cell",
       "GSE1297",       "Alzheimer's disease", "hsa05010",      9L,    7L, 14769913L,                    "Hippocampal CA1",
     "GSE5281EC",       "Alzheimer's disease", "hsa05010",     13L,   10L, 17077275L,           "Brain, Entorhinal Cortex",
    "GSE5281HIP",       "Alzheimer's disease", "hsa05010",     13L,   10L, 17077275L,                 "Brain, hippocampus",
    "GSE5281VCX",       "Alzheimer's disease", "hsa05010",     12L,   19L, 17077275L,       "Brain, primary visual cortex",
       "GSE7305",        "Endometrial cancer", "hsa05213",     10L,   10L, 17640886L,         "Endometrium/Ovarian tissue",
      "GSE36389",        "Endometrial cancer", "hsa05213",      7L,   13L,        NA,                        "Endometrium",
       "GSE8762",      "Huntington's disease", "hsa05016",     10L,   12L, 17724341L,                         "Lymphocyte",
      "GSE73655",      "Huntington's disease", "hsa05016",      7L,   13L, 26756592L,               "Subcutaneous adipose",
      "GSE37517",      "Huntington's disease", "hsa05016",      5L,    8L, 22748968L,                   "Neural stem cell",
      "GSE15471",         "Pancreatic cancer", "hsa05212",     35L,   35L, 19260470L,                           "Pancreas",
      "GSE16515",         "Pancreatic cancer", "hsa05212",     15L,   15L, 19732725L,                           "Pancreas",
      "GSE28735",         "Pancreatic cancer", "hsa05212",     45L,   45L, 23918603L,                           "Pancreas",
      "GSE20153",       "Parkinson's disease", "hsa05012",      8L,    8L, 20926834L, "Blymphocytes from peripheral blood",
      "GSE19587",       "Parkinson's disease", "hsa05012",     10L,   12L, 20837543L,                              "Brain",
      "GSE55945",           "Prostate cancer", "hsa05215",      7L,   12L, 19737960L,                           "Prostate",
      "GSE26910",           "Prostate cancer", "hsa05215",      6L,    6L, 21611158L,                           "Prostate",
      "GSE14762",      "Renal cell carcinoma", "hsa05211",     12L,    9L, 19252501L,                             "Kidney",
       "GSE6344",      "Renal cell carcinoma", "hsa05211",     10L,   10L, 17699851L,                     "Clear cell RCC",
      "GSE65144",            "Thyroid cancer", "hsa05216",     13L,   12L, 25675381L,                            "Thyroid",
      "GSE58545",            "Thyroid cancer", "hsa05216",     18L,   27L, 26625260L,                            "Thyroid",
      "GSE26887", "Type II diabetes mellitus", "hsa04930",      5L,    7L, 22427379L,                     "Left ventricle",
      "GSE39825", "Type II diabetes mellitus", "hsa04930",      6L,    4L, 23919306L,         "Fibroblasts (cell culture)"
  )%>%as.data.frame()


geo_list3=tibble::tribble(
           ~GEO, ~KO_gene, ~Impacted_path_number, ~Normal, ~Case,   ~Pubmed,                                           ~Tissue,
     "GSE22873",  "Myd88",                24L,     11L,    8L, 22075646L,                                              "Liver",
    "GSE70302a",   "Il1a",                21L,      4L,    4L, 26224856L,                                        "Spinal cord",
    "GSE70302b",   "Il1b",                41L,      4L,    4L, 26224856L,                                        "Spinal cord",
     "GSE58120",    "Il2",                20L,      6L,    6L, 25652593L,                            "Myeloid dendritic cells",
     "GSE46211", "Tgfbr2",                23L,     12L,    6L, 24496627L,                            "Anterior palatal tissue",
    "GSE138957",   "Akt1",                97L,      3L,    3L, 32937140L,                                         "Cell lines",
     "GSE85754",  "Fgfr1",                16L,      4L,    4L, 28433771L,                                             "Breast",
     "GSE88799", "Crebbp",                28L,      4L,    5L, 28069569L,                                             "B cell",
      "GSE4451",    "Met",                21L,      6L,    6L, 16710476L,                                              "Liver"
    )%>%as.data.frame()

# "Bhlhe40", "Id3", "Dusp5", "Onecut1", "Neurod1"
pos_res=lapply(c("Myd88",  "Il1a", "Il1b", "Il2", "Tgfbr2","Akt1","Fgfr1","Crebbp","Met"),
       \(i)filter(mmu_kegg_pathway$all_org_gene,gene_symbol==i)%>%select(pathway_id)%>%mutate(mutate=i))%>%
    do.call(rbind,.)

get_sig_from_all=function(comparison_res_m){
    enrich_sig_m=list(
        #RS_m=filter(comparison_res_m$RS_m,(ReporterScore)>1.64,p.adjust<0.01)%>%pull(ID),
        GRSA=filter(comparison_res_m$RS_d,abs(ReporterScore)>1.64,p.value<0.05)%>%pull(ID),
        Fisher=filter(comparison_res_m$fisher_res,p.value<0.05)%>%pull(ID),
        CP=filter(comparison_res_m$enrich_res,p.value<0.05)%>%pull(ID),
        GSEA=filter(comparison_res_m$gsea_res@result,pvalue<0.05)%>%pull(ID)
        #GSA=filter(comparison_res_m$gsa_res,p.value<0.05)%>%pull(ID)
    )
    enrich_sig_m
}

get_sig_from_all2=function(comparison_res_m){
    enrich_sig_m=list(
        #RS_m=filter(comparison_res_m$RS_m,(ReporterScore)>1.64,p.adjust<0.01)%>%pull(ID),
        GRSA=filter(comparison_res_m$RS_d,abs(ReporterScore)>1.64,p.adjust<0.05)%>%pull(ID),
        Fisher=filter(comparison_res_m$fisher_res,p.adjust<0.05)%>%pull(ID),
        CP=filter(comparison_res_m$enrich_res,p.adjust<0.05)%>%pull(ID),
        GSEA=filter(comparison_res_m$gsea_res@result,p.adjust<0.05)%>%pull(ID)
        #GSA=filter(comparison_res_m$gsa_res,p.adjust<0.05)%>%pull(ID)
    )
    enrich_sig_m
}

#10.29更新，比较富集方法，选取p-value最小，logFC最大的前10%的gene做富集
hsa_gene=unique(hsa_kegg_pathway$all_org_gene$gene_symbol)
mmu_gene=unique(mmu_kegg_pathway$all_org_gene$gene_symbol)

test_tests4=function(KO_abundance,group,metadata,name,method="wilcox.test",p.adjust=T){
    if(file.exists(paste0(name,"_methods_res_d.RDS"))) {
        methods_res_d=readRDS(paste0(name,"_methods_res_d.RDS"))
    }
    else {
        if(p.adjust){
            methods_res_d=reporter_score(KO_abundance,group = group,
                                         metadata = metadata,mode="directed",method = method,feature = "gene",type = "mmu",threads = 1)
        }
        else {
            methods_res_d=reporter_score(KO_abundance,group = group,p.adjust.method1 = "none",
                                         metadata = metadata,mode="directed",method = method,feature = "gene",type = "mmu",threads = 1)
        }
        saveRDS(methods_res_d,file = paste0(name,"_methods_res_d.RDS"))
    }

    #reporter_score
    #fisher
    res.dt=methods_res_d$ko_stat
    # if(attributes(methods_res_d$reporter_s)$type=="hsa")res.dt=filter(res.dt,KO_id%in%hsa_gene)
    # if(attributes(methods_res_d$reporter_s)$type=="mmu")res.dt=filter(res.dt,KO_id%in%mmu_gene)

    add_mini=NULL
    if(!"logFC"%in%colnames(res.dt)){
        message("No logFC in the data.frame, calculate.")
        vs_group=grep("average",colnames(res.dt),value = T)
        if(length(vs_group)!=2)stop("logFC only available for two groups")
        tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])

        if (is.null(add_mini))
            add_mini = min(tmp[tmp > 0]) * 0.05
        res.dt$logFC=log2((res.dt[,vs_group[2]]+add_mini)/(res.dt[,vs_group[1]]+add_mini))
    }

    sig_gene=filter(res.dt,p.value<0.05)%>%top_n(400,abs(logFC))%>%pull(KO_id)
    res.dt=mutate(res.dt,origin_p.adjust=ifelse(KO_id%in%sig_gene,0.04,0.1))
    fisher_res=KO_fisher(res.dt,padj_threshold = 0.05,modulelist = methods_res_d$modulelist)
    #enricher
    enrich_res=KO_enrich(res.dt,padj_threshold = 0.05,modulelist = methods_res_d$modulelist)

    #GESA
    set.seed(1234)
    gsea_res=KO_gsea(methods_res_d,padj_threshold = 1.1,weight = "Z_score")
    #GSA
    gsa_res=KO_gsa(methods_res_d)

    comparison_res_m=c(list(RS_d=methods_res_d$reporter_s),list(fisher_res=fisher_res),
                       list(enrich_res=enrich_res),gsea_res=gsea_res,list(gsa_res=gsa_res))

    saveRDS(comparison_res_m,file = paste0(name,"_comparison5.RDS"))
}

similar=\(a,b){tmp=two_set(a,b);tmp[2]/sum(tmp)}
similar_each=function(ls){
    n=length(ls)
    sim_mat=matrix(0,nrow = n,ncol = n)
    rownames(sim_mat)=colnames(sim_mat)=names(ls)
    for (i in 1:n) {
        for (j in i:n) {
            sim_mat[i,j]= sim_mat[j,i]=similar(ls[[i]],ls[[j]])
        }
    }
    as.data.frame(sim_mat)
}

#2024.1.15,添加更多的EA方法比较：
compare_EA=function(KO_abundance,group,metadata,name,method="t.test",type = "mmu",modulelist=NULL){
  if(file.exists(paste0(name,"_methods_res_d.RDS"))) {
    GRSA_res=readRDS(paste0(name,"_GRSA_res.RDS"))
  }
  else {
    if(type!="id")GRSA_res=GRSA(KO_abundance,group = group,metadata = metadata,
                                mode="directed",method = method,feature = "gene",type = type,threads = 1,perm = 999)
    else {
      GRSA_res=GRSA(KO_abundance,group = group,metadata = metadata,
                    mode="directed",method = method,modulelist = modulelist,threads = 2,perm = 999)
    }
    saveRDS(GRSA_res,file = paste0(name,"_GRSA_res.RDS"))
  }

  #reporter_score
  #fisher
  res.dt=GRSA_res$ko_stat
  # if(attributes(GRSA_res$reporter_s)$type=="hsa")res.dt=filter(res.dt,KO_id%in%hsa_gene)
  # if(attributes(GRSA_res$reporter_s)$type=="mmu")res.dt=filter(res.dt,KO_id%in%mmu_gene)

  add_mini=NULL
  if(!"logFC"%in%colnames(res.dt)){
    message("No logFC in the data.frame, calculate.")
    vs_group=grep("average",colnames(res.dt),value = T)
    if(length(vs_group)!=2)stop("logFC only available for two groups")
    tmp=c(res.dt[,vs_group[1]],res.dt[,vs_group[2]])

    if (is.null(add_mini))
      add_mini = min(tmp[tmp > 0]) * 0.05
    res.dt$logFC=log2((res.dt[,vs_group[2]]+add_mini)/(res.dt[,vs_group[1]]+add_mini))
  }
  sig_gene=filter(res.dt,p.value<0.05)%>%top_n(400,abs(logFC))%>%pull(KO_id)
  res.dt=mutate(res.dt,origin_p.adjust=ifelse(KO_id%in%sig_gene,0.04,0.1))

  #fisher_res=KO_fisher(res.dt,padj_threshold = 0.05,modulelist = GRSA_res$modulelist)

  #enricher
  enrich_res=KO_enrich(res.dt,padj_threshold = 0.05,modulelist = GRSA_res$modulelist)
  #GESA
  set.seed(1234)
  gsea_res=KO_gsea(GRSA_res,weight = "Z_score")
  #GSA
  gsa_res=KO_gsa(GRSA_res)
  #SEA
  sea_res=KO_sea(GRSA_res)
  #SAFE
  safe_res=KO_safe(GRSA_res)
  #PADOG
  padog_res=KO_padog(GRSA_res)
  #GSVA
  gsva_res=KO_gsva(GRSA_res)

  comparison_res_m=c(list(RS_d=GRSA_res$reporter_s),
                     list(enrich_res=enrich_res),
                     gsea_res=gsea_res,
                     list(gsa_res=gsa_res),
                     list(sea_res=sea_res),
                     list(safe_res=safe_res),
                     list(padog_res=padog_res),
                     list(gsva_res=gsva_res))
  comparison_res_m$gsea_res=data.frame(comparison_res_m$gsea_res)
  saveRDS(comparison_res_m,file = paste0(name,"_compare_EA.RDS"))
}
