#################################
#CANCER VARIANT EXPLORER V1.1
#WGCNAmelanoma extension
#Server script
#Andreas Mock
#University of Cambridge
#################################

#######################
#PART I: PREFIX########
#######################

#protein-changing categories of SNVs according to a functional effect score below 5,
#suggested by http://gatkforums.broadinstitute.org/gatk/discussion/4220/what-is-the-difference-between-tx-mode-best-effect-vs-canonical
pc_cat <- c("De_novo_Start_OutOfFrame", "Nonsense_Mutation",
            "Nonstop_Mutation", "Missense_Mutation",
            "De_novo_Start_InFrame","In_Frame_Del","In_Frame_Ins",
            "Frame_Shift_Del","Frame_Shift_Ins","Frame_Shift_Sub",
            "Start_Codon_SNP","Start_Codon_Del","Start_Codon_Ins",
            "Start_Codon_DNP","Start_Codon_TNP","Start_Codon_ONP",
            "Stop_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins",
            "Stop_Codon_DNP","Stop_Codon_TNP","Stop_Codon_ONP",
            "Splice_Site","Splice_Site_SNP","Splice_Site_Del",
            "Splice_Site_Ins","Splice_Site_DNP","Splice_Site_TNP",
            "Splice_Site_ONP","Splice_Site","miRNA")


#create list only containing protein changing variants (pcv)
pcv <- vector("list",length(v))
for (i in 1:length(v)){
  pcv[[i]] <- v[[i]][v[[i]]$variant_classification %in% pc_cat & !v[[i]]$protein_change=="",]
}
names(pcv) <- names(v)
if(0 %in% sapply(pcv,nrow)){
  stop("At least one oncotator input file has no protein-changing variants.")
}

#retreive rankscores (rs) from variant information
#rs_names <- colnames(pcv[[1]])[grep("rankscore",colnames(pcv[[1]]))]
#rs_names <- gsub("dbNSFP_",replacement = "",rs_names)
#rs_names = gsub("_rankscore",replacement = "",rs_names)
rs_names = c("CADD_raw","FATHMM","GERP.._RS","LRT_converted" ,"LR",
             "MutationAssessor","MutationTaster_converted","Polyphen2_HDIV",
             "Polyphen2_HVAR","RadialSVM","SIFT_converted","SiPhy_29way_logOdds" ,
             "phastCons100way_vertebrate","phastCons46way_placental", "phastCons46way_primate",
             "phyloP100way_vertebrate","phyloP46way_placental","phyloP46way_primate" )

suppressWarnings(
  for(n in 1:length(pcv)){
  rs_matrix = matrix(NA,nrow(pcv[[n]]),length(rs_names))
  for (i in 1:length(rs_names)){
    rs = strsplit(as.vector(pcv[[n]][,grep("rankscore",colnames(v[[1]]))[i]]), "[|]")
    rs = as.numeric(sapply(rs,"[", i=1))
    rs_matrix[,i] = rs
  }
  pcv[[n]]$rs_matrix = rs_matrix
  pcv[[n]]$sum_is_na = apply(apply(pcv[[n]]$rs_matrix,1,is.na),2,sum)
  colnames(pcv[[n]]$rs_matrix) = rs_names
  }
)

#define colours for plotting
hmcol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))
palette_breaks <- seq(0,1,length.out = 257)

#####################################################
#PART II: OUTPUTS####################################
#####################################################

shinyServer(function(input, output, session) {
  ####ANNOTATION####
  output$nVariants <- renderText({
    paste("Total number of variants:", nrow(v[[input$sample]]))
  })

  output$algorithm_choice <- renderText({
    paste(input$algorithm, "cutoff")
  })

  output$varianttype <- renderPlot({
    ggplot(v[[input$sample]], aes(x=variant_type)) +
      geom_bar(stat="count") + xlab("") +
      theme_bw() + ylab("count")
  })

  output$variantclass <- renderPlot({
    var_col = as.vector(v[[input$sample]]$variant_classification)
    var_col[var_col %in% pc_cat] = "protein-changing"
    var_col[!var_col =="protein-changing"] = "other"
    df = as.data.frame(cbind(v[[input$sample]], var_col))
    ggplot(df, aes(x=variant_classification, fill=var_col)) +
      geom_bar(stat="count") + coord_flip() +
      xlab("") + theme_bw() +
      guides(fill=guide_legend(reverse=TRUE)) +
      scale_fill_manual(values = c("gold1","blue")) +
      labs(fill="Variant classification")
  })

  output$nCoding_SNVs <- renderText({
    paste("Protein changing single-nucleotide variants:",
          nrow(pcv[[input$sample]]))
  })
  #consensus clustering
  ConsClustInput <- reactive({
    if((length(which(pcv[[input$sample]]$sum_is_na==18))==length(pcv[[input$sample]]$sum_is_na))==FALSE){
    ConsClust = vector("list",length(v))
    for (n in 1:length(v)){
      M_temp = pcv[[n]]$rs_matrix[!sapply(1:nrow(pcv[[n]]$rs_matrix),
                                          function(x) TRUE %in% is.na(pcv[[n]]$rs_matrix)[x,]),]
      ConsClust[[n]] = ConsensusClusterPlus(
        sweep(M_temp,1,apply(M_temp,1,median,na.rm=T)),maxK=6,reps=input$nperm,
        pItem=0.8,pFeature=0.5, clusterAlg="hc",distance="pearson",
        seed=1262118388.71279,plot=NULL, verbose=FALSE)
    }
    names(ConsClust) = names(v)
    for(n in 1:length(v)){
      for(k in 2:6){
        rownames(ConsClust[[n]][[k]]$ml) = rs_names
        colnames(ConsClust[[n]][[k]]$ml) = rs_names
      }
    }
    return(ConsClust)
    }
  })

  #combination score of rankscores
  CombScoreInput <- reactive({
    if((length(which(pcv[[input$sample]]$sum_is_na==18))==length(pcv[[input$sample]]$sum_is_na))==TRUE){
    rep(0,nrow(pcv[[input$sample]]))
    } else {
    ConsClust = ConsClustInput()
    df = as.data.frame(t(pcv[[input$sample]]$rs_matrix),stringsAsFactors = FALSE)
    centroid = t(ddply(df, .(ConsClust[[input$sample]][[as.numeric(
      input$pred_modules)]]$consensusClass), colwise(mean, na.rm = TRUE)))
    colnames(centroid) = paste0("module",1:ncol(centroid))
    centroid = centroid[2:nrow(centroid),]
    centroid = apply(centroid,2,as.numeric)
    cutoff = 0.75
    comb_score = as.vector(apply(centroid*(centroid>cutoff*1),1,
                                 function(x) sum(x,na.rm=TRUE)))
    comb_score
    }
  })

  #update slider according to algorithm
  observe({
    if(!input$algorithm=="dbNSFP combination score"){
      updateSliderInput(session, inputId = "comb_score",
                        value = 0, min = 0, max=1,step = 0.05)
    }else{
      x = as.numeric(input$pred_modules)
      updateSliderInput(session, inputId = "comb_score",
                        value = 0, min = 0, max =x,step = 0.1)
    }
  })

  output$ConsHM <- renderPlot({
    if((length(which(pcv[[input$sample]]$sum_is_na==18))==length(pcv[[input$sample]]$sum_is_na))==TRUE){
      plot.new()
      text(0.5,1,"Heatmap cannot be displayed as no variant has dbNSFP annotation",cex=1.3)
    } else {
    ConsClust = ConsClustInput()
    score_col = c("orange","blueviolet","yellow","blueviolet","orange",
                  rep("blueviolet",4),"orange","blueviolet",rep("yellow",7))
    plot.new()
    heatmap.2(ConsClust[[input$sample]][[as.numeric(
      input$pred_modules)]]$ml,margins=c(17,17),col=hmcol,
      ColSideColors = ConsClust[[input$sample]][[as.numeric(
        input$pred_modules)]]$clrs[[1]], RowSideColors = score_col,trace="n",
      scale = "n", breaks=palette_breaks,key.xlab="consensus index",
      key.title="")
    legend("bottomleft",title="category",legend = c(
      "conservation","functional prediction","ensemble"),
      fill=c("yellow","orange", "blueviolet"),
      border =c("yellow", "orange", "blueviolet"),xpd=TRUE, cex=0.9)
    legend("topright",legend = paste0("cluster ",1:as.numeric(input$pred_modules)),
           fill=ConsClust[[input$sample]][[as.numeric(input$pred_modules)]]$clrs[[3]],
           border=ConsClust[[input$sample]][[as.numeric(
             input$pred_modules)]]$clrs[[3]], xpd=TRUE, cex=0.9)
    }
  })


  ####PRIORITISATION####
  output$Centroid <- renderPlot({
    if(input$algorithm=="dbNSFP combination score"){
      comb_score = CombScoreInput()
    } else {
      comb_score = rs_matrix[,which(rs_names==input$algorithm)]
      comb_score[is.na(comb_score)] = 0
    }
    cutoff_boolean = comb_score[order(comb_score,decreasing = TRUE)]>(
      as.numeric(input$comb_score)-0.000001)
    colors = rep("cornflowerblue",length(comb_score))
    colors[cutoff_boolean==FALSE] = "brown"
    plot(comb_score[order(comb_score,decreasing = TRUE)],
         xlab=c("variant (decreasing order)"), ylab="score",
         col=colors, main=paste(input$algorithm), pch=20)
    abline(h = as.numeric(input$comb_score),lty=2, lwd=2, col="brown")
    legend("topright",legend=c("included","excluded"),col=c("cornflowerblue","brown"),pch=20)
  })

  output$ThousandG <- renderPlot({
    ThousandG = as.vector(droplevels(as.factor(pcv[[input$sample]]$X1000Genome_AF)))
    ThousandG[ThousandG==""] = NA
    if(input$algorithm=="dbNSFP combination score"){
      comb_score = CombScoreInput()
    } else {
      comb_score = rs_matrix[,which(rs_names==input$algorithm)]
      comb_score[is.na(comb_score)] = 0
    }
    cutoff_boolean = comb_score<as.numeric(input$comb_score)
    col1 = rep("cornflowerblue", sum(!is.na(ThousandG)))
    col1[cutoff_boolean[!is.na(ThousandG)]] = "brown"
    if(1 %in% input$db){col1 = rep("brown", sum(!is.na(ThousandG)))}
    if(sum(is.na(ThousandG))==length(ThousandG)){
      plot.new()
      text(0.5,1,"No germline variant",cex=1.3)
    } else {
      plot(ThousandG[!is.na(ThousandG)][
      order(ThousandG[!is.na(ThousandG)],decreasing=TRUE)],xlab=c("variant"),
      main=paste0(sum(!is.na(ThousandG))," variants in 1000 Genomes Project"),
      col=col1, ylab="AF",pch=19)
      legend("topright",legend=c("included","excluded"),
             col=c("cornflowerblue","brown"),pch=20)}
  })

  VarFilterInput <- reactive({
    if(input$algorithm=="dbNSFP combination score"){
      comb_score = CombScoreInput()
    } else {
      comb_score = rs_matrix[,which(rs_names==input$algorithm)]
      comb_score[is.na(comb_score)] = 0
    }
    f = (comb_score>=input$comb_score)*1 #keep if higher than comb.score cutoff
    f[pcv[[input$sample]][[3]]==18] = 0 #exclude variants with missing rankscores
    if (1 %in% input$db){ #if selection 1 - exclude all 1000 Genomes variants
      ThousandG = as.vector(droplevels(as.factor(pcv[[input$sample]]$X1000Genome_AF)))
      ThousandG[ThousandG==""] = NA
      f[ThousandG>0] = 0
    }
    if (2 %in% input$db){ #keep all overlapping mutations in COSMIC if selected
      f[!as.vector((droplevels(as.factor(pcv[[input$sample]]$COSMIC_overlapping_mutation_AAs))))==""] =1
    }
    if (3 %in% input$db){ #keep all non-SNVs
      nonSNVs = pcv[[input$sample]]$variant_type
      f[!nonSNVs=="SNP"] =1
    }
    if (4 %in% input$db){ #keep all DNA repair gene variants
      repair = pcv[[input$sample]]$HumanDNARepairGenes_Role
      repair[repair==""] = NA
      is.repair = !is.na(repair)
      f[is.repair] =1
    }
    sel_genes = unique(pcv[[input$sample]]$gene)
    gene_rescue = unlist(strsplit(input$rescue,", "))
    gene_rescue_intersect = sel_genes[sel_genes %in% gene_rescue]
    f[pcv[[input$sample]]$gene %in% gene_rescue_intersect] =1
    f = f==1
    f
  })

  output$COSMIC <- renderPrint({
    f = VarFilterInput()
    sum(!as.vector((droplevels(as.factor(pcv[[input$sample]]$COSMIC_overlapping_mutation_AAs)[f])))=="")
  })

  output$DNArepair <- renderTable({
    f = VarFilterInput()
    repair = pcv[[input$sample]]$HumanDNARepairGenes_Role
    repair[repair==""] = NA
    is.repair = (((!is.na(repair))*1) *(f*1)) == 1
    if(sum(is.repair)>0){
      data.frame(gene=pcv[[input$sample]]$gene[is.repair],role=repair[is.repair])
    }
  })

  ####TOP TABLE####
  output$top = renderText({
    f = VarFilterInput()
    paste0(sum(f)," variants prioritised for further exploration")
  })

  output$toptable = renderDataTable({
    f = VarFilterInput()
    if(input$algorithm=="dbNSFP combination score"){
      comb_score = CombScoreInput()
    } else {
      comb_score = rs_matrix[,which(rs_names==input$algorithm)]
      comb_score[is.na(comb_score)] = 0
    }
    genes = pcv[[input$sample]]$gene[f]
    data.frame(gene=genes,
               protein_change=pcv[[input$sample]]$protein_change[f],
               type=pcv[[input$sample]]$variant_type[f],
               classification=pcv[[input$sample]]$variant_classification[f],
               score=comb_score[f],
             COSMIC_entity = pcv[[input$sample]]$COSMIC_overlapping_primary_sites[f]

    )
  })

  output$downloadTop <- downloadHandler(
    filename = function() {
      paste('CVE_top_table_', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      f = VarFilterInput()
      if(input$algorithm=="dbNSFP combination score"){
        comb_score = CombScoreInput()
      } else {
        comb_score = rs_matrix[,which(rs_names==input$algorithm)]
        comb_score[is.na(comb_score)] = 0
      }
      genes = pcv[[input$sample]]$gene[f]
      toptable=data.frame(gene=genes,
                 protein_change=pcv[[input$sample]]$protein_change[f],
                 type=pcv[[input$sample]]$variant_type[f],
                 classification=pcv[[input$sample]]$variant_classification[f],
                 score=comb_score[f],
                 COSMIC_entity = pcv[[input$sample]]$COSMIC_overlapping_primary_sites[f])
      writeLines(paste("#Cancer Variant Explorer top table output |", Sys.Date(),"| sample file:",
                 names(v),"| number of algorithm clusters:",
                 input$pred_modules,"| algorithm chosen for prioritisation:",
                 input$algorithm,
                 "| algorithm cutoff:",
                 input$comb_score,"| filters: exclude 1000 Genomes:",
                 1 %in% input$db, ", include all COSMIC:", 2 %in% input$db,
                 ", include non-SNVs:", 3 %in% input$db, ", include DNA repair genes:",
                 4 %in% input$db, " | rescued genes:", input$rescue),
                 file)
      write.table(toptable, file,append=TRUE,sep=",")
    }
  )

  output$rescue <- renderPrint({ input$rescue })

  output$DGIdb_input = renderPrint({
    f = VarFilterInput()
    f_which = which(f==1)
    paste(as.vector(pcv[[input$sample]]$gene[f_which]),collapse=", ")
  })

  JSON = reactive({
    if(input$DGIdb_anno){
      f = VarFilterInput()
      f_which = which(f==1)
      genes = pcv[[input$sample]]$gene[f_which]
      url = paste0("http://dgidb.genome.wustl.edu/api/v1/interactions.json?genes=",
                   paste(genes, collapse=","),
                   "&drug_types=antineoplastic&interaction_sources=TEND,MyCancerGenome")
      JSON = fromJSON(url)
    }
  })

  output$DGIdb_status = renderPrint({
    if(input$DGIdb_anno){
      JSON = JSON()
      if(class(JSON)=="list"){
        cat("annotation retrieved")
      }
    } else {
    cat("click on filter to retreive annotation")
    }
  })

  output$DGIdb_output = renderPrint({
    if(input$DGIdb_anno){
      JSON = JSON()
      if(length(JSON$matchedTerms)>0){
        JSON$matchedTerms[[1]]
      } else {
     "no interactions found"
      }
    }
  })

  output$DGIdb_table = renderDataTable({
    if(input$DGIdb_anno){
    JSON = JSON()
      ngenes_JSON = nrow(JSON$matchedTerms)
      if (ngenes_JSON>0){
        df_gene = c()
        for (i in 1:ngenes_JSON){
          df_gene = c(df_gene, rep(JSON$matchedTerms[[1]][[i]],
                                   length(JSON$matchedTerms[[5]][[i]]$drugName)))
        }
        df_druginfo = matrix(NA,0,0)
        for (i in 1:ngenes_JSON){
          df_druginfo = rbind(df_druginfo, JSON$matchedTerms[[5]][[i]])
        }
        data.frame(geneName=as.vector(df_gene), drugName=df_druginfo$drugName,
                   interactionType=df_druginfo$interactionType,
                   source = df_druginfo$source)
      }
    }
  })

  output$downloadDrugs <- downloadHandler(
    filename = function() {
      paste('CVE_drug_interactions_', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      JSON = JSON()
      ngenes_JSON = nrow(JSON$matchedTerms)
      if (ngenes_JSON>0){
        df_gene = c()
        for (i in 1:ngenes_JSON){
          df_gene = c(df_gene, rep(JSON$matchedTerms[[1]][[i]],
                                   length(JSON$matchedTerms[[5]][[i]]$drugName)))
        }
        df_druginfo = matrix(NA,0,0)
        for (i in 1:ngenes_JSON){
          df_druginfo = rbind(df_druginfo, JSON$matchedTerms[[5]][[i]])
        }
        drugtable = data.frame(geneName=as.vector(df_gene), drugName=df_druginfo$drugName,
                   interactionType=df_druginfo$interactionType,
                   source = df_druginfo$source)
      }
      writeLines(paste("#Cancer Variant Explorer top table output |", Sys.Date(),"| sample file:",
                       names(v),"| number of algorithm clusters:",
                       input$pred_modules,"| algorithm chosen for prioritisation:",
                       input$algorithm,
                       "| algorithm cutoff:",
                       input$comb_score,"| filters: exclude 1000 Genomes:",
                       1 %in% input$db, ", include all COSMIC:", 2 %in% input$db,
                       ", include non-SNVs:", 3 %in% input$db, ", include DNA repair genes:",
                       4 %in% input$db, " | rescued genes:", input$rescue),
                 file)
      write.table(drugtable, file,append=TRUE,sep=",")
    }
  )

  ####MELANOMA CO-EXPRESSION NETWORKS####
  output$topgenes_in_modules <- renderText({
    f = VarFilterInput()
    if(input$algorithm=="dbNSFP combination score"){
      comb_score = CombScoreInput()
    } else {
      comb_score = rs_matrix[,which(rs_names==input$algorithm)]
      comb_score[is.na(comb_score)] = 0
    }
    topgenes=pcv[[input$sample]]$gene[f]
    paste("Number of prioritized variants in top 5000 mad genes:",
          sum(topgenes %in% genes_WGCNA), "out of ",length(topgenes %in% genes_WGCNA))
  })


  output$MS <- renderPlot({
    f = VarFilterInput()
    topgenes=pcv[[input$sample]]$gene[f]
    selection = genes_WGCNA %in% topgenes
    topgenes_per_module = as.data.frame(table(modules,selection))$Freq[43:84]
    METree_GO$labels[24] = "regulation of transmembrane receptor serine"
    METree_GO$labels[27] = "homophilic cell adhesion"
    METree_GO$labels[30] = "positive reg. of cytosolic calcium concentration"

    par(mfrow=c(1,2))
    par(mar=c(8,0,0,0))
    a = plot.phylo(as.phylo(METree_GO), font=1, cex=0.90)
    plot(0.5,1, xlim=c(0,10),ylim=c(1,42),axes=FALSE, xlab="")
    text(x=0, y=c(1:42), labels  = c(0:41)[label_order][METree_GO$order], xpd=TRUE)
    mtext("module number", las=2, side=1, at=0)
    points(rep(0.5,42),c(1:42),xpd=TRUE,pch=15, cex=1.8,
           col=c(labels2colors(c(0:41)))[label_order][METree_GO$order])
    mtext("module colour", las=2, side=1, at=0.5)
    text(x=1, y=c(1:42), labels  = topgenes_per_module[label_order][METree_GO$order], xpd=TRUE)
    mtext("# variants", las=2, side=1, at=1)
    lscore_rectangles = cbind(as.numeric(c(MS_lscore_bar[label_order][METree_GO$order]/5)),rep(0.7,42))
    symbols(1.5+as.numeric(c(MS_lscore_bar[label_order][METree_GO$order]/5))/2,c(1:42),
            rectangles = lscore_rectangles,add=TRUE, xpd=TRUE,inches = FALSE, adj=1, fg="green3", bg="green3",lwd=2)
    mtext("lymphocyte score",las=2,side = 1, at = 1.5)
    pmet_rectangles = cbind(as.numeric(c(MS_pmet_bar[label_order][METree_GO$order]/5)),rep(0.7,42))
    symbols(3+as.numeric(c(MS_pmet_bar[label_order][METree_GO$order]/5))/2,c(1:42),
            rectangles = pmet_rectangles,add=TRUE, xpd=TRUE,inches = FALSE, adj=1, fg="red3", bg="red3",lwd=2)
    mtext("primary vs met",las=2,side = 1, at = 3)
    UV_rectangles = cbind(as.numeric(c(MS_UV_bar[label_order][METree_GO$order]/5)),rep(0.7,42))
    symbols(4.5+as.numeric(c(MS_UV_bar[label_order][METree_GO$order]/5))/2,c(1:42),
            rectangles = UV_rectangles,add=TRUE, xpd=TRUE,inches = FALSE, adj=1, fg="magenta", bg="magenta",lwd=2)
    mtext("UV signature",las=2,side = 1, at = 4.5)
    survival_rectangles = cbind(as.numeric(c(MS_survival_bar[label_order][METree_GO$order]/5)),rep(0.7,42))
    symbols(6+as.numeric(c(MS_survival_bar[label_order][METree_GO$order]/5))/2,c(1:42),
            rectangles = survival_rectangles,add=TRUE, xpd=TRUE,inches = FALSE, adj=1, fg="orange", bg="orange",lwd=2)
    mtext("survival",las=2,side = 1, at = 6)
    Vem_rectangles = cbind(as.numeric(c(MS_Vem_bar[label_order][METree_GO$order]/5)),rep(0.7,42))
    symbols(7.5+as.numeric(c(MS_Vem_bar[label_order][METree_GO$order]/5))/2,c(1:42),
            rectangles = Vem_rectangles,add=TRUE, xpd=TRUE,inches = FALSE, adj=1, fg="blue", bg="blue",lwd=2)
    mtext("vemurafinib res",las=2,side = 1, at = 7.5)
    rect(9, 5, 10, 5.7, xpd=TRUE, col = "grey", border="grey")
    text(labels="0",x=9,y=6.4); text(labels="5",x=10,y=6.4); text(labels = "p-value",x = 9.5,y =  4); text(labels = "(-log10)", x=9.5, y=3)
    rect(8.75, 2, 10.25, 8)
  })

  output$MExploration <- renderPlot({
    f = VarFilterInput()
    topgenes=pcv[[input$sample]]$gene[f]
    module = as.numeric(input$module)

    if(input$measure==1){
      measure_col = "green3"
      GS = GS_lscore
    }
    if(input$measure==2){
      measure_col = "red3"
      GS = GS_pmet
    }
    if(input$measure==3){
      measure_col = "magenta"
      GS = GS_UV
    }
    if(input$measure==4){
      measure_col = "orange"
      GS = GS_survival
    }
    if(input$measure==5){
      measure_col = "blue"
      GS = GS_Vem
    }

    par(mar=c(5,4,1,1))
    plot(-log10(GS[modules==module,1]), MM[modules==module,module+1], xlab="p-value (-log10)",
         ylab="module membership", pch=20, bty="n",
         cex=(GS[modules==module,4]/max(GS[modules==module,4],na.rm=TRUE))*4,
         xlim=c(0,max(-log10(GS[modules==module,1]),na.rm=TRUE)+1.2),
         ylim=c(0,9.758229e-01), col="darkgrey")
    box(lty=1, col=measure_col, lwd=4)
    if(sum(genes_WGCNA[modules==module] %in% topgenes)>0){
      text(-log10(GS[modules==module,1])[genes_WGCNA[modules==module] %in% topgenes],
           MM[modules==module,module+1][genes_WGCNA[modules==module] %in% topgenes],
           labels=genes_WGCNA[modules==module][genes_WGCNA[modules==module] %in% topgenes],
           cex=1.3, pos=4, col=measure_col)
    }
    points(-log10(GS[modules==module,1][genes_WGCNA[modules==module] %in% topgenes]),
           MM[modules==module,module+1][genes_WGCNA[modules==module] %in% topgenes],
           pch=20, bty="n",
           cex=((GS[modules==module,4]/max(GS[modules==module,4],na.rm=TRUE))*4)[
             genes_WGCNA[modules==module] %in% topgenes],col=measure_col)
    if(input$names==1){
      text(-log10(GS[modules==module,1]), MM[modules==module,module+1],
           labels=genes_WGCNA[modules==module],cex=.8, pos=4)
    }
    abline(v=-log10(0.05),lty=2,col="darkgrey",lwd=2)
  })


})
