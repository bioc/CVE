#################################
#CANCER VARIANT EXPLORER V1.0
#Server script
#Andreas Mock
#University of Cambridge
#################################

#######################
#PART I: PREFIX########
#######################

#protein changing categories of SNVs
pc_cat <- c("De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
            "Frame_Shift_Del", "Frame_Shift_Ins",
            "Missense_Mutation","Nonsense_Mutation","Splice_Site",
            "Start_Codon_SNP")

#create list only containing protein changing variants (pcv)
pcv <- vector("list",length(v))
for (i in 1:length(v)){
  pcv[[i]] <- v[[i]][v[[i]]$Variant_Classification %in% pc_cat,]
}
names(pcv) <- names(v)
if(0 %in% sapply(pcv,nrow)){
  stop("At least one oncotator input file has no protein coding variants.")
}

#retreive rankscores (rs) from variant information
rs_names <- colnames(pcv[[1]])[grep("rankscore",colnames(pcv[[1]]))]
rs_names <- gsub("dbNSFP_",replacement = "",rs_names)
rs_names = gsub("_rankscore",replacement = "",rs_names)
suppressWarnings(
  for(n in 1:length(pcv)){
  rs_matrix = matrix(NA,nrow(pcv[[n]]),length(rs_names))
  for (i in 1:length(rs_names)){
    rs = strsplit(as.vector(pcv[[n]][,grep("rankscore",colnames(v[[1]]))[i]]), "[|]")
    rs = strsplit(pcv[[n]][,grep("rankscore",colnames(v[[1]]))[i]], "[|]")
    rs = as.numeric(sapply(rs,"[", i=1))
    rs_matrix[,i] = rs
  }
  pcv[[n]][[2]] = rs_matrix
  pcv[[n]][[3]] = apply(apply(pcv[[n]][[2]],1,is.na),2,sum)
  colnames(pcv[[n]][[2]]) = rs_names
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

  output$varianttype <- renderPlot({
    ggplot(v[[input$sample]], aes(x=Variant_Type)) +
      geom_bar(stat="count") + xlab("") +
      theme_bw() + ylab("count")
  })

  output$variantclass <- renderPlot({
    var_col = as.vector(v[[input$sample]]$Variant_Classification)
    var_col[var_col %in% pc_cat] = "protein changing"
    var_col[!var_col =="protein changing"] = "other"
    df = as.data.frame(cbind(v[[input$sample]], var_col))
    ggplot(df, aes(x=Variant_Classification, fill=var_col)) +
      geom_bar(stat="count") + coord_flip() +
      xlab("") + theme_bw() +
      guides(fill=guide_legend(reverse=TRUE)) +
      scale_fill_manual(values = c("gold1","blue")) +
      labs(fill="Variant classification")
  })

  output$nCoding_SNVs <- renderText({
    paste("Protein changing single-nucleotide variants:",
          sum(pcv[[input$sample]][[3]]<18))
  })
  #consensus clustering
  ConsClustInput <- reactive({
    ConsClust = vector("list",length(v))
    for (n in 1:length(v)){
      M_temp = pcv[[n]][[2]][!sapply(
        1:nrow(pcv[[n]][[2]]),function(x) TRUE %in% is.na(pcv[[n]][[2]])[x,]),]
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
  })

  #combination score of rankscores
  CombScoreInput <- reactive({
    ConsClust = ConsClustInput()
    df = as.data.frame(t(pcv[[input$sample]][[2]]),stringsAsFactors = FALSE)
    centroid = t(ddply(df, .(ConsClust[[input$sample]][[as.numeric(
      input$pred_modules)]]$consensusClass), colwise(mean, na.rm = TRUE)))
    colnames(centroid) = paste0("module",1:ncol(centroid))
    centroid = centroid[2:nrow(centroid),]
    centroid = apply(centroid,2,as.numeric)
    cutoff = 0.75
    comb_score = as.vector(apply(centroid*(centroid>cutoff*1),1,
                                 function(x) sum(x,na.rm=TRUE)))
    comb_score
  })

  observe({
    x = as.numeric(input$pred_modules)
    updateSliderInput(session, inputId = "comb_score",
                      value = 0, min = 0, max =x,step = 0.1)
  })

  output$ConsHM <- renderPlot({
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
    legend("bottomleft",legend = c(
      "conservation score","functional prediction score","ensemble score"),
      fill=c("yellow","orange", "blueviolet"),
      border =c("yellow", "orange", "blueviolet"),xpd=TRUE, cex=0.9)
    legend("topright",legend = paste0("cluster ",1:as.numeric(input$pred_modules)),
           fill=ConsClust[[input$sample]][[as.numeric(input$pred_modules)]]$clrs[[3]],
           border=ConsClust[[input$sample]][[as.numeric(
             input$pred_modules)]]$clrs[[3]], xpd=TRUE, cex=0.9)
  })


  ####PRIORITISATION####
  output$Centroid <- renderPlot({
    comb_score = CombScoreInput()
    cutoff_boolean = comb_score[order(comb_score,decreasing = TRUE)]>(
      as.numeric(input$comb_score)-0.000001)
    colors = rep("cornflowerblue",length(comb_score))
    colors[cutoff_boolean==FALSE] = "brown"
    plot(comb_score[order(comb_score,decreasing = TRUE)],
         xlab=c("variant (decreasing order)"), ylab=c("dbNSFP combination score"),
         col=colors, main="dbNSFP combination score", pch=20)
    abline(h = as.numeric(input$comb_score),lty=2, lwd=2, col="brown")
    legend("topright",legend=c("included","excluded"),col=c("cornflowerblue","brown"),pch=20)
  })

  output$ThousandG <- renderPlot({
    ThousandG = pcv[[input$sample]]$X1000gp3_AF
    comb_score = CombScoreInput()
    cutoff_boolean = comb_score<as.numeric(input$comb_score)
    col1 = rep("cornflowerblue", sum(!is.na(ThousandG)))
    col1[cutoff_boolean[!is.na(ThousandG)]] = "brown"
    if(1 %in% input$db){col1 = rep("brown", sum(!is.na(ThousandG)))}
    plot(ThousandG[!is.na(ThousandG)][
      order(ThousandG[!is.na(ThousandG)],decreasing=TRUE)],xlab=c("variant"),
      main=paste0(sum(!is.na(ThousandG))," variants in 1000 Genomes Project"),
      col=col1, ylab="AF",pch=19)
    legend("topright",legend=c("included","excluded"),
           col=c("cornflowerblue","brown"),pch=20)
  })

  VarFilterInput <- reactive({
    comb_score = CombScoreInput()
    f = (comb_score>=input$comb_score)*1 #keep if higher than comb.score cutoff
    f[pcv[[input$sample]][[3]]==18] = 0 #exclude variants with missing rankscores
    if (1 %in% input$db){ #if selection 1 - exclude all 1000 Genomes variants
      f[pcv[[input$sample]]$X1000gp3_AF>0] = 0
    }
    if (2 %in% input$db){ #keep all overlapping mutations in COSMIC if selected
      f[as.numeric(pcv[[input$sample]]$COSMIC_n_overlapping_mutations)>0] =1
    }
    if (3 %in% input$db){ #keep all non-SNVs
      nonSNVs = pcv[[input$sample]]$Variant_Type
      f[!nonSNVs=="SNP"] =1
    }
    if (4 %in% input$db){ #keep all DNA repair gene variants
      repair = pcv[[input$sample]]$DNARepairGenes_Role
      repair[repair==""] = NA
      is.repair = !is.na(repair)
      f[is.repair] =1
    }
    sel_genes = unique(pcv[[input$sample]]$Hugo_Symbol)
    gene_rescue = unlist(strsplit(input$rescue,", "))
    gene_rescue_intersect = sel_genes[sel_genes %in% gene_rescue]
    f[pcv[[input$sample]]$Hugo_Symbol %in% gene_rescue_intersect] =1
    f = f==1
    f
  })

  output$COSMIC <- renderPrint({
    f = VarFilterInput()
    sum((as.numeric(pcv[[input$sample]]$COSMIC_n_overlapping_mutations[f])>0)*1)
  })

  output$DNArepair <- renderTable({
    f = VarFilterInput()
    repair = pcv[[input$sample]]$DNARepairGenes_Role
    repair[repair==""] = NA
    is.repair = (((!is.na(repair))*1) *(f*1)) == 1
    if(sum(is.repair)>0){
      data.frame(gene=pcv[[input$sample]]$Hugo_Symbol[is.repair],role=repair[is.repair])
    }
  })

  ####TOP TABLE####
  output$top = renderText({
    f = VarFilterInput()
    paste0(sum(f)," variants prioritised for further exploration")
  })

  output$toptable = renderDataTable({
    f = VarFilterInput()
    comb_score = CombScoreInput()
    topgenes=pcv[[input$sample]]$Hugo_Symbol[f]
    genes = pcv[[input$sample]]$Hugo_Symbol[f]
    data.frame(gene=genes,
               variant=pcv[[input$sample]]$Protein_Change[f],
               type=pcv[[input$sample]]$Variant_Type[f],
               dbNSFP_score=comb_score[f],
               COSMIC_variant = pcv[[input$sample]]$COSMIC_overlapping_mutations[f]

    )
  })

  output$downloadTop <- downloadHandler(
    filename = function() {
      paste('CVE_top_table_', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      f = VarFilterInput()
      comb_score = CombScoreInput()
      topgenes=pcv[[input$sample]]$Hugo_Symbol[f]
      genes = pcv[[input$sample]]$Hugo_Symbol[f]
      toptable = data.frame(gene=genes,
                 variant=pcv[[input$sample]]$Protein_Change[f],
                 type=pcv[[input$sample]]$Variant_Type[f],
                 dbNSFP_score=comb_score[f],
                 COSMIC_variant = pcv[[input$sample]]$COSMIC_overlapping_mutations[f])
      writeLines(paste("#Cancer Variant Explorer top table output |", Sys.Date(),"| sample file:",
                 names(v),"| number of dbNSFP modules",
                 input$pred_modules,"| dbNSFP combination score cutoff",
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
    paste(as.vector(pcv[[input$sample]]$Hugo_Symbol[f_which]),collapse=", ")
  })

  JSON = reactive({
    if(input$DGIdb_anno){
      f = VarFilterInput()
      f_which = which(f==1)
      genes = pcv[[input$sample]]$Hugo_Symbol[f_which]
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
                       names(v),"| number of dbNSFP modules",
                       input$pred_modules,"| dbNSFP combination score cutoff",
                       input$comb_score,"| filters: exclude 1000 Genomes:",
                       1 %in% input$db, ", include all COSMIC:", 2 %in% input$db,
                       ", include non-SNVs:", 3 %in% input$db, ", include DNA repair genes:",
                       4 %in% input$db, " | rescued genes:", input$rescue),
                 file)
      write.table(drugtable, file,append=TRUE,sep=",")
    }
  )

})
