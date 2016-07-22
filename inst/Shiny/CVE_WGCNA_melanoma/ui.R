#######################################################
#CANCER VARIANT EXPLORER V1.0 - WGCNAmelanoma extension
#Graphical user interface script
#Andreas Mock
#University of Cambridge
#######################################################

shinyUI(fluidPage(
  titlePanel("Cancer Variant Explorer"),

  #########################
  ####Layout of sidebar####
  #########################
  sidebarLayout(
    sidebarPanel(
      h3("WGCNAmelanoma extension"),
      selectInput("sample", label = h5("Sample ID"), choices=names(v)),
      textOutput("nVariants"),
      tags$hr(),
      h4("ANNOTATION"),
      textOutput("nCoding_SNVs"),
      tags$hr(),
      h4("PRIORITISATION"),
      sliderInput("comb_score", label = "dbNSFP combination score cutoff",
                  min = 0, step = 0.1, max = 4, value = 0),
      checkboxGroupInput("db", label = "Filters",
                         choices = list("exclude 1000 Genomes Project variants"=1,
                                        "include all COSMIC variants"=2,
                                        "include non-SNVs"=3,
                                        "include all DNA repair gene variants"=4)
      ),
      textInput("rescue", label = h5("Rescue genes of interest (comma-seperated HUGO symbols)"),
                value = "Enter gene symbol here ..."),
      tags$hr(),
      h4("TOP TABLE"),
      span(textOutput("top"), style="color:cornflowerblue"),
      downloadButton('downloadTop', 'download top table'),
      tags$hr(),
      h4("Melanoma WGNCA"),
         selectInput("module", label=h5("Co-expression module"),selected = 3,
                     choices = list(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                                    20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42)),
         selectInput("measure", label=h5("Gene significance measure"),selected = 1,
                     choices = list("lymphocyte score" = 1, "primary vs met" = 2, "UV signature" = 3, "survival"=4, "vemurafenib resistance"=5)),
         checkboxInput("names",p("show all gene names")),
      tags$hr(),
      h4("DRUGGABILITY"),
      p("Step 1: go to Druggability tab in main panel"),
      p("Step 2: get DGIdb annotation for variants of interest"),
      p("Step 3: download drug-gene interaction table"),
      downloadButton('downloadDrugs','download drug interactions')),

    ############################
    ####Layout of main panel####
    ############################
    mainPanel(
      tabsetPanel(
        ####ANNOTATION####
        tabPanel("Annotation",
                 h4("Functional consequence annotation from",
                    a("GENCODE", href = "http://www.gencodegenes.org")),
                 fluidRow(

                   column(3,
                          plotOutput("varianttype", height=300, width=150)
                   ),

                   column(9,
                          plotOutput("variantclass", height=300, width=450)
                   )
                 ),
                 tags$hr(),
                 h4("Identification of",
                    a("dbNSFP", href = "https://sites.google.com/site/jpopgen/dbNSFP"),
                    "prediction modules"),
                 fluidRow(

                   column(3,
                          numericInput("nperm", label = h5("number of permutations"),
                                       value = 20)
                   ),
                   column(3,
                          sliderInput("pred_modules", label=h5("number of modules"),
                                      min = 2, step = 1, max = 6, value = 4))
                 ),
                 plotOutput("ConsHM", height=550, width=700)
        ),
        ####PRIORITISATION####
        tabPanel("Prioritisation",
                 fluidRow(
                   column(6,
                          plotOutput("Centroid", height=300, width=300)
                   ),

                   column(6,
                          plotOutput("ThousandG", height=300, width=300)
                   )),
                 tags$hr(),
                 fluidRow(
                   column(6,
                          h5("Overlapping COSMIC variants",
                             verbatimTextOutput("COSMIC"))
                   ),
                   column(6,
                          h5("DNA repair genes",
                             tableOutput("DNArepair"))
                   ))),
        ####TOP TABLE####
        tabPanel("Top table",
                 dataTableOutput("toptable")
        ),
        ####Melanoma WGCNA####
        tabPanel("Melanoma WGCNA",
                 h4("Distribution of variants over melanoma co-expression modules"),
                 textOutput("topgenes_in_modules"),
                 plotOutput("MS",height=650, width=720),
                 plotOutput("MExploration", width=705, height=500)
                 ),
        ####DRUGGABILITY####
        tabPanel("Druggability",
                 h4("Drug-gene interactions from",
                    a("DGIdb",href = "http://dgidb.genome.wustl.edu")),
                 p("Input genes"), verbatimTextOutput("DGIdb_input"),
                 verbatimTextOutput("DGIdb_status"),
                 checkboxInput("DGIdb_anno",p("Get DGIdb annotation")),
                 p("Gene-drug interactions in",
                   a("TEND", href = "http://www.nature.com/nrd/journal/v10/n8/full/nrd3478.html"),
                   "and", a("My Cancer Genome",href = "http://www.mycancergenome.org")),
                 verbatimTextOutput("DGIdb_output"),
                 dataTableOutput("DGIdb_table")
        )
      )
    )
  )
))
