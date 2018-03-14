# besPLOT
# Written by Laura Carroll, Cornell University

# load required packages
library(shiny)
library(ggplot2)
library(vegan)
library(plyr)
library(dplyr)
library(cluster)
library(ggrepel)


# begin ui
ui <- fluidPage(
  
  titlePanel(h1("besPLOT")),
  
  sidebarLayout(
    sidebarPanel(
      
      conditionalPanel(condition = "input.tabs=='Home'",
                       fileInput("file1","Matrix file (n x m matrix)",
                                 multiple=FALSE),
                       helpText("File should have n rows and m columns: 
                                row names in column 1, and column names in row 1"),
                       selectInput("filesep",
                                   label="Matrix delimiter",
                                   choices=c(
                                     'Whitespace ("")'="",
                                     'Tab ("\t")'="\t",
                                     'Space (" ")'=" ",
                                     'Comma (",")'=",",
                                     'Semicolon (";")'=";",
                                     'Colon (":")'=":"), selected = 'Select matrix delimiter'),
                       helpText("Character that separates the columns in your n x m matrix file."),
                       fileInput("file1m","Metadata file (n x 2 matrix)",
                                 multiple=FALSE),
                       helpText("File should have n rows, 2 columns, and no header: column 1 should contain
                                names identical to those in column 1 of your matrix file,
                                and column 2 should contain the corresponding metadata."),
                       selectInput("metasep",
                                   label="Metadata delimiter",
                                   choices=c(
                                     'Whitespace ("")'="",
                                     'Tab ("\t")'="\t",
                                     'Space (" ")'=" ",
                                     'Comma (",")'=",",
                                     'Semicolon (";")'=";",
                                     'Colon (":")'=":"), selected = 'Select matrix delimiter'),
                       helpText("Character that separates the columns in your n x 2 metadata file.")),# end conditional panel
      conditionalPanel(condition = "input.tabs=='Plot'",
                       selectInput("vplot",
                                   label="besPLOT Analyses",
                                   choices=c( 
                                     "Principal Component Analysis (PCA)"=1,
                                     "Non-Metric Multidimensional Scaling (NMDS)"=2))),
      conditionalPanel(condition = "input.vplot=='2' && input.tabs=='Plot'",
                       checkboxInput("nmdsMeta",
                                     label="Overlay Metadata",value = FALSE),
                       checkboxInput("species",
                                     label="Show Column Names", value = FALSE),
                       radioButtons("dmetric",
                                    label="Dissimilarity Index",
                                    choices=list("Binomial"="binomial",
                                                 "Bray-Curtis"="bray",
                                                 "Canberra"="canberra",
                                                 "Cao"="cao",
                                                 "Chao"="chao",
                                                 "Euclidean"="euclidean",
                                                 "Gower"="gower",
                                                 "Horn-Morisita"="horn",
                                                 "Jaccard"="jaccard",
                                                 "Kulczynski"="kulczynski",
                                                 "Mahalanobis"="mahalanobis",
                                                 "Manhattan"="manhattan",
                                                 "Morisita"="morisita",
                                                 "Mountford"="mountford",
                                                 "Raup-Crick"="raup"),
                                    selected="bray")),
      conditionalPanel(condition = "input.vplot=='1' && input.tabs=='Plot'",# && output.listPC",
                       checkboxInput("pcaMeta",
                                     label="Overlay Metadata", value=FALSE),
                       uiOutput("listPC1"),
                       uiOutput("listPC2"),
                       uiOutput("listPC3")),
      
      
      conditionalPanel(condition = "input.tabs!='Home' && input.tabs != 'Metadata'",
                       downloadButton(outputId = "downloadPlot",label="Download this plot"))
                       ),#end sidebarpanel
    
    mainPanel(tabsetPanel(id = "tabs", 
                          tabPanel("Home",htmlOutput("hometext1")),
                          tabPanel("Plot",
                                   #textOutput("yay"),
                          plotOutput("virbarchart",click = "plotclick",dblclick="plotdb",brush=brushOpts(id="plotbrush",resetOnNew = TRUE)),
                          conditionalPanel(condition = "input.tabs=='Plot'",
                          fluidRow(column(width=12,h3("Selected point(s)"),
                          helpText("Click on a point to view name(s) and the associated coordinates."),
                          helpText("The 'dist_' column refers to the distance of a point from coordinates of a mouse click"),
                          helpText("Drag mouse and double-click on plot to zoom in; double-click to reset plot."))),
                          tableOutput("clickinfo"))),
                          tabPanel("Metadata",tableOutput("metadata"))
    )) # end mainPanel
)) # end sidebar layout+ui


# server
server <- function(input, output) {
  
  # raise maximum file upload size
  options(shiny.maxRequestSize=10000*1024^2)
  
  # set xy ranges for plot clicking
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  
  # f'n to create empty plot if data does not exist
  empty.plot<-function(text){
    ggplot()+
      annotate("text", x = 4, y = 25, size=12, label = text) + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank())}
  
  # f'n to run nmds
  run.nmds<-function(vmat,metadata,dmetric,species_scores){
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    vmat<-vmat[rowSums(vmat[, -1])>0,]
    # if metadata is input
    if(!(is.null(metadata))){
      metadata<-data.frame(metadata$V1,as.factor(metadata$V2))
      ordered_meta<-metadata[match(rownames(vmat),metadata[,1]),]
      mymat<-data.frame(ordered_meta[,2],vmat)
      mymds<-metaMDS(comm=mymat[,2:ncol(mymat)],distance=dmetric,k=2,trymax = 10000)
      mymds.scores <- as.data.frame(scores(mymds)) 
      mymds.scores$site <- rownames(mymds.scores)  
      mymds.scores$metadata <- mymat[,1]
      species.scores <- as.data.frame(scores(mymds, "species"))  
      species.scores$species <- rownames(species.scores) 
      hull.data<-NULL
      for(i in 1:length(levels(mymds.scores$metadata))){
        mylevel<-levels(mymds.scores$metadata)[i]
        print("mymds.scores")
        print(mymds.scores)
        mygroup<-mymds.scores[mymds.scores$metadata == mylevel, ][chull(mymds.scores[mymds.scores$metadata==mylevel,c("NMDS1","NMDS2")]),]
        hull.data<-rbind(hull.data,mygroup)}
      nmdsplot<-ggplot()
      nmdsplot<- nmdsplot+ 
        geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,group=metadata,colour=metadata,fill=metadata),alpha=0.30)+
        #geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5,size=3) +  
        geom_point(data=mymds.scores,aes(x=NMDS1,y=NMDS2,colour=metadata,shape="a"),size=2) +
        guides(shape=FALSE) +
        coord_cartesian(xlim = ranges$x, ylim = ranges$y)
      nmds.final<-mymds.scores[,1:2]
      return(list(nmdsplot,nmds.final))}
    # if there is no metadata
    else{
      mymds<-metaMDS(comm=vmat,distance=dmetric,k=2,trymax = 10000)
      mymds.scores <- as.data.frame(scores(mymds)) 
      mymds.scores$site <- rownames(mymds.scores)  
      species.scores <- as.data.frame(scores(mymds, "species"))  
      species.scores$species <- rownames(species.scores) 
      nmdsplot<-ggplot()
      nmdsplot<- nmdsplot+ 
        #geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5,size=3) +  
        geom_point(data=mymds.scores,aes(x=NMDS1,y=NMDS2,colour=factor(mymds.scores$site),shape="a"),size=2) + 
        coord_cartesian(xlim = ranges$x, ylim = ranges$y) +
        theme(legend.position = "none")
      if (!(is.null(species_scores))){
        nmdsplot<-nmdsplot+
        geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5,size=3)}
      nmds.final<-mymds.scores[,1:2]
      return(list(nmdsplot,nmds.final))
    }}
  
  
  # run pca to get total # of PCs
  test.pca<-function(vmat){
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    pca<-prcomp(vmat,scale = TRUE,center = TRUE)
    nPC<-ncol(pca$x)
    return(nPC)}
  
  # run PCA for plotting
  run.pca<-function(vmat,metadata,pc1,pc2,pc3){
    pc1<-as.integer(pc1)
    pc2<-as.integer(pc2)
    pc3<-as.integer(pc3)
    vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
    pca<-prcomp(vmat,scale = TRUE,center = TRUE)
    pca3col <- as.data.frame(cbind(pca$x[,pc1],pca$x[,pc2],pca$x[,pc3]))
    colnames(pca3col) <- c(paste("PC",pc1,sep=""),paste("PC",pc2,sep=""),paste("PC",pc3,sep=""))
    pc1name<-paste("PC",pc1,sep="")
    pc2name<-paste("PC",pc2,sep="")
    pc3name<-paste("PC",pc3,sep="")
    g <- ggplot(pca3col, aes_string(x=pc1name, y=pc2name, size=pc3name)) 
    if (!is.null(metadata)){
      covar<-metadata[match(rownames(pca3col),metadata[,1]),]
      metadata<-covar[,2]
      # take this line out if you ever want to make it continuous
      metadata<-as.factor(metadata)
      g <- g + geom_point(aes(color=metadata)) + xlab(label = paste("PC",pc1,sep = "")) + ylab(label = paste("PC",pc2,sep="")) + scale_size_continuous(name=paste("PC",pc3,sep="")) + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    }
    else{
      g <- g + geom_point(aes(color=rownames(pca3col))) + guides(colour=FALSE) + xlab(label = paste("PC",pc1,sep = "")) + ylab(label = paste("PC",pc2,sep="")) + scale_size_continuous(name=paste("PC",pc3,sep="")) + coord_cartesian(xlim = ranges$x, ylim = ranges$y)
    }
    return(list(g,pca3col))}
  
  # home text
  output$hometext1 <- renderText({
    hometext1<-"<h1>Thanks for using besPLOT!</h1><br><br><h3>Use the panel to your left to upload a n x m matrix and/or any associated categorical metadata.
    After uploading your data, click on the tabs to view, analyze, and interact with your data.</h3>"
    return(hometext1)
  })
  
  # metadata
  output$metadata <- renderTable({
    infile<-input$file1m
    meta.delim <- input$metasep
    if (is.null(infile)){
      return(NULL)
    } 
    if (is.null(input$file1)){
      return(NULL)
    }else{
      mtable<-read.table(infile$datapath,header=FALSE,sep = meta.delim)
      validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (matrix row names in column 1, corresponding metadata in column 2)."))
      finalfiles<-input$file1
      finalsep<-input$filesep
      vmat<-read.table(finalfiles$datapath, header=T, sep=finalsep)
      rownames(vmat)<-vmat[,1]
      vmat<-vmat[,2:ncol(vmat)]
      target<-rownames(vmat)
      validate(need(nrow(mtable)==length(target),"The number of rows in your metadata doesn't match the number of rows in your matrix. Please correct your metadata table and try uploading it again."))
      validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"Names in column 1 of your metadata table should match names in column 1 of your matrix. Please correct your metadata table and try uploading it again."))
      ordered_table<-mtable[match(target,mtable$V1),]
      return(ordered_table)
    }
  })# end metadata
  
  
  output$yay<- renderText({
    #print(input$tabs)
    #print(input$vplot)
    print(input$file1$datapath)
  })
  
  
  # get x-axis PC for PCA (virulence)
  output$listPC1 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    matrix.delim <-input$filesep
    if (is.null(infile)){
      return(NULL)
    } else{
      vmat<-read.table(infile$datapath, header = T, sep=matrix.delim)
      rownames(vmat)<-vmat[,1]
      vmat<-vmat[,2:ncol(vmat)]
      if(vselect=='1'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc1",
                    label="X-axis Principal Component (PC)",
                    choices=nPC,
                    selected=1)
        
      }
    }})# end renderTable
  
  # get y axis PC for PCA (virulence)
  output$listPC2 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    matrix.delim <-input$filesep
    if (is.null(infile)){
      return(NULL)
    } else{
      vmat<-read.table(infile$datapath, header = T, sep=matrix.delim)
      rownames(vmat)<-vmat[,1]
      vmat<-vmat[,2:ncol(vmat)]
      if(vselect=='1'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc2",
                    label="Y-axis Principal Component (PC)",
                    choices=nPC,
                    selected=2)
        
      }
    }})# end renderTable
  
  
  # get z-axis (size) PC for PCA (virulence)
  output$listPC3 <- renderUI({
    infile <- input$file1
    vselect <- input$vplot
    matrix.delim <-input$filesep
    if (is.null(infile)){
      return(NULL)
    } else{
      vmat<-read.table(infile$datapath, header = T, sep=matrix.delim)
      rownames(vmat)<-vmat[,1]
      vmat<-vmat[,2:ncol(vmat)]
      if(vselect=='1'){
        totalPC<-test.pca(vmat = vmat)
        nPC<-c(1:totalPC)
        print(nPC)
        selectInput("pc3",
                    label="Z-axis Principal Component (PC)",
                    choices=nPC,
                    selected=3)
        
      }
    }})# end renderTable
  
  
  # main plot
  besplot<-function(){
    #output$virbarchart <- renderPlot({
      infile <- input$file1
      vselect <- input$vplot
      matrix.delim <-input$filesep
    if (is.null(infile)){
      return(NULL)
    } 
    # split file by new line character
    vmat<-read.table(infile$datapath, header=T, sep=matrix.delim)
    rownames(vmat)<-vmat[,1]
    vmat<-vmat[,2:ncol(vmat)]
    # if NMDS selected
    if (vselect=="2"){
      meta.in<-input$file1m
      dmetric<-input$dmetric
      species_scores<-input$species
      if (is.null(meta.in)){
        vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
        vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric,species_scores=species_scores)
        return(vnmds[[1]])}
      else{
          overlay<-input$nmdsMeta
        if (overlay==FALSE){
          vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
          vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric,species_scores=species_scores)
          return(vnmds[[1]])}
        if (overlay==TRUE){
          meta.delim <- input$metasep
          mtable<-read.table(meta.in$datapath,header=FALSE,sep=meta.delim)
          validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (matrix row names in column 1, corresponding metadata in column 2)."))
          finalfiles<-input$file1
          finalsep<-input$filesep
          vmat<-read.table(finalfiles$datapath, header=T, sep=finalsep)
          rownames(vmat)<-vmat[,1]
          vmat<-vmat[,2:ncol(vmat)]
          target<-rownames(vmat)
          validate(need(nrow(mtable)==length(target),"It looks like the number of rows listed in your metadata doesn't match the number of rows in your matrix file. Please correct your metadata table and try uploading it again."))
          validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"Names in column 1 of your metadata table should match names in column 1 of your matrix. Please correct your metadata table and try uploading it again."))
          ordered_table<-mtable[match(target,mtable$V1),]
          vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
          vnmds<-run.nmds(vmat=vmat,metadata=ordered_table,dmetric = dmetric,species_scores=species_scores)
          return(vnmds[[1]])}
      }}
    # if PCA
    else if (vselect=="1"){
      meta.in<-input$file1m
      meta.delim<-input$metasep
        if(is.null(input$pc1)){
          pc1<-1
        }
        else{
          pc1<-input$pc1
        }
        if(is.null(input$pc2)){
          pc2<-2
        }
        else{
          pc2<-input$pc2
        }
        if(is.null(input$pc3)){
          pc3<-3
        }else{
          pc3<-input$pc3
        }
      if (is.null(meta.in)){
        print(c(pc1,pc2,pc3))
        pcaplot<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]
      }
      else{
        if (input$pcaMeta==FALSE && input$tabs=="Plot"){
          pcaplot<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]}
        else{
          mtable<-read.table(meta.in$datapath,header=FALSE,sep=meta.delim)
          validate(need(ncol(mtable)==2,"Your metadata file needs to have exactly 2 columns (matrix row names in column 1, corresponding metadata in column 2)."))
          finalfiles<-input$file1
          finalsep<-input$filesep
          vmat<-read.table(finalfiles$datapath, header=T, sep=finalsep)
          rownames(vmat)<-vmat[,1]
          vmat<-vmat[,2:ncol(vmat)]
          target<-rownames(vmat)
          validate(need(nrow(mtable)==length(target),"It looks like the number of rows listed in your metadata doesn't match the number of rows in your matrix file. Please correct your metadata table and try uploading it again."))
          validate(need(any(is.na(match(target,mtable$V1)))==FALSE,"Names in column 1 of your metadata table should match names in column 1 of your matrix. Please correct your metadata table and try uploading it again."))
          ordered_table<-mtable[match(target,mtable$V1),]
          ordered_table<-data.frame(ordered_table$V1,ordered_table$V2)
          pcaplot<-run.pca(vmat = vmat,metadata = ordered_table,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[1]]
        }}
      
      return(pcaplot)}
  }# end virbarchart


output$virbarchart <- renderPlot({
  vmat<-read.table(input$file1$datapath, header=T, sep=input$filesep)
  validate(need(ncol(vmat)>1,'Your n x m matrix only has 1 column; please go to the "Home" tab and select the appropriate column delimiter in the "Matrix delimiter" drop-down menu, or reformat your matrix file.'))
  myplot<-besplot()
  print(input$file1)
  print(myplot)
  myplot
  #return(myplot)
})



# get plot coordinates on click
getclicks <- function(){
  if(input$tabs=="Plot"){
    infile <- input$file1
    vselect <- input$vplot
    matrix.delim <-input$filesep}
  if (is.null(infile)){
    return(NULL)
  } 
  vmat<-read.table(infile$datapath, header=T, sep=matrix.delim)
  rownames(vmat)<-vmat[,1]
  vmat<-vmat[,2:ncol(vmat)]
  if (input$tabs=="makeplot"){
    dmetric<-input$dmetric}
  vmat<-vmat[,apply(vmat, 2, var, na.rm=TRUE) != 0]
  if(vselect=='2'){
    dmetric<-input$dmetric
    vnmds<-run.nmds(vmat=vmat,metadata = NULL,dmetric = dmetric,species_scores=NULL)
    mypoints<-vnmds[[2]]}
  else if(vselect=='1'){
      if(is.null(input$pc1)){
        pc1<-1
      }
      else{
        pc1<-input$pc1
      }
      if(is.null(input$pc2)){
        pc2<-2
      }
      else{
        pc2<-input$pc2
      }
      if(is.null(input$pc3)){
        pc3<-3
      }else{
        pc3<-input$pc3
      }
    mypoints<-run.pca(vmat = vmat,metadata = NULL,pc1 = pc1,pc2 = pc2,pc3 = pc3)[[2]]}
  return(mypoints)}# end getclick


# output plot clicks
output$clickinfo <- renderTable({
  mypoints<-getclicks()
  mypoints$Sample<-rownames(mypoints)
  print(head(mypoints))
  print(class(mypoints))
  print(input$plotclick)
  clickinfo<-nearPoints(mypoints,input$plotclick, addDist=TRUE)
  print("clickinfo")
  print(clickinfo)
  return(clickinfo)
})




# observe double click and zoom
observeEvent(input$plotdb, {
  brush <- input$plotbrush
  if (!is.null(brush)) {
    ranges$x <- c(brush$xmin, brush$xmax)
    ranges$y <- c(brush$ymin, brush$ymax)
    
  } else {
    ranges$x <- NULL
    ranges$y <- NULL
  }
})# end double-click

# download plot
output$downloadPlot<-downloadHandler(
  filename<-function(){
    timext<-gsub(" ","_",Sys.time())
    paste("besplot_",timext,".pdf",sep="")},
  content<-function(file){
    pdf(file,onefile = TRUE,paper = "a4r",width = "10",height = "7")
    if(input$tabs=="Plot"){
      print(besplot())}
    dev.off()}
)



}

# Run the application 
shinyApp(ui = ui, server = server)
