library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(plotly)
library(tidyr)
library(DT)
library(shinycssloaders)
library(ontologyIndex)
library(ontologyPlot)
library(visNetwork)
library(data.table)

# Load keys
cell_key <- read.table("../celltype_key.txt", header = TRUE)
rownames(cell_key) <- cell_key$id

gene_key <- read.table("../gene_key.txt", header = TRUE)
rownames(gene_key) <- gene_key$feature_id

all_tissues <- read.table("../all_tissues.txt", sep = "\n", stringsAsFactors = FALSE)
unique_tissues <- c("All", all_tissues$V1)

unique_genes <- unique(gene_key$feature_name)

# Load the ontology file
obo_file <- "../cl-basic.obo"
ob <- get_ontology(obo_file, propagate_relationships = c("is_a"), extract_tags = "everything")
all_cellxgene_terms <- names(ob$subset)[grep("cellxgene", ob$subset)]

# Function to get top correlated genes using Spearman correlation
get_top_correlated_genes <- function(expression_data, target_gene, top_n) {
  wide_data <- expression_data %>%
    select(celltype, gene_symbol, average_expression) %>%
    pivot_wider(names_from = gene_symbol, values_from = average_expression, values_fill = list(average_expression = 0))
  
  target_data <- wide_data[, target_gene, drop = FALSE]
  other_genes <- setdiff(names(wide_data)[-1], target_gene)
  other_genes <- other_genes[which(colSums(wide_data[, other_genes]) > 0)]
  correlations <- sapply(other_genes, function(gene) {
    cor(target_data, wide_data[, gene, drop = FALSE], method = "spearman")
  })
  
  top_correlated <- names(sort(correlations, decreasing = TRUE))[1:top_n]
  top_genes <- c(target_gene, top_correlated)
  
  return(list(top_genes = top_genes, correlations = correlations))
}

# Function to get correlations between genes
get_gene_correlations <- function(expression_data, genes) {
  wide_data <- expression_data %>%
    select(celltype, gene_symbol, average_expression) %>%
    pivot_wider(names_from = gene_symbol, values_from = average_expression, values_fill = list(average_expression = 0))
  
  cor_matrix <- cor(wide_data[, genes, drop = FALSE], method = "spearman", use = "everything")
  return(cor_matrix)
}

# Function to create a heatmap
create_heatmap <- function(cor_matrix) {
  cor_df <- as.data.frame(as.table(cor_matrix))
  colnames(cor_df) <- c("Gene1", "Gene2", "Correlation")
  
  p <- ggplot(cor_df, aes(Gene1, Gene2, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Correlation Heatmap", x = "Gene", y = "Gene")
  
  ggplotly(p)
}

# Function to get ontology plot
get_ontology_tissue <- function(ob, tissue, gene_symbol) {
  rownames(gene_key) <- gene_key$feature_name
  gene <- gene_key[paste(gene_symbol), 1]
  rownames(gene_key) <- gene_key$feature_id
  tissue <- gsub(" ", "_", tissue)
  
  cl_files <- list.files(".", "CL", full.names = F)
  cl_files <- cl_files[grep(tissue, cl_files)]
  
  all_cl <- list()
  for (f in cl_files) {
    cur <- fread(f, sep = ",")
    cur <- as.data.frame(cur)
    cl <- gsub(".txt", "", f)
    cl <- gsub("_", ":", cl)
    cl <- gsub("-", "", cl)
    cl <- gsub(gsub("_", ":", tissue), "", cl)
    
    cur$cl <- cl
    all_cl[[cl]] <- cur
  }
  
  all <- do.call(rbind, all_cl)
  
  terms <- unique(all$cl)
  terms <- terms[which(terms != "unknown")]
  
  all_terms <- ontologyIndex::get_ancestors(terms, ontology = ob)
  all_terms <- all_terms[which(all_terms %in% all_cellxgene_terms)]
  
  all_terms <- unique(c(terms, all_terms))
  
  av_exp_matrix <- all %>%
    select(gene, cl, average_expression) %>%
    pivot_wider(names_from = cl, values_from = average_expression, values_fill = list(average_expression = 0))
  
  av_exp_matrix <- as.data.frame(av_exp_matrix)
  rownames(av_exp_matrix) <- av_exp_matrix$gene
  av_exp_matrix <- av_exp_matrix[,-1]
  
  pct_matrix <- all %>%
    select(gene, cl, pct.expressed) %>%
    pivot_wider(names_from = cl, values_from = pct.expressed, values_fill = list(pct.expressed = 0))
  
  pct_matrix <- as.data.frame(pct_matrix)
  rownames(pct_matrix) <- pct_matrix$gene
  pct_matrix <- pct_matrix[,-1]
  
  ex_ex <- av_exp_matrix[gene, terms]
  scaled_ex_ex <- scale(as.numeric(ex_ex))
  rownames(scaled_ex_ex) <- terms
  
  express_table <- cbind(scaled_ex_ex, t(ex_ex), t(pct_matrix[gene, terms]))
  colnames(express_table) <- c("Average Expression (scaled)", "Average Expression", "% expressed")
  
  cols <- c("lightgrey", "blue")
  color_palette <- colorRampPalette(cols)
  
  scaled_colors <- color_palette(length(scaled_ex_ex))[as.numeric(cut(scaled_ex_ex, breaks = length(scaled_ex_ex)))]
  scaled_colors <- c(scaled_colors, rep('white', length(all_terms) - length(terms)))
  
  term_colors <- setNames(scaled_colors, all_terms)
  
  op <- onto_plot(ob, terms = all_terms, fillcolor = term_colors, color = "grey", title = gene, main = gene)
  
  adj_matrix <- op$adjacency_matrix
  nodes <- as.data.frame(op$node_attributes)
  edges <- as.data.frame(op$edge_attributes)
  nodes$color <- term_colors[rownames(adj_matrix)]
  nodes$id <- paste(ob$name[paste(rownames(nodes))])
  rownames(adj_matrix) <- paste(ob$name[paste(rownames(adj_matrix))])
  colnames(adj_matrix) <- paste(ob$name[paste(colnames(adj_matrix))])
  
  edges <- data.frame()
  for (term in rownames(adj_matrix)) {
    connected_terms <- names(which(adj_matrix[term, ] == 1))
    for (connected_term in connected_terms) {
      edges <- rbind(edges, data.frame(from = term, to = connected_term, stringsAsFactors = FALSE))
    }
  }
  
  nodes$color <- nodes$fillcolor
  
  nodes$value <- 25
  nodes$pct.exp <- NA
  nodes$pct.exp[which(rownames(nodes) %in% rownames(express_table))] <- 
    express_table[rownames(nodes)[which(rownames(nodes) %in% rownames(express_table))], 3]
  
  nodes$shape <- "dot"
  nodes$shape[which(nodes$color == "white")] <- "square"
  nodes$size <- nodes$pct.exp
  
  v <- visNetwork(nodes, edges, main = gene_symbol) %>%
    visNodes(opacity = .8, borderWidth = 1) %>%
    visEdges(arrows = "to", color = list(color = "gray", highlight = "blue")) %>%
    visIgraphLayout(layout = 'layout_as_tree') %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visLayout(randomSeed = 34, improvedLayout = TRUE)
  
  list(op = op, v = v)
}

# UI
ui <- fluidPage(
  titlePanel("Gene Expression Across Cell Types"),
  
  tabsetPanel(
    tabPanel("Multiple Gene Input",
             sidebarLayout(
               sidebarPanel(
                 selectInput("tissue_multiple", "Tissue:", choices = unique_tissues, selected = "eye"),
                 selectizeInput("genes_multiple", "Genes:", choices = NULL, selected = c('APOE','MAL'), multiple = TRUE, options = list(maxOptions = 100)),
                 radioButtons("sort_celltypes_multiple", "Sort Cell Types By:", choices = c("Hierarchical Clustering" = "hc")),
                 radioButtons("sort_genes_multiple", "Sort Genes By:", choices = c("Input Order" = "input", "Hierarchical Clustering" = "hc")),
                 actionButton("submit_multiple", "Submit")
               ),
               mainPanel(
                 plotlyOutput("expressionPlot_multiple", width = 900, height = 500) %>% withSpinner(),
                 plotlyOutput("correlationHeatmap_multiple",height = 500, width = 500) %>% withSpinner()
               )
             )
    ),
    tabPanel("Single Gene Input",
             sidebarLayout(
               sidebarPanel(
                 selectInput("tissue_single", "Tissue:", choices = unique_tissues, selected = "eye"),
                 selectizeInput("gene_single", "Gene:", choices = NULL, selected = 'APOE', options = list(maxOptions = 100)),
                 actionButton("submit_single", "Submit")
               ),
               mainPanel(
                 DTOutput("correlationTable_single") %>% withSpinner(),
                 downloadButton("downloadData", "Download Correlation Table"),
                 plotlyOutput("expressionPlot_single", width = 900, height = 500) %>% withSpinner(),
                 visNetworkOutput("ontologyPlot_single") %>% withSpinner(),
                 downloadButton("downloadOntologyPlot", "Download Ontology Plot"),
                 plotOutput("ontologyPlot_op_single",width = 1000, height =1200 ) %>% withSpinner()
               )
             )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Update the selectize input for genes with server-side processing
  updateSelectizeInput(session, 'genes_multiple', choices = unique_genes, server = TRUE)
  updateSelectizeInput(session, 'gene_single', choices = unique_genes, server = TRUE)
  
  observeEvent(input$submit_multiple, {
    req(input$tissue_multiple, input$genes_multiple)
    
    tissue <- input$tissue_multiple
    genes <- input$genes_multiple
    
    tissue <- gsub(" ", "_", tissue)
    
    print(paste("Selected tissue:", tissue))
    print(paste("Selected genes:", paste(genes, collapse = ", ")))
    print("Loading expression data...")
    
    # Load the expression data
    expression_files <- list.files(".", full.names = TRUE)
    if (tissue != "All") {
      expression_files <- expression_files[grep(tissue, expression_files)]
    } else {
      expression_files <- expression_files[grep("-All", expression_files)]
    }
    
    exp_list <- lapply(expression_files, function(f) {
      cur <- read.csv(f)
      cur$gene_symbol <- gene_key[paste(cur$gene), "feature_name"]
      ct <- sapply(strsplit(as.character(f), "-"), "[[", 1)
      ct <- gsub("./", "", ct)
      ct <- gsub("_", ":", ct)
      cur$celltype <- cell_key[paste(ct), 1]
      cur
    })
    
    expression_data <- do.call(rbind, exp_list)
    expression_data$average_expression <- log(expression_data$average_expression + 1)
    
    print("Expression data loaded.")
    print("Filtering data for selected genes...")
    
    # Filter the data based on the user input
    filtered_data <- expression_data[which(expression_data$gene_symbol %in% genes),]
    
    print("Data filtered.")
    print("Sorting cell types and genes...")
    
    # Sort cell types and genes
    if (input$sort_celltypes_multiple == "hc" && length(genes) > 2) {
      # Reshape data to have cell types as rows and genes as columns
      reshaped_data <- filtered_data %>%
        select(celltype, gene_symbol, average_expression) %>%
        pivot_wider(names_from = gene_symbol, values_from = average_expression, values_fill = list(average_expression = 0))
      
      # Perform hierarchical clustering on the reshaped data
      dist_matrix <- dist(reshaped_data[, -1])  # Calculate distance matrix
      hc <- hclust(dist_matrix)  # Perform hierarchical clustering
      sorted_celltypes <- reshaped_data$celltype[hc$order]  # Get ordered cell types
    } else {
      sorted_celltypes <- unique(filtered_data$celltype)
    }
    
    if (input$sort_genes_multiple == "hc" && length(genes) > 2) {
      # Reshape data to have genes as rows and cell types as columns
      reshaped_data <- filtered_data %>%
        select(gene_symbol, celltype, average_expression) %>%
        pivot_wider(names_from = celltype, values_from = average_expression, values_fill = list(average_expression = 0))
      
      # Perform hierarchical clustering on the reshaped data
      dist_matrix <- dist(reshaped_data[, -1])  # Calculate distance matrix
      hc <- hclust(dist_matrix)  # Perform hierarchical clustering
      sorted_genes <- reshaped_data$gene_symbol[hc$order]  # Get ordered genes
    } else {
      sorted_genes <- genes
    }
    
    filtered_data$celltype <- factor(filtered_data$celltype, levels = sorted_celltypes)
    filtered_data$gene_symbol <- factor(filtered_data$gene_symbol, levels = sorted_genes)
    
    print("Cell types and genes sorted.")
    
    output$expressionPlot_multiple <- renderPlotly({
      print("Rendering expression plot...")
      p <- ggplot(filtered_data, aes(x = celltype, y = gene_symbol)) +
        geom_point(aes(size = pct.expressed, color = average_expression)) +
        scale_color_gradient(low = "lightgrey", high = "blue") +
        labs(title = paste("Expression in", tissue),
             x = "Cell Type", y = "Gene", size = "% Expressed", color = "Average Expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "right",   # Position legends at the top
              legend.box = "vertical")  # Arrange legends horizontally
      ggplotly(p) %>% layout(yaxis = list(automargin = TRUE))
    })
    
 
    output$correlationHeatmap_multiple <- renderPlotly({
      print("Rendering correlation heatmap...")
      correlations <- get_gene_correlations(expression_data, genes)
      create_heatmap(correlations)
    })
  })
  
  observeEvent(input$submit_single, {
    req(input$tissue_single, input$gene_single)
    
    tissue <- input$tissue_single
    tissue <- gsub(" ", "_", tissue)
    
    gene <- input$gene_single
    top_n <- 8
    
    print(paste("Selected tissue:", tissue))
    print(paste("Selected gene:", gene))
    print("Loading expression data...")
    
    # Load the expression data
    expression_files <- list.files(".", full.names = TRUE)
    if (tissue != "All") {
      expression_files <- expression_files[grep(tissue, expression_files)]
    } else {
      expression_files <- expression_files[grep("-All", expression_files)]
    }
    
    exp_list <- lapply(expression_files, function(f) {
      cur <- read.csv(f)
      cur$gene_symbol <- gene_key[paste(cur$gene), "feature_name"]
      ct <- sapply(strsplit(as.character(f), "-"), "[[", 1)
      ct <- gsub("./", "", ct)
      ct <- gsub("_", ":", ct)
      cur$celltype <- cell_key[paste(ct), 1]
      cur
    })
    
    expression_data <- do.call(rbind, exp_list)
    expression_data$average_expression <- log(expression_data$average_expression + 1)
    
    print("Expression data loaded.")

    
    # order the celltypes by input gene expression
    gene_ex = expression_data[which(expression_data$gene_symbol==gene),]
    gene_ex = gene_ex[order(gene_ex$average_expression,decreasing = T),]
    
    sorted_celltypes = paste(unique(gene_ex$celltype))
    
    expression_data$celltype <- factor(expression_data$celltype, levels = sorted_celltypes)
    
    print("Cell types sorted.")
    
    output$correlationTable_single <- renderDT({
      print("Rendering correlation table...")
      top_genes_data <- get_top_correlated_genes(expression_data, gene, top_n)
      correlations <- top_genes_data$correlations
      correlation_table <- data.frame(Gene = names(correlations), Correlation = correlations)
      correlation_table = correlation_table[order(correlation_table$Correlation,decreasing = T),]
      datatable(correlation_table)
    })
    
    
    output$expressionPlot_single <- renderPlotly({
      print("Rendering expression plot...")
      p <- ggplot(expression_data[which(expression_data$gene_symbol %in% top_cor_genes),], aes(x = celltype, y = gene_symbol)) +
        geom_point(aes(size = pct.expressed, color = average_expression)) +
        scale_color_gradient(low = "lightgrey", high = "blue") +
        labs(title = paste("Expression of", gene, "and top correlated genes in", tissue),
             x = "Cell Type", y = "Gene", size = "% Expressed", color = "Average Expression") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = "right",   # Position legends at the top
              legend.box = "vertical")  # Arrange legends horizontally
      ggplotly(p) %>% layout(yaxis = list(automargin = TRUE))
    })
    
    ontology_plots <- get_ontology_tissue(ob, tissue, gene)
    
    output$ontologyPlot_op_single <- renderPlot({
      print(ontology_plots$op)
    })
    
    output$ontologyPlot_single <- renderVisNetwork({
      print("Rendering ontology plot...")
      ontology_plots$v
    })
    
    output$downloadOntologyPlot <- downloadHandler(
      filename = function() { paste("ontology_plot_",tissue,"_", gene, ".pdf", sep = "") },
      content = function(file) {
        pdf(file, width = 10, height = 10)
        print(ontology_plots$op)
        dev.off()
      }
    )
    
    output$downloadData <- downloadHandler(
      filename = function() { paste("correlations_", gene, ".csv", sep = "") },
      content = function(file) {
        write.csv(data.frame(Gene = names(correlations), Correlation = correlations), file, row.names = FALSE)
      }
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)

