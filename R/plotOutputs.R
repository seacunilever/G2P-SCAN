# Copyright (C) 2022 as Unilever Global IP Limited
# This file is part of G2P-SCAN.
# G2P-SCAN is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# G2P-SCAN is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
# You should have received a copy of the GNU Lesser General Public License along with G2P-SCAN. If not, see https://www.gnu.org/licenses/.

#' Return analysis inputs and pathway summary
#'
#' @description Take data of first sheet of G2P-SCAN excel outputs (input summary) and output the following:
#' \itemize{
#' \item{1. n_gene_searched: the number on input genes run in the analysis}
#' \item{2. n_genes_mapped: the number of input genes mapped to Reactome pathways }
#' \item{3. n_genes_not_mapped: the number of input genes not mapped to Reactome pathways}
#' \item{4. n_genes_not_found: number of input genes not found in Reactome}
#' \item{5. species_run: the selected species run in the analysis}
#' \item{6. ortho: the selected orthologue filter applied in the analysis}
#' \item{7. pathway: the selected pathways filters applied in the analysis}
#' }
#' @param input_summary data.frame of read in 'input summary' sheet from G2P-SCAN output
#' @importFrom logger log_fatal log_warn
#' @importFrom stringi stri_extract_first_regex
#' @return a list of summarise input and analysis run data
#' @export
get_analysis_summary <- function(input_summary) {
  input_gene_row <- which(grepl("Input Gene[(]s[)] Searched", input_summary[, 1]))
  mapped_gene_row <- which(grepl("Gene[(]s[)] Mapped To Pathways", input_summary[, 1]))
  genes_not_mapped_row <-  which(grepl("Genes[(]s[)] Not Mapped To Pathways", input_summary[,1]))
  genes_not_found_row <-  which(grepl("Genes[(]s[)] Not Found", input_summary[, 1])) 
  
  n_genes_searched <- as.numeric(stringi::stri_extract_first_regex(input_summary[input_gene_row, 1], "[0-9]+"))
  n_genes_mapped <- as.numeric(stringi::stri_extract_first_regex(input_summary[mapped_gene_row, 1], "[0-9]+"))
  n_genes_not_mapped <-as.numeric(stringi::stri_extract_first_regex(input_summary[genes_not_mapped_row, 1], "[0-9]+"))
  n_genes_not_found <- as.numeric(stringi::stri_extract_first_regex(input_summary[genes_not_found_row, 1], "[0-9]+"))
  
  species_row <- which(grepl("Species Run", input_summary[, 3]))
  ortho_row <- which(grepl("Orthologue filter", input_summary[, 3]))
  species_run <- input_summary[species_row+1:(ortho_row-2), 3]
  ortho <-  input_summary[ortho_row+1, 3]
  pathway_row <- which(grepl("Pathway Levels", input_summary[, 3]))
  pathway <- na.omit(input_summary[pathway_row+1:length(input_summary[, 3]), 3])[1]
  return(list("n_genes_searched" = n_genes_searched,
              "n_genes_mapped" = n_genes_mapped,
              "n_genes_not_mapped" = n_genes_not_mapped,
              "n_genes_not_found" = n_genes_not_found,
              "species_run" = species_run,
              "ortho" = ortho,
              "pathway" = pathway))
}


#' Convert counts to percentages
#'
#' @description takes in data.frame of count_summary from G2P-SCAN counts file and calculates gene, entity and reaction percentages relative to human for each pathway and species.
#' @param count_summary data.frame of read in 'count summary' sheet from G2P-SCAN counts output
#' @importFrom stringr str_replace_all
#' @return data.frame of count_summary with additional species-percent columns
#' @export
convert_to_percent <- function(count_summary) {
  #copy 
  count_percent <- count_summary
  #genes
  gene_columns <- colnames(count_summary)[grepl("Gene.Count", colnames(count_summary))]
  count_percent[, gsub("Count", "Gene_Percent", gene_columns)] <- count_summary[, gene_columns] / count_summary$Human_Total.Gene.Count *100
  
  # Entities 
  ent_columns <- colnames(count_summary)[grepl("Entity", colnames(count_summary))]
  count_percent[, gsub("Entity", "Entity_Percent", ent_columns)] <- count_summary[, ent_columns] / count_summary$Human_Entity.Count *100
  
  # Reactions 
  react_columns <- colnames(count_summary)[grepl("Reaction", colnames(count_summary))]
  count_percent[, gsub("Reaction", "Reaction_Percent", react_columns)] <- count_summary[, react_columns] / count_summary$Human_Reaction.Count *100
  colnames(count_percent) <- stringr::str_replace_all(colnames(count_percent), c("R.norvegicus" = "Rat",
                                                                        "M.musculus" = "Mouse",
                                                                        "D.rerio"  = "Zebrafish",
                                                                        "C.elegans" = "Nematode",
                                                                        "D.melanogaster" = "Fruitfly",
                                                                        "S.cerevisiae" = "Yeast"))
  return(count_percent)
}

#' Plot bubble plots for gene, entities or reactions
#'
#' @description Takes count percentages (produced by convert_to_percent function) and creates a plot displaying percents by bubble size per pathways and species. Bar charts to of pathway sizes can be appended to bubble plot or created separately.
#' This function can only plot gene, entity or reaction percentages one at a time. 
#' Due to varying number of pathways outputted from a single G2P-SCAN analysis, number of pathways plotted in one plot is customisable - use 'n_pathways_per_plot' to set this size. If n_pathways_per_plot is lower than total number of pathways, multiple plots will be created and saved, iterating through all pathways.
#' 
#' @param count_percent data.frame output of convert_to_percent function
#' @param count_type 'genes', 'entities' or 'reactions' indicating what percent type to plot
#' @param n_genes_searched number of genes run in the analysis (output from get_analysis_summary)
#' @param data_name string prefixing output png (recommended to match G2P-SCAN file output names)
#' @param join_plots TRUE for concatenating bubble and bar plot together, FALSE to produce 2 separate plots
#' @param n_pathways_per_plot integer of number of pathways to display in one plot.
#' @param height height of output png
#' @param width width of output png
#' @inheritParams run_all_plots
#' @importFrom stringr str_replace_all
#' @importFrom logger log_fatal
#' @import ggplot2
#' @importFrom ggpubr ggarrange
#' @importFrom stringr str_wrap
#' @return ggplot objects of combined plot or individual plots given in a list
#' @export
plot_percent_bubbles <- function(count_percent, count_type = "genes", plot_output_dir = "./", n_genes_searched, data_name, join_plots = TRUE, n_pathways_per_plot = 22, height = 12, width = 15) {
  dir.create(plot_output_dir, showWarnings = F, recursive = T)
  if(tolower(count_type) == "genes"){
    plot_data <- count_percent[, grepl("Pathway", colnames(count_percent)) |
                                 grepl("Gene_Percent", colnames(count_percent)) |
                                 grepl("Input.Found.Count", colnames(count_percent))]
    human_column <- "Human_Total.Gene.Count"
  } else if(tolower(count_type) == "entities") {
    plot_data <- count_percent[, grepl("Pathway", colnames(count_percent)) |
                                 grepl("Entity_Percent", colnames(count_percent)) |
                                 grepl("Input.Found.Count", colnames(count_percent))]
    human_column <- "Human_Total.Gene.Count"
  } else if(tolower(count_type) == "reactions") {
    plot_data <- count_percent[, grepl("Pathway", colnames(count_percent)) |
                                 grepl("Reaction_Percent", colnames(count_percent)) |
                                 grepl("Input.Found.Count", colnames(count_percent))]
    human_column <- "Human_Total.Gene.Count"
  } else {
    logger::log_fatal("count_type must be one of 'genes', 'entities' or 'reactions')")
    stop()
  }
  
  species <- unique(lapply(strsplit(colnames(plot_data), "_"), "[[", 1))
  species <- species[!grepl("Pathway", species)]
  pathways <- unique(plot_data$Pathway.Name)
  table_for_plot <- data.frame()
  for(s in species) {
    if(s != "Human") {
      for(p in pathways) {
        s_data <- plot_data[grepl(s, colnames(plot_data))]
        input_found <- max(s_data[,grepl("Input.Found.Count", colnames(s_data)) ])
        s_p_data <- plot_data[plot_data$Pathway.Name == p, grepl(s, colnames(plot_data))]
        human_genes <- count_percent[count_percent$Pathway.Name == p, human_column]
        
        row <- data.frame("Pathway" = p,
                          "Species" = paste0(s, " (", input_found, "/", n_genes_searched, ")"),
                          "Human_path_size" = human_genes,
                          "Log_path_size" = log10(human_genes),
                          "Values" = s_p_data[,grepl("Percent", colnames(s_p_data)) ])
        table_for_plot <- rbind(table_for_plot, row)
      }
    }
  }
  table_for_plot$Species<-factor(table_for_plot$Species, levels = unique(table_for_plot$Species))
  table_for_plot$Values <- round(table_for_plot$Values, 0)
  np <- length(pathways)
  n_plots <- ceiling(np/n_pathways_per_plot) 
  xn <- 1
  for(n in 1:n_plots) {
    pathways_sub <- pathways[xn:(n_pathways_per_plot*n)]
    xn <- xn +n_pathways_per_plot
    table_for_plot_sub <- table_for_plot[table_for_plot$Pathway %in% pathways_sub,]
    p <- ggplot(table_for_plot_sub, aes(x = Species, y = Pathway, colour = Species, size = Values**2)) +
      geom_point(aes(x = Species, y = Pathway)) +
      geom_text(aes(label = Values), colour="black", size = 4) +
      scale_x_discrete(position="top",labels = function(x) stringr::str_wrap(x, width = 1)) +
      scale_y_discrete(labels = function(y) stringr::str_wrap(y, width = 50)) +
      scale_size_continuous(range = c(1,10)) +
      scale_color_brewer(palette = "Set2") +
      labs(x = NULL, y = NULL) +
      theme(legend.position = "none",
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x=element_text(size=10,angle=55,vjust=0.5,hjust=0.1),
            axis.text.y=element_text(size=10),
            text = element_text(size=10))
    
    p2 <- ggplot(unique(table_for_plot_sub[grepl(s, table_for_plot_sub$Species), 1:4]), aes(x = Pathway, y = Log_path_size)) +
      geom_bar(stat="identity", width = 0.3) + 
      geom_text(aes(label = Human_path_size, y = Log_path_size), size=3, hjust = -0.5, position = position_dodge(1)) +
      coord_flip() + 
      theme_classic() + labs(x=NULL, y="Log Scale Pathway Size (genes)") +
      theme(text = element_text(size=15))
    
    if(join_plots) {
      p_combine <- ggpubr::ggarrange(p, p2, ncol = 1, nrow= 2, heights = c(0.7, 0.3))
      print(p_combine)
      ggsave(file.path(plot_output_dir, paste0(data_name,"_", count_type, "_bubbles_bar_plt", n, ".png")),width = width, height = height, limitsize = FALSE)
      return(p_combine)
    } else {
      print(p)
      ggsave(file.path(plot_output_dir, paste0(data_name,"_", count_type, "_bubbles_plt", n,".png")), width = width, height = height*.7, limitsize = FALSE)
      print(p2)
      ggsave(file.path(plot_output_dir, paste0(data_name,"_", count_type, "_bar_plt", n,".png")), width = width, height = height *.3, limitsize = FALSE)
      return(list("bubble_plt" = p, "bar_plot" = p2))
      
    }
  }
}
#plot_percent_bubbles(count_percent = count_percent, count_type = "genes", n_genes_searched = n_genes_searched,data_name = data_name, join_plots = TRUE, height = 6)


#' Plot proteins and family bar chart
#'
#' @description Takes count percentages (produced by convert_to_percent function) and creates a plot with a barchart per pathway showing number of proteins, number or families and number of unmapped proteins per species per pathway.
#' Due to varying number of pathways outputted from a single G2P-SCAN analysis, number of pathways plotted in one plot is customisable - use 'n_pathways_per_plot' to set this size. If n_pathways_per_plot is lower than total number of pathways, multiple plots will be created and saved, iterating through all pathways.
#' @inheritParams plot_percent_bubbles 
#' @inheritParams run_all_plots
#' @import ggplot2
#' @importFrom forcats fct_relabel
#' @return list of ggplot objects produced by the function
#' @export
plot_protein_bars <- function(count_percent, data_name, plot_output_dir = "./", n_pathways_per_plot = 25, width = NULL, height = NULL) {
  dir.create(plot_output_dir, showWarnings = F, recursive = T)
  species <- unlist(unique(lapply(strsplit(colnames(count_percent), "_"), "[[", 1)))
  species <- species[!grepl("Pathway", species)]
  pathways <- unique(count_percent$Pathway.Name)
  table_for_plot <- data.frame()
  for(s in species) {
    for(p in pathways) {
      s_p_data <- count_percent[count_percent$Pathway.Name == p, grepl(s, colnames(count_percent))]
      rows <- data.frame("Pathway" = rep(p, 3),
                         "Species" = rep(s, 3),
                         "Status" = c("Number of Total proteins","Number of Proteins unassigned to a family", "Number of unique families assigned"),
                         "Values" = c(s_p_data[,grepl("Protein.Count", colnames(s_p_data)) ],
                                      s_p_data[,grepl("Proteins.Unmapped.to.families.Count", colnames(s_p_data)) ],
                                      s_p_data[,grepl("Family.Count", colnames(s_p_data)) ]))
      
      table_for_plot <- rbind(table_for_plot, rows)
    }
  }
  table_for_plot$Species <- factor(table_for_plot$Species, levels = c("Human",
                                                                      "Rat",
                                                                      "Mouse",
                                                                      "Zebrafish",
                                                                      "Nematode",
                                                                      "Fruitfly",
                                                                      "Yeast"))
  table_for_plot$Status <-factor(table_for_plot$Status, levels = c("Number of Total proteins","Number of unique families assigned","Number of Proteins unassigned to a family"))
  table_for_plot$Status <- forcats::fct_relabel(table_for_plot$Status, str_wrap,22)
  table_for_plot$Values <-as.numeric(table_for_plot$Values)
  

  np <- length(pathways)
  n_plots <- ceiling(np/n_pathways_per_plot) 
  xn <- 1
  list_of_plots <- list()
  for(n in 1:n_plots) {
    pathways_sub <- pathways[xn:(n_pathways_per_plot*n)]
    xn <- xn +n_pathways_per_plot
    table_for_plot_sub <- table_for_plot[table_for_plot$Pathway %in% pathways_sub,]
    new_np <- length(unique(table_for_plot_sub$Pathway))
    p <- ggplot(table_for_plot_sub, aes(x = Species, y = Values,fill = Status)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      ggtitle("Proteins & families count per pathway") +
      theme(strip.background = element_rect(colour="black",fill="white"),
            panel.background = element_rect("white"),
            axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0.5),
            legend.key.height=unit(1, "cm"),
            legend.position="top",
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =15),
            text = element_text(size = 15)) + 
      scale_fill_brewer(direction = -1) +
      labs(fill="Legend")
    if(n_pathways_per_plot <= 5) {
      w <- 18
      h <- 4
      split <- FALSE
    } else if(n_pathways_per_plot >5 & n_pathways_per_plot <= 10) {
      w <- 24
      h <- 8
      split <- TRUE
      row <- 2
      col <- 5
    } else if(n_pathways_per_plot >5 & n_pathways_per_plot <= 15) {
      w <- 35
      h <- 12
      split <- TRUE
      row <- 3
      col <- 5
    }else if(n_pathways_per_plot >15 & n_pathways_per_plot <= 20) {
      w <- 30
      h <- 16
      split <- TRUE
      row <- 4
      col <- 5
    }else if(n_pathways_per_plot >20 & n_pathways_per_plot <= 25) {
      w <- 35
      h <- 20
      split <- TRUE
      row <- 5
      col <- 5
    } else if(n_pathways_per_plot > 25) {
      col <- 6
      row <- ceiling(n_pathways_per_plot/6)
      w <- 40
      h <- 30
      split <- TRUE
    }
    if(is.null(height)){
      height = h
    } 
    if(is.null(width)) {
      width = w
    }
    
    if(split) {
      p <- p + facet_wrap(.~Pathway, nrow = row, ncol=col, scales = "free" )
    }
    print(p)
    ggsave(file.path(plot_output_dir, paste0(data_name, "_protein_bars_plt", n, ".png")),width = width, height = height)
    list_of_plots <- append(list_of_plots, p)
  }
  return(list_of_plots)
}


#' Plot entities and Reaction Box Plots
#' 
#' @description Takes count percentages (produced by convert_to_percent function) and creates a box plot of entity counts across pathways per species and of reaction counts across pathways per species.
#' @inheritParams plot_percent_bubbles
#' @inheritParams run_all_plots
#' @import ggplot2
#' @return list of entity and reaction ggplot objects 
#' @export
plot_entities_reaction_boxplots <- function(count_percent, data_name, plot_output_dir = "./", width = 14, height = 6) {
  dir.create(plot_output_dir, showWarnings = F, recursive = T)
  plot_data <- count_percent[, grepl("Pathway", colnames(count_percent)) |
                               grepl("Entity_Percent", colnames(count_percent)) |
                               grepl("Reaction_Percent", colnames(count_percent))]
  species <- unlist(unique(lapply(strsplit(colnames(plot_data), "_"), "[[", 1)))
  species <- species[!grepl("Pathway", species)]
  table_for_plot <- data.frame()
  e_r_data <- data.frame()
  for(s in species) {
    if(s != "Human") {
      s_data <- plot_data[, grepl(s, colnames(plot_data)) |grepl("Pathway", colnames(plot_data)) ]
      s_data$Species <- s
      colnames(s_data) <- c("Pathway.ID", "Pathway.Name", "Entity_Percent", "Reaction_Percent", "Species")
      
      e_r_data <- rbind(e_r_data, s_data)
    }
  }
  e_r_data$Species <- factor(e_r_data$Species, levels = c("Rat",
                                                          "Mouse",
                                                          "Zebrafish",
                                                          "Nematode",
                                                          "Fruitfly",
                                                          "Yeast"))
  
  e <- ggplot(e_r_data, aes(x=Species, y=Entity_Percent, fill=Species)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x = "", y = "Conservation (%)") +
    ggtitle("Entities count (% across pathways relative to human)") + 
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10),
          legend.position = "none")
  print(e)
  ggsave(file.path(plot_output_dir, paste0(data_name, "_entity_box_plt.png")), width = width, height = height)
  
  r <- ggplot(e_r_data, aes(x=Species, y=Reaction_Percent, fill=Species)) + 
    geom_boxplot() +
    theme_bw() +
    labs(x = "", y = "Conservation (%)") +
    ggtitle("Reaction count (% across pathways relative to human)") + 
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=10),
          legend.position = "none")
  print(r)
  ggsave(file.path(plot_output_dir, paste0(data_name, "_reaction_box_plt.png")), width = width, height = height)
  return(list("entities" = e, "reactions" = r))
}


#' Wrapper function to run all G2P-SCAN output plots
#' Given a directory, this function will read each file with suffix "_counts.xlsx" assuming a G2P-SCAN output. For each file the following plots are created:
#' \itemize{
#' \item{1. Bubble plots for gene percentage relative to human for each species and pathway}
#' \item{2. Protein and family bar charts}
#' \item{3. Entity box plots}
#' \item{4. Reaction box plots}
#'}
#'
#' @param g2p_output_dir the directory path for '_counts.xlsx' files to have plots created.
#' @param plot_output_dir output directory to write png files to. If output-dir does not exist, the directory will be created.
#' @param data_name output file prefix. If left NULL, input file name will be used. 
#' 
#' @importFrom openxlsx read.xlsx
#' @export
run_all_plots <- function(g2p_output_dir = "../", plot_output_dir = "./", data_name = NULL) {
  g2p_output_file <- list.files(path = g2p_output_dir, pattern = "_counts.xlsx", full.names = T)
  g2p_output_file_name <- list.files(path = g2p_output_dir, pattern = "_counts.xlsx", full.names = F)
  dir.create(plot_output_dir, showWarnings = F, recursive = T)
  for(i in 1:length(g2p_output_file)) {
    if(is.null(data_name)){
      data_name <- gsub("_counts.*","",g2p_output_file_name[i])
    }
    
    count_summary <- openxlsx::read.xlsx(xlsxFile = g2p_output_file[i], sheet = "count summary")
    input_summary <- openxlsx::read.xlsx(xlsxFile = g2p_output_file[i], sheet = "Input Summary")
    
    analysis_summary <- get_analysis_summary(input_summary)
    n_genes_searched <- analysis_summary["n_genes_searched"]
    count_percent <- convert_to_percent(count_summary)
    plot_percent_bubbles(count_percent = count_percent, count_type = "genes", plot_output_dir = plot_output_dir, n_genes_searched = n_genes_searched, data_name = data_name, join_plots = FALSE)
    plot_protein_bars(count_percent, data_name, plot_output_dir = plot_output_dir)
    plot_entities_reaction_boxplots(count_percent, data_name, plot_output_dir = plot_output_dir)
  }
}
