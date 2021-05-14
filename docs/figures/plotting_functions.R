#' Plot sample demographics
#'
#' @param sample_info Dataframe with columns: "sample_id", "Disease_Group",
#'   "Sex", "AoO", "AoD", "DD", "PMI", "aSN", "TAU", "thal AB", "aCG aSN score",
#'   "CSF_pH", "brain_weight_(g)", "RIN" , "ng_per_ul", "ng", "sent_to_bulk_seq"
#'
#' @return Plot of sample demogrpahics
#'

plot_sample_info <- function(sample_info){
  
  # Create named list with correct column names
  measures <- 
    c("Sex (count)", "Age of onset (years)", "Age of death (years)", "Disease duration (years)", "Post-mortem interval (hours)", "Braak alpha-syn staging (count)", 
      "Braak neuritic tau stage (count)", "Thal Abeta phase estimate (count)", "Alafuzoff alpha-syn staging in ACG (count)", "CSF pH", "Brain weight (g)", "RIN", "Sent to bulk-tissue RNA-seq (count)") %>% 
    str_wrap(., width = 25)
  names(measures) <- c("Sex", "AoO", "AoD", "DD", "PMI", "aSN", "TAU", "thal AB", "aCG aSN score", "CSF_pH", "brain_weight_(g)", "RIN", "sent_to_bulk_seq")
  
  # Split into numeric and categorical data
  numeric_data <- 
    sample_info %>% 
    dplyr::select(-Sex, -ng_per_ul, -ng, -sent_to_bulk_seq, -aSN, -TAU, -`thal AB`, -`aCG aSN score`) %>% 
    tidyr::gather(key = sample_characteristic, value = value, -sample_id, -Disease_Group)
  
  categorical_data <- 
    sample_info %>% 
    dplyr::select(sample_id, Disease_Group, aSN, TAU, `thal AB`, `aCG aSN score`, Sex, sent_to_bulk_seq) %>% 
    tidyr::gather(key = sample_characteristic, value = value, -sample_id, -Disease_Group)
  
  # Plot of numerical values
  numeric_plot <- 
    numeric_data %>% 
    dplyr::mutate(sample_characteristic = fct_relevel(sample_characteristic, 
                                                      unique(numeric_data$sample_characteristic))) %>% 
    ggplot(aes(x = Disease_Group,
               y = value,
               fill = Disease_Group)) + 
    geom_boxplot(alpha = 1, outlier.shape = NA, na.rm = T) + 
    # geom_point(binaxis = 'y', size = 0.8,   
    #            method = 'dotdensity', shape = 21, binpositions = 'all', binwidth = NULL, alpha = 0.6) +
    ggbeeswarm::geom_beeswarm(priority = "density", shape = 21, alpha = 0.6, size = 0.8) +
    facet_wrap(vars(sample_characteristic),
               scales = "free_y", 
               labeller = labeller(sample_characteristic = measures),
               ncol = 3) +
    expand_limits(y = 0) +
    labs(x = "", y = "") +
    scale_fill_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    scale_colour_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    theme_rhr +
    theme(legend.position = "none")
  
  # Plot of categorical values
  plot_list <- vector(mode = "list", length = length(unique(categorical_data$sample_characteristic)))
  
  for(i in 1:length(unique(categorical_data$sample_characteristic))){
    
    # If plot isn't on bottom, get rid of x-axis text, x-axis ticks and change margins so the plots appear to be overlapping
    if(i %in% 1:c(length(unique(categorical_data$sample_characteristic)) -3)){
      plot_list[[i]] <- 
        categorical_data %>% 
        dplyr::filter(sample_characteristic == unique(categorical_data$sample_characteristic)[i]) %>% 
        ggplot(aes(x = Disease_Group,
                   # y = ..prop..,
                   fill = fct_rev(value))) + 
        geom_bar() +
        facet_wrap(vars(sample_characteristic),
                   scales = "free_y", labeller = labeller(sample_characteristic = measures)) +
        labs(x = "", y = "", legend = "") +
        scale_fill_grey() +
        theme_rhr +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "right",
              legend.title = element_blank(),
              legend.key.size = unit(0.3, "cm"),
              plot.margin = margin(t = 0, r = 0.1, b= 0, l = 0.1, "cm"))
      
    } else{
      
      plot_list[[i]] <- 
        categorical_data %>% 
        dplyr::filter(sample_characteristic == unique(categorical_data$sample_characteristic)[i]) %>% 
        ggplot(aes(x = Disease_Group,
                   # y = ..prop..,
                   fill = fct_rev(value))) + 
        geom_bar() +
        facet_wrap(vars(sample_characteristic),
                   scales = "free_y", labeller = labeller(sample_characteristic = measures)) +
        labs(x = "", y = "", legend = "") +
        scale_fill_grey() +
        theme_rhr +
        theme(legend.position = "right",
              legend.title = element_blank(),
              legend.key.size = unit(0.3, "cm"),
              plot.margin = margin(t = 0, r = 0.1, b = 0.1, l = 0.1, "cm")) 
      
    }
    
  }
  
  categorical_plot <- 
    cowplot::plot_grid(plotlist = plot_list,
                       ncol = 3,
                       align = "hv")
  
  # Kruskal-wallis test for numeric data
  kw_test <- 
    numeric_data %>% 
    dplyr::mutate(sample_characteristic = fct_relevel(sample_characteristic, 
                                                      unique(numeric_data$sample_characteristic))) %>% 
    dplyr::group_by(sample_characteristic) %>% 
    dplyr::do(kruskal.test(x= .$value, g = .$Disease_Group) %>% 
                broom::tidy()) %>% 
    dplyr::inner_join(tibble(sample_characteristic = names(measures),
                             label = measures)) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(sample_characteristic = label %>% str_replace_all("\n", " ")) %>% 
    dplyr::select(sample_characteristic, p.value, method)
  
  # Chi-squared tests for categorical data
  # NAs omitted
  chisq_test <- 
    categorical_data %>% 
    dplyr::filter(sample_characteristic != "sent_to_bulk_seq" | is.na(value)) %>% 
    dplyr::group_by(sample_characteristic) %>% 
    summarise(p.value= chisq.test(Disease_Group, value)$p.value, 
              method = chisq.test(Disease_Group, value)$method) %>% 
    dplyr::inner_join(tibble(sample_characteristic = names(measures),
                             label = measures)) %>% 
    dplyr::mutate(sample_characteristic = label %>% str_replace_all("\n", " ")) %>% 
    dplyr::select(sample_characteristic, p.value, method)
  
  # Summary table of testing of sample demographics
  stable <- 
    kw_test %>% 
    dplyr::bind_rows(chisq_test) %>% 
    dplyr::mutate(sample_characteristic = str_wrap(sample_characteristic, width = 20),
                  p.value = round(p.value, digits = 5),
                  method = case_when(method == "Kruskal-Wallis rank sum test" ~ "Kruskal-Wallis test",
                                     method == "Pearson's Chi-squared test" ~ "Chi-squared test")) %>% 
    dplyr::arrange(desc(method), p.value) %>% 
    dplyr::rename(`Sample measure` = sample_characteristic,
                  `p-value` = p.value,
                  `Statistical test`= method)
  
  # Ggplot table
  stable.p <- 
    ggpubr::ggtexttable(stable, rows = NULL, theme = ttheme(base_style = "default", base_size = 6)) %>% 
    ggpubr::table_cell_font(row = c(which(stable$`p-value` < 0.05) + 1), column = 2, face = "bold", size = 6)
  
  
  # Arrange together
  ggarrange(ggarrange(numeric_plot,
                      stable.p,
                      ncol = 2,
                      labels = c("a", "c"),
                      widths = c(1.25,1)),
            categorical_plot,
            nrow = 2,
            labels = c("", "b"),
            heights = c(1.33, 1))
  
  
  
  
}

#' Load and process MultiQC statistics.
#'
#' @param filepath_to_multiqc File path to data folder for MultiQC report.
#'
#' @return Dataframe with processed MultiQC statistics.
#'

process_multiqc_metrics <- function(filepath_to_multiqc){
  
  # MultiQC reported metrics
  file_df <- 
    tibble(file_paths = list.files(filepath_to_multiqc, full.names = T, pattern = "txt"),
           file_name = list.files(filepath_to_multiqc, full.names = F, pattern = "txt") %>% 
             str_remove("multiqc_") %>% 
             str_remove(".txt"))
  
  for(i in 1:nrow(file_df)){
    
    file <- read_delim(file_df$file_paths[i] %>% as.character(), delim = "\t")
    
    if(i == 1){
      
      multiqc <- list(file)
      
    }else{
      
      multiqc <- c(multiqc, list(file))
      
    }
    
  }
  
  names(multiqc) <- file_df$file_name
  
  # Fastp stats
  fastp <- 
    multiqc$general_stats %>% 
    dplyr::filter(str_detect(Sample, "_fastp")) %>% 
    dplyr::mutate(sample_id = Sample %>% 
                    str_remove("_fastp")) %>% 
    dplyr::select(sample_id, pct_reads_passing_filter = `fastp_mqc-generalstats-fastp-pct_surviving`, gc_content = `fastp_mqc-generalstats-fastp-after_filtering_gc_content`) %>% 
    dplyr::mutate(gc_content = gc_content * 100)
  
  .scale <- function(x) (x/1e6)
  
  # STAR stats
  star <-
    multiqc$star %>% 
    dplyr::mutate(sample_id = Sample %>% 
                    str_remove("_.*")) %>% 
    dplyr::select(sample_id, contains("mapped"), -contains("percent")) %>% 
    dplyr::mutate_at(vars(uniquely_mapped:unmapped_tooshort), list(millions = .scale)) %>% 
    dplyr::select(sample_id, avg_mapped_read_length, contains("millions"))
  
  # Salmon stats
  salmon <-
    multiqc$salmon %>% 
    dplyr::mutate(sample_id = Sample %>% 
                    str_remove("_.*")) %>% 
    dplyr::mutate(num_mapped_millions = .scale(num_mapped),
                  num_processed_millions = .scale(num_processed)) %>% 
    dplyr::select(sample_id, num_mapped_millions, num_processed_millions, percent_mapped)
  
  
  # RSeQC read distribution
  rseqc <- 
    multiqc$rseqc_read_distribution %>% 
    dplyr::mutate(sample_id = Sample %>% 
                    str_remove("_.*")) %>% 
    dplyr::select(sample_id, contains("tag_pct"))
  
  multiqc_processed <- 
    setNames(list(fastp, star, salmon, rseqc),
             c("Fastp", "STAR", "Salmon", "RSeQC")) %>% 
    lapply(., function(df) df %>% tidyr::gather(key = "metric", value = "value", -sample_id)) %>% 
    qdapTools::list_df2df(col1 = "tool") %>% 
    dplyr::filter(sample_id != "Undetermined")
  
  return(multiqc_processed)
  
}



#' Plot multiqc metrics.
#'
#' @param multiqc_processed Dataframe. Output of \code{process_multiqc_metrics}.
#' @param sample_info Dataframe with columns: "sample_id", "Disease_Group",
#'   "Sex", "AoO", "AoD", "DD", "PMI", "aSN", "TAU", "thal AB", "aCG aSN score",
#'   "CSF_pH", "brain_weight_(g)", "RIN" , "ng_per_ul", "ng", "sent_to_bulk_seq".
#'
#' @return Plot of multiqc metrics.
#'


plot_multiqc_metrics <- function(multiqc_processed, sample_info){
  
  # Facet labels
  facet_labels <- c("Reads passing filter (%)", "GC content (%)", "Average mapped read length", "Uniquely mapped reads (millions)", "Multimapped reads (millions) - too many", "Unmapped mismatched reads (millions",
                    "Unmapped reads (millions) - other", "Multimapped reads (millions)", "Unmapped reads (millions) - too short", "Reads mapped (millions)", "Reads processed (millions)", "Reads mapped (%)",
                    "TSS_up_1kb", "3'UTR_exons", "TES_down_5kb", "TES_down_10kb", "CDS_exons", "TSS_up_10kb", "5'UTR_exons", "TSS_up_5kb", "Other_intergenic", "TES_down_1kb", "Introns")
  names(facet_labels) <- multiqc_processed$metric %>% unique()
  
  # Joint sample info
  multiqc_joint <- 
    multiqc_processed %>% 
    dplyr::inner_join(sample_info)
  
  boxplots <- 
    multiqc_joint %>% 
    dplyr::filter(
      tool == "Fastp" | 
        tool == "STAR" & metric == "avg_mapped_read_length" | 
        tool == "Salmon" & metric %in% c("num_mapped_millions", "percent_mapped")
    ) %>% 
    ggplot(
      aes(
        x = Disease_Group,
        y = value,
        fill = Disease_Group)
    ) + 
    geom_boxplot(alpha = 1, outlier.shape = NA, na.rm = T) + 
    # geom_point(binaxis = 'y', size = 0.8,   
    #            method = 'dotdensity', shape = 21, binpositions = 'all', binwidth = NULL, alpha = 0.6) +
    ggbeeswarm::geom_beeswarm(priority = "density", shape = 21, alpha = 0.6, size = 0.8) +
    facet_wrap(vars(tool, metric),
               scales = "free_y",
               labeller = labeller(metric = facet_labels),
               ncol = 1) +
    labs(x = "", y = "") +
    scale_fill_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    scale_colour_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    theme_rhr +
    theme(legend.position = "none")
  
  star_plot <- 
    multiqc_joint %>% 
    dplyr::filter(tool == "STAR", str_detect(metric, "millions")) %>% 
    dplyr::filter(!metric %in% c("multimapped_millions")) %>% 
    dplyr::mutate(metric = metric %>% 
                    str_remove("_millions") %>% 
                    str_replace_all("_", " ") %>% 
                    str_to_sentence() %>% 
                    str_wrap(width = 15)) %>% 
    ggplot(
      aes(
        x = metric, 
        y = value, 
        colour = Disease_Group, 
        fill = Disease_Group)
    ) +
    geom_boxplot(outlier.size = 0.8) +
    facet_wrap(vars(tool)) +
    labs(x = "", y = "Reads (millions)") +
    scale_fill_manual(values = alpha(as.character(colour_schemes_list$disease_greens$colour_hex), 0.6)) +
    scale_colour_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    theme_rhr +
    theme(legend.position = "none")
  
  rseqc_plot <- 
    multiqc_joint %>% 
    dplyr::filter(tool == "RSeQC") %>% 
    dplyr::mutate(metric = metric %>% 
                    str_remove("_tag_pct") %>% 
                    str_replace_all("_", " ") %>%
                    as.factor() %>% 
                    fct_relevel(.,
                                rev(c("cds exons", "5 utr exons", "3 utr exons", "introns", 
                                      "tss up 1kb", "tss up 5kb", "tss up 10kb", 
                                      "tes down 1kb", "tes down 5kb", "tes down 10kb", 
                                      "other intergenic")))) %>% 
    ggplot(aes(
      x = fct_relevel(sample_id,
                      multiqc_joint %>% 
                        dplyr::filter(tool == "RSeQC", metric == "cds_exons_tag_pct") %>% 
                        dplyr::arrange(Disease_Group, -value) %>% 
                        .[["sample_id"]]), 
      y = value, 
      fill = metric)
    ) +
    geom_col() +
    facet_grid(cols = vars(tool, Disease_Group), scales = "free_x", space = "free_x") +
    labs(x = "Sample ID", y = "% reads distributed over genomic features") +
    scale_fill_manual(values = RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"),
                      guide = guide_legend(nrow = 3)) +
    theme_rhr +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, units = "cm"))
  
  plot <- 
    ggarrange(boxplots,
              ggarrange(star_plot, rseqc_plot,
                        nrow = 2,
                        heights = c(1,1.5)),
              widths = c(1,2))
  
  return(plot)
  
}



#' Plot single-nucleus sequencing metrics.
#'
#' @param snrna_qc Dataframe of single-nucleus metrics outputted by Cell Ranger
#'   pipeline. Columns should include: "sample_id", "disease_group",
#'   "Total.Viable.Nuclei", "Fraction.Reads.in.Nuclei",
#'   "Mean.Reads.per.Nucleus", "Median.Genes.per.Nucleus",
#'   "Median.UMI.Counts.per.Nucleus", "Total.Genes.Detected",
#'   "Q30.Bases.in.Barcode", "Q30.Bases.in.RNA.Read",
#'   "Q30.Bases.in.Sample.Index", "Q30.Bases.in.UMI",
#'   "Reads.Mapped.Confidently.to.Genome",
#'   "Reads.Mapped.Confidently.to.Intergenic.Regions",
#'   "Reads.Mapped.Confidently.to.Exonic.and.Intronic.Regions",
#'   "Reads.Mapped.Confidently.to.Transcriptome".
#'
#' @return Plot of single-nucleus sequencing metrics grouped by type of metric
#'   they represent.
#'   

plot_snrna_metrics <- function(snrna_qc){
  
  # Tidy data for plotting
  data_to_plot <-
    snrna_qc %>% 
    dplyr::select(sample_id, 
                  disease_group, 
                  contains("Nuclei"),
                  contains("Nucleus"),
                  contains("Genes"), 
                  contains("Q30"), 
                  contains("Confidently"), 
                  contains("UMI"),
                  -Total.Genes.Detected2) %>%
    dplyr::rename_with(
      ~stringr::str_replace_all(.x, "\\.", " ") %>% 
        stringr::str_to_sentence() %>% 
        stringr::str_replace("umi", "UMI") %>% 
        stringr::str_replace("rna", "RNA"),
      -c("sample_id", "disease_group")
    ) %>% 
    tidyr::pivot_longer(
      cols = -c("sample_id", "disease_group"), 
      names_to = "metric",
      values_to = "value"
    ) %>% 
    dplyr::mutate(
      metric_group =
        case_when(
          str_detect(metric, "nucleus|nuclei|Total genes") ~ "Other QC",
          str_detect(metric, "Q30") ~ "Base call QC",
          str_detect(metric, "mapped") ~ "Alignment QC"
        ) %>% 
        fct_relevel(.,
                    c("Base call QC", "Alignment QC", "Other QC")),
      metric =
        case_when(str_detect(metric,
                             "Q30|Reads") ~ str_c(metric, " (%)"),
                  TRUE ~ metric) %>% 
        stringr::str_wrap(., width = 32),
      # For those metrics that are meant to be % transform
      value =
        case_when(str_detect(metric,
                             "Q30|Reads") ~ value * 100,
                  TRUE ~ value)
    ) 
  
  data_to_plot %>% 
    ggplot(
      aes(
        x = disease_group,
        y = value,
        fill = disease_group
      )
    ) +
    geom_boxplot(alpha = 1, outlier.shape = NA, na.rm = T) + 
    ggbeeswarm::geom_beeswarm(
      groupOnX = TRUE,
      priority = "density", 
      shape = 21, 
      alpha = 0.6, 
      size = 0.8) +
    facet_wrap(
      vars(metric_group, metric),
      ncol = 4,
      scales = "free_y") +
    labs(x = "", y = "") +
    scale_fill_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    scale_colour_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    theme_rhr +
    theme(legend.position = "none") +
    guides(x=guide_axis(angle = 45))
  
}



#' Plot snRNA-seq normalised gene expression of marker genes
#'
#' @param marker_gene_data Dataframe containing normalised expression of genes
#'   to plot, derived using \code{Seurat::FetchData} together with seurat object
#'   containing normalised expression across all nuclei. Additional columns that
#'   should be included are: "cell_type"
#' @param cols chr. Vector of columns that correspond to genes to plot.
#' @param marker_genes_df Dataframe containing marker genes and which cell type
#'   they are supposed to "mark". Must include columns: "gene", "marker".
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("Class") and wrapped labels ("Class_wrapped").
#' @param palette chr. Character vector of colour palette to be used for cell
#'   types.
#'
#' @return Violin plot of normalised marker gene expression across cell types.

plot_snrna_marker_genes <- function(marker_gene_data, cols, marker_genes_df, ct_class, palette){
  
  marker_genes_df <- 
    marker_genes_df %>% 
    dplyr::mutate(
      marker = 
        case_when(
          marker == "Astrocytes" ~ "Astro",
          marker == "Oligodendrocytes" ~ "Oligo",
          TRUE ~ marker
        )
    ) %>% 
    dplyr::inner_join(
      ct_class, 
      by = c("marker" = "CellType")
    ) %>% 
    dplyr::select(
      gene, marker_of = class
    )
      
  data_to_plot <-
    marker_gene_data %>% 
    tidyr::pivot_longer(
      cols = {{ cols }},
      names_to = "gene",
      values_to = "norm_expr"
    ) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(
      cell_type = 
        case_when(
          cell_type == "Astrocytes" ~ "Astro",
          cell_type == "Oligodendrocytes" ~ "Oligo",
          TRUE ~ cell_type
        ),
      # Add noise to replicate Seurat VlnPlot (https://github.com/satijalab/seurat/issues/3322)
      n = n(),
      noise = rnorm(n) / 100000,
      norm_expr_w_noise = norm_expr + noise
    ) %>% 
    dplyr::inner_join(
      ct_class,
      by = c("cell_type" = "CellType")
    ) %>% 
    dplyr::inner_join(
      marker_genes_df,
      by = c("gene")
    ) 
  
  data_to_plot %>% 
    ggplot(
      aes(
        x = class,
        y = norm_expr_w_noise,
        fill = class
      )
    ) +
    geom_violin(scale = "width", trim = T, lwd = 0.25) +
    ggh4x::facet_nested(marker_of + gene ~ ., 
                        scales = "free_x") +
    labs(
      x = "",
      y = "Expression level",
      fill = ""
    ) +
    scale_fill_manual(
      values = palette
    ) +
    theme_rhr +
    theme(legend.position = "none")  +
    guides(x=guide_axis(angle = 45))
  
}

#' Plot expression of XIST and DDX3Y by sample.
#'
#' @param dds DESeq2 object. Unfiltered DESeq2 object.
#' @param sample_info Dataframe with sample id and sex.
#'
#' @return Plot of XIST and DDX3Y gene expression coloured by reported sex.

plot_sex_genes <- function(dds, sample_info){
  
  sex_gene_expr <- 
    dds %>% 
    assay() %>% # Access count data
    as_tibble(rownames = "gene_id") %>%
    dplyr::mutate(gene_id = gene_id %>% 
                    str_replace("\\..*", "")) %>% 
    dplyr::filter(gene_id %in% c("ENSG00000229807", "ENSG00000067048")) %>% 
    tidyr::gather(key = sample_id, value = counts, -gene_id) %>% 
    dplyr::mutate(gene_name = case_when(gene_id == "ENSG00000229807" ~ "XIST",
                                        gene_id == "ENSG00000067048" ~ "DDX3Y")) %>% 
    inner_join(sample_info) 
  
  ggplot(sex_gene_expr, aes(x = sample_id, y = as.numeric(counts), fill = Sex)) +
    geom_col(colour = "black") + 
    facet_grid(gene_name~., scales = "free_y") +
    labs(x = "Sample ID", y = "Counts") +
    scale_fill_manual(limits = c("F", "M"), values = rev(gray.colors(2)), drop = FALSE) +
    theme_rhr
  
}


#' Plot significance of correlations between covariates and PC axes.
#'
#' @param correlations Output of \code{correlatePCs}.
#' @param prcomp_obj A \code{prcomp} object.
#' @param total_pc int. Total number of PC axes.
#' @param num_pc_to_plot int. Number of PC axes to plot.
#' @param fdr_threshold int. FDR threshold for significance. This will determine
#'   which tiles have text relating to the value of FDR. Default = 0.05.
#' @param sample_info Dataframe of sample information including only
#'   demographics. This is used to categorise covariates into either "sample
#'   demographics" or "cell-type proportions".
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("Class") and wrapped labels ("Class_wrapped").
#'
#' @return Plot of significantly correlated covariates and PC axes.

plot_covar_PC_corr <- function(correlations, prcomp_obj, total_pc, num_pc_to_plot = NULL, fdr_threshold = 0.05, sample_info, ct_class){
  
  # Calculate variation explained by each axis
  pc_var <- 
    tibble(PC_order = c(1:length(prcomp_obj$sdev)),
           sd = prcomp_obj$sdev,
           var = prcomp_obj$sdev ^ 2,
           percent_var = var/sum(var)*100,
           PC_label = str_c("PC", PC_order, " (", round(percent_var, digits = 0), "%)")) %>% 
    dplyr::select(PC_order, PC_label)
  
  
  # FDR correct correlations
  pc_corr_fdr <- 
    correlations$pvalue %>%
    as_tibble() %>%
    dplyr::mutate(PC = str_c("PC ", c(1:total_pc)),
                  PC_order = c(1:total_pc)) %>%
    dplyr::select(-sample_id) %>%
    tidyr::gather(key = covariate, value = p, -PC, -PC_order) %>%
    dplyr::mutate(FDR = p.adjust(p, method = "fdr"),
                  covariate_type = case_when(covariate %in% colnames(sample_info) ~ "Sample demographics",
                                             str_detect(covariate, "_proportions") ~ "Cell-type proportions") %>% 
                    fct_rev(),
                  covariate = covariate %>%
                    str_remove("_proportions") %>% 
                    str_replace_all("_", " "),
                  fdr_text = case_when(FDR < fdr_threshold ~ scientific(FDR, digits = 2))) %>%
    dplyr::left_join(ct_class, by = c("covariate" = "CellType")) %>%
    dplyr::inner_join(pc_var) %>% 
    dplyr::mutate(covariate = case_when(!is.na(as.character(class)) ~ as.character(class),
                                        TRUE ~ as.character(covariate)))
  
  if(!is.null(num_pc_to_plot)){
    
    pc_corr_fdr <- 
      pc_corr_fdr %>%
      dplyr::filter(PC_order %in% c(1:num_pc_to_plot))
    
  }
  
  
  
  pc_corr_fdr %>%
    ggplot(aes(x = fct_relevel(PC_label, 
                               pc_corr_fdr %>% 
                                 dplyr::arrange(PC_order) %>% 
                                 .[["PC_label"]] %>% 
                                 unique()), 
               y = fct_relevel(covariate,
                               pc_corr_fdr %>% 
                                 .[["covariate"]] %>% 
                                 unique() %>% 
                                 rev()), 
               fill = -log10(FDR))) +
    geom_tile() +
    # coord_equal() +
    geom_text(aes(label = fdr_text), size = 2) +
    labs(x = "", y = "") +
    facet_grid(rows = vars(covariate_type), scales = "free_y", space = "free") +
    scale_fill_viridis_c(option = "cividis") +
    theme_rhr +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0))
  
}

#' Plot genes detected pre- and post-filtering.
#'
#' @param genes_detected Dataframe with columns: "sample_id",
#'   "n_genes_detected", "n_genes_detected_after_filtering".
#' @param sample_info Dataframe with columns: "sample_id", "Disease_Group",
#'   "Sex", "AoO", "AoD", "DD", "PMI", "aSN", "TAU", "thal AB", "aCG aSN score",
#'   "CSF_pH", "brain_weight_(g)", "RIN" , "ng_per_ul", "ng",
#'   "sent_to_bulk_seq".
#'
#' @return Plot of genes detected pre- and post-filtering.
#' 

plot_genes_detected <- function(genes_detected, sample_info){
  
  genes_detected <- 
    genes_detected %>% 
    dplyr::inner_join(sample_info, by = "sample_id")
  
  supp_fig <- vector(mode = "list", length = 2)
  
  facet_labels <- c("Genes detected", "Genes detected post-filtering") %>% str_wrap(width = 15)
  names(facet_labels) <- c("n_genes_detected", "n_genes_detected_after_filtering")
  
  supp_fig[[1]] <-
    genes_detected %>% 
    dplyr::select(sample_id, contains("genes"), Disease_Group) %>% 
    tidyr::gather(key = "gene_category", value = "n", -sample_id, -Disease_Group) %>% 
    ggplot(aes(x = Disease_Group, y = n, fill = Disease_Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 1) +
    # geom_point(binaxis = 'y', size = 0.8,   
    #            method = 'dotdensity', shape = 21, binpositions = 'all', binwidth = NULL, alpha = 0.6) +
    ggbeeswarm::geom_beeswarm(priority = "density", shape = 21, alpha = 0.6, size = 0.8) +
    facet_wrap(vars(gene_category), labeller = labeller(gene_category = facet_labels)) +
    labs(x = "", y = "Number of genes") +
    scale_fill_manual(values = as.character(colour_schemes_list$disease_greens$colour_hex)) +
    theme_rhr +
    theme(legend.position = "none")
  
  supp_fig[[2]] <- 
    genes_detected %>% 
    dplyr::group_by(Disease_Group) %>% 
    dplyr::summarise(`Median genes detected` = median(n_genes_detected),
                     `Median genes detected post-filtering` = median(n_genes_detected_after_filtering)) %>% 
    dplyr::rename(`Disease group` = Disease_Group) %>% 
    ggpubr::ggtexttable(rows = NULL, theme = ttheme(base_style = "default", base_size = 6))
  
  return(supp_fig)
  
}

#' Plot heatmap-style plot with number of DEG per cell and comparison.
#'
#' @param DEG_df_filtered df. Dataframe of differentially expressed genes. Must
#'   contain columns: "comparison", "cell_type", "gene".
#' @param fct_levels chr. Vector of order in which comparisons should be
#'   factored. Bear in mind any "_" in current comparison names will be replaced
#'   with " ", so this should also be reflected in the fct_levels vector.
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("class") and wrapped labels ("class_wrapped").
#'
#' @return Heatmap-style plot of number of DEG per cell type and comparison.

plot_snrna_summary_table <- function(DEG_df_filtered, fct_levels, ct_class){
  
  # summary table of n genes
  stable <-
    DEG_df_filtered %>% 
    dplyr::distinct(comparison, cell_type, gene) %>% 
    dplyr::group_by(comparison, cell_type) %>% 
    dplyr::summarise(n_genes = n(),
                     n_genes_text = n() %>% 
                       scales::comma()) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(comparison = comparison %>% 
                    str_replace_all("_", " ") %>% 
                    fct_relevel(.,
                                fct_levels),
                  n_genes_text_col = case_when(n_genes > 2650 ~ "high",
                                               TRUE ~ "low")) %>% 
    dplyr::inner_join(ct_class, by = c("cell_type" = "CellType"))
  
  # Plot
  stable.p <- 
    stable %>% 
    ggplot(aes(x = comparison, y = fct_rev(class), fill = n_genes)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = n_genes_text, colour = n_genes_text_col), size = 1.5) +
    coord_equal() +
    labs(x = "", y = "", title = "Number of DEGs") +
    scale_fill_distiller(palette = "Greys", direction = 1) +
    scale_colour_manual(values = c("white", "black")) +
    guides(colour = F, fill = F) +
    theme_minimal(base_size = 6) + 
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          plot.title = element_text(size = 6, hjust = 0.5))
  
  return(stable.p)
  
}

#' Generate dummy df with unique DEG gene ids from each comparison combined with
#' and cell types.
#'
#' This function will first extract all unique differentially expressed gene ids
#' from each comparison across all cell types and then generate a dummy
#' dataframe, wherein the unique set of DEG ids from a particular comparison is
#' combined with all cell types. This dummy dataframe can be joined (either left
#' of right depending on order of dataframes) with the dataframe of DEG results,
#' which will ensure that for those genes that are not detected as expressed in
#' a particular cell type, these will now appear as a row in the dataframe of
#' DEG results (with values of NA).
#'
#' @param DEG_df df. Dataframe of DEG results. Must contain the columns
#'   "comparison","gene", "coef", "fdr", "direction_of_effect". In the column
#'   "direction_of_effect", directions should be coded as "up" or "down".
#' @param fc_threshold num. Fold change threshold used to call a gene as
#'   differentially expressed. **Remember to apply log2 to the fold-change
#'   threshold desired as MAST outputs log2(FC)**. E.g. if you want DEGs that
#'   have an expression 1.5 times higher/lower between two groups, the
#'   equivalent threshold would be log2(1.5).
#' @param fdr_threshold num. FDR threshold used to call a gene as differentially
#'   expressed.
#' @param DEG_direction chr. Direction of effect to filter results for.
#'
#' @return Dataframe with gene ids of differentially expressed genes from each
#'   comparison combined with all cell types.

generate_DEG_dummy_df <- function(DEG_df, fc_threshold, fdr_threshold, DEG_direction = c("up", "down")){
  
  ids <- 
    DEG_df %>% 
    dplyr::filter(direction_of_effect == DEG_direction) %>% 
    dplyr::filter(fdr < fdr_threshold, abs(coef) > fc_threshold) %>%
    dplyr::distinct(comparison, gene)
  
  cell_types <- unique(DEG_df$class)
  
  ids_list <- vector(mode = "list", length = length(cell_types)) 
  
  for(i in 1:length(cell_types)){
    
    ids_list[[i]] <- 
      ids %>% 
      dplyr::mutate(class = cell_types[i])
    
  }
  
  ids <- 
    ids_list %>% 
    qdapTools::list_df2df(col1 = "list_name") %>% 
    dplyr::select(-list_name)
  
  return(ids)
  
}

#' Plot DEG genes as a yes/no overlapping bar plot.
#'
#' Function will plot DEG genes by comparison and cell type. Within a
#' comparison, DEG genes will be plotted as a binary 1/0 indicating whether they
#' are DEG or not within a cell type.
#'
#' @param data_to_plot df. Dataframe of DEGs. Must contain columns:
#'   "comparison", "gene", "DEG" (binary 1/0 of whether gene is DEG or not),
#'   and "class" (factored cell types).
#' @param ylab chr. Character vector for labelling y-axis.
#' @param palette chr. Character vector of colour palette to be used for cell
#'   types.
#'
#' @return Plot of DEG overlaps between cell types within a comparison.

plot_snrna_gene_overlaps <- function(data_to_plot, ylab, palette){
  
  theme_rhr_no_axis_text <-  
    theme_bw(base_family = "Helvetica") + 
    theme(panel.grid.major.x = element_blank(),
          legend.position = "top",
          strip.text = element_text(size = 5),
          strip.text.y = element_text(size = 5, angle = 0),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_text(vjust = 0.6),
          axis.title = element_text(size = 5),
          panel.spacing = unit(0.1, "lines"))
  
  plot_components <- setNames(vector(mode = "list", length = length(unique(sort(data_to_plot$comparison)))),
                              unique(sort(data_to_plot$comparison)))
  
  for(i in 1:length(plot_components)){
    
    plot <- 
      data_to_plot %>%   
      dplyr::filter(comparison == names(plot_components)[i])
    
    genes_ordered <- 
      plot %>% 
      dplyr::arrange(-DEG, class) %>% 
      .[["gene"]] %>% 
      unique()
    
    n_genes <- str_c("Genes (n = ",length(unique(plot$gene)), ")") 
    
    plot_components[[i]] <- setNames(list(plot, genes_ordered, n_genes),
                                     c("data", "x_order", "x_lab"))
    
  }
  
  .plot_function <- function(data, gene_order, xlab, ylab){
    
    data %>% 
      ggplot(aes(x = fct_relevel(gene, gene_order), 
                 y = DEG,
                 fill = class)) +
      geom_col(width = 1) +
      facet_grid(cols = vars(comparison), rows = vars(class), scales = "free_x") +
      labs(x = xlab, y = ylab) +
      scale_fill_manual(values = palette) +
      theme_rhr_no_axis_text + 
      theme(legend.position = "none") 
    
  }
  
  cowplot::plot_grid(.plot_function(data = plot_components[[1]]$data, gene_order = plot_components[[1]]$x_order, xlab = plot_components[[1]]$x_lab, ylab = ylab) +
                       theme(strip.text.y = element_blank()),
                     .plot_function(data = plot_components[[2]]$data, gene_order = plot_components[[2]]$x_order, xlab = plot_components[[2]]$x_lab, ylab = ylab) +
                       theme(strip.text.y = element_blank(),
                             axis.title.y = element_blank()),
                     .plot_function(data = plot_components[[3]]$data, gene_order = plot_components[[3]]$x_order, xlab = plot_components[[3]]$x_lab, ylab = ylab) +
                       theme(axis.title.y = element_blank()), rel_widths = c(1.05,1,1.55),
                     align = "h", axis = "bt",
                     nrow = 1) 
  
}

#' Generate a dummy dataframe with all comparisons, cell types, pathways and DEG
#' directions of effect.
#'
#' This function will generate a dummy dataframe with all combinations of
#' comparison, cell types, pathways, and DEG directions of effect. This can be
#' joined (either left of right depending on order of dataframes) with real
#' results to ensure that if a pathway was not enriched in some combination of
#' comparison/cell type/direction of effect a value of NA will be assigned to
#' the FDR column. This is important for plotting with
#' \code{ggplot2::geom_tile()}, which will not plot an outline of a tile if
#' there is no value within the FDR column (which is used as the
#' \code{ggplot2::geom_tile()} fill aesthetic) of that particular
#' comparison/cell type/direction of effect combination.
#'
#' @param df df. Dataframe of pathway enrichment results.
#' @param pathway_col chr. Vector with the name of the column denoting pathway
#'   names.
#' @param cell_type_col chr. Vector with the name of the column denoting cell
#'   types.
#' @param direction_col chr. Vector with the name of the column denoting DEG
#'   direction of effect.
#' @param comparison_col chr. Vector with the name of the column denoting
#'   comparisons.
#'
#' @return Dummy dataframe with all combinations of pathway/comparison/cell
#'   type/direction of effect.

generate_dummy_pathway_df <- function(df, pathway_col, cell_type_col, direction_col, comparison_col){
  
  pathways <- 
    df %>% 
    .[[pathway_col]] %>% 
    unique()
  
  cell_types <- 
    df %>% 
    .[[cell_type_col]] %>% 
    unique()
  
  direction <- 
    df %>% 
    .[[direction_col]] %>% 
    unique()
  
  comparisons <- 
    df %>% 
    .[[comparison_col]] %>% 
    unique()
  
  dummy_df <- 
    tibble(pathway = rep(pathways, c(length(comparisons) * length(cell_types) * length(direction)))) %>% 
    dplyr::arrange(pathway) %>% 
    dplyr::mutate(comparison = rep(comparisons, c(length(pathways) * length(cell_types) * length(direction)))) %>% 
    dplyr::arrange(comparison, pathway) %>% 
    dplyr::mutate(cell_type = rep(cell_types, c(length(pathways) * length(comparisons) * length(direction)))) %>% 
    dplyr::arrange(comparison, cell_type, pathway) %>%
    dplyr::mutate(direction_of_effect = rep(direction, c(length(cell_types) * length(comparisons) * length(pathways))))
  
  return(dummy_df)
  
}

#' Plot GSEA results from snRNA.
#'
#' Function will generate a heatmap style plot of GSEA enrichment.
#'
#' @param pathway_GO_sim_df df. Dataframe of GSEA pathway results that have been
#'   processed through \code{rutils::go_reduce()}.
#' @param gsea_dummy_df df. Dataframe of GSEA pathway results processed through
#'   \code{generate_dummy_pathway_df()}.
#' @param fct_levels chr. Vector of order in which comparisons should be
#'   factored. Bear in mind any "_" in current comparison names will be replaced
#'   with " ", so this should also be reflected in the fct_levels vector.
#'
#' @return Heatmap-style plot of GSEA pathway results. On y-axis, GSEA terms are
#'   represented by their reduced parent term. The number inside the brackets
#'   indicates the number of child GO terms assigned to each parent term. Parent
#'   terms are represented in a comparison/cell-type by the minimum FDR assigned
#'   to a child term that is associated with the parent term (within that
#'   comparison/cell-type). Parent terms will be ordered by the times they
#'   appear across comparsions and cell types.

plot_snrna_gsea_pathways <- function(pathway_GO_sim_df, gsea_dummy_df, fct_levels){
  
  child_parent_df <- 
    pathway_GO_sim_df %>% 
    dplyr::distinct(description, parent_term) %>% 
    dplyr::group_by(parent_term) %>% 
    dplyr::summarise(n_child_terms = n()) %>% 
    dplyr::mutate(parent_labels = str_c(parent_term, " (", n_child_terms, ")"))
  
  data_to_plot <- 
    pathway_GO_sim_df %>% 
    dplyr::group_by(comparison, class, direction_of_effect, parent_id, parent_term) %>% 
    dplyr::summarise(n_enriched_pathways_per_celltype_comparison = n(), min_FDR = min(FDR)) %>% 
    dplyr::ungroup() %>% 
    dplyr::right_join(gsea_dummy_df) %>% 
    dplyr::inner_join(child_parent_df) %>% 
    dplyr::mutate(parent_term = str_to_sentence(parent_term),
                  min_FDR = case_when(min_FDR == 0 ~ 2.22e-16,
                                      TRUE ~ min_FDR),
                  direction_of_effect = case_when(direction_of_effect == "down" ~ "Down-regulated DEG",
                                                  direction_of_effect == "up" ~ "Up-regulated DEG"),
                  comparison = comparison %>% 
                    str_replace_all("_", " ") %>% 
                    fct_relevel(.,
                                fct_levels))
  
  
  
  # Re-factor pathways to appear in order of most frequent first
  pathway_fct <- 
    data_to_plot %>% 
    dplyr::filter(!is.na(min_FDR)) %>% 
    dplyr::group_by(parent_labels) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::arrange(n, desc(parent_labels)) %>% 
    .[["parent_labels"]] %>% 
    unique() %>% 
    as.character()
  
  plot <- data_to_plot %>% 
    ggplot(aes(x = class, 
               y = fct_relevel(parent_labels,
                               pathway_fct), 
               fill = -log10(min_FDR))) +
    geom_tile(colour = "black") +
    # coord_equal() +
    ggh4x::facet_nested(. ~ direction_of_effect + comparison) + 
    labs(x = "", y = "Reduced GO Term") +
    scale_fill_viridis_c(name = "-log"[1][0]~"(FDR)",
                         breaks = c(0,5,10,15),
                         na.value = "white", 
                         option = "cividis") +
    theme_rhr +
    theme(legend.key.size = unit(0.35, "cm"))
  
  return(plot)
  
}

#' Plot DEG fold change
#'
#' Function will plot fold change of DEGs as heatmap, ordered such that DEGs
#' that are DE across the most conditions are at the top of plot.
#'
#' @param DEG_df df. Dataframe of DEGs to be plotted.
#' @param fct_levels chr. Vector of order in which comparisons should be
#'   factored. Bear in mind any "_" in current comparison names will be replaced
#'   with " ", so this should also be reflected in the fct_levels vector.
#'
#' @return Heatmap-style plot of DEG fold change.

plot_deg_fc <- function(DEG_df, fct_levels){
  
  dummy_df <- 
    tibble(gene = rep(unique(DEG_df$gene), length(unique(DEG_df$comparison)) * length(unique(DEG_df$class)))) %>% 
    dplyr::arrange(gene) %>% 
    dplyr::mutate(class = rep(unique(DEG_df$class), length(unique(DEG_df$comparison)) * length(unique(DEG_df$gene)))) %>% 
    dplyr::arrange(class, gene) %>% 
    dplyr::mutate(comparison = rep(unique(DEG_df$comparison), length(unique(DEG_df$class)) * length(unique(DEG_df$gene)))) 
  
  gene_fct <- 
    DEG_df %>% 
    dplyr::group_by(gene) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::arrange(n) %>% 
    .[["gene"]]
  
  # get max absolute fold change to set colour scale
  limit <- DEG_df %>% 
    dplyr::top_n(n = 1, wt = abs(coef)) %>% 
    .[["coef"]] %>%
    abs() %>% 
    round()
  
  plot <- 
    DEG_df %>% 
    dplyr::right_join(dummy_df) %>% 
    dplyr::mutate(comparison = comparison %>% 
                    str_replace_all(., "_", " ") %>% 
                    fct_relevel(.,
                                fct_levels)) %>% 
    ggplot(aes(x = class, 
               y = fct_relevel(gene,
                               gene_fct), 
               fill = coef)) +
    geom_tile(colour = "black") +
    # coord_equal() +
    facet_grid(cols = vars(comparison), scales = "free", space = "free",
               labeller = labeller(class = setNames(ct_class$class_wrapped, ct_class$class))) +
    scale_fill_distiller(palette = "RdBu", direction = -1, na.value = "#cccccc", limits = c(-limit, limit)) +
    labs(x = "", y = "Gene", fill = "log"[2]~"(fold change)") +
    theme_rhr +
    theme(strip.text.y = element_text(angle = 0),
          legend.key.size = unit(0.35, "cm"))
  
  return(plot)
  
}

#' Plot cumulative frequency of normalised single-nucleus gene expression.
#'
#' Function will plot the cumulative distribution function of single-nucleus
#' gene expression across cell types. Plots are facetted by disease group and
#' gene.
#'
#' @param expr_data Dataframe containing normalised expression of genes
#'   to plot, derived using \code{Seurat::FetchData} together with seurat object
#'   containing normalised expression across all nuclei. Additional columns that
#'   should be included are: "cell_type", "disease_group"
#' @param cols chr. Vector of columns that correspond to genes to plot.
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("Class") and wrapped labels ("Class_wrapped").
#' @param palette chr. Character vector of colour palette to be used for cell
#'   types.
#'
#' @return Cumulative frequency plot of normalised gene expression across cell
#'   types in each disease group.

plot_snrna_cdf_expr <- function(expr_data, cols, ct_class, palette){
  
  data_to_plot <-
    expr_data %>% 
    tidyr::pivot_longer(
      cols = {{ cols }},
      names_to = "gene",
      values_to = "norm_expr"
    ) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(
      cell_type = 
        case_when(
          cell_type == "Astrocytes" ~ "Astro",
          cell_type == "Oligodendrocytes" ~ "Oligo",
          TRUE ~ cell_type
        ),
      # Add noise to replicate Seurat VlnPlot (https://github.com/satijalab/seurat/issues/3322)
      n = n(),
      noise = rnorm(n) / 100000,
      norm_expr_w_noise = norm_expr + noise
    )  %>% 
    dplyr::inner_join(
      ct_class,
      by = c("cell_type" = "CellType")
    )
  
  data_to_plot %>%
    ggplot(
      aes(
        x = norm_expr_w_noise,
        colour = class
      )
    ) +
    stat_ecdf(geom = "step") + 
    geom_hline(yintercept = 0.5, linetype = "dotted") +
    labs(
      x = "Expression level", 
      y = "Cumulative frequency",
      colour = ""
    ) +
    facet_grid(cols = vars(gene),
               rows = vars(disease_group)) +
    expand_limits(x = 0) +
    scale_colour_manual(values = palette) +
    theme_rhr +
    theme(
      axis.text.x = element_text(angle = 0),
      legend.position = "right"
    )
  
}

#' Draws a ridge plot of normalised single-nucleus gene expression.
#'
#' Function will draw a ridge plot function of single-nucleus gene expression
#' across each disease group within a cell type. Ridge plots include quantiles
#' displaying 0-50% of the data, 50-90% and 90-100%.
#'
#' @param expr_data Dataframe containing normalised expression of genes to plot,
#'   derived using \code{Seurat::FetchData} together with seurat object
#'   containing normalised expression across all nuclei. Additional columns that
#'   should be included are: "cell_type", "disease_group"
#' @param cols chr. Vector of columns that correspond to genes to plot.
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("Class") and wrapped labels ("Class_wrapped").
#'
#' @return Ridgeplot of normalised single-nucleus gene expression, split by
#'   quantiles.

plot_snrna_ridge_expr <- function(expr_data, cols, ct_class){
  
  data_to_plot <-
    expr_data %>% 
    tidyr::pivot_longer(
      cols = {{ cols }},
      names_to = "gene",
      values_to = "norm_expr"
    ) %>% 
    dplyr::group_by(gene) %>% 
    dplyr::mutate(
      cell_type = 
        case_when(
          cell_type == "Astrocytes" ~ "Astro",
          cell_type == "Oligodendrocytes" ~ "Oligo",
          TRUE ~ cell_type
        ),
      disease_group =
        fct_relevel(
          disease_group,
          c("Control", "PD", "PDD", "DLB")
        ),
      # Add noise to replicate Seurat VlnPlot/Ridgeplot (https://github.com/satijalab/seurat/issues/3322)
      n = n(),
      noise = rnorm(n) / 100000,
      norm_expr_w_noise = norm_expr + noise
    )  %>% 
    dplyr::inner_join(
      ct_class,
      by = c("cell_type" = "CellType")
    )
  
  data_to_plot %>% 
    ggplot(
      aes(
        x = norm_expr_w_noise,
        y = disease_group,
        fill = factor(stat(quantile))
      )
    ) +
    ggridges::stat_density_ridges(
      geom = "density_ridges_gradient",
      quantile_lines = TRUE,
      quantiles = c(0.5, 0.90), 
      calc_ecdf = TRUE,
      alpha = 0.5
    ) +
    labs(x = "Expression level", y = "Disease group") +
    facet_wrap(vars(class, gene)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    scale_fill_manual(
      name = "Cumulative frequency", 
      values = c("#74ADD1B2", "#A0A0A0A0", "#F46D43B2"),
      labels = c("0-50%", "50-90%", "90-100%")
    ) +
    theme_rhr +
    theme(
      axis.text.x = element_text(size = 5, angle = 0),
      legend.position = "top",
      legend.key.size = unit(0.35, "cm")
    )
  
}

#' Plot number of enriched pathways by comparison and cell type.
#'
#' @param pd_pathways_df df. Dataframe of pathway enrichment results using
#'   PD-specific pathways. Must contain columns: "class" (factored cell types),
#'   "comparison", "direction_of_effect" (DEG direction of effect) and "pathway_enriched"
#'   (binary 1/0 denoting whether pathway is enriched).
#'
#' @return Bar plot showing number of pathways enriched per cell type and
#'   comparison. Fill denotes DEG direction of effect.


plot_n_snrna_pd_pathways <- function(pd_pathways_df){
  
  plot <- 
    pd_pathways_df %>% 
    dplyr::group_by(class, comparison, direction_of_effect) %>% 
    dplyr::summarise(n_enriched_pathways = n()) %>% 
    ggplot(aes(x = comparison, y = n_enriched_pathways, fill = direction_of_effect)) +
    geom_col(colour = "black", size = 0.25) +
    facet_wrap(vars(class), ncol = 4) +
    labs(x = "",y = "Number of enriched pathways") +
    scale_fill_manual(name = "",
                      values = alpha(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")[c(9,3)], alpha = 0.7)) +
    theme_rhr +
    theme(legend.key.size = unit(0.35, "cm"))
  
  return(plot)
  
}



#' Plot log10(FDR) of PD-specific pathway enrichments.
#'
#' @param pd_pathways_df df. Dataframe of pathway enrichment results using
#'   PD-specific pathways. Must contain columns: "class" (factored cell types),
#'   "comparison", "direction_of_effect" (DEG direction of effect), "pathway"
#'   (pathway name) and "pathway_group" (group to which pathway belongs).
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("class") and wrapped labels ("class_wrapped").
#' @param pathway_group logical. Whether plot should include a facet for pathway
#'   groups. Default is TRUE.
#'
#' @return Heatmap-style plot of GSEA pathway results using only PD-specific
#'   pathways.

plot_snrna_pd_pathways <- function(pd_pathways_df, ct_class, pathway_group = TRUE){
  
  data_to_plot <- 
    pd_pathways_df %>% 
    dplyr::mutate(FDR = case_when(FDR == 0 ~ 2.22e-16,
                                  TRUE ~ FDR),
                  log_fdr = -log10(FDR),
                  pathway = str_wrap(pathway, width = 58))
  
  # Re-factor pathways to appear in order of most frequent first within each pathway group
  pathway_fct <-
    data_to_plot %>% 
    dplyr::filter(!is.na(log_fdr)) %>% 
    dplyr::group_by(pathway_group, pathway) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::arrange(pathway_group, n, desc(pathway)) %>% 
    .[["pathway"]] %>% 
    unique() %>% 
    as.character()
  
  
  if(pathway_group == TRUE){
    
    plot <- 
      data_to_plot %>% 
      ggplot(aes(x = comparison, 
                 y = pathway %>% fct_relevel(pathway_fct), 
                 fill = log_fdr)) +
      geom_tile(colour = "black") +
      # coord_equal() +
      ggh4x::facet_nested(pathway_group ~ direction_of_effect + class, 
                          scales = "free", 
                          space = "free",
                          labeller = labeller(pathway_group = setNames(unique(pd_pathways_df$pathway_group_wrapped), unique(pd_pathways_df$pathway_group)),
                                              class = setNames(ct_class$class_wrapped, ct_class$class))) +
      scale_fill_viridis_c(na.value = "white", 
                           option = "cividis") +
      labs(x = "", y = "Pathway", fill = "-log"[1][0]~"(FDR)") +
      theme_rhr +
      theme(strip.text.x = element_text(angle = 90),
            strip.text.y = element_text(angle = 0),
            legend.key.size = unit(0.35, "cm"))
    
  } else {
    
    plot <- 
      data_to_plot %>% 
      ggplot(aes(x = comparison, 
                 y = pathway %>% fct_relevel(pathway_fct), 
                 fill = log_fdr)) +
      geom_tile(colour = "black") +
      # coord_equal() +
      ggh4x::facet_nested(. ~ direction_of_effect + class, 
                          scales = "free", 
                          space = "free",
                          labeller = labeller(pathway_group = setNames(unique(pd_pathways_df$pathway_group_wrapped), unique(pd_pathways_df$pathway_group)),
                                              class = setNames(ct_class$class_wrapped, ct_class$class))) +
      scale_fill_viridis_c(na.value = "white", 
                           option = "cividis") +
      labs(x = "", y = "Pathway", fill = "-log"[1][0]~"(FDR)") +
      theme_rhr +
      theme(strip.text.x = element_text(angle = 90),
            strip.text.y = element_text(angle = 0),
            legend.key.size = unit(0.35, "cm"))
  }
  
  return(plot)
  
} 



#' Plot cell-type proportions coloured by fill variable.
#'
#' Either plot of cell-type proportions by fill variable OR if fill_var == NULL,
#' cell-type proportions will be plotted in order of median cell-type proportion
#' across all groups.
#'
#' @param ct_proportions Dataframe with columns: "sample_id", "Celltype",
#'   "Celltype_class", "Disease_Group", "cell_type_proportion", "class",
#'   "class_wrapped"
#' @param fill_var chr. Variable to fill by.
#' @param fill_palette chr. Colours to use for fill. Only necessary if fill_var
#'   = T.
#' @param fill_var_title chr. Legend title for fill variable.
#' @param include_x_axes logical. Set to TRUE if x axis ticks and labels to be
#'   included. Default = FALSE.
#' @param add_points logical. Add points to boxplot. Default = NULL.
#' @param y_lab Character. Y-axis label.
#' @param y_lim Numeric. Y-axis limits.
#'
#' @return Plot cell-type proportions coloured by fill variable.
#'   

plot_cell_prop <- function(ct_proportions, fill_var = NULL, fill_palette = NULL, fill_var_title = "", add_points = NULL, include_x_axes = FALSE, y_lab, y_lim){
  
  plot <- 
    ct_proportions %>% 
    ggplot(aes(x = fct_reorder(class,
                               .x = cell_type_proportion,
                               .fun = median,
                               .desc = T), 
               y = cell_type_proportion)) + 
    geom_boxplot(fill = "grey", outlier.alpha = 0) +
    coord_cartesian(ylim = y_lim) +
    labs(x = NULL, y = y_lab) + 
    theme_rhr + 
    theme(legend.key.size = unit(0.5, "cm"),
          legend.position = "right")
  
  if(!is.null(add_points)){
    
    plot <- 
      plot +
      ggbeeswarm::geom_beeswarm(fill = "grey", priority = "density", shape = 21, alpha = 0.5, size = 0.5)
    
  }
  
  if(!is.null(fill_var)){
    
    if(is.null(fill_palette)) stop("If fill_var = TRUE, must supply a vector of hex codes to fill_palette.")
    
    plot <- 
      ct_proportions %>% 
      ggplot(aes(x = class, 
                 y = cell_type_proportion, 
                 fill = .data[[fill_var]])) + 
      geom_boxplot(outlier.size = 0.5, outlier.alpha = 0) +
      scale_fill_manual(values = fill_palette) +
      coord_cartesian(ylim = y_lim) +
      labs(x = NULL, y = y_lab, fill = fill_var_title) + 
      theme_rhr + 
      theme(legend.key.size = unit(0.5, "cm"),
            legend.position = "right")
    
    if(!is.null(add_points)){
      
      plot <- 
        plot +
        ggbeeswarm::geom_beeswarm(priority = "density", shape = 21, alpha = 0.5, dodge.width=0.75, size = 0.5)
      
    }
    
    if(include_x_axes == FALSE){
      
      plot <- 
        plot +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = margin(t = 0.5, r= 0.5, b = 0, l = 0.5, "cm"))
      
    }
    
  }
  
  return(plot)
  
}

#' Plot cell-type proportions derived from snRNA-seq data.
#'
#' @param snrna_proportions Dataframe with columns: "sample_id", "Celltype",
#'   "cell_type_proportion.
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("class") and wrapped labels ("class_wrapped").
#' @param sample_info Dataframe with columns: "sample_id", "Disease_Group",
#'   "Sex", "AoO", "AoD", "DD", "PMI", "aSN", "TAU", "thal AB", "aCG aSN score",
#'   "CSF_pH", "brain_weight_(g)", "RIN" , "ng_per_ul", "ng",
#'   "sent_to_bulk_seq".
#' @param palette Chracter. Vector with hex codes for colour palette.
#'
#' @return Plot of snRNA-seq cell-type proportions grouped by disease and by individual.
#' 

plot_snrna_cell_prop <- function(snrna_proportions, ct_class, sample_info, palette){
  
  data_to_plot <- 
    snrna_proportions %>% 
    dplyr::inner_join(sample_info) %>% 
    dplyr::inner_join(ct_class, by = c("Celltype" = "CellType")) 
  
  sample_plot <- 
    data_to_plot %>% 
    ggplot(aes(x = sample_id, 
               y = cell_type_proportion, 
               fill = class)) +
    geom_col() +
    facet_grid(cols = vars(Disease_Group), space = "free_x", scales = "free_x") +
    labs(x = "", y = "Cell-type proportion") +
    scale_fill_manual(values = palette) +
    theme_rhr +
    theme(legend.position = "none")
  
  grouped_plot <-
    data_to_plot %>% 
    dplyr::group_by(Disease_Group) %>% 
    dplyr::mutate(total = sum(cell_type_proportion),
                  normalised_cell_type_proportion = cell_type_proportion/total) %>% 
    ggplot(aes(x = Disease_Group, 
               y = normalised_cell_type_proportion, 
               fill = class)) +
    geom_col() +
    labs(x = "", y = "Cell-type proportion") +
    scale_fill_manual(values = palette) +
    theme_rhr +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, "cm"))
  
  plots <- list(grouped_plot, sample_plot)
  
}

#' Plot results of Wilcoxon test on cell-type proportions.
#'
#' @param ct_proportions_test Dataframe with factored cell type ("class"),
#'   wrapped cell type labels ("class_wrapped"), name of pairwise comparison
#'   tested ("comparison"), p-value of test ("p.value") and FDR-corrected
#'   p-value ("fdr").
#'
#' @return Heatmap style plot of FDR values.

test_cell_prop <- function(ct_proportions_test){
  
  ct_proportions_test %>%
    dplyr::mutate(log_fdr = case_when(fdr <= 0.1 ~ -log10(fdr)),
                  fdr_text = case_when(fdr < 0.05 ~ "**",
                                       fdr <= 0.1 ~ "*")) %>%
    ggplot(aes(x = class, y = fct_rev(comparison), fill = log_fdr)) +
    geom_tile(colour = "black") + 
    # coord_equal() +
    geom_text(aes(label = fdr_text)) +
    labs(x = "Cell type", y = "", fill = "-log"[1][0]~"(FDR)") +
    scale_fill_gradient(low = "#dddddd", high = "#777777", na.value = "white") +
    theme_rhr + 
    theme(legend.key.size = unit(0.5, "cm"),
          legend.position = "right", 
          plot.margin = margin(t = 0, r = 0.5, b = 0.5, l = 0.5, "cm"))
  
}

#' Plot scatter plot of snRNA-seq and deconvolution proportions.
#'
#' @param gathered_results Dataframe with columns: "sample_id", "Celltype",
#'   "cell_type_proportion" and "source".
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("class") and wrapped labels ("class_wrapped").
#' @param ncol int. Number of columns in facet. Default = 3.
#'
#' @return Scatterplot of snRNA-seq and deconvolution proportions.

plot_cell_prop_scatterplot <- function(gathered_results, ct_class, ncol = 3){
  
  gathered_results %>% 
    dplyr::select(sample_id, Celltype, cell_type_proportion, source) %>% 
    tidyr::spread(key = source, value = cell_type_proportion)  %>% 
    dplyr::inner_join(ct_class, by = c("Celltype" = "CellType")) %>%  
    dplyr::filter(Celltype != "Unclassified") %>% 
    ggplot(aes(x = `Scaden deconvolution`, y = `snRNA-seq`)) +
    geom_point(alpha = 0.5, size = 0.5) +
    ggpubr::stat_cor(method = "spearman", 
                     cor.coef.name = "rho", 
                     label.x.npc = "left", 
                     label.y.npc = "top", 
                     size = 2) +
    facet_wrap(~ class, 
               ncol = ncol,
               # scales = "free",
               labeller = labeller(class = setNames(ct_class$class_wrapped, ct_class$class))) +
    labs(x = "Cell-type proportion (%) - Scaden deconvolution", y = "Cell-type proportion (%) - single nucleus labelling") +
    scale_x_continuous(limits = c(0,100)) +
    scale_y_continuous(limits = c(0,100)) +
    coord_fixed(ratio = 1) +
    theme_rhr +
    theme(axis.text.x = element_text(angle = 0))
  
}

#' Extract DEGs from pairwise comparisons in a DESeqDataSet.
#'
#' Function to extract DEGs from multiple pairwise comparisons in a DESeqDataSet. Uses \code{DESeq2::results()} and \code{DESeq2::lfcShrink()}. Provides both original fold change estimates as well as those shrunk using  using the normal prior distribution in DESeq2.
#'
#' @param dds DESeqDataSet.
#' @param sample_info Sample info as extract from dds using \code{DESeq2::colData()}.
#' @param group_column_name chr. Name of column that groups samples.
#'
#' @return Dataframe with: gene; baseMean (mean of normalized counts for all samples); log2FoldChange; lfcSE (log2 fold change standard error); pvalue; stat (test statistic); padj (FDR-adjusted p-value); shrunk_log2FC (log2FoldChange shrunk with normal prior distribution); shrunk_lfcSE; comparison (where group on left-hand-side of "vs" is reference).

extract_results_dds <- function(dds, sample_info, group_column_name){
  
  # Create vector of groups
  groups <- sample_info %>%
    .[[group_column_name]] %>%
    unique() %>% 
    sort() %>% 
    as.character()
  
  # Create dataframe of comparisons
  comparisons <-
    combn(x = groups,
          m = 2) %>%
    t() %>%
    as_tibble()
  
  # Loop through comparisons and for each extract results
  for(i in 1:nrow(comparisons)){
    
    # Filter for comparisons
    comparison <- comparisons[i,] %>%
      purrr::as_vector()
    
    # Contrast in this order will result in logFC(comparison[2]/comparison[1])
    contrast <- c(group_column_name, comparison[2], comparison[1])
    
    print(contrast)
    
    # Unshrunken results
    results <- DESeq2::results(dds, contrast = contrast, pAdjustMethod = "BH") 
    
    # Shrink estimates
    shrunk <- DESeq2::lfcShrink(dds, contrast = contrast, res = results, type = "normal") 
    
    # Combine
    results <- results %>% 
      as_tibble(rownames = "gene") %>%
      dplyr::inner_join(shrunk %>% 
                          as_tibble(rownames = "gene") %>% 
                          dplyr::select(gene, shrunk_log2FC = log2FoldChange, shrunk_lfcSE = lfcSE))
    
    if(i == 1){
      
      all_results <- results %>% 
        dplyr::mutate(comparison = str_c(comparison[1], "_vs_", comparison[2]))
      
    } else {
      
      all_results <- all_results %>% 
        dplyr::bind_rows(results %>% 
                           dplyr::mutate(comparison = str_c(comparison[1], "_vs_", comparison[2])))
      
      
    }
    
  }
  
  return(all_results)
  
}

#' Plot pathway results from ClusterProfiler.
#'
#' @param pathway_results df. Dataframe derived from clusterProfiler object
#'   using \code{fortify}.
#' @param facet_var var. Variable(s) to facet by. Provide in the format vars(x,
#'   y). Default = NULL.
#' @param cividis_palette logical. Use "cividis" in viridis palette. Default =
#'   FALSE.
#'
#' @return Plot where size of dot indicates gene ratio, which is the proportion
#'   of genes in the input list that are annotated to the pathway. Fill of dot
#'   indicates FDR-adjusted p-value.
#'   

plot_pathways <- function(pathway_results, facet_var = NULL, cividis_palette = FALSE){
  
  cluster_fct <-
    pathway_results %>% 
    dplyr::distinct(GO_type, comparison, Cluster) %>% 
    dplyr::arrange(GO_type, comparison) %>% 
    .[["Cluster"]] %>% 
    as.character() %>% 
    unique()
  
  pathway_fct <- 
    pathway_results %>% 
    dplyr::distinct(GO_type, comparison, Cluster, Description) %>% 
    dplyr::arrange(GO_type, comparison, desc(Description)) %>% 
    .[["Description"]] %>% 
    as.character() %>% 
    unique()
  
  p <- 
    pathway_results %>% 
    ggplot(aes(x = Cluster %>% 
                 fct_relevel(., cluster_fct), 
               y = Description %>% 
                 fct_relevel(., pathway_fct), 
               size = GeneRatio, colour = -log10(p.adjust))) +
    geom_point() +
    geom_point(shape = 1,colour = "black") +
    facet_wrap(vars(GO_type), scales = "free") +
    scale_size_continuous(name = "Gene Ratio") +
    scale_colour_viridis_c(name = "-log"[1][0]~"(FDR)") +
    labs(x = "", y = "Pathway") +
    theme_rhr +
    theme(legend.position = "right",
          legend.key.size = unit(0.35, "cm"))
  
  
  if(cividis_palette == TRUE){
    p <- 
      p +
      scale_colour_viridis_c(name = "-log"[1][0]~"(FDR)", option = "cividis")
    
  }
  
  return(p)
}

#' Plot reduced pathway results from ClusterProfiler run with multiple input
#' gene lists (e.g. differential splicing results from multiple pairwise
#' comparisons)
#'
#' ClusterProfiler results should have been reduced with
#' \code{rutils::go_reduce()}.
#'
#' @param pathway_results df. Dataframe of ClusterProfiler pathway results that
#'   have been processed through \code{rutils::go_reduce()}.
#' @param fct_levels chr. Vector of order in which input gene lists should be
#'   factored. Bear in mind any "_" in current input gene list names will be
#'   replaced with " ", so this should also be reflected in the fct_levels
#'   vector.
#'
#' @return Dot plot. On y-axis, GO terms are represented by their reduced parent
#'   term. The number inside the brackets indicates the number of child GO terms
#'   assigned to each parent term. On x-axis, the input gene lists run, with
#'   size of input list in brackets. Size of dot indicates the number of child
#'   terms annotated to the parent term within a particular pairwise input gene
#'   list. Fill of dot indicates the minimum FDR assigned to a child term that
#'   is associated with the parent term (within that input gene list). Parent
#'   terms are ordered first by the factoring of the input gene list and then by
#'   the number of associated child terms.
#'   

plot_cprofiler_reduced <- function(pathway_results, fct_levels){
  
  child_parent_df <-
    pathway_results %>%
    dplyr::distinct(Description, parent_term) %>%
    dplyr::group_by(parent_term) %>%
    dplyr::summarise(n_child_terms = n()) %>%
    dplyr::mutate(parent_labels = str_c(parent_term, " (", n_child_terms, ")"))
  
  data_to_plot <-
    pathway_results %>%
    dplyr::group_by(comparison, parent_id, parent_term) %>%
    dplyr::summarise(n_enriched_pathways_per_comparison = n(), 
                     min_FDR = min(p.adjust),
                     # Number of genes in input list should be same for each comparison, thus mean() to summarise
                     n_genes = mean(n_genes)) %>%
    dplyr::ungroup() %>%
    # dplyr::right_join(gsea_dummy_df) %>%
    dplyr::inner_join(child_parent_df) %>%
    dplyr::mutate(prop_enriched_terms_per_comparison = n_enriched_pathways_per_comparison/n_child_terms,
                  parent_term = str_to_sentence(parent_term),
                  min_FDR = case_when(min_FDR == 0 ~ 2.22e-16,
                                      TRUE ~ min_FDR),
                  comparison = comparison %>%
                    str_replace_all("_", " ") %>% 
                    fct_relevel(.,
                                fct_levels),
                  comparison_labels = str_c(comparison, "\n(", n_genes, ")"))
  
  
  
  # Re-factor pathways to appear in order of comparisons on x-axis and thereafter n child terms in parent term
  pathway_fct <-
    data_to_plot %>% 
    dplyr::distinct(comparison, parent_labels, n_child_terms) %>% 
    dplyr::arrange(comparison, -n_child_terms, parent_labels) %>% 
    .[["parent_labels"]] %>% 
    unique() %>% 
    as.character()
  
  # Re-factor comparison labels
  comparison_fct <- 
    data_to_plot %>% 
    dplyr::distinct(comparison, comparison_labels) %>% 
    dplyr::arrange(comparison) %>% 
    .[["comparison_labels"]] %>% 
    unique() %>% 
    as.character()
  
  plot <- data_to_plot %>%
    ggplot(aes(x = fct_relevel(comparison_labels,
                               comparison_fct),
               y = fct_relevel(parent_labels,
                               pathway_fct),
               colour = -log10(min_FDR),
               size = prop_enriched_terms_per_comparison)) +
    geom_point() +
    geom_point(shape = 1,colour = "black") +
    labs(x = "", y = "Reduced GO Term") +
    scale_size_continuous(name = "Prop. enriched terms", breaks = c(0.2, 0.4, 0.6, 0.8, 1)) +
    scale_colour_viridis_c(name = "-log"[1][0]~"(FDR)",
                           # breaks = c(0,5,10,15),
                           na.value = "white",
                           option = "cividis") +
    theme_rhr +
    theme(legend.key.size = unit(0.35, "cm"))
  
  return(plot)
  
}

#' Plot reduced pathway results from ClusterProfiler run with one input gene
#' list.
#'
#' ClusterProfiler results should have been reduced with
#' \code{rutils::go_reduce()}.
#'
#' @param pathway_results df. Dataframe of ClusterProfiler pathway results that
#'   have been processed through \code{rutils::go_reduce()}.
#'
#' @return Bar plot. On y-axis, GO terms are represented by their reduced parent
#'   term. The number inside the brackets indicates the number of child GO terms
#'   assigned to each parent term. On x-axis, the -log10(FDR) is shown. For each
#'   parent term, the minimum FDR assigned to a child term is used to represent
#'   the parent term.. Parent terms are ordered by significance.
#'   


plot_cprofiler_reduced_barplot <- function(pathway_results){
  
  child_parent_df <-
    pathway_results %>%
    dplyr::distinct(Description, parent_term) %>%
    dplyr::group_by(parent_term) %>%
    dplyr::summarise(n_child_terms = n()) %>%
    dplyr::mutate(parent_labels = str_c(parent_term, " (", n_child_terms, ")"))
  
  data_to_plot <-
    pathway_results %>%
    dplyr::group_by(parent_id, parent_term) %>%
    dplyr::summarise(n_enriched_pathways_per_comparison = n(), 
                     min_FDR = min(p.adjust)) %>%
    dplyr::ungroup() %>%
    # dplyr::right_join(gsea_dummy_df) %>%
    dplyr::inner_join(child_parent_df) %>%
    dplyr::mutate(parent_term = str_to_sentence(parent_term),
                  min_FDR = case_when(min_FDR == 0 ~ 2.22e-16,
                                      TRUE ~ min_FDR))
  
  
  
  # Re-factor pathways to appear by significance
  pathway_fct <-
    data_to_plot %>% 
    dplyr::arrange(-min_FDR) %>% 
    .[["parent_labels"]] %>% 
    unique() %>% 
    as.character()
  
  plot <-  
    data_to_plot %>%
    ggplot(aes(x = -log10(min_FDR),
               y = fct_relevel(parent_labels,
                               pathway_fct))) +
    geom_col() +
    labs(x = "-log"[1][0]~"(FDR)", y = "Reduced GO Term") +
    theme_rhr +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  
  return(plot)
  
}

#' Plot EWCE results.
#'
#' @param ewce_results Output of \code{MarkerGenes::run_ewce_controlled}, with
#'   added columns: FDR.p, FDR-correction of p-value; comparison_type, whether
#'   comparison is with reference to the control or a disease group.
#' @param ct_class Dataframe with cell type ("CellType"), factored cell types
#'   ("class") and wrapped labels ("class_wrapped").
#' @param rotate logical. If true, plot is rotated. Default = FALSE.
#' @param cividis_palette logical. Use "cividis" in viridis palette. Default =
#'   FALSE.
#'
#' @return Geom_tile plot of EWCE results.

plot_ewce <- function(ewce_results, ct_class, rotate = FALSE, cividis_palette = FALSE){
  
  data_to_plot <- 
    ewce_results %>% 
    dplyr::inner_join(ct_class)
    
  p <- 
    data_to_plot %>% 
    ggplot(aes(x = GeneSet, y = forcats::fct_rev(Study))) +
    geom_tile(aes(fill = sd_from_mean), colour = "black") + 
    # coord_equal() +
    facet_grid(cols = vars(comparison_type), 
               rows = vars(class),
               scales = "free", 
               space = "free_y", 
               labeller = labeller(class = setNames(ct_class$class_wrapped, ct_class$class))) +
    scale_fill_viridis_c(name = "s.d. from\nthe mean", na.value = "white") +
    labs(x = "", y = "Disease status of cell type") +
    theme_rhr + 
    theme(panel.grid = element_blank(),
          strip.text.y = element_text(angle = 0),
          legend.key.size = unit(0.35, "cm"),
          legend.position = "right")
  
  if(rotate == TRUE){
    
    p <- 
      data_to_plot %>% 
      ggplot(aes(x = Study, y = forcats::fct_rev(GeneSet))) +
      geom_tile(aes(fill = sd_from_mean), colour = "black") + 
      # coord_equal() +
      facet_grid(rows = vars(comparison_type), 
                 cols = vars(class),
                 scales = "free", 
                 space = "free_x", 
                 labeller = labeller(class = setNames(ct_class$class_wrapped, ct_class$class))) +
      scale_fill_viridis_c(name = "s.d. from\nthe mean", na.value = "white") +
      labs(x = "Disease status of cell type", y = "") +
      theme_rhr + 
      theme(panel.grid = element_blank(),
            strip.text.y = element_text(angle = 0),
            legend.key.size = unit(0.35, "cm"),
            legend.position = "top")
    
    
  }
  
  if(cividis_palette == TRUE){
    p <- 
      p +
      scale_fill_viridis_c(name = "s.d. from\nthe mean", na.value = "white", option = "cividis")
    
  }
  
  return(p)
  
}


#' Plot proportions of each junction category.
#'
#' @param intron_list Output of \code{dasper::annotate_junc_ref} applied to
#'   leafcutter intron usage data.
#' @param legend_position chr. Where legend should be positioned.
#'
#' @return Plot of junction category proportions across disease comparisons.

plot_intron_annot_prop <- function(intron_annotated, legend_position = c("top", "bottom", "left", "right"), novel_combo = FALSE){
  
  
  count <- 
    intron_annotated %>% 
    dplyr::filter(!junc_cat %in% c("ambig_gene", "novel_combo")) %>% 
    dplyr::group_by(comparison, junc_cat) %>% 
    dplyr::summarise(n = n()) %>% 
    dplyr::mutate(prop = n/sum(n) * 100) %>% 
    dplyr::ungroup()
  
  if(novel_combo == TRUE){
    
    count <- 
      intron_annotated %>% 
      dplyr::filter(!junc_cat %in% c("ambig_gene")) %>% 
      dplyr::group_by(comparison, junc_cat) %>% 
      dplyr::summarise(n = n()) %>% 
      dplyr::mutate(prop = n/sum(n) * 100) %>% 
      dplyr::ungroup()
    
  }
  
  total <- 
    count %>% 
    dplyr::group_by(comparison) %>% 
    dplyr::summarise(total = sum(n))
  
  labels <- 
    count %>% 
    dplyr::inner_join(total) %>% 
    dplyr::distinct(comparison, total) %>% 
    dplyr::mutate(comparison_fct = fct_relevel(comparison, 
                                               c("Control_vs_PD", "Control_vs_PDD", "Control_vs_DLB", "PD_vs_PDD", "PD_vs_DLB", "PDD_vs_DLB")),
                  comparison_label = str_replace_all(comparison, "_", " ") %>% 
                    str_c(., "\n(", total, ")")) %>% 
    dplyr::arrange(comparison_fct)
  
  count_to_plot <- 
    count %>% 
    dplyr::mutate(junc_cat = junc_cat %>% 
                    str_replace_all("_", " ") %>% 
                    str_to_sentence() %>% 
                    fct_relevel(levels = c("Annotated", "Novel exon skip", "Novel donor", "Novel acceptor", "None")),
                  comparison_type = case_when(str_detect(comparison, "Control") ~ "Ref: control",
                                              TRUE ~ "Ref: disease")) %>% 
    dplyr::inner_join(labels)
  
  fill_colours <- c("#3C5488", "#4DBBD5", "#E64B35", "#00A087", "grey")
  
  if(novel_combo == TRUE){
    
    count_to_plot <- 
      count %>% 
      dplyr::mutate(junc_cat = case_when(junc_cat == "novel_combo" ~ "novel_combination",
                                         TRUE ~ junc_cat), 
                    junc_cat = junc_cat %>% 
                      str_replace_all("_", " ") %>% 
                      str_to_sentence() %>% 
                      fct_relevel(levels = c("Annotated", "Novel exon skip", "Novel combination", "Novel donor", "Novel acceptor", "None")),
                    comparison_type = case_when(str_detect(comparison, "Control") ~ "Ref: control",
                                                TRUE ~ "Ref: disease")) %>% 
      dplyr::inner_join(labels)
    
    fill_colours <- c("#3C5488", "#4DBBD5", "#F3E361", "#E64B35", "#00A087", "grey")
    
  }
  
  plot <- 
    count_to_plot %>% 
    ggplot(aes(x = fct_relevel(comparison_label,
                               count_to_plot %>% 
                                 arrange(comparison_fct) %>% 
                                 .[["comparison_label"]] %>% 
                                 unique()), 
               y = prop, 
               fill = fct_rev(junc_cat)), colour = "black") +
    geom_col() +
    facet_grid(~ comparison_type, scales = "free_x") +
    labs(x = "", y = "Proportion of introns (%)") +
    scale_fill_manual(name = "Splicing event",
                      values = rev(fill_colours)) +
    theme_rhr +
    theme(
      legend.position = legend_position,
      legend.key.size = unit(0.6, "cm"))
  
  return(plot)
  
}



#' Plot pathway enrichments from single-nucleus DEG, bulk-tissue DS and bulk-tissue PC.
#' 
#' This function will plot pathways that were (i) enriched in one of the bulk-tissue analyses (DS or PC) and (ii) enriched across at least two of the analyses (bulk DS/PC or snRNA DEG).
#'
#' @param pathway_summary_list Named list ("snrna_deg", "bulk_pc", "bulk_deg") of dataframes, containing the columns: "comparison", "class", "go_type", "go_id" and "description". Ensure that descriptions are standardised across dataframes.
#' @param comparison_pattern chr. Vector of patterns to use in \code{str_detect()} in order to filter which comparisons to display.
#' @param output_dir chr. File path to output directory where png will be saved.
#' @param output_name chr. Name given to saved png.
#' @param column_title_rotation int. The rotation that should be used for column labels. Default is 0.
#' @param figure_width int. Width in mm that figure should be saved as. Default is 150.
#' @param figure_height int. Height in mm that figure should be saved as. Default is 100.
#'
#' @return
#'

plot_pathway_summary_heatmap <- 
  function(
    pathway_summary_list,
    comparison_pattern,
    output_dir,
    output_name,
    column_title_rotation = 0,
    figure_width = 150,
    figure_height = 100
  ){
    
    # Find terms that are found in one of bulk sets and overlap at least one other set
    # Visualisation will be restricted to these terms
    pathway_summary_count <-
      pathway_summary_list %>%
      qdapTools::list_df2df(col1 = "method") %>%
      dplyr::group_by(go_id, description) %>%
      dplyr::summarise(
        n = n()
      ) %>%
      dplyr::mutate(
        present_in_bulk_ds =
          case_when(go_id %in% c(pathway_summary_list$bulk_ds %>%
                                   dplyr::filter(str_detect(comparison, comparison_pattern)) %>%
                                   .[["go_id"]] %>%
                                   unique()) ~ TRUE,
                    TRUE ~ FALSE),
        present_in_bulk_pc =
          case_when(go_id %in% pathway_summary_list$bulk_pc$go_id ~ TRUE,
                    TRUE ~ FALSE),
        present_in_snrna =
          case_when(go_id %in% c(pathway_summary_list$snrna %>%
                                   dplyr::filter(str_detect(comparison, comparison_pattern)) %>%
                                   .[["go_id"]] %>%
                                   unique()) ~ TRUE,
                    TRUE ~ FALSE),
        sum_present = sum(present_in_bulk_ds, present_in_bulk_pc, present_in_snrna)
      ) %>%
      dplyr::filter(
        # sum_present >= 1,
        sum_present > 1,
        present_in_bulk_ds == TRUE | present_in_bulk_pc == TRUE
      )
    
    # Filter pathway summary for comparison of interest and go ids
    data_to_plot <- 
      pathway_summary_list %>% 
      qdapTools::list_df2df(col1 = "method") %>% 
      dplyr::mutate(
        method = 
          case_when(method == "bulk_ds" ~ "Bulk DS",
                    method == "bulk_pc" ~ "Bulk PC",
                    method == "snrna_deg" ~ "snRNA DEG")
      ) %>% 
      tidyr::unite(
        col = "gene_list", 
        c("method", "comparison", "class"),
        remove = F,
        sep = ":"
      ) %>% 
      dplyr::filter(str_detect(comparison, comparison_pattern),
                    go_id %in% pathway_summary_count$go_id) %>% 
      dplyr::arrange(method, gene_list)
    
    
    # Create matrix of 1 for present and 0 for absent
    mat_df <-
      data_to_plot %>%
      dplyr::mutate(
        present =
          case_when(go_id %in% pathway_summary_count$go_id ~ 1,
                    TRUE ~ 0)
      ) %>%
      dplyr::select(gene_list, description, present) %>%
      tidyr::pivot_wider(names_from = gene_list,
                         values_from = present,
                         values_fill = 0) %>%
      # Use row totals to reorder go_ids
      dplyr::mutate(
        sum = rowSums(.[-1])
      ) %>%
      dplyr::arrange(-sum) %>%
      dplyr::select(-sum)
    
    mat <- as.matrix(mat_df[,-1])
    row.names(mat) <- mat_df$description
    
    # Set graphic parameters of heatmaps globally
    ht_opt(
      heatmap_column_names_gp = gpar(fontsize = 5), 
      heatmap_column_title_gp = gpar(fontsize = 5, fontface = "bold"),
      heatmap_row_names_gp = gpar(fontsize = 5), 
      heatmap_row_title_gp = gpar(fontsize = 5, fontface = "bold"),
      legend_title_gp = gpar(fontsize = 5, fontface = "bold"),
      legend_title_position = "topleft", 
      legend_labels_gp = gpar(fontsize = 5),
      legend_grid_width = unit(2, "mm"),
      legend_grid_height = unit(2, "mm"),
      legend_border = "black",
      heatmap_border = TRUE,
      annotation_border = TRUE,
      simple_anno_size = unit(2, "mm")
    )
    
    # Colour palettes
    method_pal <- 
      setNames(
        object = brewer_pal(palette = "Greys", direction = -1)(9)[c(2,4,7)],
        nm = c(data_to_plot$method %>%  
                 unique() %>% 
                 sort())
      )
    
    cell_type_pal <- 
      setNames(
        object = as.character(colour_schemes_list$cell_type$colour_hex),
        nm = colour_schemes_list$cell_type$group
      )
    
    # Count number of lists terms appear in
    term_count <-
      matrix(
        data = 
          c(mat[,str_detect(colnames(mat), pattern = "Bulk DS")] %>% 
              rowSums(),
            mat[,str_detect(colnames(mat), pattern = "Bulk PC")],
            mat[,str_detect(colnames(mat), pattern = "snRNA")] %>% 
              rowSums()),
        nrow = length(rownames(mat)), 
        ncol = 3,
        byrow = F,
        dimnames = list(
          rownames(mat),
          c("Bulk DS", "Bulk PC", "snRNA")
        )
      )
    
    term_count_ha <- 
      ComplexHeatmap::HeatmapAnnotation(
        count = 
          anno_barplot(term_count,
                       axis_param = 
                         list(
                           gp = gpar(fontsize = 5),
                           labels_rot = 0
                         ),
                       gp = gpar(fill = method_pal)
          ), 
        which = "row",
        width = unit(1.5, "cm"),
        annotation_label = "Count by analysis",
        annotation_name_gp = gpar(fontsize = 5, fontface = "bold"),
        annotation_name_rot = 0
      )
    
    simple_anno_df <- 
      data_to_plot %>% 
      dplyr::select(method, comparison, class) %>% 
      dplyr::distinct()
    
    comparison_split <- 
      simple_anno_df$comparison %>% 
      str_replace_all("_", " ")
    
    simple_ha <- 
      ComplexHeatmap::HeatmapAnnotation(
        df = 
          simple_anno_df %>% 
          dplyr::select(-comparison),
        annotation_legend_param =
          list(method = 
                 list(
                   title = "Analysis"
                 ),
               class = 
                 list(
                   title = "Cell type"
                 )
          ),
        col = 
          list(
            method = method_pal,
            class = cell_type_pal
          ),
        annotation_label = 
          c("Analysis", "Cell type"),
        annotation_name_gp = gpar(fontsize = 5, fontface = "bold"),
        na_col = "white"
      )
    
    ht <- 
      ComplexHeatmap::Heatmap(
        matrix = mat, 
        width = unit(5.5, "cm"),
        height = unit(8, "cm"),
        clustering_distance_columns = "pearson",
        column_split = comparison_split,
        column_dend_height = unit(1.5, "cm"),
        col = c("white", "black"),
        column_title_side = "bottom",
        column_title_rot = column_title_rotation,
        show_column_names = F,
        top_annotation = simple_ha,
        cluster_rows = F,
        row_names_side = "left",
        right_annotation = term_count_ha,
        rect_gp = gpar(col = "grey", lwd = 0.5),
        heatmap_legend_param = 
          list(
            title = "GO term", 
            labels = c("Present", "Absent"),
            direction = "horizontal"
          )
      )
    
    png(file.path(output_dir, str_c(output_name, ".png")),
        width = figure_width,
        height = figure_height,
        units = "mm",
        res = 300
    )
    draw(ht, merge_legend = T)
    dev.off()
    
  }

#' Process H-MAGMA and LDSC results.
#'
#' Function will process H-MAGMA and LDSC results that have been run using the
#' same gene sets. The function will categorise H-MAGMA and LDSC outputs by
#' p-value and FDR significance (multiple test correction by the number of cell
#' types). Gene sets that are significant using one method (FDR < 0.05) are
#' marked with **, while gene sets that were nominal with one method (p < 0.05)
#' are marked with *. If a gene set is significant (FDR < 0.05) using both H-MAGMA
#' and LDSC, it is assigned the category "Both**". If a gene set is nominally
#' significant (p < 0.05) using both H-MAGMA and LDSC, it is assigned the category
#' "Both*"
#'
#' @param ldsc_results df. LDSC output. Must have columns: "GWAS", "group",
#'   "direction_of_effect", "cell_type", "Z_score_P", "class" (factored cell
#'   types) and "class wrapped" (wrapped labels).
#' @param magma_results df. H-MAGMA output. Must have columns: "GWAS", "group",
#'   "direction_of_effect", "VARIABLE", "P", "class", "class_wrapped".
#' @param fct_levels chr. Vector with gene sets from "group" column ordered how
#'   user wishes groups to appear. Bear in mind any "_" in current comparison
#'   names will be replaced with " ", so this should also be reflected in the
#'   fct_levels vector.
#'
#' @return Dataframe with H-MAGMA and LDSC outputs categorised by significance.

process_magma_ldsc <- function(ldsc_results, magma_results, fct_levels){
  
  joint_results <- 
    magma_results %>% 
    dplyr::select(GWAS, group, direction_of_effect, cell_type = VARIABLE, magma_P = P, class, class_wrapped) %>%
    dplyr::left_join(ldsc_results %>% 
                       dplyr::select(GWAS, group, direction_of_effect, cell_type, ldsc_P = Z_score_P, class, class_wrapped)) %>% 
    dplyr::mutate(group = group %>% 
                    str_replace_all(., "_", " ") %>% 
                    fct_relevel(., 
                                fct_levels)) %>% 
    dplyr::group_by(GWAS, group, direction_of_effect) %>% 
    dplyr::mutate(ldsc_FDR = p.adjust(ldsc_P, method = "BH"),
                  magma_FDR = p.adjust(magma_P, method = "BH"),
                  p_category = case_when(ldsc_FDR < 0.05 & magma_FDR < 0.05 ~ "Both**",
                                         ldsc_FDR < 0.05 & magma_P < 0.05 ~ "LDSC** & H-MAGMA*",
                                         ldsc_P < 0.05 & magma_FDR < 0.05 ~ "H-MAGMA** & LDSC*",
                                         ldsc_P < 0.05 & magma_P < 0.05 ~ "Both*",
                                         magma_FDR < 0.05 ~ "H-MAGMA**",
                                         ldsc_FDR < 0.05 ~ "LDSC**",
                                         magma_P < 0.05 ~ "H-MAGMA*",
                                         ldsc_P < 0.05 ~ "LDSC*",
                                         TRUE ~ "None")) 
  
  
  
  return(joint_results)
  
}

#' Plot joint MAGMA/LDSC results.
#'
#' @param joint_results df. Output of \code{process_magma_ldsc}.
#' @param facet_by_direction logical. Facet results by direction of effect.
#'   Default = F.
#' @param facet_labels chr. Named vector with how GWAS facets should be
#'   labelled.
#'
#' @return Figure of joint MAGMA/LDSC results.

plot_magma_ldsc <- function(joint_results, facet_by_direction = F, facet_labels){
  
  p_colours <- 
    tibble(p_category = c("Both**", "LDSC** & H-MAGMA*", "H-MAGMA** & LDSC*", "Both*", "LDSC**", "H-MAGMA**", "LDSC*", "H-MAGMA*", "None"), 
           p_category_fct = p_category %>%
             fct_relevel(.,
                         c("Both**", "Both*", "LDSC** & H-MAGMA*", "LDSC**", "LDSC*", "H-MAGMA** & LDSC*", "H-MAGMA**", "H-MAGMA*", "None"))) %>% 
    dplyr::arrange(p_category_fct) %>% 
    dplyr::mutate(col_hex = c("#00A087FF", "#91D1C2FF", "#E64B35FF", "#f58b7aff", "#f6c6bfff", "#3c5488ff", "#7288beff", "#a9c1f6ff",  "#C1C1C1"))
  
  joint_results <- 
    joint_results %>% 
    dplyr::inner_join(p_colours)
  
  plot <- 
    joint_results %>% 
    ggplot(aes(x = group, y = fct_rev(class), fill = p_category_fct)) +
    geom_tile(colour = "black") +
    # coord_equal() +
    labs(x = "", y = "") +
    scale_fill_manual(limits = p_colours$p_category, values = p_colours$col_hex, drop = FALSE) +
    theme_rhr +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.key.size = unit(0.3, units = "cm"))
  
  if(facet_by_direction == T){
    
    plot <- 
      plot +
      facet_grid(
        cols = vars(GWAS), 
        rows = vars(direction_of_effect),
        labeller = labeller(GWAS = facet_labels)
      )
    
  } else{
    
    plot <- 
      plot +
      facet_grid(
        cols = vars(GWAS), 
        labeller = labeller(GWAS = facet_labels)
      )
    
  }
  
  return(plot)
  
}
