library(shiny)
library(tidyverse)
library(ggpubr)
theme_set(theme_bw())



calculate_outcome <- function(row, stat_col, sig_thres) {
    pass_thres <- as.numeric(row[stat_col]) < sig_thres
    is_de <- row["type"] == "DE"
    if (pass_thres && is_de) "True pos."
    else if (pass_thres && !is_de) "False pos."
    else if (!pass_thres && is_de) "False neg."
    else "True neg."
}

#raw_data %>% ggplot(aes(x=pval)) + geom_histogram(bins=100) + scale_fill_manual(values=c("steelblue", "gray"))

color_levels <- c("orange", "darkred", "steelblue", "gray")
names(color_levels) <- c("False neg.", "False pos.", "True pos.", "True neg.")

calc_stats <- function(outcome_col, top, bottom) {
    outcome_table <- table(outcome_col)
    if (!(bottom %in% names(outcome_table))) 1
    else round(outcome_table[top] / (outcome_table[top] + outcome_table[bottom]), 2)
}

# features <- 1000
# sig_features <- 50
# base_level <- 20
# diff <- 5
# noise <- 1


ui <- fluidPage(

    titlePanel("Inspection of cutoffs for P-values and FDR corrections"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("diff", "Regulated features difference", min=0, max=10, value=5, step=0.1),
            sliderInput("noise", "Standard deviation", min=0, max=20, value=1, step=0.1),
            sliderInput("features", "Background features", min=0, max=10000, value=1000, step=50),
            sliderInput("sig_features", "Regulated features", min=0, max=1000, value=50, step=5),
            checkboxInput("advanced_settings", "Show advanced settings", value=FALSE),
            conditionalPanel(
                "input.advanced_settings == 1",
                sliderInput("p_thres", "P-value threshold", min=0, max=1, value=0.05, step=0.01),
                sliderInput("fdr_thres", "FDR Threshold", min=0, max=1, value=0.1, step=0.01),
                sliderInput("base_level", "Base level", min=0, max=50, value=20, step=1),
                sliderInput("bins", "Number of bins:", min = 1, max = 200, value = 100, step=1),
                numericInput("seed", "Random seed", min=0, step=1, value=37)
            )
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("source_distributions", height = 200),
           plotOutput("pvalue_histogram_pcut", height = 200),
           plotOutput("pvalue_histogram_bonfcut", height = 200),
           plotOutput("pvalue_histogram_bhcut", height = 200)
        )
    )
)

make_ggplot <- function(target_data, outcome_col, title, color_levels, number_bins=100) {
    
    target_data %>% 
        ggplot(aes_string(x="pval", fill=outcome_col)) + 
        geom_histogram(bins=number_bins) + 
        scale_fill_manual(values=color_levels) + 
        ggtitle(title) + 
        xlab("P-value") + 
        ylab("Count") + 
        labs(fill="Outcome")
}

# raw_data %>% ggplot(aes(x=pval, fill=fdr_bh_outcome)) + geom_histogram(bins=100) + scale_fill_manual(values=color_levels) + ggtitle(fdr_bh_title) + xlab("") + ylab("Count"),
# raw_data %>% ggplot(aes(x=pval, fill=fdr_bonf_outcome)) + geom_histogram(bins=100) + scale_fill_manual(values=color_levels) + ggtitle(fdr_bonf_title) + xlab("P-value") + ylab(""), 

server <- function(input, output) {

    simulated_clean_data <- reactive({
        set.seed(input$seed)
        
        get_conda <- function() {
            rnorm(input$features+input$sig_features, mean=input$base_level, sd=input$noise)
        }
        
        get_condb <- function() {
            c(
                rnorm(input$features, mean=input$base_level, sd=input$noise), 
                rnorm(input$sig_features, mean=input$base_level+input$diff, sd=input$noise)
            )
        }
        
        data.frame(
            type = c(rep("background", input$features), rep("DE", input$sig_features)),
            s1 = get_conda(),
            s2 = get_conda(),
            s3 = get_conda(),
            s4 = get_condb(),
            s5 = get_condb(),
            s6 = get_condb()
        )
    })
    
    results_df <- reactive({
        
        target_data <- simulated_clean_data()
        
        target_data$pval <- apply(target_data[, 2:7], 1, function(row) {
            t.test(row[1:3], row[4:6])$p.value
        })
        target_data$fdr_bh <- p.adjust(target_data$pval, method = "BH")
        target_data$fdr_bonf <- p.adjust(target_data$pval, method = "bonferroni")
        
        target_data$pval_outcome <- apply(target_data, 1, calculate_outcome, stat_col="pval", sig_thres=input$p_thres)
        target_data$fdr_bh_outcome <- apply(target_data, 1, calculate_outcome, stat_col="fdr_bh", sig_thres=input$fdr_thres)
        target_data$fdr_bonf_outcome <- apply(target_data, 1, calculate_outcome, stat_col="fdr_bonf", sig_thres=input$fdr_thres)
        
        target_data
    })
    
    output$source_distributions <- renderPlot({
        source_background <- rnorm(1000, mean=input$base_level, sd=input$noise)
        source_diff <- rnorm(1000, mean=input$base_level+input$diff, sd=input$noise)
        plot_df <- rbind(
            data.frame(type="Background", value=source_background),
            data.frame(type="Diff", value=source_diff)
        )
        plot_df %>%
            ggplot(aes(x=value, color=type)) + geom_density()
    })
    
    output$pvalue_histogram_pcut <- renderPlot({
        
        target_data <- results_df()
        pcut_title <- sprintf("P-value < 0.05 (Sensitivity: %s, Precision: %s)", 
                         calc_stats(target_data$pval_outcome, "True pos.", "False neg."), 
                         calc_stats(target_data$pval_outcome, "True pos.", "False pos."))
        make_ggplot(target_data, "pval_outcome", pcut_title, color_levels, number_bins=input$bins)
    })
    
    output$pvalue_histogram_bonfcut <- renderPlot({
        target_data <- results_df()
        fdr_bh_title <- sprintf("FDR (BH) < 0.1 (Sensitivity: %s, Precision: %s)", 
                                calc_stats(target_data$fdr_bh_outcome, "True pos.", "False neg."), 
                                calc_stats(target_data$fdr_bh_outcome, "True pos.", "False pos."))
        make_ggplot(target_data, "fdr_bh_outcome", fdr_bh_title, color_levels, number_bins=input$bins)
    })
    
    output$pvalue_histogram_bhcut <- renderPlot({
        target_data <- results_df()
        fdr_bonf_title <- sprintf("FDR (Bonf.) < 0.1 (Sensitivity: %s, Precision: %s)", 
                                  calc_stats(target_data$fdr_bonf_outcome, "True pos.", "False neg."), 
                                  calc_stats(target_data$fdr_bonf_outcome, "True pos.", "False pos."))
        make_ggplot(target_data, "fdr_bonf_outcome", fdr_bonf_title, color_levels, number_bins=input$bins)
    })
}

shinyApp(ui = ui, server = server)
