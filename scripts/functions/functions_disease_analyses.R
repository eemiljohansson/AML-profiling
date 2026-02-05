
# Run differential expression using limma
de_limma_disease <-
  function(data_wide, 
           metadata,
           disease,
           correct = T) {
    
    # Filter data (cohort, variable is 0 or 1)
    dat <-
      data_wide %>% 
      inner_join(metadata %>% 
                   select(DAid, Sex, Age, Disease), by = "DAid") %>% 
      rename(Group = Disease) %>% 
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control")) 
    
    # Design a model - add Group, and Sex, Age, BMI (removed for HIV)
    if(correct == T) {
      design <- model.matrix(~0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age) 
      colnames(design) <- c("control", "case",  "Sex", "Age")
    } else {
      design <- model.matrix(~0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <- makeContrasts(Diff = case - control, levels = design)
    
    # Fit linear model to each protein assay
    dat_fit <- 
      dat %>% 
      select(-Sex, -Age, -Group)  %>% 
      column_to_rownames("DAid") %>% 
      t()
    
    fit <- lmFit(dat_fit, design=design,  method="robust", maxit=10000)
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit)
    
    # Extract DE results
    DE_results <-
      topTable(ebays_fit,
               n=nrow(ebays_fit$p.value), 
               adjust.method="fdr", 
               confint=TRUE)
    
    DE_res <- 
      DE_results %>% 
      as_tibble(rownames = "Assay") %>% 
      mutate(Disease = disease)
    
    return(DE_res)
  }

# Generate balanced controls from the other cancers
generate_balanced_controls <- function(cancer, metadata, set_data, n_control_groups) {
  
  set_data <-
    set_data %>% 
    filter(Disease != "CAD healthy")
  
  if(cancer %in% c("BRC","CVX","ENDC","OVC")) {
    
    female_samples <- 
      metadata %>% 
      filter(Sex == "Female") %>% 
      pull(DAid)
    
    set_data <- 
      set_data %>% 
      filter(DAid %in% female_samples)
    
  } else if (cancer == "PRC") {
    
    male_samples <- 
      metadata %>% 
      filter(Sex == "Male") %>% 
      pull(DAid)
    
    set_data <- 
      set_data %>% 
      filter(DAid %in% male_samples)
    
  } else {
    set_data <- 
      set_data
  }
  
  control_pool <- 
    set_data %>% 
    filter(Disease != cancer)
  
  cancer_samples_train <- 
    set_data %>% 
    filter(Disease == cancer,
           set == "train") %>% 
    pull(DAid) %>% 
    length()
  
  n_train <- 
    cancer_samples_train/n_control_groups 
  
  set.seed(213)
  train_set <- 
    control_pool %>%
    filter(set == "train") %>% 
    group_by(Disease) %>% 
    sample_n(size = min(ceiling(n_train), n())) %>% 
    bind_rows(cancers_split %>% 
                filter(Disease == cancer,
                       set == "train")) %>% 
    mutate(Disease = case_when(Disease == cancer ~ Disease,
                               T ~ "Other"))
  
  cancer_samples_test <- 
    set_data %>% 
    filter(Disease == cancer,
           set == "test") %>% 
    pull(DAid) %>% 
    length()
  
  n_test <- 
    cancer_samples_test/n_control_groups
  
  set.seed(213)
  test_set <- 
    control_pool %>%
    filter(set == "test") %>% 
    group_by(Disease) %>% 
    sample_n(size = min(ceiling(n_test), n())) %>% 
    bind_rows(cancers_split %>% 
                filter(Disease == cancer,
                       set == "test")) %>% 
    mutate(Disease = case_when(Disease == cancer ~ Disease,
                               T ~ "Other"))
  
  train_set %>%
    bind_rows(test_set) 
  
}

# Generate balanced controls from the IGT cohort
generate_balanced_controls_igt <- function(cancer, set_data, metadata) {
  
  # Number of samples in train and test
  n_train <-
    set_data %>% 
    filter(Disease == cancer,
           set == "train") %>% 
    pull(DAid) %>% 
    unique() %>% 
    length()
  
  n_test <-
    set_data %>% 
    filter(Disease == cancer,
           set == "test") %>% 
    pull(DAid) %>% 
    unique() %>% 
    length()
  
  # Female / male cancers 
  if(cancer %in% c("BRC","CVX","ENDC","OVC")) {
    
    female_samples <- 
      metadata %>% 
      filter(Sex == "Female") %>% 
      pull(DAid)
    
    set_data <- 
      set_data %>% 
      filter(DAid %in% female_samples)
    
  } else if (cancer == "PRC") {
    
    male_samples <- 
      metadata %>% 
      filter(Sex == "Male") %>% 
      pull(DAid)
    
    set_data <- 
      set_data %>% 
      filter(DAid %in% male_samples)
    
  } else {
    set_data <- 
      set_data
  }
  
  # Random subset IGT
  set.seed(213)
  set_data %>% 
    filter(Disease == "CAD healthy",
           set == "train") %>% 
    sample_n(size = min(ceiling(n_train), n())) %>% 
    bind_rows(set_data %>% 
                filter(Disease == "CAD healthy",
                       set == "test") %>% 
                sample_n(size = min(ceiling(n_test), n()))  
    ) %>% 
    bind_rows(set_data %>% 
                filter(Disease == cancer))
  
}
