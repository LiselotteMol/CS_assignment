TVs_merged <- read_excel("Downloads/data_excel.xlsx")

# Step 1: Data Cleaning
TV_list <- data.frame(TVs_merged)
titles <- TV_list %>% dplyr::select(title, Name)

model_words <- list()
for(i in 1:nrow(titles)) { 
  
  titles[i,1] <- tolower(titles[i,1])
  titles[i,1] <- gsub('Hertz|hertz|Hz|HZ| hz|-hz|hz', 'hz', titles[i,1])
  titles[i,1]  <- gsub('Inch|inches|"|-"|-inch| inch|inch|”', 'inch', titles[i,1])
  titles[i,1] <- gsub('-|(|)|(-|-)|[()]|[[]]|,|series|diagonal|diag|/|of|class|newegg|neweggcom|amazon|amazoncom|bestbuycom|bestbuy|best|buy|the|com|tv', "",titles[i,1])
  word_ended <- FALSE
  j <- 1
  
  while(word_ended == FALSE){
    word <- word(titles[i,1], start = j)
    word <- gsub('\\[|\\]', "", word)
    j <- j + 1
    if(is.na(word)){
      word_ended <- TRUE 
    }
    else {
      if(word %in% model_words)  {
        
      } else if (grepl("([a-z0-9]*(([0-9]+[ˆ0-9, ])|([ˆ0-9, ]+[0-9]+))[a-z0-9]*)", word)) {
        model_words <- append(model_words, word)
      }
    }
  }
}

# Step 2: Minhashing
# Step 2a: Creating a binary input matrix
names = TV_list[,1]
bin_matrix = matrix(0, length(model_words), nrow(titles))
input_matrix = data.frame(bin_matrix, row.names = model_words)
colnames(input_matrix) = names

for(i in 1:nrow(titles)){
  for(j in 1:length(model_words)) {
    if (grepl(model_words[[j]],titles[i,1])) {
      input_matrix[j,i] = 1
    }
  }
}

#Step 2b: Minhashing
hashes <- 600 
signatures = data.frame(matrix(Inf, hashes, nrow(titles)))
colnames(signatures) = names

a <- list(sample(1 : 600, size = 600, replace = T))
b <- list(sample(1 : 600, size = 600, replace = T))
a <- as.vector(unlist(a))
b <- as.vector(unlist(b))
for(i in 1:nrow(input_matrix)) {
  for(j in 1:ncol(input_matrix)) {
    if(input_matrix[i,j] == 1) {
      for(k in 1:hashes) {
        value = (a[k] + b[k] * i) %% 617
        if(value < signatures[k,j]) {
          signatures[k,j] = value
        }
      }
    }
  }
}
        

#Split train en test data
set.seed(3)
boot <- sample(1:1624,replace = TRUE)
boot <- sort(unique(boot))
test_index = setdiff(c(1:1624), boot)


bootstrap_train <- TV_list[boot,]
bootstrap_test  <-TV_list[test_index,]

# Step 3: Local-Sensitive Hashing
# examine scores for different values of r
signature_train = signatures[,boot]

signature_test = signatures[,test_index]
r = seq(1, nrow(signature_test), 5)

## Find actual duplicates
real_dupl_test <- real_dupl(bootstrap_test)
real_dupl_train <- real_dupl(bootstrap_train)

input_test = input_matrix[,test_index]
input_train = input_matrix[,boot]

index = 1
f1_star_results_4 = rep(NA, length(r))
pq_results_4 = rep(NA, length(r))
pc_results_4 = rep(NA, length(r))
candidate_test = func_candidates(signature_test, 20)
f1_star(candidate_test, real_dupl_test, bootstrap_test$modelID)


for(i in seq_along(r)) {
  candidate_test = candidates(i, signature_test)
  results_4 = f1_star(candidate_test, real_dupl_test, bootstrap_test$modelID)
  f1_star_results_4[index]= results_4[3]
  pq_results_4[index] = results_4[1]
  pc_results_4[index] = results_4[2]
  index = index + 1
}
f1_star_results_4 <- as.matrix(f1_star_results_4)


max_h_4 = rep(NA, length(r))
index = 1
for(i in seq_along(r)) {
  candidate_train = candidates(i, signature_train)
  dis_matrix_train = dis_matrix(candidate_train, bootstrap_train, input_train)
  
  hc_train <- agnes(dis_matrix_train, diss = TRUE, method = "complete")
  h_seq = seq(0,1000, 100)
  f1_hr <- list()
  for(m in 1:length(h_seq)) {
    h = h_seq[m]
    cut_cluster_train <- cutree(hc_train, h = h)
    duplicates_found_train <- cut_cluster_train[duplicated(cut_cluster_train)]
    clusters_train <- list()
    for(j in 1:length(duplicates_found_train)) {
      clusters_train[[j]] = which(cut_cluster_train %in% duplicates_found_train[[j]])
    }
    clusters_train <- unique(clusters_train)
    f1_hr[m] <- f1(clusters_train, real_dupl_train)
  }
  
  maximum <- max(unlist(f1_hr), na.rm = TRUE)
  best_h <- which(f1_hr == maximum)
  max_h_4[index] <- h_seq[best_h[1]]
  index = index + 1
}

possible_comparisons_4 = rep(NA, length(r))
foc_4 = rep(NA, length(r))
f1_test_4 = rep(NA, length(r))

index = 1
for(i in seq_along(r)) {
  candidate_test = candidates(i, signature_test)
  dis_matrix_test = dis_matrix(candidate_test, bootstrap_test, input_test)
  
  hc_test <- agnes(dis_matrix_test, diss = TRUE, method = "complete")
  cut_test <- cutree(hc_test, h = max_h_4[index])
  
  duplicates_test <- cut_test[duplicated(cut_test)]
  clusters_test <- list()
  for(k in 1:length(duplicates_test)) {
    clusters_test[[k]] = which(cut_test %in% duplicates_test[[k]])
  }
  clusters_test <- unique(clusters_test)
  
  possible_comparisons_4 = (nrow(bootstrap_test) *(nrow(bootstrap_test)-1)) / 2
  foc_4[index] = length(which(candidate_test==1)) / possible_comparisons_4
  
  f1_test_4[index] = f1(clusters_test, real_dupl_test)
  index = index + 1
}
