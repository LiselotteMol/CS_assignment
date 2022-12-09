install.packages("sets")
install.packages("factoextra")
install.packages("dendextend")
library(sets)
library(readxl)
library(tidyverse)
library(numbers) 
library(textreuse)
library(knitr)
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms

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

#Step 3: Local-Sensitive Hashing
# examine scores for different values of r

r = seq(1, nrow(signatures), 5)

## Find actual duplicates
real_duplicates <- real_dupl(TV_list)

f1_star_results = rep(0,length(r))
index = 1
for(i in seq_along(r)) {
  candidate = candidates(i, signatures)
  f1_star_results[index]= f1_star(candidate, real_duplicates, TV_list$modelID)[3]
  index = index + 1
}
f1_star_results



# Step 4: Calculating dissimilarity matrix
#candidate <- candidates(20, signatures)
#dissim_matrix = dis_matrix(candidate, TV_list, input_matrix)

for(i in seq_along(r)) {
  candidate <- candidates(20, signatures)
  dissim_matrix = dis_matrix(candidate, TV_list, input_matrix)
  
  hc1 <- agnes(dissim_matrix, diss = TRUE, method = "complete")
  h_seq = seq(0,1000, 100)
  f1_results <- rep(0, length(h_seq))
  for(h in seq_along(h_seq)) {
    cut_cluster <- cutree(hc1, h = h)
    duplicates_found <- cut_cluster[duplicated(cut_cluster)]
    clusters <- list()
    for(j in 1:length(duplicates_found)) {
      clusters[[j]] = which(cut_cluster %in% duplicates_found[[j]])
    }
    clusters <- unique(clusters)
    print(f1(clusters, real_duplicates = real_duplicates))
    f1_results[i] <- f1(clusters, real_duplicates = real_duplicates)
  }
}

############## Bootstrap results
f1_star_combi = cbind(f1_star_results, f1_star_results_2, f1_star_results_3, f1_star_results_4, f1_star_results_5)
f1_star_ave = rowMeans(f1_star_combi)

pc_results_combi = cbind(pc_results, pc_results_2, pc_results_3, pc_results_4, pc_results_5)
pq_results_combi = cbind(pq_results, pq_results_2, pq_results_3, pq_results_4, pq_results_5)

pc_results_ave = rowMeans(pc_results_combi)
pq_results_ave = rowMeans(pq_results_combi)

f1_test_combi = cbind(f1_test, f1_test_2, f1_test_3, f1_test_4, f1_test_5)
f1_test_ave = rowMeans(f1_test_combi)


max_h_combi = cbind(max_h, max_h_2, max_h_3, max_h_4, max_h_5)
max_h_ave = rowMeans(max_h_combi)


plot(foc_3, f1_star_ave, type = "l",  xlab = "Fraction of comparisons", ylab = "F1*-measure")
plot(foc_3, pc_results_3, type = "l", xlab = "Fraction of comparisons", ylab = "Pair completeness")
plot(foc_3, pq_results_ave, type = "l", xlab = "Fraction of comparisons", ylab = "Pair quality")
plot(foc_3, f1_test_ave, type = "l", xlab = "Fraction of comparisons", ylab = "F1-measure")

################## Functions #########################
#Jaccard Similarity
jaccard <- function(x,y) {
  return(length(intersect(x,y)) / length(union(x,y)))
}

#Function to make buckets in LSH 
buckets <- function(vector, r) { 
  b = floor(length(vector)/r)
  r_sequence = seq(0, length(vector), r)
  buckets <- list()
  for(i in 1:b) {
    begin = r_sequence[i] + 1
    end = r_sequence[i] + r
    bucket = paste(vector[begin:end], collapse = "_")
    buckets <- append(buckets, bucket)
  }
  return(buckets)
}

#Function to find candidate pairs
candidates <- function(r, signature_matrix) {
  candidate_matrix = matrix(0, ncol(signature_matrix), ncol(signature_matrix))
  all_buckets <- list()
  for(i in 1:ncol(signature_matrix)) {
    all_buckets[[i]] = buckets(signature_matrix[,i], r)
  }
  
  for(i in 1:ncol(signature_matrix)) {
    bucket_i = all_buckets[[i]]
    for(j in 1:ncol(signature_matrix)) {
      if(i < j) {
        bucket_j = all_buckets[[j]]
        for(k in 1:length(bucket_i)) {
          if(bucket_i[[k]] == bucket_j[[k]]) {
            candidate_matrix[i,j] = 1
            break
          }
        }
      }
    }
  }
  return(candidate_matrix)
}

f1_star <- function(candidate_matrix, real_duplicates, names) {
  coordinates <- which(candidate_matrix == 1, arr.ind=TRUE)
  duplicates_found <- 0 
  real <- 0
  for(i in 1:nrow(coordinates)) {
    if(names[coordinates[i,1]] == names[coordinates[i,2]]) {
      duplicates_found <- duplicates_found + 1
    }
  }
  
  for(j in 1:length(real_duplicates)) {
    real <- real + choose(length(real_duplicates[[j]]),2)
  }
  
  compl <- duplicates_found/real
  qual <- duplicates_found/nrow(coordinates)
  F1_star <- (2*qual*compl)/(qual + compl)
  return(c(qual, compl, F1_star))
}

real_dupl <- function(df) {
  real_duplicates = list()
  for(i in 1:nrow(df)) {
    for(j in 1:nrow(df)) {
      if(df$modelID[i] == df$modelID[j] && i != j) {
        if(length(real_duplicates) > 0) {
          index1 <- grep(df$modelID[i], names(real_duplicates))
          index2 <- grep(df$modelID[j], names(real_duplicates))
          if(length(index1) > 0) {
            real_duplicates[[index1]]<- append(real_duplicates[[index1]], c(i,j))
          } else if(length(index2) > 0) {
            real_duplicates[[index2]]<- append(real_duplicates[[index2]], c(i,j))
          } else {
            list_names <- c(names(real_duplicates), df$modelID[i])
            real_duplicates[[length(real_duplicates)+1]] = c(i,j)
            names(real_duplicates) = list_names
          }
        } else {
          list_names <- c(names(real_duplicates), df$modelID[i])
          real_duplicates[[length(real_duplicates)+1]] = c(i,j)
          names(real_duplicates) = list_names
        }
      }
    }
  }
  for(i in 1:length(real_duplicates)) {
    real_duplicates[[i]] <- real_duplicates[[i]][!duplicated(real_duplicates[[i]])]
  }
  return(real_duplicates)
}


dis_matrix <- function(candidate, df, input) {
  dissim_matrix = matrix(0, nrow(candidate), ncol(candidate))
  for(i in 1:nrow(candidate)) {
    for(j in 1:ncol(candidate)) {
      if(i<= j) {
        if(df$shop[i] == df$shop[j]) {
          dissim_matrix[i,j] = 1000
          dissim_matrix[j,i] = 1000
        } else if (!is.na(df$Brand[i]) && !is.na(df$Brand[j]) && df$Brand[i] != df$Brand[j]) {
          dissim_matrix[i,j] = 1000
          dissim_matrix[j,i] = 1000
        } else if(candidate[i,j] == 0) {
          dissim_matrix[i,j] = 1000
          dissim_matrix[j,i] = 1000
        } else {
          dissim_matrix[i,j] = 1 - jaccard(input[,i], input[,j])
          dissim_matrix[j,i] = dissim_matrix[i,j]
        }
      }
    }
  }
  return(dissim_matrix)
}

f1 <- function(clusters, real_duplicates) {
  real <- 0
  truepos <- 0
  positives <- 0
  for(j in 1:length(real_duplicates)) {
    real <- real + choose(length(real_duplicates[[j]]),2)
  }
  
  for(i in 1:length(clusters)) {
    positives <- positives + choose(length(clusters[[i]]), 2) 
    for(j in 1:length(real_duplicates)) {
      number_of_products = length(intersect(as.set(clusters[[i]]), real_duplicates[[j]]))
      if(!is_empty(intersect(as.set(clusters[[i]]), real_duplicates[[j]]))) {
        truepos <- truepos + choose(number_of_products,2)
      }
    }
  }

  falsepos <- positives - truepos
  falseneg <- real - truepos
  
  recall <- truepos / (truepos + falsepos)
  precision <- truepos / (truepos + falseneg)
  f1_score <- 2 * precision * recall / (precision + recall)
  return(f1_score)
} 


