#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################

# In all functions the following parameters are the same:
# seqA: the first sequence to align
# seqB: the second sequence to align
# score_gap: score for a gap
# score_match: score for a character match
# score_mismatch: score for a character mismatch
# local: True if alignment is local, False otherwise

init_score_matrix = function(nrow, ncol, local, score_gap) {
    # Initialize the score matrix with zeros.
    # If the alignment is global, the leftmost column and the top row will have incremental gap scores,
    # i.e. if the gap score is -2 and the number of columns is 4, the top row will be [0, -2, -4, -6].
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix

    score_matrix <- matrix(0, nrow, ncol)
    
    if (local==FALSE){
      score_matrix[-1,1]<- c(seq(-2, -(2*(nrow-1)), by=-2))
      score_matrix[1,-1]<- c(seq(-2, -(2*(ncol-1)), by=-2))
    }
    
    # Return the initialized empty score matrix
    # score_matrix: nrow by ncol matrix
    return(score_matrix)
}

init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values, except the top row and the leftmost column if global alignment.
    # If global alignment, the top row has "left" on all positions except 1st.
    # Similarly, leftmost column has "up" on all positions except 1st.
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix

    path_matrix <- matrix('', nrow, ncol)
    
    if (local==FALSE){
      path_matrix[-1,1]<- 'up'
      path_matrix[1,-1]<- 'left'
    }

    # Return the initialized empty path matrix
    # path_matrix: nrow by ncol matrix
    return(path_matrix)
}

get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix
    # score_matrix: the score_matrix that is being filled out
    
    dia <- score_matrix[row - 1, col - 1] + score_mismatch
    if (nucA==nucB){
      dia <- score_matrix[row - 1, col - 1] + score_match
    }
  
    gapB <- score_matrix[row-1, col]+score_gap
    gapA <- score_matrix[row, col-1]+score_gap
    
    #local alignment
    
    if (local){
      score <- max(dia, gapB, gapA, 0)
      path <- '-'
      if (score==score_matrix[row-1, col]+score_gap){
        path <- 'up'
      } else if (score==score_matrix[row, col-1]+score_gap){
        path <- 'left'
      } else if(score==score_matrix[row - 1, col - 1] + score_match | score==score_matrix[row - 1, col - 1] + score_mismatch){
          path <- 'diag'
        }
      } else { #global alignment
      score <- max(dia, gapB, gapA)
      if (score==score_matrix[row-1, col]+score_gap){
        path <- 'up'
      } else if (score==score_matrix[row, col-1]+score_gap){
        path <- 'left'
      } else if (score == score_matrix[row - 1, col - 1] + score_match | score==score_matrix[row - 1, col - 1] + score_mismatch){
          path <- 'diag'
        } 
      }
      
    # Return the best score for the particular position in the score matrix
    # In the case that there are several equally good paths available, return any one of them.
    # score: best score at this position
    # path: path corresponding to the best score, one of ["diag", "up", "left"] in the global case and of ["diag", "up", "left", "-"] in the local case
    return(list("score"=score, "path"=path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    # score_matrix: initial matrix of the scores
    # path_matrix: initial matrix of paths

    for (i in 2:nrow(score_matrix)){
      for (j in 2:ncol(score_matrix)){
        nucA <- substr(seqA, i-1, i-1)
        nucB <- substr(seqB, j-1, j-1)
        sp_list <- get_best_score_and_path(i, j, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local)
        score_matrix[i,j] <- sp_list$score
        path_matrix[i,j] <- sp_list$path
      }
    }

    # Return the full score and path matrices
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}

get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    # i.e. if the path is diagonal both the characters in seqA and seqB should be added,
    # if the path is up or left, there is a gap in one of the sequences.
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # path: best path pre-computed for the given position
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix

    newrow <- row-1
    newcol <- col-1
    char1 <- nucA
    char2 <- nucB
  
    if (path=='up'){
      newcol <- col
      char2 <- '-'
    } else if (path=='left'){
      newrow <- row
      char1 <- '-'
    }
  
    
    # Return the new row and column and the aligned characters
    # newrow: row if gap in seqA, row - 1 otherwise
    # newcol: col if gap in seqB, col - 1 otherwise
    # char1: '-' if gap in seqA, appropriate character if a match
    # char2: '-' if gap in seqB, appropriate character if a match
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}

get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths

    #initialize empty allignments
    aliA <- ""
    aliB <- ""
    
    #initialize current position and get score
    if (local){
      pi <- which(score_matrix == max(score_matrix), arr.ind = TRUE)[[1,1]]
      pj <- which(score_matrix == max(score_matrix), arr.ind = TRUE)[[1,2]]
    } else {
      pi <- nrow(score_matrix)
      pj <- ncol(score_matrix)
    }
    
    score <- score_matrix[pi,pj]
    
    #get the alignment
    if (local){ #for Smith-Waterman
      while (score_matrix[pi,pj]>0){
        path <- path_matrix[pi,pj]
        move <- get_best_move(substr(seqA, pi-1, pi-1), substr(seqB, pj-1, pj-1), path, pi, pj)
        pi <- move$newrow
        pj <- move$newcol
        aliA <- paste0(move$char1, aliA)
        aliB <- paste0(move$char2, aliB)
      }
      alignment <- c(aliA, aliB)
    } else { #for Needleman_Wunsch
      while (pi>1 || pj>1){
        path <- path_matrix[pi,pj]
        move <- get_best_move(substr(seqA, pi-1, pi-1), substr(seqB, pj-1, pj-1), path, pi, pj)
        pi <- move$newrow
        pj <- move$newcol
        aliA <- paste0(move$char1, aliA)
        aliB <- paste0(move$char2, aliB)
      }
      alignment <- c(aliA, aliB) 
    }
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # score: score of the best alignment
    # alignment: the actual alignment in the form of a vector of two strings
    return(list("score"=score, "alignment"=alignment))
}

align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # For testing purposes, use seqA for the rows and seqB for the columns of the matrices
  
    # Initialize score and path matrices
    score_matrix <- init_score_matrix(nchar(seqA)+1, nchar(seqB)+1, local, score_gap)
    path_matrix <- init_path_matrix(nchar(seqA)+1, nchar(seqB)+1, local)
  
    # Fill in the matrices with scores and paths using dynamic programming
    matrices <- fill_matrices(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)
    score_matrix <- matrices$score_matrix
    path_matrix <- matrices$path_matrix
    
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    result <- get_best_alignment(seqA, seqB, score_matrix, path_matrix, local)
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # Returns the same value types as get_best_alignment
    return(result)
}

test_align = function() {
    seqA = "TCACACTA"
    seqB = "AGCACACA"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = T
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()
