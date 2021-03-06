#########################################################
#######        COMPUTATIONAL BIOLOGY         ############
#######             HOMEWORK 2               ############
#########################################################
#                                                       #
# Simulate the evolution of sequences on a fixed tree   #
# under the TN93 sequence evolution model               #
#                                                       #
#########################################################
#########################################################

#########################################################
######    Code below this should not be changed   #######
#########################################################

library(Matrix)
library(ape)

nucleotides = c("T", "C", "A", "G")

get_nucleotide_from_number = function(i) {
    # Transform a nucleotide number into the appropriate letter
    return(nucleotides[i])
}

transform_to_nucleotides = function(sequence) {
    # Transform a sequence of nucleotide numbers into the appropriate character sequence
    nucl_sequence = paste0(lapply(sequence, get_nucleotide_from_number), collapse = "")
    return(nucl_sequence)
}

#########################################################
######    Code above this should not be changed   #######
#########################################################

# In all functions the following parameters are the same:
# pi: the stationary frequencies of nucleotides
# alpha1: rate coefficient for the C <-> T transition
# alpha2: rate coefficient for the A <-> G transition
# beta: rate coefficient for transversions
# N: number of sites in the simulated alignment

create_TN93_Q_matrix = function(pi, alpha1, alpha2, beta) {
    # Create the TN93 transition rate matrix Q as specified in the assignment.
    
    Q <- matrix(c(0,alpha1*pi[2], beta*pi[3], beta*pi[4],
                  alpha1*pi[1], 0, beta*pi[3], beta*pi[4],
                  beta*pi[1], beta*pi[2], 0, alpha2*pi[4],
                  beta*pi[1], beta*pi[2], alpha2*pi[3], 0), nrow=4, byrow=TRUE)
    for (i in 1:4){
      Q[i,i] <- -(rowSums(Q)[i])
    }
    # Return the transition rate matrix
    # Q: 4 by 4 matrix of rates
    return(Q)
}

get_starting_nucleotide = function(pi) {
    # Sample a starting nucleotide from the stationary distribution
    
    nucleotide <- sample(1:4, size = 1,
         prob = pi)

    # Return the sampled nucleotide
    # nucleotide: integer nucleotide value
    return(nucleotide)
}

get_starting_sequence = function(pi, N) {
    # Sample a starting sequence of length N
  
    starting_sequence = c()
    for (i in 1:N){
      starting_sequence = c(starting_sequence,get_starting_nucleotide(pi))
    }
    
    # Return the sampled sequence
    # starting_sequence: vector of integer nucleotide values
    return(starting_sequence)
}

get_evolved_sequence = function(sequence, branch_length, Q) {
    # Evolve a given nucleotide sequence along a branch of specified length.
    # sequence: nucleotide sequence at the beginning of the branch
    # branch_length: the length of the branch along which evolution happens
    # Q: the transition rate matrix
 
    P_tb <- expm(branch_length*Q)
    evolved_sequence <- c()
    for (i in sequence){
      evolved_sequence <- c(evolved_sequence, get_starting_nucleotide(P_tb[i,]))
    }
    
    # Return the nucleotide sequence after all positions have evolved along the given branch.
    # evolved_sequence: the vector of new integer nucleotide values at the end of the branch
    return(evolved_sequence)
}


simulate_evolution = function(newick_tree, pi, alpha1, alpha2, beta, N) {
    # Simulate evolution along the given tree.
    # newick_tree: the tree in newick text format
    
    # Transfrom the tree from text format to an object of the phylo class which represents the tree in R
    tree = read.tree(text = newick_tree)
    # Reorder the tree for easier traversing
    tree = reorder(tree, order = "cladewise")
    
    # Set up the Q matrix
    Q <- create_TN93_Q_matrix(pi, alpha1, alpha2, beta)
    
    # Set up the starting sequence @ the root of the tree
    starting_sequence <- get_starting_sequence(pi, N)
    
    # Prepare a list to store evolved sequences at each node
    sequence_per_node = list()
    sequence_per_node[[tree$edge[1,1]]] = starting_sequence
    
    # Walk the tree while evolving sequences along appropriate branches
    for (i in 1:length(tree$edge.length)) {
        node_parent = tree$edge[i, 1]
        node_child = tree$edge[i, 2]
        branch_length = tree$edge.length[i]
        parent_sequence = sequence_per_node[[node_parent]]
        
        child_sequence <- get_evolved_sequence(parent_sequence, branch_length, Q)

        sequence_per_node[[node_child]] = child_sequence
    }

    # Transform the alignment from nucleotide indices to nucleotide characters
    # and filter out the sequences at the tips
    alignment = list()
    for (i in 1:length(tree$tip.label)) {
        alignment[[tree$tip.label[i]]] = transform_to_nucleotides(sequence_per_node[[i]])
    }
    
    # Return the simulated alignment.
    # The alignment should be in the form of a list where the tip label corresponds to the
    # appropriate simulated sequence, e.g. alignment$human = ACTG
    return(alignment)
}

test_simulation = function() {
  library(ape)
  newick_tree = "(orangutan:13,(gorilla:10.25,(human:5.5,chimp:5.5):4.75):2.75);"
  N = 40
  beta = 0.035
  alpha1 = 0.044229
  alpha2 = 0.021781
  pi = c(0.22, 0.26, 0.33, 0.19)
  result = simulate_evolution(newick_tree, pi, alpha1, alpha2, beta, N)
  print(result$orangutan)
  print(result$gorilla)
  print(result$human)
  print(result$chimp)
  result
}

minimal_t <- function(Q) {
  res <- expm(1000*Q)
  t=1
  while (sum(expm(t*Q)==res)!=16){
    t <- t+1
  }
  return(t)
}

minimal_t <- function(Q) {
  res <- expm(1000*Q)
  t=1
  while (sum(abs(expm(t*Q)-res))>1e-10){
    t <- t+1
  }
  return(t)
}
