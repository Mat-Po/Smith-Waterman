
import numpy as np

#### input

error = True
while error:
    seqA = input("please insert the first sequence to align: ")
    seqA = seqA.upper()
    lseqA = len(seqA)
    seqB = input("please insert the second sequence to align: ")
    seqB = seqB.upper()
    lseqB = len(seqB)
    match = input("please insert the score you want to give to a match: ")
    match = +(abs(int(match)))
    mismatch = input("please insert the score you want to give to a mismatch: ")
    mismatch = -(abs(int(mismatch)))
    gapPenalty = input("please insert the score you want to give to a gap: ")
    gapPenalty = -(abs(int(gapPenalty)))
    error = False

####1) defining the substitution matrix 

def substiMatrix(seqA, seqB, match, mismatch):
    
    """ 
    this function creates the substitution matrix for the two given
    sequences, given a match and mismatch score
    """
    lseqA = len(seqA)
    lseqB = len(seqB)

    substiMatrix = np.zeros((lseqA,lseqB),np.int) #initialisation of the matrix
    
    i = -1 #initializing index
    j = -1
    for baseA in seqA: #each letter in seqA
        i = i + 1 #updating position of the current letter
        j = -1
        for baseB in seqB:
                j = j+1
                if baseA == baseB: #for each letter in seqA, each in seqB
                    substiMatrix[i,j] = match #modifying the matrix with 
                                            #appropiate match-mismatch values
                else:
                    substiMatrix [i,j] = mismatch
    return substiMatrix
  
####2) defining the scoring matrix and a way of traceback the alignment

def ScoreTracebackMatrix (substiMatrix, gapPenalty):
    
    """
    This function create a Scoring matrix and a trace matrix needed for the traceback process,
    given a substitution matrix and a gap penalty score.
    This function return, as tuple,
    0)the scoring matrix
    1)the traceback matrix
    2)the highest score
    3)the row index of the highest score
    4)the column index of the highest score
    """
    #initialising the two matrices 
    scoringMatrix = np.zeros((lseqA+1,lseqB+1),np.int)
    traceMatrix = np.zeros((lseqA+1,lseqB+1), np.int)

    highestScore = 0 #initializing the best score

    for i in range (1, lseqA+1): #skipping the first zeros row
        for j in range (1, lseqB+1):
            #now i have to compute the score considering the tree different
            #possible direction:
            
            sx = scoringMatrix[i,j-1] + gapPenalty #value of the cell + penalty
            up = scoringMatrix[i-1,j] + gapPenalty #""""""
            diagonal = scoringMatrix [i-1,j-1] + substiMatrix[i-1,j-1] #values of the cell + match or mismatch
            scoringMatrix[i,j] = max(sx,up,diagonal,0) #inserting the max value
            
            #updaing the traceback matrix
            if scoringMatrix [i,j] == sx:
                traceMatrix [i,j] = 1
            if scoringMatrix [i,j] == up:
                traceMatrix [i,j] = 2
            if scoringMatrix [i,j] == diagonal:
                traceMatrix [i,j] = 3
            if scoringMatrix [i,j] == 0:
                traceMatrix [i,j] = 0
                
            #updatating the highest score if needed:
            if scoringMatrix[i][j] >= highestScore:
                highestScore = scoringMatrix[i][j]
                #coordinates of the top score
                letter_i = i #row
                letter_j = j #column

    return scoringMatrix, traceMatrix, highestScore, letter_i, letter_j 

#3)doing the actual alignment

def traceback (traceMatrix,letter_i,letter_j):
    """
    This function does the real alignment of the two sequences taking as input
    the trace back matrix and starting from the highest score cell.
    It returns the two, locally aligned, sequences, number of match, mismatches, gaps.
    """
    i = letter_i
    j = letter_j
    align1= ""
    align2= ""

    match = 0
    mismatch = 0 
    gap = 0

    while traceMatrix[i,j] != 0: #stop when i get zero!
        if traceMatrix[i,j] == 3: #-->from diagonal
            pos1 = seqA[i-1] #getting the corrispondent letter
            pos2 = seqB[j-1]
            i = i-1
            j = j-1
            if pos1 == pos2:
                match = match +1    #updating num of matches
            else:
                mismatch = mismatch +1   #updating num of mismatches

        elif traceMatrix [i,j] == 2: #-->from above
            pos1 = seqA[i-1]
            pos2 = "-"        #gap insertion
            i -= 1
            gap = gap +1      #updating num of gaps
        elif traceMatrix [i,j] == 1: #--> from left
            pos1 = "-"         #gap
            pos2 = seqB [j-1]
            j -= 1
            gap = gap +1       #updating num of gaps
        
        align1 += pos1
        align2 += pos2
    align1 = align1[::-1] #reversing the sequences
    align2 = align2[::-1]

    return align1, align2, match, mismatch, gap



####output:


#print(seqA,seqB)

substitutionMatrix = substiMatrix(seqA,seqB,match,mismatch)
score_traceback_matrix = ScoreTracebackMatrix(substitutionMatrix,gapPenalty)
scoringMatrix = score_traceback_matrix[0]
alignment = traceback(score_traceback_matrix[1],score_traceback_matrix[3],score_traceback_matrix[4])
max_score = score_traceback_matrix[2]
#print("\nbest score:{} ".format(score_traceback_matrix[2]))
#print("\nthe alignment:")
#print(alignment[0])
#print(alignment[1])
#print("\nlenght of the sub-alignment: {}".format(len(alignment[0])))
#print("number of matches: {}".format(alignment[2]))
#print("number of mismatches: {}".format(alignment[3]))
#print("number of gaps: {}".format(alignment[4]))

print("\n")
print("trace matrix:\n",score_traceback_matrix[1])    #trace matrix
print("\nscoring matrix:\n",score_traceback_matrix[0])     #score matrix
#print("\nsubstitution matrix:",substitutionMatrix)
#print("best score indexes:",score_traceback_matrix[3],score_traceback_matrix[4])


print("\nhieghst score:{}".format(max_score))
print("\n")

#the idea is to search in scoring matrix where score is greater the 80% of the best score
#once i founde that cells, from that i can do the traceback and save the infos in 
#a list of dictionaries. Once i have it i can use the lambda function to sort it on the 
#key that i prefer, in this case the number of gaps.
TOTalignment = []
for i in range (lseqA+1):
    for j in range (lseqB+1):
        if scoringMatrix [i,j] > round((max_score*80)/100):
            al1,al2,m,mm,g = traceback(score_traceback_matrix[1],i,j)
            TOTalignment.append({"alignment1": al1, "alignment2": al2, "matches": m, 
            "mismatches": mm, "gaps": g, "length": len(al1), "score": scoringMatrix[i,j]})


#print(TOTalignment)
sortedTOTal = sorted(TOTalignment, key = lambda x: x["gaps"], reverse = False)
for i in sortedTOTal:
    print("gaps:{}, match:{}, mismatch:{}, length:{}, score:{}".format(
        i["gaps"],i["matches"],i["mismatches"],i["length"],i["score"]
    ))
    
    print(i["alignment1"])
    print(i["alignment2"])
    print("#"*100)
    print("\n")






