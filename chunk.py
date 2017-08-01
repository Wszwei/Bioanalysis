def chunk_pair():
	dirn = "/workspace/Python/chunk_ls/Chunks/"
	omit = .1
	chunk_to_pair, chunk_paired = chunks()
	to_pair_seq = []
	paired_seq = []
	pass



def chunks():
       chunk_to_pair = ["AA2", 
									 "FF1", 
									 "FF3", 
									 "FF5", 
									 "HH4", 
									 "II1", 
									 "L2", 
									 "M2", 
									 "M5", 
									 "Q1", 
									 "Q5", 
									 "R3", 
									 "X1", 
									 "X1", 
									 "X4", 
									 "X5", 
									 "Z1", 
									 "Z5", 
									 "Z6"]
       chunk_be_paired = ["AA4", 
                     "AA4", 
                     "AA5", 
                     "AA5", 
                     "AA6", 
                     "AA6", 
                     "BB1", 
                     "BB1", 
                     "BB1", 
                     "BB2", 
                     "BB2", 
                     "BB3", 
                     "BB3", 
                     "CC5", 
                     "CC5", 
                     "DD1", 
                     "DD1", 
                     "DD2", 
                     "DD2", 
                     "DD3", 
                     "DD3", 
                     "EE4", 
                     "EE4", 
                     "FF1", 
                     "FF1", 
                     "GG1", 
                     "GG1", 
                     "GG2", 
                     "GG3", 
                     "GG3", 
                     "GG3", 
                     "GG3", 
                     "GG4", 
                     "GG4", 
                     "II5", 
                     "JJ2", 
                     "JJ2", 
                     "JJ3", 
                     "JJ3", 
                     "JJ5", 
                     "K3", 
                     "K3", 
                     "K4", 
                     "K4", 
                     "L1", 
                     "L1", 
                     "L3", 
                     "L3", 
                     "L4", 
                     "L4", 
                     "L5", 
                     "L5", 
                     "M1", 
                     "M4", 
                     "M4", 
                     "N1", 
                     "N1", 
                     "N3", 
                     "N3", 
                     "O2", 
                     "O2", 
                     "O3", 
                     "O3", 
                     "O4", 
                     "O4", 
                     "O5", 
                     "P1", 
                     "P1", 
                     "P4", 
                     "P4", 
                     "Q2", 
                     "Q2", 
                     "Q4", 
                     "Q4", 
                     "Q5", 
                     "Q5", 
                     "R2", 
                     "R2", 
                     "S5", 
                     "S5", 
                     "T4", 
                     "T4", 
                     "T5", 
                     "T5", 
                     "U1", 
                     "V2", 
                     "V2", 
                     "W3", 
                     "W3", 
                     "W5", 
                     "W5", 
                     "X2", 
                     "X2", 
                     "X4", 
                     "Y3", 
                     "Y3"]
       chunk_p = []
       for c in chunk_be_paired:
       	if c not in chunk_p:
       		chunk_p.append(c)
       return chunk_to_pair, chunk_p



if __name__ == "__main__":
       op = open("out.txt", "w")
       cp, c = chunks()
       for ch in c:
              op.write("%s\n" %ch)
              print ch
       op.close()