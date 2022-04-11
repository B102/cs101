#Project Milestone 2
#More info: https://relate.cs.illinois.edu/course/cs101-fa21/page/ms2/
#Please go here and check your code against all possible cases: https://relate.cs.illinois.edu/course/cs101-fa21/file-version/6e82d86863d1e012a33175992d2ca1455f1d13ec/resources/ms-tests.pdf
# def find_splice(dna_motif, dna)
def find_splice(dna_motif, dna):
  p = [-1] #The list of positions. Start with -1 to check for the position 0
  loc_dna = 0 #We are using this so we don't have to check the entire dna string every time. We only need to look for the part after the previous nucleotide.
  for motif in dna_motif:
    for i in range(len(dna[loc_dna:])):
      if motif == dna[i+loc_dna] and i+loc_dna > p[-1]:
        #Check to see if the motif is in the dna string and the position is after the previous one. This is why we are affing the -1.
        p += [i+loc_dna]
        loc_dna = i+loc_dna
        break #We don't need to check the rest if we find one
  if len(p)-1 < len(dna_motif):
    #If some motif is not in the dna string, the list of position should be less than the last of motifs. Then we should return the empty list
    return []
  return p[1:] #cutting off the -1 in front.

# def shared_motif(dna_list)
def shared_motif(dna_list):
  first = dna_list[0] #first element
  first_len = len(first) #length of the first element
  for i in range(first_len): 
    #This for-loop controls the length of the substring. As we loop back, i increases, so the substring will be one shorter.
    start = 0 #we start at the first letter
    step = first_len - i #the range of the substring will then be from start to start+step
    while start+step <= first_len:
      check = True #this serves as a way of remembering if the substring is in the other strings
      for dna in dna_list: #check other strings to see if the substring is in them
        if first[start:start+step] not in dna:
          #if it is not in one of the strings, then we can go all the way back and start with another substring
          check = False #to remember this is indeed not in other strings
          break
      if check: #This step is necessary because break only terminates the current loop, so we need something to remember if the substring is in all other strings. If so, return the current substring as all other substrings will be shorter than or the same length as it.
        return first[start:start+step]
      start +=1 #if it is not, then we start at the next letter in the shortest string and try that substring.
  return '' #return the empty string if we cannot find any common substring

# def get_edges(dna_dict)
def get_edges(dna_dict):
  result = []
  #We essentially compare every two distinct ROSALIND identifiers and see if they are adjacent
  for first in dna_dict.keys(): #Our first identifier
    for second in dna_dict.keys(): #second identifier
      if first != second and dna_dict[first][-3:] == dna_dict[second][:3]: 
        result += [(first, second)] #If they are adjacent, add tuple to result
  return result  

# def assemble_genome(dna_list)
def assemble_genome(dna_list):
  r = two_assem(dna_list)
  while(len(ec(dna_list,r))==0):
    if(two_assem(r)==[]):
      break
    r=two_assem(r)
  fl=ec(dna_list,r)
  return lmf(fl)

#ec的作用：list B里所有包含”list A的所有项“的项，例如ec(["a","c","b"],["ab","abc"])会返回['abc']。
#How to use ec: check all elements in listB that contain all of the listA
def ec(lista,listb):
  l=[]
  for eb in listb:
    for i in range(0,len(lista)):
      if lista[i] not in eb:
        break
      if (i==(len(lista)-1)):
        l.append(eb)
        break
  return l

#two_assem的作用：把一个list里的所有项目两两组和
#How to use two_assem: add elements in a list(in pairs) 
def two_assem(dna_list):
  r = []
  for n in range(0,len(dna_list)):
    if n+1 == len(dna_list):break
    for m in range(n+1,len(dna_list)):
      r.append(merge(dna_list[n],dna_list[m],0))
  return r+dna_list

#merge的作用：text A与text B按先后顺序合并，如果无法合并就调换顺序再次尝试，如果还是不行就直接返回两个text的总和
#How to use merge: combine textA and textB. If they can't be merged, just add them up directly
def merge(textA,textB,pn):
  r=''
  lo = None
  if len(textA) > len(textB):
    lo = len(textA) 
  else: 
    lo = len(textB)
  for le in range(-lo,0):
    if textA[le:] == textB[0:-le]:
      r=textA+textB.replace(textA[le:],'',1)
      break
  if (r=='' and pn==0):
    return merge(textB,textA,1)
  if (r=='' and pn==1):
    r=textB+textA
  return r

#lmf的作用：一个list里最短项目
def lmf(x):
  x.sort(key = lambda y: len(y))
  return x[0]

# def perfect_match(rna)
def perfect_match(rna):
  rna = rna.upper()
  A = rna.count('A')
  U = rna.count('U')
  C = rna.count('C')
  G = rna.count('G')
  AU = 1 
  CG = 1 
  if A != U or C != G: 
    return 0 #No perfect match if the AU and CG doesnt have the same number
  while A>1:
    AU *= A #Number of combinations for AU
    A -= 1 
  while C>1:
    CG *= C #Number of combinations for CG
    C -= 1
  return AU*CG #total number of combinations

# def random_genome(dna, gc_content)
def random_genome_for_one(dna, gc_content):
  #Function for calculating one probability given one gc content rate
  from math import log
  dna = dna.upper()
  #Dictionary for probability of each type of nuleotide
  P_match = {'G':gc_content/2, 'C':gc_content/2, 'A':(1-gc_content)/2,'T':(1-gc_content)/2}
  prob = 1
  for na in dna:
    prob *= P_match[na] #multiply all probability together
  return log(prob,10) #Log base 10 result
def random_genome(dna, gc_content):
  result = []
  for gc_rate in gc_content: #Use the previous functio to iterate over all gc content rates
    result += [random_genome_for_one(dna, gc_rate)]
  return result

# def rev_palindrome(dna)
def rev_palindrome(dna):
  tuplist = []  #create empty list
  for i in range(len(dna)):  #iterates through all indices of dna
    for j in range(4,13):  #iterates through all possible substring lengths from 4 to 12
      if i+j >len(dna):   #if the substring length goes out of bounds we break from the loop
        break
      sub = dna[i:i+j]
      revcomp = reverse_complement(sub)
      if revcomp == sub:   #if the substring created is equal to its reverse complement we record its start position (i) and its length (j)
        tuplist = tuplist + [(i, j)]
  return tuplist
  
def reverse_complement(dna):
    DNA=dna.upper() #Make sure it's upper case
    DNA_reverse_comp = '' #Make a string to which we add the complement letters in reverse
    switch={'A':'T', 'G':'C', 'T':'A', 'C':'G'} #Dictionary for switching the letter into the complement
    for letter in DNA[::-1]:
      DNA_reverse_comp+=switch[letter] #Add the complement into our string in reverse
    return DNA_reverse_comp
