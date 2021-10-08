#Shared file for group project.
#Information for the first milestone:
#https://relate.cs.illinois.edu/course/cs101-fa21/page/ms1/

#1. s(dna)
def s(dna):
  DNA=dna.upper() #make sure it's upper case
  count={'A':0, 'T':0, 'G':0, 'C':0} #create a dictionary for counting letters
  for letter in DNA:
    count[letter] += 1 #count and add to the according value
  return count


#2. dna2rna(dna)
def dna2rna(dna):
  li=list(dna) #make dna into a list of its nucleotides so we can edit them
  for n in range(0,len(dna)):
    if li[n] == 'T': #check and replace T with U
      li[n] = 'U'
  return ''.join(li) #join it back together into a string


#3. reverse_complement(dna)
def reverse_complement(dna):
  DNA=dna.upper() #Make sure it's upper case
  DNA_reverse_comp = '' #Make a string to which we add the complement letters in reverse
  switch={'A':'T', 'G':'C', 'T':'A', 'C':'G'} #Dictionary for switching the letter into the complement
  for letter in DNA[::-1]:
    DNA_reverse_comp+=switch[letter] #Add the complement into our string in reverse
  return DNA_reverse_comp


#4. mendels_law(hom, het, rec)
def mendels_law(hom, het, rec):
  total=hom+het+rec #total number of organisms
  if total == 1:
    return False
  p_het=het/total #probability of getting a heterozygous organism 
  p_rec=rec/total #probability of getting a homozygous recessive organism
  #We are calculating the probability of getting a double recessive genome. The probability of getting a dominant allele is then 1 minus that.
  p_AaAa=1/4* p_het*((het-1)/(total-1)) #The probability of getting two heterozygous and that of getting an 'aa' resulting genome if we get two heterozygous.
  p_aaaa=(p_rec)*(rec-1)/(total-1) #The case for two homozygous recessive
  p_Aaaa=1/2*(rec)/(total-1)*(p_het) 
  p_aaAa=1/2*p_rec*(het/(total-1)) #The case for one 'Aa' and one 'aa'. 
  p_dom=1-(p_AaAa+p_aaaa+p_Aaaa+p_aaAa) #Probability of getting a dominant allele is everything else.
  return p_dom


#5. fibonacci_rabbits(n, k)
def fibonacci_rabbits(n,k):
  # k represents the number of baby rabbit pairs produced by one adult pair per month
  # n represents the number of months
  baby_pairs = 1 # we begin with one pair of baby rabbits
  adult_pairs = 0
  
  for i in range(1, n): 
    # every month, there are k new baby pairs produced for each adult pair
    new_pairs = adult_pairs*k 
    # in each month the babies grow to become adults
    adult_pairs += baby_pairs
    # the new babies produced in a month remain babies through the next month
    baby_pairs = new_pairs
  return baby_pairs + adult_pairs


#6. gc_content(dna_list)
def GC_content(dna_list):
  place=0 #the index of the string w/ the highest gc content. Starts at 0
  largest=0 #Ratio of the largest gc content
  for index in range(len(dna_list)): 
    dna_string=dna_list[index].upper() #check each string
    g=dna_string.count('G') #number of g
    c=dna_string.count('C') #number of c
    total=len(dna_string) #total nucleotides
    ratio=(g+c)/total #ratio
    if ratio>largest: #renew the location and ratio if it's larger
      place = index
      largest = ratio
  return (place, largest*100)


#7. rna2codon(rna)
def triplet2codon(triplet): #Turning one triplet to one codons
    triplet = triplet.upper() #make sure it's upper case
    #A dictionary for the corresponding amino acids and codons
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    allowed_codons = set('ACGU')
    for symbol in triplet:
        if symbol not in allowed_codons:
            return "Invalid" 
    amino = genetic_code.get(triplet) #get the corresponding amino acide
    return amino


def rna2codon(triplets): #Turn an rna sequence into codons
  #into upper case
  triplets = triplets.upper()
  #create an empty string
  amino=''
  #look at the rna sequence by codons
  for i in range( 0,int( len( triplets ) / 3 ) ): 
    triplet= triplets[ 3*i:3*i+3 ]
    #if we have the stop codon, we end the program
    if triplet2codon(triplet)=='*':
      return amino
    #accumulate string
    amino += triplet2codon(triplet)
  return amino


#8. locate_substring(dna_snippet, dna)
def locate_substring(d,ds):
  #create an empty list
  li=[]
  #get the length of the original string 
  l=len(ds)
  #len(ds) is a big number. Large enough to account for the situation where every single element of your string is the substring you're looking for
  for i in range(0,len(ds)):
    #if the substring DNE in the big string, break
    if ds.find(d)==-1:break
    #put the index of the substring into the list li
    li.append((ds.find(d))+(l-len(ds)))
    #cut the original string
    ds=ds[ds.find(d)+1:l]
  return li


#9. hamming_dist(dna1, dna2)
def hamming_dist(dna1, dna2):
  d=0
  for l in range(0,len(dna1)):
    if dna1[l]!=dna2[l]: #check each nucleotide to see if they're different
      d+=1
  return d


#10. count_dom_phenotype(genotypes)
def count_dom_phenotype(genotypes):
  AAAA=2*genotypes[0] #AAAA f1 generation will always have dominant allele
  AAAa=2*genotypes[1] #AAAa f1 always dominant
  AAaa=2*genotypes[2] #AAaa f1 always dominant
  AaAa=2*3/4*genotypes[3] #AaAa will have 3/4 dominant
  Aaaa=2*1/2*genotypes[4] #Aaaa will have 1/2 dominant
  #aaaa will have no dominant
  return AAAA+AAAa+AAaa+AaAa+Aaaa #add them together


#11. source_rna(protein)
def source_rna(protein):
  #We are using the amino acide dict to make a new dict that uses amino acids as keys to get how many different codons wil allow for it
  genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',

        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',

        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',

        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
  protein += '*'
  a_to_num = {} #Dictionary for how many unique codons an amino acid represents
  #set up the dictionary
  for amino in genetic_code.values():
    a_to_num[amino] = 0 #First set the keys in place
  for amino in genetic_code.values():
    a_to_num[amino] += 1 #Count the number of ways
  count = 1
  for amino in protein:
    count *= a_to_num[amino] #use the dict to get all the possibility for each amino acid
  return count


#12. splice_rna(dna, intron_list)
def splice_rna(dna, intron_list):
  for intro in intron_list:
    dna=dna.replace(intro, '') #replace introns with empty string, essentially removing them
  rna = dna2rna(dna) #turn into rna sequence
  codon = rna2codon(rna) #turn into amino acid sequence
  return codon

