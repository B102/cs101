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
  li=list(dna)
  for n in range(0,len(dna)):
    if li[n] == 'T':
      li[n] = 'U'
  return li


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
  p_het=het/total #probability of getting a heterozygous organism 
  p_rec=rec/total #probability of getting a homozygous recessive organism
  #We are calculating the probability of getting a double recessive genome. The probability of getting a dominant allele is then 1 minus that.
  p_AaAa=1/4* (p_het)**2 #The probability of getting two heterozygous and that of getting an 'aa' resulting genome if we get two heterozygous.
  p_aaaa=(p_rec)**2 #The case for two homozygous recessive
  p_Aaaa=(p_rec)*(p_het) #The case for one 'Aa' and one 'aa'. 
  p_dom=1-(p_AaAa+p_aaaa+p_Aaaa) #Probability of getting a dominant allele is everything else.
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
def gc_content(dna_list):
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
    triplet = triplet.upper()
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
    amino = genetic_code.get(triplet)
    return amino
def rna2codon(triplets): #Turn an rna sequence into codons
    triplets = triplets.upper()
    amino=''
    for i in range( 0,int( len( triplets ) / 3 ) ):
        triplet= triplets[ 3*i:3*i+3 ]
        amino = amino + triplet2codon(triplet)
    return amino


#8. locate_substring(dna_snippet, dna)
def locate_substring(ds,d):
  li=[]
  l=len(ds)
  i=0
  while True:
    if ds.find(d)==-1:
      break
    li.append(ds.find(d))
    ds=ds[ds.find(d)+1:l]
    i=ds.find(d)
  return li


#9. hamming_dist(dna1, dna2)



#10. count_dom_phenotype(genotypes)



#11. source_rna(protein)



#12. splice_rna(dna, intron_list)


