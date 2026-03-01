# Schumacher-Advanced-Computational-Biology-Python-BISC-4503
This is the portfolio that I have created for BISC 4503 Advanced Computational Biology Winter 2026


#Sequence Objects pt 1,2,3,4

```python




```python
from Bio.Seq import Seq 
```


```python
my_seq = Seq("GATCG")
```


```python
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```

    0 G
    1 A
    2 T
    3 C
    4 G



```python
# We can also print the length of each sequence
print(len(my_seq))
```

    5



```python
# We can also print the letter using the length 
print (my_seq[0])
```

    G



```python
print (my_seq[4])
```

    G



```python
print (my_seq[2])
```

    T



```python
#We can also count patterns within a sequence 
Seq("AAAA").count("AA")
```




    2




```python
# We can count the number of "G" within a seq
my_seq = Seq("GATTTGAGAAATCCCGCGCTTTA")
```


```python
len(my_seq)
```




    23




```python
# We can also count the number of repeating letters within a sequence 
my_seq.count("G")
```




    5




```python
# We can find the percentage of G and C content within our sequences 
100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
```




    43.47826086956522




```python
from Bio.SeqUtils import gc_fraction
```


```python
my_seq = Seq("GATTTGAGAAATCCCGCGCTTTA")
```


```python
gc_fraction(my_seq)
```




    0.43478260869565216




```python
# We can choose sections of our sequence 
my_seq[4:12]
```




    Seq('TGAGAAAT')




```python
my_seq[0::3]
```




    Seq('GTAACGCT')




```python
my_seq[1::3]
```




    Seq('ATGACCTA')




```python
my_seq[2:3]
```




    Seq('T')




```python
my_seq[::-1]
```




    Seq('ATTTCGCGCCCTAAAGAGTTTAG')




```python
str(my_seq)
```




    'GATTTGAGAAATCCCGCGCTTTA'




```python
fasta_format_string = ">Name\n%s\n" % my_seq
```


```python
print(fasta_format_string)
```

    >Name
    GATTTGAGAAATCCCGCGCTTTA
    



```python
seq1 = Seq("ACGT")
seq2 = Seq("AACCGG")
```


```python
seq1 + seq2
```




    Seq('ACGTAACCGG')




```python
seq2+seq1
```




    Seq('AACCGGACGT')




```python
contigs = [Seq("ATG"), Seq("ATCCCG"), Seq("TTGCA")]
```


```python
spacer = Seq("N" *10)
```


```python
spacer.join(contigs)
```




    Seq('ATGNNNNNNNNNNATCCCGNNNNNNNNNNTTGCA')




```python
dna_seq = Seq("acgtACGT")
```


```python
dna_seq 
```




    Seq('acgtACGT')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
dna_seq.lower()
```




    Seq('acgtacgt')




```python
dna_seq.upper()
```




    Seq('ACGTACGT')




```python
"gtac" in dna_seq
```




    False




```python
"GTAC" in dna_seq
```




    False




```python
dna_seq
```




    Seq('acgtACGT')




```python
dna_seq = dna_seq.upper()
```


```python
"GTAC" in dna_seq
```




    True




```python
my_seq = Seq("GTAGCTAGGTACAATCGATTAGTAGCAGCTAGCTATA")
```


```python
my_seq.complement()
```




    Seq('CATCGATCCATGTTAGCTAATCATCGTCGATCGATAT')




```python
my_seq.reverse_complement()
```




    Seq('TATAGCTAGCTGCTACTAATCGATTGTACCTAGCTAC')




```python
protein_seq = Seq("EVRNAK")
protein_seq.complement()
```




    Seq('EBYNTM')




```python
coding_dna = Seq("ATCGATCGTTAGCTTAACGCGCTAGCTATCGA")
```


```python
coding_dna
```




    Seq('ATCGATCGTTAGCTTAACGCGCTAGCTATCGA')




```python
template_dna = coding_dna.reverse_complement()
```


```python
template_dna
```




    Seq('TCGATAGCTAGCGCGTTAAGCTAACGATCGAT')




```python
coding_dna
```




    Seq('ATCGATCGTTAGCTTAACGCGCTAGCTATCGA')




```python
# We can transcribe dna sequences 
messenger_rna = coding_dna.transcribe()
```


```python
messenger_rna
```




    Seq('AUCGAUCGUUAGCUUAACGCGCUAGCUAUCGA')




```python
template_dna.reverse_complement().transcribe()
```




    Seq('AUCGAUCGUUAGCUUAACGCGCUAGCUAUCGA')




```python
messenger_rna.back_transcribe()
```




    Seq('ATCGATCGTTAGCTTAACGCGCTAGCTATCGA')




```python
messenger_rna
```




    Seq('AUCGAUCGUUAGCUUAACGCGCUAGCUAUCGA')




```python
# We can translate rna sequence to proteins 
messenger_rna.translate()
```

    /home/student/anaconda3/lib/python3.7/site-packages/Bio/Seq.py:2808: BiopythonWarning: Partial codon, len(sequence) not a multiple of three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
      BiopythonWarning,





    Seq('IDR*LNALAI')




```python
coding_dna.translate(table="Vertebrate Mitochondrial")
```




    Seq('IDR*LNALAI')




```python
coding_dna.translate(table =2)
```




    Seq('IDR*LNALAI')




```python
coding_dna.translate(to_stop = True)
```




    Seq('IDR')




```python
coding_dna.translate(table = 2, to_stop=True)
```




    Seq('IDR')




```python
coding_dna.translate(table = 2, stop_symbol = "!")
```




    Seq('IDR!LNALAI')




```python
gene =Seq("ATGCGATGCCCTTCCTTATGCTATTCGATCGCAGCTACGATCGATCGATGTTATCGATCATCCTATCGCGCGATATCGATGATCGATCATCGATCGATCGATCGATCGUGA")
```


```python
#We can translate our gene to bacteria 
gene.translate(table = "Bacterial")
```




    Seq('MRCPSLCYSIAATIDRCYRSSYRAISMIDHRSIDRS*')




```python
gene.translate(table="Bacterial", to_stop = True)
```




    Seq('MRCPSLCYSIAATIDRCYRSSYRAISMIDHRSIDRS')




```python
gene.translate(table= "Bacterial", cds = True)
```




    Seq('MRCPSLCYSIAATIDRCYRSSYRAISMIDHRSIDRS')




```python
#We can generate Codon Table 
from Bio.Data import CodonTable 
```


```python
standard_table = CodonTable. unambiguous_dna_by_name["Standard"]
```


```python
mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
```


```python
print(standard_table)
```

    Table 1 Standard, SGC0
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
    T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I   | ACT T   | AAT N   | AGT S   | T
    A | ATC I   | ACC T   | AAC N   | AGC S   | C
    A | ATA I   | ACA T   | AAA K   | AGA R   | A
    A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V   | GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
print(mito_table)
```

    Table 2 Vertebrate Mitochondrial, SGC1
    
      |  T      |  C      |  A      |  G      |
    --+---------+---------+---------+---------+--
    T | TTT F   | TCT S   | TAT Y   | TGT C   | T
    T | TTC F   | TCC S   | TAC Y   | TGC C   | C
    T | TTA L   | TCA S   | TAA Stop| TGA W   | A
    T | TTG L   | TCG S   | TAG Stop| TGG W   | G
    --+---------+---------+---------+---------+--
    C | CTT L   | CCT P   | CAT H   | CGT R   | T
    C | CTC L   | CCC P   | CAC H   | CGC R   | C
    C | CTA L   | CCA P   | CAA Q   | CGA R   | A
    C | CTG L   | CCG P   | CAG Q   | CGG R   | G
    --+---------+---------+---------+---------+--
    A | ATT I(s)| ACT T   | AAT N   | AGT S   | T
    A | ATC I(s)| ACC T   | AAC N   | AGC S   | C
    A | ATA M(s)| ACA T   | AAA K   | AGA Stop| A
    A | ATG M(s)| ACG T   | AAG K   | AGG Stop| G
    --+---------+---------+---------+---------+--
    G | GTT V   | GCT A   | GAT D   | GGT G   | T
    G | GTC V   | GCC A   | GAC D   | GGC G   | C
    G | GTA V   | GCA A   | GAA E   | GGA G   | A
    G | GTG V(s)| GCG A   | GAG E   | GGG G   | G
    --+---------+---------+---------+---------+--



```python
mito_table.stop_codons
```




    ['TAA', 'TAG', 'AGA', 'AGG']




```python
mito_table.start_codons
```




    ['ATT', 'ATC', 'ATA', 'ATG', 'GTG']




```python
seq = Seq("ACGT")
```


```python
"ACGT" == seq1
```




    True




```python
seq1 == "ACGT"
```




    True




```python
unknown_seq = Seq(None, 10)
```


```python
unknown_seq
```




    Seq(None, length=10)




```python
len(unknown_seq)
```




    10




```python
seq = Seq({117512683: "TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT"}, length = 159345973)
```


```python
seq[1000:1020]
```




    Seq(None, length=20)




```python
seq[117512690:117512700]
```




    Seq('CCTGAATGTG')




```python
seq[117512670:]
```




    Seq({13: 'TTGAAAACCTGAATGTGAGAGTCAGTCAAGGATAGT'}, length=41833303)




```python
seq = Seq("ACGT")
```


```python
undefined_seq = Seq(None, length =10)
```


```python
seq + undefined_seq + seq 
```




    Seq({0: 'ACGT', 14: 'ACGT'}, length=18)




```python
my_seq = Seq("GCCATTGTAATGGGCCGTGAAAGGGTGCCCGA")
```


```python
from Bio.Seq import MutableSeq
```


```python
mutable_seq = MutableSeq(my_seq)
```


```python
mutable_seq
```




    MutableSeq('GCCATTGTAATGGGCCGTGAAAGGGTGCCCGA')




```python
mutable_seq[5] = "C"
```


```python
mutable_seq
```




    MutableSeq('GCCATCGTAATGGGCCGTGAAAGGGTGCCCGA')




```python
mutable_seq.remove("T")
```


```python
mutable_seq
```




    MutableSeq('GCCACGTAATGGGCCGTGAAAGGGTGCCCGA')




```python
mutable_seq.reverse()
```


```python
mutable_seq
```




    MutableSeq('AGCCCGTGGGAAAGTGCCGGGTAATGCACCG')




```python
new_seq = Seq(mutable_seq)
```


```python
new_seq
```




    Seq('AGCCCGTGGGAAAGTGCCGGGTAATGCACCG')




```python
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
```


```python
my_string = "GCGATCGATCGGAGGCTTAGCTAGCTAGGCTACGATGAGTGATA"
```


```python
reverse_complement(my_string)
```




    'TATCACTCATCGTAGCCTAGCTAGCTAAGCCTCCGATCGATCGC'




```python
transcribe(my_string)
```




    'GCGAUCGAUCGGAGGCUUAGCUAGCUAGGCUACGAUGAGUGAUA'




```python
back_transcribe(my_string)
```




    'GCGATCGATCGGAGGCTTAGCTAGCTAGGCTACGATGAGTGATA'




```python
translate(my_string)
```




    'AIDRRLS*LGYDE*'


