# Schumacher-Advanced-Computational-Biology-Python-BISC-4503
This is the portfolio that I have created for BISC 4503 Advanced Computational Biology Winter 2026


# Sequence Objects pt 1,2,3,4

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





# Sequence Annotations pt 1,2,3,4



    ```python
from Bio.SeqRecord import SeqRecord 
```


```python
help(SeqRecord)
```

    Help on class SeqRecord in module Bio.SeqRecord:
    
    class SeqRecord(builtins.object)
     |  SeqRecord(seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |  
     |  A SeqRecord object holds a sequence and information about it.
     |  
     |  Main attributes:
     |   - id          - Identifier such as a locus tag (string)
     |   - seq         - The sequence itself (Seq object or similar)
     |  
     |  Additional attributes:
     |   - name        - Sequence name, e.g. gene name (string)
     |   - description - Additional text (string)
     |   - dbxrefs     - List of database cross references (list of strings)
     |   - features    - Any (sub)features defined (list of SeqFeature objects)
     |   - annotations - Further information about the whole sequence (dictionary).
     |     Most entries are strings, or lists of strings.
     |   - letter_annotations - Per letter/symbol annotation (restricted
     |     dictionary). This holds Python sequences (lists, strings
     |     or tuples) whose length matches that of the sequence.
     |     A typical use would be to hold a list of integers
     |     representing sequencing quality scores, or a string
     |     representing the secondary structure.
     |  
     |  You will typically use Bio.SeqIO to read in sequences from files as
     |  SeqRecord objects.  However, you may want to create your own SeqRecord
     |  objects directly (see the __init__ method for further details):
     |  
     |  >>> from Bio.Seq import Seq
     |  >>> from Bio.SeqRecord import SeqRecord
     |  >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |  ...                    id="YP_025292.1", name="HokC",
     |  ...                    description="toxic membrane protein")
     |  >>> print(record)
     |  ID: YP_025292.1
     |  Name: HokC
     |  Description: toxic membrane protein
     |  Number of features: 0
     |  Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |  
     |  If you want to save SeqRecord objects to a sequence file, use Bio.SeqIO
     |  for this.  For the special case where you want the SeqRecord turned into
     |  a string in a particular file format there is a format method which uses
     |  Bio.SeqIO internally:
     |  
     |  >>> print(record.format("fasta"))
     |  >YP_025292.1 toxic membrane protein
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  <BLANKLINE>
     |  
     |  You can also do things like slicing a SeqRecord, checking its length, etc
     |  
     |  >>> len(record)
     |  44
     |  >>> edited = record[:10] + record[11:]
     |  >>> print(edited.seq)
     |  MKQHKAMIVAIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  >>> print(record.seq)
     |  MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |  
     |  Methods defined here:
     |  
     |  __add__(self, other)
     |      Add another sequence or string to this sequence.
     |      
     |      The other sequence can be a SeqRecord object, a Seq object (or
     |      similar, e.g. a MutableSeq) or a plain Python string. If you add
     |      a plain string or a Seq (like) object, the new SeqRecord will simply
     |      have this appended to the existing data. However, any per letter
     |      annotation will be lost:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = record + "ACT"
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNNACT
     |      >>> print(list(new.letter_annotations))
     |      []
     |      
     |      The new record will attempt to combine the annotation, but for any
     |      ambiguities (e.g. different names) it defaults to omitting that
     |      annotation.
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      
     |      Now let's cut the plasmid into two pieces, and join them back up the
     |      other way round (i.e. shift the starting point on this plasmid, have
     |      a look at the annotated features in the original file to see why this
     |      particular split point might make sense):
     |      
     |      >>> left = plasmid[:3765]
     |      >>> right = plasmid[3765:]
     |      >>> new = right + left
     |      >>> print("%s %i" % (new.id, len(new)))
     |      pBAD30 4923
     |      >>> str(new.seq) == str(right.seq + left.seq)
     |      True
     |      >>> len(new.features) == len(left.features) + len(right.features)
     |      True
     |      
     |      When we add the left and right SeqRecord objects, their annotation
     |      is all consistent, so it is all conserved in the new SeqRecord:
     |      
     |      >>> new.id == left.id == right.id == plasmid.id
     |      True
     |      >>> new.name == left.name == right.name == plasmid.name
     |      True
     |      >>> new.description == plasmid.description
     |      True
     |      >>> new.annotations == left.annotations == right.annotations
     |      True
     |      >>> new.letter_annotations == plasmid.letter_annotations
     |      True
     |      >>> new.dbxrefs == left.dbxrefs == right.dbxrefs
     |      True
     |      
     |      However, we should point out that when we sliced the SeqRecord,
     |      any annotations dictionary or dbxrefs list entries were lost.
     |      You can explicitly copy them like this:
     |      
     |      >>> new.annotations = plasmid.annotations.copy()
     |      >>> new.dbxrefs = plasmid.dbxrefs[:]
     |  
     |  __bool__(self)
     |      Boolean value of an instance of this class (True).
     |      
     |      This behaviour is for backwards compatibility, since until the
     |      __len__ method was added, a SeqRecord always evaluated as True.
     |      
     |      Note that in comparison, a Seq object will evaluate to False if it
     |      has a zero length sequence.
     |      
     |      WARNING: The SeqRecord may in future evaluate to False when its
     |      sequence is of zero length (in order to better match the Seq
     |      object behaviour)!
     |  
     |  __bytes__(self)
     |  
     |  __contains__(self, char)
     |      Implement the 'in' keyword, searches the sequence.
     |      
     |      e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> "GAATTC" in record
     |      False
     |      >>> "AAA" in record
     |      True
     |      
     |      This essentially acts as a proxy for using "in" on the sequence:
     |      
     |      >>> "GAATTC" in record.seq
     |      False
     |      >>> "AAA" in record.seq
     |      True
     |      
     |      Note that you can also use Seq objects as the query,
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> Seq("AAA") in record
     |      True
     |      
     |      See also the Seq object's __contains__ method.
     |  
     |  __eq__(self, other)
     |      Define the equal-to operand (not implemented).
     |  
     |  __format__(self, format_spec)
     |      Return the record as a string in the specified file format.
     |      
     |      This method supports the Python format() function and f-strings.
     |      The format_spec should be a lower case string supported by
     |      Bio.SeqIO as a text output file format. Requesting a binary file
     |      format raises a ValueError. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      ...
     |      >>> format(record, "fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(f"Here is {record.id} in FASTA format:\n{record:fasta}")
     |      Here is YP_025292.1 in FASTA format:
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      See also the SeqRecord's format() method.
     |  
     |  __ge__(self, other)
     |      Define the greater-than-or-equal-to operand (not implemented).
     |  
     |  __getitem__(self, index)
     |      Return a sub-sequence or an individual letter.
     |      
     |      Slicing, e.g. my_record[5:10], returns a new SeqRecord for
     |      that sub-sequence with some annotation preserved as follows:
     |      
     |      * The name, id and description are kept as-is.
     |      * Any per-letter-annotations are sliced to match the requested
     |        sub-sequence.
     |      * Unless a stride is used, all those features which fall fully
     |        within the subsequence are included (with their locations
     |        adjusted accordingly). If you want to preserve any truncated
     |        features (e.g. GenBank/EMBL source features), you must
     |        explicitly add them to the new SeqRecord yourself.
     |      * With the exception of any molecule type, the annotations
     |        dictionary and the dbxrefs list are not used for the new
     |        SeqRecord, as in general they may not apply to the
     |        subsequence. If you want to preserve them, you must explicitly
     |        copy them to the new SeqRecord yourself.
     |      
     |      Using an integer index, e.g. my_record[5] is shorthand for
     |      extracting that letter from the sequence, my_record.seq[5].
     |      
     |      For example, consider this short protein and its secondary
     |      structure as encoded by the PDB (e.g. H for alpha helices),
     |      plus a simple feature for its histidine self phosphorylation
     |      site:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> from Bio.SeqFeature import SeqFeature, SimpleLocation
     |      >>> rec = SeqRecord(Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLAT"
     |      ...                     "EMMSEQDGYLAESINKDIEECNAIIEQFIDYLR"),
     |      ...                 id="1JOY", name="EnvZ",
     |      ...                 description="Homodimeric domain of EnvZ from E. coli")
     |      >>> rec.letter_annotations["secondary_structure"] = "  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  "
     |      >>> rec.features.append(SeqFeature(SimpleLocation(20, 21),
     |      ...                     type = "Site"))
     |      
     |      Now let's have a quick look at the full record,
     |      
     |      >>> print(rec)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      >>> rec.letter_annotations["secondary_structure"]
     |      '  S  SSSSSSHHHHHTTTHHHHHHHHHHHHHHHHHHHHHHTHHHHHHHHHHHHHHHHHHHHHTT  '
     |      >>> print(rec.features[0].location)
     |      [20:21]
     |      
     |      Now let's take a sub sequence, here chosen as the first (fractured)
     |      alpha helix which includes the histidine phosphorylation site:
     |      
     |      >>> sub = rec[11:41]
     |      >>> print(sub)
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('RTLLMAGVSHDLRTPLTRIRLATEMMSEQD')
     |      >>> sub.letter_annotations["secondary_structure"]
     |      'HHHHHTTTHHHHHHHHHHHHHHHHHHHHHH'
     |      >>> print(sub.features[0].location)
     |      [9:10]
     |      
     |      You can also of course omit the start or end values, for
     |      example to get the first ten letters only:
     |      
     |      >>> print(rec[:10])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLAD')
     |      
     |      Or for the last ten letters:
     |      
     |      >>> print(rec[-10:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 0
     |      Per letter annotation for: secondary_structure
     |      Seq('IIEQFIDYLR')
     |      
     |      If you omit both, then you get a copy of the original record (although
     |      lacking the annotations and dbxrefs):
     |      
     |      >>> print(rec[:])
     |      ID: 1JOY
     |      Name: EnvZ
     |      Description: Homodimeric domain of EnvZ from E. coli
     |      Number of features: 1
     |      Per letter annotation for: secondary_structure
     |      Seq('MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE...YLR')
     |      
     |      Finally, indexing with a simple integer is shorthand for pulling out
     |      that letter from the sequence directly:
     |      
     |      >>> rec[5]
     |      'K'
     |      >>> rec.seq[5]
     |      'K'
     |  
     |  __gt__(self, other)
     |      Define the greater-than operand (not implemented).
     |  
     |  __init__(self, seq, id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=None, features=None, annotations=None, letter_annotations=None)
     |      Create a SeqRecord.
     |      
     |      Arguments:
     |       - seq         - Sequence, required (Seq or MutableSeq)
     |       - id          - Sequence identifier, recommended (string)
     |       - name        - Sequence name, optional (string)
     |       - description - Sequence description, optional (string)
     |       - dbxrefs     - Database cross references, optional (list of strings)
     |       - features    - Any (sub)features, optional (list of SeqFeature objects)
     |       - annotations - Dictionary of annotations for the whole sequence
     |       - letter_annotations - Dictionary of per-letter-annotations, values
     |         should be strings, list or tuples of the same length as the full
     |         sequence.
     |      
     |      You will typically use Bio.SeqIO to read in sequences from files as
     |      SeqRecord objects.  However, you may want to create your own SeqRecord
     |      objects directly.
     |      
     |      Note that while an id is optional, we strongly recommend you supply a
     |      unique id string for each record.  This is especially important
     |      if you wish to write your sequences to a file.
     |      
     |      You can create a 'blank' SeqRecord object, and then populate the
     |      attributes later.
     |  
     |  __iter__(self)
     |      Iterate over the letters in the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a protein FASTA file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/loveliesbleeding.pro", "fasta")
     |      >>> for amino in record:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      This is just a shortcut for iterating over the sequence directly:
     |      
     |      >>> for amino in record.seq:
     |      ...     print(amino)
     |      ...     if amino == "L": break
     |      X
     |      A
     |      G
     |      L
     |      >>> print(record.seq[3])
     |      L
     |      
     |      Note that this does not facilitate iteration together with any
     |      per-letter-annotation.  However, you can achieve that using the
     |      python zip function on the record (or its sequence) and the relevant
     |      per-letter-annotation:
     |      
     |      >>> from Bio import SeqIO
     |      >>> rec = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(rec.letter_annotations))
     |      ['solexa_quality']
     |      >>> for nuc, qual in zip(rec, rec.letter_annotations["solexa_quality"]):
     |      ...     if qual > 35:
     |      ...         print("%s %i" % (nuc, qual))
     |      A 40
     |      C 39
     |      G 38
     |      T 37
     |      A 36
     |      
     |      You may agree that using zip(rec.seq, ...) is more explicit than using
     |      zip(rec, ...) as shown above.
     |  
     |  __le__(self, other)
     |      Define the less-than-or-equal-to operand (not implemented).
     |  
     |  __len__(self)
     |      Return the length of the sequence.
     |      
     |      For example, using Bio.SeqIO to read in a FASTA nucleotide file:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> len(record)
     |      309
     |      >>> len(record.seq)
     |      309
     |  
     |  __lt__(self, other)
     |      Define the less-than operand (not implemented).
     |  
     |  __ne__(self, other)
     |      Define the not-equal-to operand (not implemented).
     |  
     |  __radd__(self, other)
     |      Add another sequence or string to this sequence (from the left).
     |      
     |      This method handles adding a Seq object (or similar, e.g. MutableSeq)
     |      or a plain Python string (on the left) to a SeqRecord (on the right).
     |      See the __add__ method for more details, but for example:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      
     |      >>> new = "ACT" + record
     |      >>> print("%s %s" % (new.id, new.seq))
     |      slxa_0001_1_0001_01 ACTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(new.letter_annotations))
     |      []
     |  
     |  __repr__(self)
     |      Return a concise summary of the record for debugging (string).
     |      
     |      The python built in function repr works by calling the object's __repr__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(Seq("MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKAT"
     |      ...                     "GEMKEQTEWHRVVLFGKLAEVASEYLRKGSQVYIEGQLRTRKWTDQ"
     |      ...                     "SGQDRYTTEVVVNVGGTMQMLGGRQGGGAPAGGNIGGGQPQGGWGQ"
     |      ...                     "PQQPQGGNQFSGGAQSRPQQSAPAAPSNEPPMDFDDDIPF"),
     |      ...                 id="NP_418483.1", name="b4059",
     |      ...                 description="ssDNA-binding protein",
     |      ...                 dbxrefs=["ASAP:13298", "GI:16131885", "GeneID:948570"])
     |      >>> print(repr(rec))
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      At the python prompt you can also use this shorthand:
     |      
     |      >>> rec
     |      SeqRecord(seq=Seq('MASRGVNKVILVGNLGQDPEVRYMPNGGAVANITLATSESWRDKATGEMKEQTE...IPF'), id='NP_418483.1', name='b4059', description='ssDNA-binding protein', dbxrefs=['ASAP:13298', 'GI:16131885', 'GeneID:948570'])
     |      
     |      Note that long sequences are shown truncated. Also note that any
     |      annotations, letter_annotations and features are not shown (as they
     |      would lead to a very long string).
     |  
     |  __str__(self)
     |      Return a human readable summary of the record and its annotation (string).
     |      
     |      The python built in function str works by calling the object's __str__
     |      method.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein, small")
     |      >>> print(str(record))
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      In this example you don't actually need to call str explicitly, as the
     |      print command does this automatically:
     |      
     |      >>> print(record)
     |      ID: YP_025292.1
     |      Name: HokC
     |      Description: toxic membrane protein, small
     |      Number of features: 0
     |      Seq('MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF')
     |      
     |      Note that long sequences are shown truncated.
     |  
     |  count(self, sub, start=None, end=None)
     |      Return the number of non-overlapping occurrences of sub in seq[start:end].
     |      
     |      Optional arguments start and end are interpreted as in slice notation.
     |      This method behaves as the count method of Python strings.
     |  
     |  format(self, format)
     |      Return the record as a string in the specified file format.
     |      
     |      The format should be a lower case string supported as an output
     |      format by Bio.SeqIO, which is used to turn the SeqRecord into a
     |      string.  e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF"),
     |      ...                    id="YP_025292.1", name="HokC",
     |      ...                    description="toxic membrane protein")
     |      >>> record.format("fasta")
     |      '>YP_025292.1 toxic membrane protein\nMKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF\n'
     |      >>> print(record.format("fasta"))
     |      >YP_025292.1 toxic membrane protein
     |      MKQHKAMIVALIVICITAVVAALVTRKDLCEVHIRTGQTEVAVF
     |      <BLANKLINE>
     |      
     |      The Python print function automatically appends a new line, meaning
     |      in this example a blank line is shown.  If you look at the string
     |      representation you can see there is a trailing new line (shown as
     |      slash n) which is important when writing to a file or if
     |      concatenating multiple sequence strings together.
     |      
     |      Note that this method will NOT work on every possible file format
     |      supported by Bio.SeqIO (e.g. some are for multiple sequences only,
     |      and binary formats are not supported).
     |  
     |  islower(self)
     |      Return True if all ASCII characters in the record's sequence are lowercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  isupper(self)
     |      Return True if all ASCII characters in the record's sequence are uppercase.
     |      
     |      If there are no cased characters, the method returns False.
     |  
     |  lower(self)
     |      Return a copy of the record with a lower case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Fasta/aster.pro", "fasta")
     |      >>> print(record.format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      GGHVNPAVTFGAFVGGNITLLRGIVYIIAQLLGSTVACLLLKFVTNDMAVGVFSLSAGVG
     |      VTNALVFEIVMTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI
     |      <BLANKLINE>
     |      >>> print(record.lower().format("fasta"))
     |      >gi|3298468|dbj|BAA31520.1| SAMIPF
     |      gghvnpavtfgafvggnitllrgivyiiaqllgstvaclllkfvtndmavgvfslsagvg
     |      vtnalvfeivmtfglvytvyataidpkkgslgtiapiaigfivgani
     |      <BLANKLINE>
     |      
     |      To take a more annotation rich example,
     |      
     |      >>> from Bio import SeqIO
     |      >>> old = SeqIO.read("EMBL/TRBG361.embl", "embl")
     |      >>> len(old.features)
     |      3
     |      >>> new = old.lower()
     |      >>> len(old.features) == len(new.features)
     |      True
     |      >>> old.annotations["organism"] == new.annotations["organism"]
     |      True
     |      >>> old.dbxrefs == new.dbxrefs
     |      True
     |  
     |  reverse_complement(self, id=False, name=False, description=False, features=True, annotations=False, letter_annotations=True, dbxrefs=False)
     |      Return new SeqRecord with reverse complement sequence.
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the reversed sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or True to keep that of the parent, or False to
     |      omit them. The default is to keep the original features (with the
     |      strand and locations adjusted).
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent,
     |      or False to omit them. The default is to keep the original
     |      annotations (with the letter annotations reversed).
     |      
     |      To show what happens to the pre-letter annotations, consider an
     |      example Solexa variant FASTQ file with a single entry, which we'll
     |      read in as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Now take the reverse complement, here we explicitly give a new
     |      identifier (the old identifier with a suffix):
     |      
     |      >>> rc_record = record.reverse_complement(id=record.id + "_rc")
     |      >>> print("%s %s" % (rc_record.id, rc_record.seq))
     |      slxa_0001_1_0001_01_rc NNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
     |      
     |      Notice that the per-letter-annotations have also been reversed,
     |      although this may not be appropriate for all cases.
     |      
     |      >>> print(rc_record.letter_annotations["solexa_quality"])
     |      [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40]
     |      
     |      Now for the features, we need a different example. Parsing a GenBank
     |      file is probably the easiest way to get an nice example with features
     |      in it...
     |      
     |      >>> from Bio import SeqIO
     |      >>> with open("GenBank/pBAD30.gb") as handle:
     |      ...     plasmid = SeqIO.read(handle, "gb")
     |      >>> print("%s %i" % (plasmid.id, len(plasmid)))
     |      pBAD30 4923
     |      >>> plasmid.seq
     |      Seq('GCTAGCGGAGTGTATACTGGCTTACTATGTTGGCACTGATGAGGGTGTCAGTGA...ATG')
     |      >>> len(plasmid.features)
     |      13
     |      
     |      Now, let's take the reverse complement of this whole plasmid:
     |      
     |      >>> rc_plasmid = plasmid.reverse_complement(id=plasmid.id+"_rc")
     |      >>> print("%s %i" % (rc_plasmid.id, len(rc_plasmid)))
     |      pBAD30_rc 4923
     |      >>> rc_plasmid.seq
     |      Seq('CATGGGCAAATATTATACGCAAGGCGACAAGGTGCTGATGCCGCTGGCGATTCA...AGC')
     |      >>> len(rc_plasmid.features)
     |      13
     |      
     |      Let's compare the first CDS feature - it has gone from being the
     |      second feature (index 1) to the second last feature (index -2), its
     |      strand has changed, and the location switched round.
     |      
     |      >>> print(plasmid.features[1])
     |      type: CDS
     |      location: [1081:1960](-)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      >>> print(rc_plasmid.features[-2])
     |      type: CDS
     |      location: [2963:3842](+)
     |      qualifiers:
     |          Key: label, Value: ['araC']
     |          Key: note, Value: ['araC regulator of the arabinose BAD promoter']
     |          Key: vntifkey, Value: ['4']
     |      <BLANKLINE>
     |      
     |      You can check this new location, based on the length of the plasmid:
     |      
     |      >>> len(plasmid) - 1081
     |      3842
     |      >>> len(plasmid) - 1960
     |      2963
     |      
     |      Note that if the SeqFeature annotation includes any strand specific
     |      information (e.g. base changes for a SNP), this information is not
     |      amended, and would need correction after the reverse complement.
     |      
     |      Note trying to reverse complement a protein SeqRecord raises an
     |      exception:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> protein_rec = SeqRecord(Seq("MAIVMGR"), id="Test",
     |      ...                         annotations={"molecule_type": "protein"})
     |      >>> protein_rec.reverse_complement()
     |      Traceback (most recent call last):
     |         ...
     |      ValueError: Proteins do not have complements!
     |      
     |      If you have RNA without any U bases, it must be annotated as RNA
     |      otherwise it will be treated as DNA by default with A mapped to T:
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rna1 = SeqRecord(Seq("ACG"), id="Test")
     |      >>> rna2 = SeqRecord(Seq("ACG"), id="Test", annotations={"molecule_type": "RNA"})
     |      >>> print(rna1.reverse_complement(id="RC", description="unk").format("fasta"))
     |      >RC unk
     |      CGT
     |      <BLANKLINE>
     |      >>> print(rna2.reverse_complement(id="RC", description="RNA").format("fasta"))
     |      >RC RNA
     |      CGU
     |      <BLANKLINE>
     |      
     |      Also note you can reverse complement a SeqRecord using a MutableSeq:
     |      
     |      >>> from Bio.Seq import MutableSeq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> rec = SeqRecord(MutableSeq("ACGT"), id="Test")
     |      >>> rec.seq[0] = "T"
     |      >>> print("%s %s" % (rec.id, rec.seq))
     |      Test TCGT
     |      >>> rc = rec.reverse_complement(id=True)
     |      >>> print("%s %s" % (rc.id, rc.seq))
     |      Test ACGA
     |  
     |  translate(self, table='Standard', stop_symbol='*', to_stop=False, cds=False, gap=None, id=False, name=False, description=False, features=False, annotations=False, letter_annotations=False, dbxrefs=False)
     |      Return new SeqRecord with translated sequence.
     |      
     |      This calls the record's .seq.translate() method (which describes
     |      the translation related arguments, like table for the genetic code),
     |      
     |      By default the new record does NOT preserve the sequence identifier,
     |      name, description, general annotation or database cross-references -
     |      these are unlikely to apply to the translated sequence.
     |      
     |      You can specify the returned record's id, name and description as
     |      strings, or True to keep that of the parent, or False for a default.
     |      
     |      You can specify the returned record's features with a list of
     |      SeqFeature objects, or False (default) to omit them.
     |      
     |      You can also specify both the returned record's annotations and
     |      letter_annotations as dictionaries, True to keep that of the parent
     |      (annotations only), or False (default) to omit them.
     |      
     |      e.g. Loading a FASTA gene and translating it,
     |      
     |      >>> from Bio import SeqIO
     |      >>> gene_record = SeqIO.read("Fasta/sweetpea.nu", "fasta")
     |      >>> print(gene_record.format("fasta"))
     |      >gi|3176602|gb|U78617.1|LOU78617 Lathyrus odoratus phytochrome A (PHYA) gene, partial cds
     |      CAGGCTGCGCGGTTTCTATTTATGAAGAACAAGGTCCGTATGATAGTTGATTGTCATGCA
     |      AAACATGTGAAGGTTCTTCAAGACGAAAAACTCCCATTTGATTTGACTCTGTGCGGTTCG
     |      ACCTTAAGAGCTCCACATAGTTGCCATTTGCAGTACATGGCTAACATGGATTCAATTGCT
     |      TCATTGGTTATGGCAGTGGTCGTCAATGACAGCGATGAAGATGGAGATAGCCGTGACGCA
     |      GTTCTACCACAAAAGAAAAAGAGACTTTGGGGTTTGGTAGTTTGTCATAACACTACTCCG
     |      AGGTTTGTT
     |      <BLANKLINE>
     |      
     |      And now translating the record, specifying the new ID and description:
     |      
     |      >>> protein_record = gene_record.translate(table=11,
     |      ...                                        id="phya",
     |      ...                                        description="translation")
     |      >>> print(protein_record.format("fasta"))
     |      >phya translation
     |      QAARFLFMKNKVRMIVDCHAKHVKVLQDEKLPFDLTLCGSTLRAPHSCHLQYMANMDSIA
     |      SLVMAVVVNDSDEDGDSRDAVLPQKKKRLWGLVVCHNTTPRFV
     |      <BLANKLINE>
     |  
     |  upper(self)
     |      Return a copy of the record with an upper case sequence.
     |      
     |      All the annotation is preserved unchanged. e.g.
     |      
     |      >>> from Bio.Seq import Seq
     |      >>> from Bio.SeqRecord import SeqRecord
     |      >>> record = SeqRecord(Seq("acgtACGT"), id="Test",
     |      ...                    description = "Made up for this example")
     |      >>> record.letter_annotations["phred_quality"] = [1, 2, 3, 4, 5, 6, 7, 8]
     |      >>> print(record.upper().format("fastq"))
     |      @Test Made up for this example
     |      ACGTACGT
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |      
     |      Naturally, there is a matching lower method:
     |      
     |      >>> print(record.lower().format("fastq"))
     |      @Test Made up for this example
     |      acgtacgt
     |      +
     |      "#$%&'()
     |      <BLANKLINE>
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors defined here:
     |  
     |  __dict__
     |      dictionary for instance variables (if defined)
     |  
     |  __weakref__
     |      list of weak references to the object (if defined)
     |  
     |  letter_annotations
     |      Dictionary of per-letter-annotation for the sequence.
     |      
     |      For example, this can hold quality scores used in FASTQ or QUAL files.
     |      Consider this example using Bio.SeqIO to read in an example Solexa
     |      variant FASTQ file as a SeqRecord:
     |      
     |      >>> from Bio import SeqIO
     |      >>> record = SeqIO.read("Quality/solexa_faked.fastq", "fastq-solexa")
     |      >>> print("%s %s" % (record.id, record.seq))
     |      slxa_0001_1_0001_01 ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTNNNNNN
     |      >>> print(list(record.letter_annotations))
     |      ['solexa_quality']
     |      >>> print(record.letter_annotations["solexa_quality"])
     |      [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      The letter_annotations get sliced automatically if you slice the
     |      parent SeqRecord, for example taking the last ten bases:
     |      
     |      >>> sub_record = record[-10:]
     |      >>> print("%s %s" % (sub_record.id, sub_record.seq))
     |      slxa_0001_1_0001_01 ACGTNNNNNN
     |      >>> print(sub_record.letter_annotations["solexa_quality"])
     |      [4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
     |      
     |      Any python sequence (i.e. list, tuple or string) can be recorded in
     |      the SeqRecord's letter_annotations dictionary as long as the length
     |      matches that of the SeqRecord's sequence.  e.g.
     |      
     |      >>> len(sub_record.letter_annotations)
     |      1
     |      >>> sub_record.letter_annotations["dummy"] = "abcdefghij"
     |      >>> len(sub_record.letter_annotations)
     |      2
     |      
     |      You can delete entries from the letter_annotations dictionary as usual:
     |      
     |      >>> del sub_record.letter_annotations["solexa_quality"]
     |      >>> sub_record.letter_annotations
     |      {'dummy': 'abcdefghij'}
     |      
     |      You can completely clear the dictionary easily as follows:
     |      
     |      >>> sub_record.letter_annotations = {}
     |      >>> sub_record.letter_annotations
     |      {}
     |      
     |      Note that if replacing the record's sequence with a sequence of a
     |      different length you must first clear the letter_annotations dict.
     |  
     |  seq
     |      The sequence itself, as a Seq or MutableSeq object.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __hash__ = None
    



```python
from Bio.Seq import Seq
```


```python
simple_seq = Seq("GATC")
```


```python
simple_seq_r = SeqRecord(simple_seq)
```


```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
simple_seq_r.id = "AC12345"
```


```python
simple_seq_r.description = "Made up sequence for the VDB Computational Biology Class" 
```


```python
print(simple_seq_r.description)
```

    Made up sequence for the VDB Computational Biology Class



```python
simple_seq_r.seq
```




    Seq('GATC')




```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.annotations["evidence"] ="None. This is just an example"
```


```python
print(simple_seq_r.annotations["evidence"])
```

    None. This is just an example



```python
simple_seq_r
```




    SeqRecord(seq=Seq('GATC'), id='AC12345', name='<unknown name>', description='Made up sequence for the VDB Computational Biology Class', dbxrefs=[])




```python
simple_seq_r.letter_annotations["phred_equality"] = [40, 40, 38, 30]
```


```python
print(simple_seq_r.letter_annotations)
```

    {'phred_equality': [40, 40, 38, 30]}



```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna
```


```python
from Bio import SeqIO
```


```python
#We can import fasta files within jupyter lab notebook and read the file 
record = SeqIO.read("NC_005816.fna.txt", "fasta")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='gi|45478711|ref|NC_005816.1|', name='gi|45478711|ref|NC_005816.1|', description='gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'gi|45478711|ref|NC_005816.1|'




```python
record.description
```




    'gi|45478711|ref|NC_005816.1| Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.dbxrefs
```




    []




```python
record.annotations
```




    {}




```python
record.features
```




    []




```python
#https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
record.seq
```




    Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG')




```python
record.id
```




    'NC_005816.1'




```python
record.name
```




    'NC_005816'




```python
# we can get a fasta file from offline and read it using python
record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
record.letter_annotations
```




    {}




```python
len (record.annotations)
```




    13




```python
record.annotations["source"]
```




    'Yersinia pestis biovar Microtus str. 91001'




```python
record.dbxrefs
```




    ['Project:58037']




```python
len(record.features)
```




    41




```python
from Bio import SeqFeature
```


```python
start_pos = SeqFeature.AfterPosition(1)
```


```python
end_pos = SeqFeature.BetweenPosition(9, left = 8, right = 9)
```


```python
my_location = SeqFeature.SimpleLocation(start_pos, end_pos)
```


```python
print(my_location)
```

    [>1:(8^9)]



```python
my_location.start
```




    AfterPosition(1)




```python
my_location.end
```




    BetweenPosition(9, left=8, right=9)




```python
int(my_location.end)
```




    9




```python
int(my_location.start)
```




    1




```python
exact_location = SeqFeature.SimpleLocation(1,9)
```


```python
print(exact_location)
```

    [1:9]



```python
exact_location.start
```




    ExactPosition(1)




```python
from Bio.SeqRecord import SeqRecord
```


```python
record = SeqRecord(Seq("MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEIIAVTFRGPSETHLDSMVGQALFGD"     
                      "GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK" 
                      "NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM"
                      "SSAC"), 
id = "gi|14150838|gb|AAK54648.1|AF376133_1", 
description = "chalcone synthase [Cucumis sativus]",)
```


```python
#we can imput a sequence into python and it will read it as a fasta file 
print(record.format("fasta"))
```

    >gi|14150838|gb|AAK54648.1|AF376133_1 chalcone synthase [Cucumis sativus]
    MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEIIAVTFRGPSETHLDSMVGQALFGD
    GAGAVIVGSDPDLSVERPLYELVWTGATLLPDSEGAIDGHLREVGLTFHLLKDVPGLISK
    NIEKSLKEAFTPLGISDWNSTFWIAHPGGPAILDQVEAKLGLKEEKMRATREVLSEYGNM
    SSAC
    



```python
print(record)
```

    ID: gi|14150838|gb|AAK54648.1|AF376133_1
    Name: <unknown name>
    Description: chalcone synthase [Cucumis sativus]
    Number of features: 0
    Seq('MMYQQGCFAGGTVLRLAKDLAENNRGARVLVVCSEIIAVTFRGPSETHLDSMVG...SAC')



```python
from Bio import SeqIO
```


```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
len(record)
```




    9609




```python
print(record.features[20])
```

    type: gene
    location: [4342:4780](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(record.features[21])
```

    type: CDS
    location: [4342:4780](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record = record[4300:4800]
```


```python
len(sub_record)
```




    500




```python
len(sub_record.features)
```




    2




```python
sub_record.features[1]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='CDS', qualifiers=...)




```python
sub_record.features[0]
```




    SeqFeature(SimpleLocation(ExactPosition(42), ExactPosition(480), strand=1), type='gene', qualifiers=...)




```python
print(sub_record.features[0])
```

    type: gene
    location: [42:480](+)
    qualifiers:
        Key: db_xref, Value: ['GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
    



```python
print(sub_record.features[1])
```

    type: CDS
    location: [42:480](+)
    qualifiers:
        Key: codon_start, Value: ['1']
        Key: db_xref, Value: ['GI:45478716', 'GeneID:2767712']
        Key: gene, Value: ['pim']
        Key: locus_tag, Value: ['YP_pPCP05']
        Key: note, Value: ['similar to many previously sequenced pesticin immunity protein entries of Yersinia pestis plasmid pPCP, e.g. gi| 16082683|,ref|NP_395230.1| (NC_003132) , gi|1200166|emb|CAA90861.1| (Z54145 ) , gi|1488655| emb|CAA63439.1| (X92856) , gi|2996219|gb|AAC62543.1| (AF053945) , and gi|5763814|emb|CAB531 67.1| (AL109969)']
        Key: product, Value: ['pesticin immunity protein']
        Key: protein_id, Value: ['NP_995571.1']
        Key: transl_table, Value: ['11']
        Key: translation, Value: ['MGGGMISKLFCLALIFLSSSGLAEKNTYTAKDILQNLELNTFGNSLSHGIYGKQTTFKQTEFTNIKSNTKKHIALINKDNSWMISLKILGIKRDEYTVCFEDFSLIRPPTYVAIHPLLIKKVKSGNFIVVKEIKKSIPGCTVYYH']
    



```python
sub_record.annotations
```




    {'molecule_type': 'DNA'}




```python
sub_record.dbxrefs
```




    []




```python
sub_record.annotations["topology"] = "linear"
```


```python
sub_record.annotations
```




    {'molecule_type': 'DNA', 'topology': 'linear'}




```python
sub_record.id
```




    'NC_005816.1'




```python
sub_record.name
```




    'NC_005816'




```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence'




```python
sub_record.description = "Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence"
```


```python
sub_record.description
```




    'Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial sequence'




```python
print(sub_record.format("genbank")[:200] + "...")
```

    LOCUS       NC_005816                500 bp    DNA     linear   UNK 01-JAN-1980
    DEFINITION  Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, partial
                sequence.
    ACCESSION   NC_00581...



```python
record = SeqIO.read("NC_005816.gb.txt", "genbank")
```


```python
record 
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
len(record)
```




    9609




```python
len(record.features)
```




    41




```python
record.dbxrefs
```




    ['Project:58037']




```python
record.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
#we can move around in the sequence to find the section we are looking for
shifted = record[2000:] + record[:2000]
```


```python
shifted
```




    SeqRecord(seq=Seq('GATACGCAGTCATATTTTTTACACAATTCTCTAATCCCGACAAGGTCGTAGGTC...GGA'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=[])




```python
len(shifted)
```




    9609




```python
len(shifted.features)
```




    40




```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type'])




```python
shifted.dbxrefs
```




    []




```python
shifted.dbxrefs = record.dbxrefs[:]
```


```python
shifted.dbxrefs
```




    ['Project:58037']




```python
shifted.annotations = record.annotations.copy()
```


```python
shifted.annotations.keys()
```




    dict_keys(['molecule_type', 'topology', 'data_file_division', 'date', 'accessions', 'sequence_version', 'gi', 'keywords', 'source', 'organism', 'taxonomy', 'references', 'comment'])




```python
record
```




    SeqRecord(seq=Seq('TGTAACGAACGGTGCAATAGTGATCCACACCCAACGCCTGAAATCAGATCCAGG...CTG'), id='NC_005816.1', name='NC_005816', description='Yersinia pestis biovar Microtus str. 91001 plasmid pPCP1, complete sequence', dbxrefs=['Project:58037'])




```python
print("%s %i %i %i %i" % (record.id, len(record), len(record.features), len(record.dbxrefs), len(record.annotations)))
```

    NC_005816.1 9609 41 1 13



```python
rc = record.reverse_complement(id ="Testing")
```


```python
rc
```




    SeqRecord(seq=Seq('CAGGGGTCGGGGTACGCATTCCCTCATGCGTCAATATTATCTGGCATTGCGATG...ACA'), id='Testing', name='<unknown name>', description='<unknown description>', dbxrefs=[])




```python
print("%s %i %i %i %i" % (rc.id, len(rc), len(rc.features), len(rc.dbxrefs), len(rc.annotations)))
```

    Testing 9609 41 0 0



```python

```



