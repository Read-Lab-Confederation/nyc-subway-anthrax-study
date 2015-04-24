---
title: "Summary of Results"
order: 3
---

Our control experiments showed that direct mapping of sequences reads from the metagenome was a sensitive way to detect plasmids,  Plasmid reads were present when as little as 0.01x *B. anthracis* genome equivalents were added.  However, some *B. cereus* reads mapped to pXO2 because the strain we used probably contained a closely related mobile element.  What differentiated *B. anthracis* from the control was that coverage was evenly distributed across the plasmids.  We found that the NYC samples had evidence of quite deep coverage of *B. cereus*-like chromosome and plasmids, which was not an unexpected result.  The pattern of coverage of the plasmids qualitatively resembled *B. cereus* rather than *B. anthracis*. There were some reads in sample SRR1748708 that mapped across the part of the pXO1 plasmid containing the lethal toxin gene cluster but none of the reads ultimately mapped to *B. anthracis* using BLAST [Table 4](#table-4). 

We showed that Kraken was a sensitive tool for detection of *B. anthracis* down to at least 0.01X genome coverage.  However, the *B. cereus* control also gave us false positive *B. anthracis* reads, albeit at a lower incidence than the true positive.  On the NYC subway samples we obtained a reasonably large number of *B. anthracis* specific calls but this could be explained by the fact that the sample contained up to 50x chromosomal coverage of *B. cereus* group organisms.

Our novel assay of 1,793 31-mers specific to the  *B. anthracis* core genome was sensitive at 0.01X *B. anthracis* genome coverage and did not give false positives on the *B. cereus* control at up to 5x coverage.  Surprisingly, we obtained 2 positive hits on on of the NYC subway samples. There are a number of possible reasons for this:

1.  *B. anthracis* is truly present in this sample at about 0.01x genome coverage along with a much greater number of *B. cereus* group organisms.  At this level of coverage we would expect to see some hits to the lethal toxin gene in pXO1, which we did. However, the reads that mapped inside the lethal toxin gene were found by BLASTN to be false positives.  It is possible that the 0.01x genome coverage the *B. anthracis* pXO1 reads stochastically missed the toxin gene and we were seeing false positive hits from the other *B. cereus* group strains present. Alternatively, we detected a *B. anthracis* variant that had lost its plasmids.
2.  The SNPs are present in a *B. cereus* strain, very closely related to *B. anthracis* that has not been sequenced (if its sequence were available in a public database, it would have resulted in the SNPs being filtered out of consideration).
3.  The SNPs were present in a *B. cereus* strain whose ancestor had at some point undergone homologous recombination in this area of the chromosome with an ancestor of *B. anthracis*.
4.  The sequences were introduced through laboratory cross-contamination.
5.  The SNPs were an artifact of random sequencing error.


At the present time, we believe the most likely explanation for the results is option 2.  There is no direct evidence for the lethal toxin.  The B. cereus close to *B. anthracis* lkely represents a small fraction of a mixture of *B. cereus* strains present.

Deeper sequence coverage of this sample preferably combined with culture and PCR-based analysis would be helpful here.
