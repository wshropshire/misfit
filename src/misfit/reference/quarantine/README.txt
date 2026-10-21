Reference alleles removed from the active panels, kept for review.

k_pneumoniae_ref_ompK35_TRUNCATED_ALLELE.fasta
  ompK35_Kp1_11846853, 435 nt / 144 aa, versus 1080 nt / 359 aa for the
  other two alleles in ref_ompK35.fasta. This is a truncated (mutant)
  OmpK35, not a wild-type reference. While it was in the panel it won hit
  selection for 210 of 236 K. pneumoniae isolates and made their intact,
  full-length ompK35 look frameshifted. Restore only if you intend it as a
  known-mutant allele AND hit selection accounts for allele length.

Paralogs and duplicate isoforms removed 2026-08-12
RefSeq gives paralogous genes the same gene name, so one reference file can end up
holding two different genes. Both then align to their own locus in the assembly and
hit selection picks whichever yields more matching bases -- which was reporting the
wrong locus. Measured against the E. coli ortholog, each pair below is separated by
>=18 identity points, so the assignment is unambiguous:
  mrdA  k_pneumoniae/citrobacter/k_oxytoca/k_aerogenes: 69-71% paralog removed,
        93-96% true MrdA kept. Affected 216/236 Kp, 24/25 Cfr, 13/14 Kae isolates.
  ftsI  enterobacter/citrobacter: 581 aa 72-74% paralog removed, 588 aa 94-96% kept.
  arnF  e_coli: b2270 YfbK removed -- it entered via an incorrect alias, not a paralog.
  mrcB  e_coli: 799 aa PBP-1Bgamma isoform of the same locus removed; 844 aa kept.

Misnamed ompC paralogs removed 2026-08-12
In Enterobacter, K. oxytoca and K. aerogenes the CDS that RefSeq *names* ompC are not
the OmpC ortholog: they sit elsewhere in the genome and score only 72-78% against
E. coli OmpC. Cohort isolates split across them (Enterobacter 37/30/17), so 'ompC' was
not a consistent locus. The true ortholog is an UNNAMED CDS whose product is 'porin
OmpC', sitting in the same neighbourhood as E. coli ompC and K. pneumoniae ompK36
(apbE / ada / rcsD-rcsB / nrdB / glpQ) and scoring 86-92%. Two concordant lines of
evidence, 13-20 identity points clear of the alternatives.

Klebsiella porin naming corrected 2026-08-12
An earlier pass filed the ada/rcsDB-syntenic porin as ompC for K. oxytoca and
K. aerogenes. In Klebsiella that locus is ompK36, and these genomes carry a
SEPARATE gene RefSeq names ompC. Corrected against the K. pneumoniae porin set:
  k_oxytoca   ompK36=WP_016809062.1 (92.2%), ompK35=WP_004100735.1 (93.5%, new),
              ompK37=WP_004101605.1 (88.6%); no ompC ortholog found.
  k_aerogenes ompK36=WP_015706141.1 (91.5%), ompC=WP_015705989.1 (89.0%),
              ompK35=WP_015367453.1 (95.5%), ompK37=WP_015705623.1 (88.3%).
The curated ompK36_Kaerogenes_QXB08681.1 is retired: WP_015706141.1 is the same
protein with a RefSeq accession.
