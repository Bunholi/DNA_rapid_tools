This readme file was generated on 2025-04-01 by Ingrid Bunholi
This repository contains the code and data to reproduce all tables and figures presented in Alvarenga & Bunholi et al. "Rapid DNA/eDNA-based Tools for improved Chondrichthyes monitoring and management"

General Information

Title of Dataset: DNA_rapid_tools
Name: Ingrid Bunholi ORCID: 0000-0001-5489-276X Email: ingridbunholi@gmail.com
Institution: University of Texas at Austin Marine Science Institute (UTMSI); Address: 750 Channel View Dr, Port Aransas, TX 78373
Details
The elements of this project include:
	1.	metadata: this folder includes two raw metadata necessary to run the scripts
	2.	scripts: this folder includes all scripts to generate the figures and summary statistics
	3.	output: this folder includes all outputs from the R-script - include figures and summary statistics

DATA & FILE OVERVIEW

Metadata

The metadata used to generate Figure 1 is “allsets_rapid_tools_0331.csv” and Figure 2 is “primers_specific_unique_cleaned_0331.csv”. These files contain every metric collected from the 70 papers: 
Allsets_rapid_tools_0331 file: paper_doi, citation (first author’s name), source (Wos Search or snowball), title, year, country_first_author_c, country_last_author_c, country_origin_samples_c, nagoya_protocol_c (yes, no, or not specified), primer_developed_c (applied, developed, both), application_primer (eDNA or Tissue DNA), Rapid_DNAtools_c (Multiplex PCR, Droplet Digital PCR, PCR RFLP, etc), on-site_identification_device, RapidDNAtools_develop_c (applied, developed, both), OBS, Rapid_tool_reference 
 
Primers_specific_unique_cleaned_0331 file: paper_doi, citation (first author’s name), purpose (rapid_tools or just primer), rapid_DNAtools, application (eDNA or Tissue DNA), primer_name, primer_type (F or R), primer_sequence, developed_OR_applied (applied, developed or both), primer_citation, test_method (lab experiment or in silico), molecular_marker (e.g., ND2, COI, CytB, etc), GC (GC content in %), annealing_temperature, Shark_Ray_Chimaera (Shark, Ray or Chimaera), Order (e.g., Carcharhiniformes, Lamniformes, etc), Family (e.g., Carcharhinidae, Lamnidae, etc), Genus (e.g., Carcharhinus, Alopias, etc), Species (e.g., Carcharhinus_leucas, Alopias_pelagicus, etc)

Scripts

All the figures were generated using functions from packages described below. Additional graphical features were added in Illustrator.

METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: the data was obtained through a systematic literature review on Web of Science, using targeted keywords to identify studies that implemented rapid DNA/eDNA tools for Chondrichthyes research. To ensure a rigorous and comprehensive collection, we adhered to theRepOrting Standards for Systematic Evidence Syntheses (ROSES; Haddaway et al. 2018) and applied the snowball method, examining reference and citation lists for additional relevant studies (Wohlin 2014).
 
Date of data collection: 2025-03-31

Methods for processing the data:  The final dataset after WoS research comprised 742 studies, which underwent manual screening based on three inclusion criteria: (1) relevance to Chondrichthyes research, (2) development or application of specific primers for taxonomic identification, and (3) development or application of Rapid DNA-based tools for taxonomic identification. Each study was independently evaluated by three reviewers and included if at least two agreed (majority rule). The reliability of this review process was assessed using inter-rater agreement (IRR), a metric for consistency among independent reviewers (Gisev et al. 2013). The screening identified 352 studies focusing on specific primers (primer dataset), achieving an IRR of 0.766. Among these, 86 studies involved DNA-based rapid tools (rapid tools dataset), with an initial IRR of 0.554. To improve reliability, two reviewers reassessed the rapid tools dataset under a full agreement criterion, confirming 58 studies and improving the IRR to 0.927. One study was inaccessible, and another, which utilized microarray tools not tested in situ was excluded. Additionally, two studies published before 2000 were removed to ensure relevance to contemporary techniques.  The snowballing method (Wohlin 2014) was applied to expand the rapid tools dataset, incorporating both forward (citing) and backward (cited) references. This iterative process was concluded when no new studies were identified, yielding 16 additional studies. The final rapid tools dataset included 68 studies, comprising 62 laboratory-based techniques and six on-site identification tools (More details on Alvarenga & Bunholi et al. (2025) – Appendix 2).

Instrument- or software-specific information needed to interpret the data: we used the statistical software R version 4.2.1 for all data analysis and visualization. Necessary packages to run the scripts: tidyverse, stringr, sf, ggplot2, packcircles, maps, treemap, patchwork, cowplot, reshape2, ggpubr, ggrepel, ggspatial, rfishbase, rredlist, fmsb, scales.

Please contact me for any issue or question (ingridbunholi@gmail.com)
