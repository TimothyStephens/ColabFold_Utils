#!/usr/bin/env python
DESCRIPTION = '''
Plot PDB file.

This script will produce a HTML file with:
 - An interactive 3D viewer of the (best) PDB file (coloring by AlphaFold2 confidence score).
 - Any additional *.png files produced by ColabFold-AlphaFold2
    (auto detected by this script).
 - And, if the PDB file contains a multimeric complex, an additional 3D viewer with each 
    monomer colored differently.

NOTE:
 - Can be given either the output directory from colabfold_batch (will pick top ranked 
   relaxed or unrelaxed structure) or a specific PDB file to analyze.
'''
import sys
import os
import shutil
import argparse
import logging
import gzip
from pathlib import Path



## Pass arguments.
def main():
    ## Pass command line arguments. 
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=DESCRIPTION)
    parser.add_argument('results',
        type=str, 
        help='One of: 1) directory with PDB files, or 2) a specific PDB file'
    )
    parser.add_argument('plots', 
        type=str, 
        help='HTML files output directory.'
    )
    parser.add_argument('--colabfold',
        required=False, action='store_true',
        help='Directory is ColabFold output and we want to plot just the best PDB file (default: %(default)s)'
    )
    parser.add_argument('--all_chains',
        required=False, action='store_true',
        help='Calculate quality stats using all chains, not just chain A (default), this makes sense if you have a heteromeric complex but not so much for homomeric complexes (default: %(default)s)'
    )
    parser.add_argument('--dont_keep_files',
        required=False, action='store_true',
        help='Dont keep images, PDB, and Rmd files used for plotting; makes "plots" dir smaller but prevents rebuilding of Rmd file (default: %(default)s)'
    )
    parser.add_argument('--debug', 
        required=False, action='store_true', 
        help='Print DEBUG info (default: %(default)s)'
    )
    args = parser.parse_args()
    
    ## Set up basic debugger
    logFormat = "[%(levelname)s]: %(message)s"
    logging.basicConfig(format=logFormat, stream=sys.stderr, level=logging.INFO)
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logging.debug('%s', args) ## DEBUG
    
    colabfold_plot(args.results, args.plots, args.colabfold, args.all_chains, args.dont_keep_files)

RMD_HEADER='''---
title: "AlphaFold2 Structure of <<<PROTEIN_ID>>>"
author: "colabfold_plot.py"
date: "`r format(Sys.time(), '%d %B, %Y')`"
output:
  html_document
---


```{r setup, echo=FALSE}
## Setup
#Setup R env. Load packages and set default image export formats, size and resolution.

knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 16, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(r3dmol)
options(scipen = 999) #Prevent scientific notation
```

'''

RMD_3D_PLOT_CONFIDENCE_LEGEND=os.path.dirname(os.path.abspath(__file__))+"/AlphaFold2-legend.png"

RMD_3D_PLOT_CONFIDENCE='''
# Protein 3D structure

Average structure confidence: <<<CONFIDENCE>>>

```{r plot_AlphaFold2_structure, echo=FALSE}
pdb.file  <- "<<<PDB>>>"

width=800
height=800

## Setup initial protein 3D viewer
# - Color by by b-chain info (i.e., AlphaFold2 prediction confidence)
viewer <- r3dmol(
    width=width,
    height=height,
  viewer_spec = m_viewer_spec(
    cartoonQuality = 10,
    lowerZoomLimit = 50,
    upperZoomLimit = 350,
  ),
) %>%
  # Add model to scene
  m_add_model(data = pdb.file, format = "pdb") %>%
  # Zoom to encompass the whole scene
  m_zoom_to() %>%
  # Set style of structures
  m_set_style(style = m_style_cartoon(
    colorfunc = "
        function(atom) {
          if (atom.b < 50) {return '#ef7c45'};
          if (atom.b < 70 && atom.b >= 50) {return '#f8db13'};
          if (atom.b < 90 && atom.b >= 70) {return '#65cbf3'};
          if (atom.b >= 90) {return '#2a54d6'};
          return 'white';
        }"
  ))

# - Color by by chain (i.e., we have a multimeric complex; color each protein in complex as a separate color)

## Create interactive results HTML
legend <- htmltools::img(src = knitr::image_uri("AlphaFold2-legend.png"), 
                         alt = "logo", 
                         style = paste('float: bottom;padding-bottom:0px;height:',40*4,'px;width:',160*4,'px', sep="")
                         )

viewer <- htmlwidgets::prependContent(viewer, legend)
viewer
```

'''


RMD_3D_PLOT_CONFIDENCE_AND_CHAIN='''
# Protein 3D structure

Average structure confidence: <<<CONFIDENCE>>>

```{r plot_AlphaFold2_structure, echo=FALSE}
pdb.file  <- "<<<PDB>>>"

width=800
height=800

## Setup initial protein 3D viewer
# - Color by by b-chain info (i.e., AlphaFold2 prediction confidence)
viewer1 <- r3dmol(
    width=width,
    height=height,
  viewer_spec = m_viewer_spec(
    cartoonQuality = 10,
    lowerZoomLimit = 50,
    upperZoomLimit = 350,
  ),
) %>%
  # Add model to scene
  m_add_model(data = pdb.file, format = "pdb") %>%
  # Zoom to encompass the whole scene
  m_zoom_to() %>%
  # Set style of structures
  m_set_style(style = m_style_cartoon(
    colorfunc = "
        function(atom) {
          if (atom.b < 50) {return '#ef7c45'};
          if (atom.b < 70 && atom.b >= 50) {return '#f8db13'};
          if (atom.b < 90 && atom.b >= 70) {return '#65cbf3'};
          if (atom.b >= 90) {return '#2a54d6'};
          return 'white';
        }"
  ))

# - Color by by chain (i.e., we have a multimeric complex; color each protein in complex as a separate color)
viewer2 <- r3dmol(
    width=width,
    height=height,
  viewer_spec = m_viewer_spec(
    cartoonQuality = 10,
    lowerZoomLimit = 50,
    upperZoomLimit = 350,
  ),
) %>%
  # Add model to scene
  m_add_model(data = pdb.file, format = "pdb") %>%
  # Zoom to encompass the whole scene
  m_zoom_to() %>%
  # Set style of structures
  m_set_style(style = m_style_cartoon(
    colorfunc = "
        function(atom) {
          if (atom.chain == 'A') {return '#1b9e77'};
          if (atom.chain == 'B') {return '#d95f02'};
          if (atom.chain == 'C') {return '#7570b3'};
          if (atom.chain == 'D') {return '#e7298a'};
          if (atom.chain == 'E') {return '#66a61e'};
          if (atom.chain == 'F') {return '#e6ab02'};
          if (atom.chain == 'G') {return '#a6761d'};
          if (atom.chain == 'H') {return '#666666'};
          if (atom.chain == 'I') {return '#e41a1c'};
          if (atom.chain == 'J') {return '#377eb8'};
          if (atom.chain == 'K') {return '#4daf4a'};
          if (atom.chain == 'L') {return '#984ea3'};
          if (atom.chain == 'M') {return '#ff7f00'};
          if (atom.chain == 'N') {return '#ffff33'};
          if (atom.chain == 'O') {return '#a65628'};
          if (atom.chain == 'P') {return '#f781bf'};
          if (atom.chain == 'Q') {return '#66c2a5'};
          if (atom.chain == 'R') {return '#fc8d62'};
          if (atom.chain == 'S') {return '#8da0cb'};
          if (atom.chain == 'T') {return '#e78ac3'};
          if (atom.chain == 'U') {return '#a6d854'};
          if (atom.chain == 'V') {return '#ffd92f'};
          if (atom.chain == 'W') {return '#e5c494'};
          if (atom.chain == 'X') {return '#b3b3b3'};
          if (atom.chain == 'Y') {return '#80b1d3'};
          if (atom.chain == 'Z') {return '#ffffb3'};
          return 'white';
        }"
  ))

## Connect two color schemes so they rotate together
viewer <- m_grid(
  viewer = list(viewer1, viewer2),
  rows = 1,
  cols = 2,
  control_all = TRUE,
  viewer_config = m_viewer_spec(
    backgroundColor = "white"
  )
)

## Create interactive results HTML
legend <- htmltools::img(src = knitr::image_uri("AlphaFold2-legend.png"), 
                         alt = "logo", 
                         style = paste('float: bottom;padding-bottom:0px;height:',40*4,'px;width:',160*4,'px', sep="")
                         )

viewer <- htmlwidgets::prependContent(viewer, legend)
viewer
```

'''


RMD_COVERAGE_PLOT='''
![Sequence Alignment Coverage across Protein ](<<<FILE>>>)
'''
RMD_plDDT_PLOT='''
![Predicted lDDT across Protein per Rank](<<<FILE>>>)
'''
RMD_PAE_PLOT='''
![Predicted Aligned Error across Protein per Rank](<<<FILE>>>)
'''


RMD_SESSION_INFO='''
# Session Info

```{r ressionInfo, echo=FALSE}
sessionInfo()
```

'''


def colabfold_plot(results, plots, colabfold, all_chains, dont_keep_files):
    # Create output directory
    os.makedirs(plots, exist_ok=True)
    
    # For each PDB file that we have found
    for prefix, in_pdb_file, in_cov_file, in_pae_file, in_plddt_file in get_files_for_plotting(results, colabfold):
        prefix = os.path.basename(prefix) # Remove directory name from prefix
        logging.info(f"Processing: {in_pdb_file}; Output Prefix: {prefix}") ## INFO
        
        # Output variables
        out_pdb   = f"{plots}/{prefix}.pdb"
        out_rmd   = f"{plots}/{prefix}.Rmd"
        out_cov   = f"{plots}/{prefix}.coverage.png"
        out_pae   = f"{plots}/{prefix}.pae.png"
        out_plddt = f"{plots}/{prefix}.plddt.png"
        out_stats = f"{plots}/{prefix}.confidence.tsv"
        
        # Get the number of copies
        copies = 0
        previous_chain = ''
        count_total = 0
        count_vh    = 0
        count_h     = 0
        count_l     = 0
        count_vl    = 0
        b_sum       = 0
        
        shutil.copyfile(in_pdb_file, out_pdb)
        current_residue_number = ''
        for line in open(out_pdb, "r").readlines():
            if line.startswith('ATOM'):
                current_chain  = line[21]
                residue_number = line[22:26]
                
                # Skip entry if we have already seen an atom from this residue before (assume all atoms have same b-factor value)
                if current_residue_number == residue_number:
                    continue
                else:
                    current_residue_number = residue_number
                
                # Check which chain we are on (normally just chain 'A', but can be higher if multimeric protein)
                if previous_chain != current_chain:
                    previous_chain = current_chain
                    copies += 1
                
                # Count b-factor values for just A chain, or all chains if all_chains==True
                #  - Using all chains will inflate length if we have homodimers but makes sense if we have heterodimers.
                if current_chain == 'A' or all_chains:
                    b_factor = float(line[60:66])
                    b_sum += b_factor
                    count_total += 1
                    if b_factor >= 90:
                        count_vh += 1
                    elif b_factor >= 70:
                        count_h  += 1
                    elif b_factor >= 50:
                        count_l  += 1
                    else: # < 50
                        count_vl += 1
        logging.info(f"Copies in complex: {copies}") ## DEBUG
        
        ## Compute confidence stats
        pLDDT_mean = (float(b_sum)/float(count_total))
        percent_vh = (float(count_vh)/float(count_total))*100
        percent_h  = (float(count_h) /float(count_total))*100
        percent_l  = (float(count_l) /float(count_total))*100
        percent_vl = (float(count_vl)/float(count_total))*100
        
        if pLDDT_mean >= 90:
            pLDDT_overall_quality = "Very_High"
        elif pLDDT_mean >= 70:
            pLDDT_overall_quality = "High"
        elif pLDDT_mean >= 50:
            pLDDT_overall_quality = "Low"
        else: # < 50
            pLDDT_overall_quality = "Very_Low"
        
        t_head = []
        t_stat = []
        
        t_head.append(f"PDB_file")
        t_stat.append(f"{out_pdb}")
        
        t_head.append(f"length")
        t_stat.append(f"{count_total}")
        
        t_head.append(f"pLDDT_mean")
        t_stat.append(f"{pLDDT_mean}")
        
        t_head.append(f"pLDDT_overall_quality")
        t_stat.append(f"{pLDDT_overall_quality}")
        
        t_head.append(f"pLDDT_VERY_LOW_residue_count")
        t_stat.append(f"{count_vl}")
        
        t_head.append(f"pLDDT_VERY_LOW_residue_percent")
        t_stat.append(f"{percent_vl}")
        
        t_head.append(f"pLDDT_LOW_residue_count")
        t_stat.append(f"{count_l}")
        
        t_head.append(f"pLDDT_LOW_residue_percent")
        t_stat.append(f"{percent_l}")
        
        t_head.append(f"pLDDT_HIGH_residue_count")
        t_stat.append(f"{count_h}")
        
        t_head.append(f"pLDDT_HIGH_residue_percent")
        t_stat.append(f"{percent_h}")
        
        t_head.append(f"pLDDT_VERY_HIGH_residue_count")
        t_stat.append(f"{count_vh}")
        
        t_head.append(f"pLDDT_VERY_HIGH_residue_percent")
        t_stat.append(f"{percent_vh}")
        
        with open(out_stats, 'w') as f:
            f.write('\t'.join(t_head)+'\n')
            f.write('\t'.join(t_stat)+'\n')
        
        ## Assemble R Markdown script
        tmp =  RMD_HEADER.replace('<<<PROTEIN_ID>>>', out_pdb)
        
        # 3D viewer
        if copies == 1:
            tmp += RMD_3D_PLOT_CONFIDENCE.replace('<<<PDB>>>', os.path.basename(out_pdb))
        else:
            tmp += RMD_3D_PLOT_CONFIDENCE_AND_CHAIN.replace('<<<PDB>>>', os.path.basename(out_pdb))
        
        # Copy legend file
        shutil.copyfile(RMD_3D_PLOT_CONFIDENCE_LEGEND, plots+"/"+os.path.basename(RMD_3D_PLOT_CONFIDENCE_LEGEND))
        
        # Add confidence value to file
        if pLDDT_mean >= 90:
            tmp = tmp.replace("<<<CONFIDENCE>>>", "<span style=\"color:#2a54d6\">Very High</span>")
        elif pLDDT_mean >= 70:
            tmp = tmp.replace("<<<CONFIDENCE>>>", "<span style=\"color:#65cbf3\">High</span>")
        elif pLDDT_mean >= 50:
            tmp = tmp.replace("<<<CONFIDENCE>>>", "<span style=\"color:#f8db13\">Low</span>")
        else: # < 50
            tmp = tmp.replace("<<<CONFIDENCE>>>", "<span style=\"color:#ef7c45\">Very Low</span>")
        
        # Coverage plot
        if not in_cov_file is None:
            shutil.copyfile(in_cov_file, out_cov)
            tmp += RMD_COVERAGE_PLOT.replace('<<<FILE>>>',os.path.basename(out_cov))
        
        # pLDDT file    
        if not in_plddt_file is None:
            shutil.copyfile(in_plddt_file, out_plddt)
            tmp += RMD_plDDT_PLOT.replace('<<<FILE>>>',os.path.basename(out_plddt))
        
        # PAE file
        if not in_pae_file is None:
            shutil.copyfile(in_pae_file, out_pae)
            tmp += RMD_PAE_PLOT.replace('<<<FILE>>>',os.path.basename(out_pae))
        
        # Session Info
        tmp += RMD_SESSION_INFO
        
        # Write R Markdown file
        with open(out_rmd, 'w') as rmd:
            rmd.write(tmp)
        
        # Run R Markdown script
        cmd = f"Rscript -e \"rmarkdown::render('{out_rmd}')\""
        logging.debug(f"Running cmd: {cmd}") ## DEBUG
        ret = os.system(cmd)
        
        # Fail if command returned non-zero error code.
        if ret != 0:
            raise Exception(f'Failed to render R Markdown file! See the above log for details.')
        
        # Print explanation of output files
        if dont_keep_files:
            if_exist_remove(out_pdb)
            if_exist_remove(out_rmd)
            if_exist_remove(out_stats)
            if_exist_remove(plots+"/"+os.path.basename(RMD_3D_PLOT_CONFIDENCE_LEGEND))
            if_exist_remove(out_cov)
            if_exist_remove(out_plddt)
            if_exist_remove(out_pae)



def get_files_for_plotting(input_path, colabfold):
    '''
    Read one of:
      1) directory with colabfold_batch PDB files (if --colabfold set);
      2) directory with PDB files to plot;
      3) a specific PDB file;
    '''
    
    # Check 'input_path' exists in some form
    input_path = Path(input_path)
    if not input_path.exists():
        raise OSError(f"{input_path} could not be found")
    
    ## Get PDB file from 'results'
    # 'results' is a file
    if input_path.is_file():
        if input_path.suffix == ".pdb":
            prefix     = str(input_path).rstrip(".pdb")
            pdb_file   = input_path
            cov_file   = check_file_exists(f"{prefix}_coverage.png")
            pae_file   = check_file_exists(f"{prefix}_pae.png")
            plddt_file = check_file_exists(f"{prefix}_plddt.png")
            yield(prefix, pdb_file, cov_file, pae_file, plddt_file)
        else:
            raise ValueError(f"Unknown file format {input_path.suffix}")
    
    # 'results' is a dir
    else:
        assert input_path.is_dir(), "Expected either an input file or a input directory"
        
        # Run in PDB or A3M modes
        if not colabfold:
            # Plot all PDB files
            for file in sorted(input_path.iterdir()):
                if not file.is_file():
                    continue
                if file.suffix.lower() not in [".pdb"]:
                    continue
                
                # Yield each set
                prefix     = str(file).rstrip(".pdb")
                pdb_file   = file
                cov_file   = check_file_exists(f"{prefix}_coverage.png")
                pae_file   = check_file_exists(f"{prefix}_pae.png")
                plddt_file = check_file_exists(f"{prefix}_plddt.png")
                yield(prefix, pdb_file, cov_file, pae_file, plddt_file)
                
        else:
            # For each A3M file in 'results' (infer the number of queries)
            prefixes = []
            for file in sorted(input_path.iterdir()):
                if not file.is_file():
                    continue
                if file.suffix.lower() not in [".a3m"]:
                    continue
                prefixes.append(str(file).rstrip(".a3m"))
            if len(prefixes) == 0:
                raise ValueError(f"No *.a3m files found in {input_path}")
            
            logging.debug(f"A3M prefixes found: {prefixes}") ## DEBUG
            
            # For each A3M prefix (i.e., each protein query in results)
            for prefix in prefixes:
                queries = []
                for file in sorted(input_path.iterdir()):
                    logging.debug(f"Checking {prefix} {file}") ## DEBUG
                    if not file.is_file():
                        continue
                    if file.suffix.lower() not in [".pdb"]:
                        continue
                    if not str(file).startswith(prefix):
                        continue
                    queries.append(file)
                
                # Get either top ranked relaxed structure, or if relaxation not performed, top ranked unrelazed structure
                for file in queries:
                    if '_relaxed_rank_001_' in str(file):
                        pdb_file = file
                        break
                    if '_unrelaxed_rank_001_' in str(file):
                        pdb_file = file
                        break
                else:
                    raise ValueError(f"No top relaxed or unrelaxed PDB file found in {queries}")
            
                # Get additional plot paths
                cov_file   = check_file_exists(f"{prefix}_coverage.png")
                pae_file   = check_file_exists(f"{prefix}_pae.png")
                plddt_file = check_file_exists(f"{prefix}_plddt.png")
                
                # Yield each set
                yield(prefix, pdb_file, cov_file, pae_file, plddt_file)




def if_exist_remove(file):
    file = Path(file)
    if file.is_file():
        os.remove(file)


def check_file_exists(file):
    file = Path(file)
    if file.is_file():
        return(file)
    else:
        return(None)



if __name__ == '__main__':
    main()
