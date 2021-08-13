
<!-- README.md is generated from README.Rmd. Please edit that file -->
# Background

This repository contains analysis code for the publication found [here](https://pubmed.ncbi.nlm.nih.gov/34309761/), which applied paired bulk-tissue and single-nucleus RNA-sequencing to anterior cingulate cortex samples derived from 28 individuals, including healthy controls and Lewy body disease cases.

Analyses were shared between two groups, thus code herein relates to:

-   Processing and analysis bulk-tissue RNA-sequencing data
-   Generation of sLDSC outputs
-   Generation of manuscript figures

Code used to process and analyse single-nucleus RNA-sequencing data and to generate H-MAGMA outputs is available at: <https://github.com/rahfel/snRNAseqProcessingSteps>.

# Citation

If you use any of the code or data from this repository, please cite our [paper](https://pubmed.ncbi.nlm.nih.gov/34309761/).

# License

The code in this repository is released under an MIT license. This repository is distributed in the hope that it will be useful to the wider community, but without any warranty of any kind. Please see the [LICENSE](LICENSE) file for more details.

# Code contents

All `.Rmd`s can be viewed as interactive `.html`s here: <https://rhreynolds.github.io/LBD-seq-bulk-analyses/>.

Within this repository you will otherwise find:

<table>
<colgroup>
<col width="11%" />
<col width="88%" />
</colgroup>
<thead>
<tr class="header">
<th>Directory</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><a href="docs" class="uri">docs</a></td>
<td>Contains all <code>.Rmd</code>s and their corresponding <code>.html</code>s describing analyses performed for this project. These can be view interactively at: <a href="https://rhreynolds.github.io/LBD-seq-bulk-analyses/" class="uri">https://rhreynolds.github.io/LBD-seq-bulk-analyses/</a></td>
</tr>
<tr class="even">
<td><a href="man" class="uri">man</a></td>
<td>Contains a figure used in <a href="docs/overviews/RNAseq_workflow_tissue.Rmd">overview.Rmd</a></td>
</tr>
<tr class="odd">
<td><a href="misc_scripts" class="uri">misc_scripts</a></td>
<td>Contains analysis scripts that were run separate from <code>.Rmd</code>s. Each script contains a one-line description and is also referenced in its corresponding <code>.Rmd</code>.</td>
</tr>
<tr class="even">
<td><a href="nohup_logs" class="uri">nohup_logs</a></td>
<td>For any scripts that were run outside of an <code>.Rmd</code> (e.g. scripts from the <a href="misc_scripts" class="uri">misc_scripts</a> directory), a log file was recorded and can be accessed here.</td>
</tr>
<tr class="odd">
<td><a href="R" class="uri">R</a></td>
<td>Various functions called in <a href="docs" class="uri">docs</a> and <a href="misc_scripts" class="uri">misc_scripts</a>.</td>
</tr>
</tbody>
</table>
