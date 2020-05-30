
# PaperWork scRNA-GSE111527-mouseSVZ

### Reanalyze the data GSE111527 to get the Figure 2/3 from this paper.

----
The parper: [Single-Cell Transcriptomics Characterizes Cell Types in the Subventricular Zone and Uncovers Molecular Defects Impairing Adult Neurogenesis](https://www.sciencedirect.com/science/article/pii/S2211124718317327)

The data: [GSE111527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527)

## Table of Contents

- [Brief Intro to the Paper](#brief-intro-to-the-paper)
- [Brief Intro to the Data](#brief-intro-to-the-data)
- [Main Analysises and Tools Used](#main-analysis-and-tools-used)
- [Usage](#usage)
- [Contact Us](#contact-us)


## Brief Intro to the Paper

Summary:
> Neural stem cells (NSCs) contribute to plasticity and repair of the adult brain. Niches harboring NSCs regulate stem cell self-renewal and differentiation. We used comprehensive and untargeted single-cell RNA profiling to generate a molecular cell atlas of the largest germinal region of the adult mouse brain, the subventricular zone (SVZ). We characterized >20 neural and non-neural cell types and gained insights into the dynamics of neurogenesis by predicting future cell states based on computational analysis of RNA kinetics. Furthermore, we applied our single-cell approach to document decreased numbers of NSCs, reduced proliferation activity of progenitors, and perturbations in Wnt and BMP signaling pathways in mice lacking LRP2, an endocytic receptor required for SVZ maintenance. Our data provide a valuable resource to study adult neurogenesis and a proof of principle for the power of single-cell RNA sequencing to elucidate neural cell-type-specific alterations in loss-of-function models.

Highlights:
>- Single-cell transcriptomics characterizes the SVZ adult neural stem cell niche
- Different transcriptional dynamics along the neurogenic lineage
- Cell-type-specific dysfunctions underlying impaired adult neurogenesis



## Brief Intro to the Data

Drop-seq sequencing reads for microdissected SVZs from three to five mice. The author collected 9,804 cells with 734 genes and 1,137 unique molecular identifiers (UMIs) quantified per cell (medians) in five replicates. We will perform analysis for these five replicates: an002, an003_F, an003_L, an008 and an009.

## Main Analysises and Tools Used

- Drop-seq reads align and count: [Drop-seq v2.3.0 tools](http://mccarrolllab.org/dropseq/)
- Clustering and visualization: [Seurat v3.1.5](https://satijalab.org/seurat/)
- RNA velocity calculation: [velocyto v0.17.17](http://velocyto.org/velocyto.py/)
- Pseudotime calculation: [SCANPY v1.5.1](https://github.com/theislab/scanpy)

## Usage
### **scripts-preprocess**<br>
  Mapping reads and get digital gene expression matrix (DGE);<br>Geting .loom files using the aligned reads (bam files) from Drop-seq tools

  1- Download data from NCBI and extract fastq data
   ```bash
   bash wget.sh
   bash fastqdump_batch.sh
   ```  
  2- Transfer fastq format to sam format
  ```bash
  bash fastq2sam_batch.sh
  ```
  3- Drop-seq reads align
  ```bash
  bash dropseq_batch_align.sh
  ```
  4- Generate (DGE)
  ```bash
  bash dropseq_batch_count.sh
  ```
  5- Getting loom files
  ```python
  velocyto run -b an002.filtered.cellbarcode.tsv -o ../velocyto/bamfiles/ -m mm10_rmsk.gtf -@ 20 ../velocyto/bamfiles/SRR6814853.bam gencode.vM25.chr_patch_hapl_scaff.annotation.gtf

  ```
  **PS**: due to the lack of detail paramters for seleceting cells (not provied in the paper), we counldn't get the exactly same DGE as the paper, so we used the DGE from  [GSE111527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111527) in the follwing analysis.


### **scripts-plot**<br>
Clustering and visualization of the cells

  1- Generate the merged Seurat object
  ```R
  Rscript step2_seurat_mergeReplicates_new.R
  ```

  2- Clustering and tSNE plot (**Figure 2A**)
  ```R
  Rscript step3_seurat_clustering_new.R
  ```

  3- Stacked violin plots,featureplot and heatmap for marker genes (**Figure 2B-2D**)

  ```
  Rscript step4_seurat_makerExpression.R
  ```

  4- Estimating RNA velocity and visualization (**Figure 3A**)

  ```
  Rscript step5_velocyto.R
  ```
  or
  ```python
  python velocyto.test.py
  ```

  5- Subclustering of neurogenesis cells (**Figure 3C**)

  ```
  Rscript step6_seurat_subclustering.R
  ```

  6- Gene set scoring for subclustering cells (**Figure 3D**)

  ```
  Rscript step7_geneSet_score.R
  ```

  7- calculation of pseudotime and visualization for subclustering cells (**Figure 3E**)

  ```
  Rscript step9_scanpy.R
  ```
  or

  ```python
  python scanpy.py
  ```


### **data**
DGE counts; Seurat obeject; loom files
### **results**
visualization results in pdf format

<a href="https://sm.ms/image/bVa2U5OnDglRmwr" target="_blank"><img   src="https://i.loli.net/2020/05/29/bVa2U5OnDglRmwr.png" title="Figure 2A paper" width="410" hspace="20" ><img src="https://i.loli.net/2020/05/30/FA7phdNIgfxcWJv.png"  title="Figure 2A our" width="420"><p>Figure 2A</p></a>


<a href="https://sm.ms/image/MujASo7LDrbRq1x" target="_blank"><img src="https://i.loli.net/2020/05/29/MujASo7LDrbRq1x.png"  title="Figure 2B paper" width=420  height="380" hspace="20"><img src="https://i.loli.net/2020/05/30/5AnPq2cuFpbIzKv.png"  title="Figure 2B our" width=420/><p>Figure 2B</p></a>

<a href="https://sm.ms/image/4uJB2PQHfnkME61" target="_blank"><img src="https://i.loli.net/2020/05/29/4uJB2PQHfnkME61.png"  title="Figure 2C paper" width="340" height="600" hspace="30"><img src="https://i.loli.net/2020/05/30/MNYH8nUuCQpGbw4.png" title="Figure 2C our" width="300" height="700"><p>Figure 2C</p></a>

<a href="https://sm.ms/image/sPW92ApDvrojGKq" target="_blank"><img src="https://i.loli.net/2020/05/29/sPW92ApDvrojGKq.png"  title="Figure 2D paper" width="480" ><img src="https://i.loli.net/2020/05/30/HtmiZjb8x5SawVd.png" title="Figure 2D our" width="400"><p>Figure 2D</p></a>

<a href="https://sm.ms/image/yRGAutcF7CaDLfJ" target="_blank"><img src="https://i.loli.net/2020/05/29/yRGAutcF7CaDLfJ.png"  title="Figure 3A paper" width="420"><img src="https://i.loli.net/2020/05/30/oNi5U1fwmnuP2Er.png" title="Figure 3A our" width="460"><p>Figure 3A</p></a>

<a href="https://sm.ms/image/4fIpOZh5CymcAUJ" target="_blank"><img src="https://i.loli.net/2020/05/29/4fIpOZh5CymcAUJ.png"  title="Figure 3C paper" width="480"><img src="https://i.loli.net/2020/05/30/lzejQiUMYwvXbfy.png" title="Figure 3C our" width="400"><p>Figure 3C</p></a>

<a href="https://sm.ms/image/OZ3muRjAUJSfvyF" target="_blank"><img src="https://i.loli.net/2020/05/29/OZ3muRjAUJSfvyF.png"  title="Figure 3D paper" width="500"><img src="https://i.loli.net/2020/05/30/TPEM4kmWCw3AIVU.png" title="Figure 3D our" width="350"><p>Figure 3D</p></a>
<a href="https://sm.ms/image/sd6Mm1T7WyBzS2n" target="_blank"><img src="https://i.loli.net/2020/05/29/sd6Mm1T7WyBzS2n.png"  title="Figure 3E paper" width="450"><img src="https://i.loli.net/2020/05/30/v5AqOTdKSfl8I3W.png" title="Figure 3E our" width="430"></a>



## Contact Us

If you have any questions about these scripts, please feel free to contact  Yuting Liu, [lyt17@pku.edu.cn](lyt17.pku.edu.cn)
