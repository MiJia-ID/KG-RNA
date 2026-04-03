# Protein–RNA Binding Site Prediction via Mixture of Experts and Knowledge Graph


## Abstract
RNA-binding proteins (RBPs) play central roles in post-transcriptional regulation, yet accurately predicting their RNA-binding sites remains challenging due to two unresolved issues: the pronounced heterogeneity between structured domains and intrinsically disordered regions (IDRs) within RBP sequences, and the absence of traceable biological evidence supporting prediction outputs. To address both challenges, we propose a binding site prediction framework integrating a mixture of experts and a knowledge graph, built around two co-equal pillars. First, a mixture-of-experts architecture comprising a Structural Expert (EGNN-based) and a Sequence Expert (Transformer-based), adaptively fused via a gating mechanism, explicitly models within-protein heterogeneity to produce residue-level binding probability distributions. Second, a protein–RNA knowledge graph constructed from 453,959 experimentally validated interactions grounds each prediction in traceable biological evidence through homology-based retrieval and multi-dimensional evidence-weighted scoring. Experiments on multiple benchmark datasets demonstrate significant improvements over representative baselines including MPBind, GPSite, and DeepPROBind, achieving AUC 0.874, AUPRC 0.692, F1 0.675, and MCC 0.511 (p < 0.05). Ablation studies confirm the complementary contributions of both experts and the gating mechanism. Knowledge graph ranking analysis further demonstrates that evidence-weighted scoring provides biologically meaningful signals beyond sequence similarity alone. Together, this framework advances protein–RNA binding site prediction in both accuracy and interpretability, with potential applications in RBP functional annotation and disease-associated variant interpretation.

<div align=center>
<img src="KG-RNA.jpg" width=75%>
</div>


## Preparation
### Environment Setup
```python 
   git clone https://github.com/MiJia-ID/KG-RNA.git
   conda env create -f environment.yml
```
You also need to install the relative packages to run ESM2, ProtTrans and AlphaFold3 protein language model. 

## Experimental Procedure
### Create Dataset
**Firstly**, you need to use ESM3 to obtain the PDB files of proteins in the tran and test datasets(DNA_train_573, DNA_129_Test and DNA_181_Test ). More details can be found here:https://github.com/evolutionaryscale/esm

Then, run the script below to create node features (PSSM, SS, AF, One-hot encoding). The file is located in the scripts folder.
```python 
python3 ./code/script/train_186_164_test_71/data_io.py 
```

**Secondly** , run the script below to create node features(ESM2 embeddings and ProtTrans embeddings). The file can be found in feature folder.</br>

```python 
python3 ./code/extract_feature/ESM2_5120.py 
```
```python 
python3 ./code/extract_feature/ProtTrans.py 
```
We choose the esm2_t48_15B_UR50D() pre-trained model of ESM-2 which has the most parameters. More details about it can be found at: https://huggingface.co/facebook/esm2_t48_15B_UR50D   </br>
We also choose the prot_t5_xl_uniref50 pre-trained model of ProtTrans, which uses a masked language modeling(MLM). More details about it can be found at: https://huggingface.co/Rostlab/prot_t5_xl_uniref50    </br>

**Thirdly**, run the script below to create edge features. The file can be found in feature folder.
```python 
python3 ./code/feature/create_edge.py 
```

**Fourth**, the training dataset used in this study is provided in the compressed file `data.zip`, which contains all necessary data for model training and evaluation.

The processed protein–RNA knowledge graph data is provided in the compressed file `protein_rna_with_fasta_clean.zip`. This dataset includes curated protein–RNA interaction information along with corresponding FASTA sequences, which are used to construct the knowledge graph and support evidence-based prediction.

Please unzip both files before running the code:

```bash
unzip data.zip
unzip protein_rna_with_fasta_clean.zip
```

### Model Training
Run the following script to train the model.
```python
python3 train_val_bestAUPR_predicted.py 
```

