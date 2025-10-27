from os import read
import anndata as ad
from pathlib import Path
import sys
import scanpy as sc
import datetime
import numpy as np
import matplotlib.pyplot as plt
basePath = Path(__file__).resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "utils"))
import logging
from utils.read import scTERead
from utils.QC import lowquality,Doublet_scrub
from utils.batch import Run_Normalization,Run_batchRemove, combine_group
# from utils.annotation import auto_omi
from utils.annotation import Show_Markers,handful_annotate,show_annotation
logging.basicConfig(
	level=logging.INFO,
	format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
	datefmt='%Y-%m-%d %H:%M:%S'
)

samplesDict = { 'Intestine': ['CKO-chang-10XSC3', 'E2-chang-10XSC3', 'TRA-chang-10XSC3', 'WT-chang-10XSC3'],
        'Lung': ['CKO-fei-10XSC3', 'E2-fei-10XSC3', 'TRA-fei-10XSC3', 'WT-fei-10XSC3'],
        'Muscle': ['CKO-jirou-10XSC3', 'E2-jirou-10XSC3', 'TRA-jirou-10XSC3', 'WT-jirou-10XSC3']}

marker_genes_Intestine = {
    "T cytotoxic cell": ["Cd8a", "Gzmb", "Trac"],  # 细胞毒性 T 细胞
    "T helper cell": ["Cd4", "Ccr4", "Il13"],  # 辅助性 T 细胞
    "Enterocyte cell": ["Slc10a2","Alpi", "Krt20"],
    "B cell": ["Pxk","Cd19", "Cd79a", "Ms4a1"],
    "Intestinal stem cell": ["Lgr5","Msi1","Ascl2", "Olfm4"],  # Intestinal stem cell
    "Goblet": ["Muc2", "Muc5ac", "Clca1"],
    "Fibroblast": ["Dcn", "Col3a1", "Sparc","Col1a2"],
    "Enteroendocrine cell": ["Chga", "Pcsk1", "Sct","Hepacam2"],
    "Tuft": ["Avil", "Siglecf", "Pou2f3"],
    "Paneth": ["Lyz1", "Defa5", "Mmp7"],
    "Macrophage": ["Fcgr1", "Cd68", "Naaa","Lyz2"],
    "Plasma cell": ["Jchain","Mzb1"],
    "Plasmacytoid dendritic cells": ["Irf8","Siglech","Bst2","Tcf4","Lin"] #https://www.rndsystems.com/resources/cell-markers/immune-cells/dendritic-cells/plasmacytoid-dendritic-cell-markers
}
marker_genes_Lung = {
    "endothelial cells": ["Pecam1", "Cdh5", "Cd31","Cldn5"],
    "Macrophage": ["Fcgr1", "Cd68", "Naaa","Lyz2"],
    "B cell": ["Pxk","Cd19", "Cd79a", "Ms4a1"],
    "Neutrophil": ["Csf3r", "S100a8", "Ly6g","Il1r2"],
    "Airway epithelial cell": ["Aqp1","Cdh1","Sec14l3"],
    "Airway goblet cells": ["Muc16","Dusp4","Ggh","Muc5b"], # "Airway epithelial cell"的一种亚型
    "Alveolar macrophages": ["Marco","Siglecf","Itgax","Mrc1"],
    "Ciliated cell": ["Foxj1","Ccdc17","Tubb4b","Appl2"],
    "Clara cells": ["Scgb1a1","Cyp4b1","Sftpa1","Ern1"],
    "Ionocytes": ["Foxi1","Cftr","Ascl3","Foxn4"],
    "Pulmonary alveolar type I cells": ["Cldn18","Rtkn2","Col4a3","Fstl3"],
    "Pulmonary alveolar type II cells":["Sftpc","Slc34a2","Sftpa1","Nkx2-1"],
    "T helper cell": ["Cd4", "Ccr4", "Il13"],
    "smooth muscle cells": ["Acta2","Myh11","Tagln"],
    "Aerocyte": ["Apln", "Car4", "Ednrb","Tbx2"],  # 肺毛细血管特殊内皮
    "Fibroblast": ["Vim", "Lum", "Pdgfra"],
    "Pericyte": ["Pdgfrb", "Rgs5", "Acta2"],
    "monocyte": ["Cfp", "Cd14", "Rhoc","Ccr2"],
    "Airway smooth muscle cells": ["Nog", "Acta2", "Foxf1","Gata5"],  # 平滑肌细胞
    "T cytotoxic cell": ["Cd8a", "Gzmb", "Trac"],  # 细胞毒性 T 细胞 (CD8+)
    "Mki67_proliferating_cells": ["Mki67", "Top2a", "Pcna"],
    "Lymphatic endothelial cell": ["Ccl21a", "Lyve1", "Flt4","Prox1"],
    "Basophil": ["Ccl4", "Npl", "Wrn","Nfil3"],
    "Secretory cell": ["5330417C22Rik","8430408G22Rik","Agr2","Ano1"],
    "type-2 innate lymphoid cell": ["Il1rl1","Gata3","Il7r","Rora"]
}
marker_genes_Muscle = {
    "FastIIB": ["Myh4"],
    "FastIIX": ["Myh1"],
    "FastIIA": ["Myh2"],
    "Type I fibers": ["Myh7"],
    "Tenocytes": ["TNMD","SCX","MKX"],
    "Schwann cells":["CDH19","MPZ"],
    "Pericyte cells":["RGS5","ACTA2","NOTCH3","NG2"],
    "Myogenic progenitors":["PAX7","MYF5"],
    "Myocytes":["MYOG","ACTC1"],
    "Myoblasts":["MYOD1","MYF5"],
    "Monocytes":["FCER1G","CD14","C1QA"],
    "Lymphocytes":["CD52","SKAP1","IKZF1","CD19","CD28"],
    "Fibroblasts":["COL3A1","COL1A1","PDGFRA"],
    "Erythroid cells":["SNCA"],
    "Endothelial cells":["PECAM1","CDH5","KDR","ESAM"],
    "Satellite cell": ["Myf5","Myod1","Pax7"]
}

Intestine_annotation = {
    "0": "T cytotoxic cell",
    "1": "T cytotoxic cell",
    "2": "T cytotoxic cell",
    "3": "T cytotoxic cell",
    "4": "B cell",
    "5": "Plasma cell",
    "6": "Enteroendocrine cell",
    "7": "T cytotoxic cell",
    "8": "Enterocyte cell",
    "9": "Intestinal stem cell",
    "10": "Plasmacytoid dendritic cells",
    "11": "Macrophage",
    "12": "Fibroblast",
    "13": "Paneth",
    "14": "T cytotoxic cell",
    "15": "Tuft",
    "16": "Goblet"
    
}
Lung_annotation = {
    "0": "Macrophage",
    "1": "endothelial cells",
    "2": "endothelial cells",
    "3": "Neutrophil",
    "4": "T cytotoxic cell(Il7r)",
    "5": "smooth muscle cells",
    "6": "Fibroblast",
    "7": "B cell",
    "8": "Airway goblet cells",
    "9": "endothelial cells",
    "10": "Pulmonary alveolar type II cells",
    "11": "Macrophage",
    "12": "Lymphatic endothelial cell",
    "13": "Macrophage",
    "14": "Macrophage(Mki67)", 
    "15": "Pulmonary alveolar type I cells",
}
# 肌肉组织是单核测序吗，需要去问一下
Muscle_annotation = {
    "0":"hybrid IIX/IIB fibers",
    "1":"hybrid IIX/IIB fibers",
    "2":"hybrid IIX/IIB fibers",
    "3":"hybrid IIX/IIB fibers",
    "4":"hybrid IIX/IIB fibers",
    "5":"hybrid IIX/IIB fibers",
}
    
### read the output of scTE, and store it in h5ad
class  Single:
    def read():
        indir = "/disk5/luosg/scRNAseq/output/scTE"
        outdir = "/disk5/luosg/scRNAseq/output/h5ad"
        for group in samplesDict.keys():
            for sample in samplesDict[group]:
                infile = f"{indir}/{group}/{sample}.csv"
                # print(infile)
                adata = scTERead(infile,sample)
                output_path = Path(f"{outdir}/{group}/{sample}.h5ad")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(output_path)

    def qc():
        indir = "/disk5/luosg/scRNAseq/output/h5ad"
        outdir = "/disk5/luosg/scRNAseq/output/h5ad_QC"
        start = datetime.datetime.now()
        for group in samplesDict.keys():
            for sample in samplesDict[group]:
                infile = f"{indir}/{group}/{sample}.h5ad"
                # print(infile)
                adata = sc.read_h5ad(infile)
                if group == 'Muscle':
                    adata = lowquality(adata,sample,mt_cutoff=20,use_mt_outlier=False,low_count_outlier=3,low_gene_outlier=3,top20_outlier=3,
                                    fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/QC")
                    adata = Doublet_scrub(adata,sample,scrub_estimated_rate = 0.06,fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/doublet/plot",
                            table_flag=True,table_dir="/disk5/luosg/scRNAseq/output/result/doublet/table")
                elif group == 'Lung':
                    continue
                    # adata = lowquality(adata,sample,use_mt_outlier=False,low_count_outlier=3,low_gene_outlier=3,top20_outlier=3,
                    # fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/QC")
                    # adata = Doublet_scrub(adata,sample,fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/doublet/plot",
                    #     table_flag=True,table_dir="/disk5/luosg/scRNAseq/output/result/doublet/table")
                else:
                    continue
                    # adata = lowquality(adata,sample,use_mt_outlier=False,low_count_outlier=3,low_gene_outlier=3,top20_outlier=3,
                    # fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/QC")
                    # adata = Doublet_scrub(adata,sample,fig_flag = True,fig_dir = "/disk5/luosg/scRNAseq/output/result/doublet/plot",
                    #         table_flag=True,table_dir="/disk5/luosg/scRNAseq/output/result/doublet/table")
                output_path = Path(f"{outdir}/{group}/{sample}_QC.h5ad")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(output_path)
        end = datetime.datetime.now()
        print("程序运行时间："+str((end-start).seconds/3600)+"h")
        

    def batch():
        indir = "/disk5/luosg/scRNAseq/output/h5ad_QC"
        outdir = "/disk5/luosg/scRNAseq/output/result/h5ad"
        start = datetime.datetime.now()
        for group in samplesDict.keys():
            for sample in samplesDict[group]:
                infile = f"{indir}/{group}/{sample}_QC.h5ad"
                # print(infile)
                adata = sc.read_h5ad(infile)
                adata = Run_Normalization(adata,sample,batch='sample',fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/Normalization")
                adata = Run_batchRemove(adata,sample,batch='sample',fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/BatchRemove")
                output_path = Path(f"{outdir}/{group}/{sample}_final.h5ad")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(output_path)

        end = datetime.datetime.now()
        print("程序运行时间："+str((end-start).seconds/3600)+"h")

class Qc_Allcombine:
    def step3Combine():
        """
        combine all adata into one:
        - first, add sample information into adata (complete when read)
        - second, add group infromation into adata (complete in  this step)
        note:
        combine after qc (because not all group are the same qc parameters)
        """
        indir = "/disk5/luosg/scRNAseq/output/h5ad_QC"
        outdir = "/disk5/luosg/scRNAseq/output/result/h5ad"
        start = datetime.datetime.now()
        allList = []
        for group in samplesDict.keys():
            for sample in samplesDict[group]:
                infile = f"{indir}/{group}/{sample}_QC.h5ad"
                # print(infile)
                adata = sc.read_h5ad(infile)
                adata.obs['group'] = np.full(adata.n_obs, group) # read: add sample name, now : add group name
                experiment = sample.split("-")[0]
                adata.obs['experiment'] = experiment # experimet
                allList.append(adata)
        adata = combine_group(allList)
        adata = Run_Normalization(adata,'all',fig_flag=True,n_pcs = 30,n_top_genes = 2000, do_scale=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/Normalization")
        adata = Run_batchRemove(adata,'all', batch=["group","experiment"] ,n_pcs = 30,fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/BatchRemove")
        output_path = Path(f"{outdir}/all_final.h5ad")
        output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
        adata.write_h5ad(output_path)
        end = datetime.datetime.now()
        print("程序运行时间："+str((end-start).seconds/3600)+"h")

class Qc_GroupCombine:
    def step3GroupCombine():
        start = datetime.datetime.now()
        indir = "/disk5/luosg/scRNAseq/output/h5ad_QC"
        outdir = "/disk5/luosg/scRNAseq/output/result/h5ad"
        for group in samplesDict.keys():
            group_list = []
            for sample in samplesDict[group]:
                infile = f"{indir}/{group}/{sample}_QC.h5ad"
                adata = sc.read_h5ad(infile)
                experiment = sample.split("-")[0]
                adata.obs['experiment'] = experiment
                group_list.append(adata)
            adata = combine_group(group_list)
            if group == 'Muscle':
                # continue
                adata = Run_Normalization(adata,group,batch="experiment",fig_flag=True,resolution=1,n_pcs = 100,n_top_genes = 2200, do_scale=False,fig_dir = "/disk5/luosg/scRNAseq/output/result/Normalization")
                adata = Run_batchRemove(adata,group, batch=["experiment"] ,resolution=1,n_pcs = 35 ,fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/BatchRemove")
            elif group == 'Lung':
                # continue
                adata = Run_Normalization(adata,group,batch="experiment",fig_flag=True,resolution=0.4,n_pcs = 100,n_top_genes = 2200, do_scale=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/Normalization")
                adata = Run_batchRemove(adata,group, batch=["experiment"] ,resolution=0.4,n_pcs = 35,fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/BatchRemove")
            elif group == 'Intestine':
                # continue
                adata = Run_Normalization(adata,group,batch="experiment",fig_flag=True,resolution=0.4,n_pcs = 100,n_top_genes = 2200, do_scale=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/Normalization")
                adata = Run_batchRemove(adata,group, batch=["experiment"] ,resolution=0.4,n_pcs = 35,fig_flag=True,fig_dir = "/disk5/luosg/scRNAseq/output/result/BatchRemove")
            else:
                raise ValueError("group is wrong name, must be Muscle,Lung or Intestine")
            output_path = Path(f"{outdir}/{group}_class.h5ad")
            output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
            adata.write_h5ad(output_path)
        end = datetime.datetime.now()
        print("程序运行时间："+str((end-start).seconds/3600)+"h")

class AmbientRNAQC_GroupCombine:
    def __init__(self,marker_genes_Intestine,marker_genes_Lung,marker_genes_Muscle,
                 Intestine_annotation,Lung_annotation,Muscle_annotation):
        self.marker_genes_Intestine = marker_genes_Intestine
        self.marker_genes_Lung = marker_genes_Lung
        self.marker_genes_Muscle = marker_genes_Muscle
        self.Intestine_annotation = Intestine_annotation
        self.Lung_annotation = Lung_annotation
        self.Muscle_annotation = Muscle_annotation
        pass
    def step3GroupRun():
        start = datetime.datetime.now()
        infiles = {
            "Intestine": "/disk5/luosg/scRNAseq/output/result/h5ad/Intestine_qc_deambiendRNA.h5ad",
            "Lung": "/disk5/luosg/scRNAseq/output/result/h5ad/Lung_qc_deambiendRNA.h5ad",
            "Muscle": "/disk5/luosg/scRNAseq/output/result/h5ad/Muscle_qc_deambiendRNA.h5ad"
        }
        outdir = "/disk5/luosg/scRNAseq/output/result1"
        for group,infile in infiles.items():
            fig_dir = f"{outdir}/{group}"
            adata = sc.read_h5ad(infile)
            if group == 'Muscle':
                # continue
                adata_normalized = Run_Normalization(adata,
                                        group,
                                        batch="sample",
                                        fig_flag=True,
                                        resolution=1,
                                        n_pcs = 100,
                                        n_top_genes = 3000, 
                                        do_scale=True,
                                        fig_dir = fig_dir)
                adata_batch = Run_batchRemove(adata_normalized,
                                        group, 
                                        batch=["sample"] ,
                                        resolution=0.4,
                                        n_pcs = 10 ,
                                        fig_flag=True,
                                        fig_dir = fig_dir)
            elif group == 'Lung':
                continue
                # adata_normalized = Run_Normalization(adata,
                #                           group,batch="sample",
                #                           fig_flag=True,
                #                           resolution=0.4,
                #                           n_pcs = 100,
                #                           n_top_genes = 2500, 
                #                           do_scale=False,
                #                           fig_dir = fig_dir)
                # adata_batch = Run_batchRemove(adata_normalized,
                #                         group, 
                #                         batch=["sample"] ,
                #                         resolution=0.4,
                #                         n_pcs = 35,
                #                         fig_flag=True,
                #                         fig_dir = fig_dir)
            elif group == 'Intestine':
                continue
                # adata_normalized = Run_Normalization(adata,
                #                           group,
                #                           batch="sample",
                #                           fig_flag=True,
                #                           resolution=0.4,
                #                           n_pcs = 100,
                #                           n_top_genes = 2200, 
                #                           do_scale=True,
                #                           fig_dir = fig_dir)
                # adata_batch = Run_batchRemove(adata_normalized,
                #                         group, 
                #                         batch=["sample"] ,
                #                         resolution=0.4,
                #                         n_pcs = 35,
                #                         fig_flag=True,
                #                         fig_dir = fig_dir)
            else:
                raise ValueError("group is wrong name, must be Muscle,Lung or Intestine")
            output_path = Path(f"{outdir}/{group}_class.h5ad")
            output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
            adata_batch.write_h5ad(output_path)
        end = datetime.datetime.now()
        logging.info("程序运行时间："+str((end-start).seconds/3600)+"h")

    def step4GroupRun(self):
        start = datetime.datetime.now()
        infiles = {
            "Intestine": "/disk5/luosg/scRNAseq/output/result1/Intestine_class.h5ad",
            "Lung": "/disk5/luosg/scRNAseq/output/result1/Lung_class.h5ad",
            "Muscle": "/disk5/luosg/scRNAseq/output/result1/Muscle_class.h5ad"
        }
        outdir = "/disk5/luosg/scRNAseq/output/result1"
        for group,infile in infiles.items():
            adata = sc.read_h5ad(infile)
            fig_dir = Path(f"{outdir}/{group}/")
            fig_dir.mkdir(parents=True, exist_ok=True)
            table_dir = Path(f"{outdir}/{group}/table/")
            table_dir.mkdir(parents=True, exist_ok=True)
            output_path = Path(f"{outdir}/{group}_annotate.h5ad")
            output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
            if group == 'Intestine':
                # continue
                # Show_Markers(adata,fig_dir=str(fig_dir),table_dir=str(table_dir),out=group,fig_flag=True)
                # handful_annotate(adata,marker_genes=marker_genes_Intestine,fig_dir=fig_dir,out=group)
                adata = show_annotation(adata,cluster2annotation=self.Intestine_annotation,fig_dir=fig_dir,out = group)
                adata.obs["experiment"] = adata.obs["sample"].str.split("-").str[0]
                adata.write_h5ad(output_path)
            elif group == 'Lung':
                # continue
                # Show_Markers(adata,fig_dir=str(fig_dir),table_dir=str(table_dir),out=group,fig_flag=True)
                # handful_annotate(adata,marker_genes=marker_genes_Lung,fig_dir=fig_dir,out=group)
                adata = show_annotation(adata,cluster2annotation=self.Lung_annotation,fig_dir=fig_dir,out = group)
                adata.obs["experiment"] = adata.obs["sample"].str.split("-").str[0]
                adata.write_h5ad(output_path)
            elif group == 'Muscle':
                # continue
                # Show_Markers(adata,fig_dir=str(fig_dir),table_dir=str(table_dir),out=group,fig_flag=True)
                # handful_annotate(adata,marker_genes=marker_genes_Muscle,fig_dir=fig_dir,out=group)
                adata = show_annotation(adata,cluster2annotation=self.Muscle_annotation,fig_dir=fig_dir,out = group)
                adata.obs["experiment"] = adata.obs["sample"].str.split("-").str[0]
                adata.write_h5ad(output_path)
            else:
                raise ValueError("condition is wrong")

        end = datetime.datetime.now()
        print("程序运行时间："+str((end-start).seconds/3600)+"h")

class readCombine:
    __adataDict = {}
    def __init__(self,samplesDict:dict,indir:str,outdir:str,flag:str | None = None):
        self.samplesDict = samplesDict
        self.indir = indir
        self.outdir = outdir
        self.flag = flag
    def step1ReadRun(self):
        logging.info("start read scTE output and store in h5ad(optional)")
        indir = self.indir
        outdir = self.outdir
        flag = self.flag
        start = datetime.datetime.now()
        for group in samplesDict.keys():
            adata_list = []
            for sample in samplesDict[group]:
                # infile = f"{indir}/{group}/{sample}.csv"
                # adata = scTERead(infile,sample)
                adata = sc.read_h5ad(f"{indir}/{group}/{sample}.h5ad")
                adata.obs['experiment'] = sample.split("-")[0]
                adata_list.append(adata)
            adata = combine_group(adata_list)
            if flag == 'w_h5ad_1':
                output_path = Path(f"{outdir}/{group}/{group}_raw.h5ad")
                output_path.parent.mkdir(parents=True, exist_ok=True)
                adata.write_h5ad(output_path)
            self.__adataDict.update({group:adata})
        end = datetime.datetime.now()
        logging.info("start read scTE output and store in h5ad(optional) run time: "+str((end-start).seconds/3600)+"h")

    def step2QCRun(self):
        logging.info("start qc on combined adata")
        outdir = self.outdir
        start = datetime.datetime.now()
        for group,adata in self.__adataDict.items():

            adata = lowquality(adata,group,use_mt_outlier=False,low_count_outlier=3,low_gene_outlier=3,top20_outlier=3,
                            fig_flag = True,fig_dir = f"{outdir}/{group}/QC")
            adata = Doublet_scrub(adata,group,scrub_estimated_rate = 0.06,fig_flag = True,fig_dir = f"{outdir}/{group}/doublet/plot",
                    table_flag=True,table_dir=f"{outdir}/{group}/doublet/table")
            output_path = Path(f"{outdir}/{group}/{group}_qc.h5ad")
            output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
            adata.write_h5ad(output_path)
            self.__adataDict.update({group:adata})
        end = datetime.datetime.now()
        logging.info("start qc on combined adata run time:"+str((end-start).seconds/3600)+"h")
    
    def step3BatchRun(self):
        logging.info("start batch normalization and batch remove on combined adata")
        outdir = self.outdir
        start = datetime.datetime.now()
        for group,adata in self.__adataDict.items():
            adata = Run_Normalization(adata,group,batch='experiment',fig_flag=True,fig_dir = f"{outdir}/{group}/Normalization")
            adata = Run_batchRemove(adata,group, batch=['experiment'] ,fig_flag=True,fig_dir = f"{outdir}/{group}/BatchRemove")
            output_path = Path(f"{outdir}/{group}/{group}_class.h5ad")
            output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
            adata.write_h5ad(output_path)
        end = datetime.datetime.now()
        logging.info("start batch normalization and batch remove on combined adata run time:"+str((end-start).seconds/3600)+"h")
    

def step2QCRun(adataDict:dict[str,ad.AnnData],
               outdir:str):
    logging.info("start qc on combined adata")
    start = datetime.datetime.now()
    for group,adata in adataDict.items():
        if group == "Intestine":
            continue
            adata = lowquality(adata,group,mt_cutoff=8,mt_outlier=3,low_count_outlier=5,low_gene_outlier=5,top20_outlier=5,
                            fig_flag = True,fig_dir = f"{outdir}/{group}/QC")
            adata = Doublet_scrub(adata,group,scrub_estimated_rate = 0.06,fig_flag = True,fig_dir = f"{outdir}/{group}/doublet/plot",
                    table_flag=True,table_dir=f"{outdir}/{group}/doublet/table")
        elif group == "Lung":
            continue
            adata = lowquality(adata,group,mt_cutoff=8,mt_outlier=3,low_count_outlier=5,low_gene_outlier=5,top20_outlier=5,
                            fig_flag = True,fig_dir = f"{outdir}/{group}/QC")
            adata = Doublet_scrub(adata,group,scrub_estimated_rate = 0.20,fig_flag = True,fig_dir = f"{outdir}/{group}/doublet/plot",
                    table_flag=True,table_dir=f"{outdir}/{group}/doublet/table")
        elif group == "Muscle":
            # continue
            adata = lowquality(adata,group,mt_cutoff=30,mt_outlier=10,low_count_outlier=5,low_gene_outlier=5,top20_outlier=5,
                            n_genes_by_counts_cutoff=(200,99999999),total_counts_cutoff=(300,99999999),fig_flag = True,fig_dir = f"{outdir}/{group}/QC")
            adata = Doublet_scrub(adata,group,scrub_estimated_rate = 0.06,fig_flag = True,fig_dir = f"{outdir}/{group}/doublet/plot",
                    table_flag=True,table_dir=f"{outdir}/{group}/doublet/table")
        else:
            raise ValueError("the group is not be supported")
        output_path = Path(f"{outdir}/{group}/{group}_qc.h5ad")
        output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
        adata.write_h5ad(output_path)
    end = datetime.datetime.now()
    logging.info("start qc on combined adata run time:"+str((end-start).seconds/3600)+"h")

def step3BatchRun(
        adataDict:dict[str,ad.AnnData],
        outdir:str
):
    logging.info("start batch normalization and batch remove on combined adata")
    start = datetime.datetime.now()
    for group,adata in adataDict.items():
        if group == 'Intestine':
            continue
            adata = Run_Normalization(adata,group,batch='experiment',n_top_genes=2200,do_scale=True,n_pcs = 100,fig_flag=True,fig_dir = f"{outdir}/{group}/Normalization")
            adata = Run_batchRemove(adata,group, batch=['experiment'],resolution=0.4,n_pcs= 25,fig_flag=True,fig_dir = f"{outdir}/{group}/BatchRemove")
        elif group == 'Lung':
            continue
            adata = Run_Normalization(adata,group,batch='experiment',n_top_genes=2200,do_scale=True,n_pcs = 100,fig_flag=True,fig_dir = f"{outdir}/{group}/Normalization")
            adata = Run_batchRemove(adata,group, batch=['experiment'],resolution=0.4,n_pcs= 25,fig_flag=True,fig_dir = f"{outdir}/{group}/BatchRemove")
        elif group == 'Muscle':
            # continue
            adata = Run_Normalization(adata,group,batch='experiment',n_top_genes=2200,do_scale=True,n_pcs = 100,fig_flag=True,fig_dir = f"{outdir}/{group}/Normalization")
            adata = Run_batchRemove(adata,group, batch=['experiment'],resolution=0.4,n_pcs= 6,fig_flag=True,fig_dir = f"{outdir}/{group}/BatchRemove")
        else:
            raise ValueError("the group is not be supported")
        output_path = Path(f"{outdir}/{group}/{group}_class.h5ad")
        output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
        adata.write_h5ad(output_path)
    end = datetime.datetime.now()
    logging.info("start batch normalization and batch remove on combined adata run time:"+str((end-start).seconds/3600)+"h")

def step4AnnotateRun(
        adataDict:dict[str,ad.AnnData],
        outdir:str
):
    logging.info("annotate begin")
    start = datetime.datetime.now()
    for group,adata in adataDict.items():
        if group == "Intestine":
            continue
            # adata = Show_Markers(adata,fig_dir=f"{outdir}/{group}/annotate/plot",table_dir=f"{outdir}/{group}/annotate/table",out=group,fig_flag=True)
            # adata = handful_annotate(adata,marker_genes=marker_genes_Intestine,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
            adata = show_annotation(adata,cluster2annotation=Intestine_annotation,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
        elif group == "Lung":
            continue
            adata = Show_Markers(adata,fig_dir=f"{outdir}/{group}/annotate/plot",table_dir=f"{outdir}/{group}/annotate/table",out=group,fig_flag=True)
            adata = handful_annotate(adata,marker_genes=marker_genes_Lung,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
            adata = show_annotation(adata,cluster2annotation=Lung_annotation,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
        elif group == "Muscle":
            # continue
            # adata = Show_Markers(adata,fig_dir=f"{outdir}/{group}/annotate/plot",table_dir=f"{outdir}/{group}/annotate/table",out=group,fig_flag=True)
            # adata = handful_annotate(adata,marker_genes=marker_genes_Muscle,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
            adata = show_annotation(adata,cluster2annotation=Muscle_annotation,fig_dir=f"{outdir}/{group}/annotate/plot",out=group)
        else:
            raise ValueError("group is not be supported")
        output_path = Path(f"{outdir}/{group}/{group}_annotate.h5ad")
        output_path.parent.mkdir(parents=True, exist_ok=True) # Create only missing directories
        adata.write_h5ad(output_path)
    end = datetime.datetime.now()
    logging.info("annotate run time:"+str((end-start).seconds/3600)+"h")       
if __name__ == '__main__':
    
    # ambientRNAQC_GroupCombine = AmbientRNAQC_GroupCombine(marker_genes_Intestine,marker_genes_Lung,marker_genes_Muscle,
    #                          Intestine_annotation,Lung_annotation,Muscle_annotation)
    # ambientRNAQC_GroupCombine.step3GroupRun()
    # ambientRNAQC_GroupCombine.step4GroupRun()
    # readCombine = readCombine(samplesDict,
    #                           indir="/disk5/luosg/scRNAseq/output/h5ad",
    #                           outdir="/disk5/luosg/scRNAseq/output/combine",flag="w_h5ad_1")
    # readCombine.step1ReadRun()
    # readCombine.step2QCRun()
    # readCombine.step3BatchRun()
    # adataDict = {
    #     "Intestine": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_deambiendRNA.h5ad"),
    #     "Lung": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_deambiendRNA.h5ad"),
    #     "Muscle": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Muscle/Muscle_deambiendRNA.h5ad")
    # }
    # step2QCRun(adataDict,outdir="/disk5/luosg/scRNAseq/output/combine")
    # adataDict = {
    #     "Intestine": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_qc.h5ad"),
    #     "Lung": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_qc.h5ad"),
    #     "Muscle": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Muscle/Muscle_qc.h5ad")
    # }
    # step3BatchRun(adataDict,
    #               outdir="/disk5/luosg/scRNAseq/output/combine")
    adataDict = {
        "Intestine": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Intestine/Intestine_class.h5ad"),
        "Lung": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Lung/Lung_class.h5ad"),
        "Muscle": sc.read_h5ad("/disk5/luosg/scRNAseq/output/combine/Muscle/Muscle_class.h5ad")
    }
    step4AnnotateRun(adataDict,outdir="/disk5/luosg/scRNAseq/output/combine")


    