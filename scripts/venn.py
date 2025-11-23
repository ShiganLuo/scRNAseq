from ast import expr
from pathlib import Path
import sys
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import logging
from upsetplot import from_contents, plot
basePath = Path(__file__).resolve()
baseDir = basePath.parent
sys.path.append(str(baseDir / "pyvenn"))
import pyvenn.venn as venn
from matplotlib_venn import venn3,venn2
# conda: scRNAseq_rpy2
logging.basicConfig(

    level=logging.INFO,

    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',

    datefmt='%Y-%m-%d %H:%M:%S'

)
def Intersection(infiles:list):
    for file in infiles:
        df = pd.read_csv(file,sep="\t")
    
def search_files_with_regex(directory, pattern):
    """
    Recursively searches a directory for files whose names match a given regex pattern.

    Args:
        directory (str): The starting directory for the search.
        pattern (str): The regular expression string used to match filenames.

    Returns:
        list: A list containing the full paths of all matching files.
    """
    matching_files = []
    # 编译正则表达式以提高效率
    regex = re.compile(pattern)

    # os.walk() 会递归遍历目录树
    for root, dirs, files in os.walk(directory):
        for file in files:
            # 使用正则表达式匹配文件名
            if regex.search(file):
                # 拼接完整路径并添加到列表中
                matching_files.append(os.path.join(root, file))
    
    return matching_files

def classify_files_by_regex(file_paths, patterns):
    """
    Classifies file paths based on a dictionary of regex patterns that match filenames.

    Args:
        file_paths (list): A list of full file paths.
        patterns (dict): A dictionary where keys are category names (str) and
                         values are corresponding regex pattern strings.

    Returns:
        dict: A dictionary where keys are category names and values are lists
              of file paths that match the corresponding pattern.
    """
    classified_files = {category: [] for category in patterns}

    for path in file_paths:
        file_name = os.path.basename(path)
        for category, pattern in patterns.items():
            if re.match(pattern, file_name):
                classified_files[category].append(path)
                break  # 匹配到一个类别后就停止，避免重复归类
    
    return classified_files

def classify_files_by_tissue_and_cell(file_paths):
    """
    Classifies file paths based on tissue and cell type extracted from the path.

    Args:
        file_paths (list): A list of full file paths.

    Returns:
        dict: A nested dictionary with the following structure:
              {tissue_name: {cell_type: [list_of_matching_file_paths]}}.
    """
    classified_data = {}
    
    for path in file_paths:
        dir_name, file_name = os.path.split(path)
        
        # 提取组织名称（倒数第三个目录）
        tissue = os.path.basename(os.path.dirname(dir_name))
        
        # 提取细胞类型名称
        # 找到第一个 '10XSC3_' 和 '-combined_group' 之间的字符串
        try:
            start_index = file_name.find('10XSC3_') + len('10XSC3_')
            end_index = file_name.find('-combined_group', start_index)
            cell_type = file_name[start_index:end_index]
            
            # 将 'T_cytotoxic_cell' 这样的名称标准化，例如，去除尾部的 '_cell'
            # cell_type = cell_type.replace('_cell', '')

            if tissue not in classified_data:
                classified_data[tissue] = {}

            if cell_type not in classified_data[tissue]:
                classified_data[tissue][cell_type] = []
            
            classified_data[tissue][cell_type].append(path)
            
        except ValueError:
            print(f"Warning: Could not parse cell type from file name: {file_name}")
            continue
            
    return classified_data

if __name__ == '__main__':
    dir = "/disk5/luosg/scRNAseq/output/result/DEG"
    pattern = r'^combined_group.*TEsubfamily\.tsv$'
    files = search_files_with_regex(dir,pattern)
    # print(files)
    group = classify_files_by_tissue_and_cell(files)
    targetGroup = ["CKO_vs_WT","E2_vs_CKO","TRA_vs_CKO"]
    for tissue,cellDcit in group.items():
        for cell,files in cellDcit.items():
            # print(cell,files)
            logging.info(f"tissue: {tissue}\ncell: {cell}\nfiles: {files}")
            geneList = []
            labels = []
            geneDict = {}
            for file in files:
                dir_name, file_name = os.path.split(file)
                start_index1 = file_name.find('combined_group') + len('combined_group')
                end_index1 = file_name.find('_',start_index1)
                start_index2 = file_name.find('-combined_group', start_index1) + len('-combined_group')
                end_index2 = file_name.find('_',start_index2)
                exp = file_name[start_index1:end_index1]
                ctl = file_name[start_index2:end_index2]
                groupName = f"{exp}_vs_{ctl}"
                if groupName not in targetGroup:
                    continue
                labels.append(groupName)
                
                df = pd.read_csv(file,sep="\t")

                geneList.append(df['subfamily'].to_list())
                geneDict[groupName] = df['subfamily'].to_list()
            logging.info(f"geneList's length: {len(geneList)}")
            if len(geneList) == 2:
                venn2([set(g) for g in geneList], set_labels=tuple(labels),set_colors=('#4e79a7', '#f28e2b'),alpha=0.7)
                # labels = venn.get_labels(geneList, fill=['number'])
                # fig, ax = venn.venn2(labels, names=["2C_GFPp_SA1","TPS"],dpi=300,figsize=(15,15),fontsize=30)
            elif len(geneList) == 3:
                venn3([set(g) for g in geneList], set_labels=tuple(labels),set_colors=('#4e79a7', '#f28e2b', "#8BFBC3"),alpha=0.7)
                # labels = venn.get_labels(geneList, fill=['number'])
                # fig, ax = venn.venn3(labels, names=["2C_GFPp_SA1","TPS","c"],dpi=300,figsize=(15,15),fontsize=30)
            elif len(geneList) == 1 or len(geneList) == 0:
                continue           
            else:
                sets = from_contents(geneDict)
                plot(geneList)
                # print(f"{tissue} {cell} failed")
                # continue
            # plt.text(0.5, -0.1, cell, ha='center', va='top', fontsize=12)
            cellName = ""
            for i in cell.split("_"):
                cellName = cellName + " " + i
            # plt.title(f"DEG gene overlap in {cell}",fontsize = 16)
            plt.tight_layout()
            outfile = f"{dir}/{tissue}/plot/venn_DEG_gene_{cellName}.png"
            plt.savefig(outfile,dpi=300)
            plt.close()
    logging.info("end")