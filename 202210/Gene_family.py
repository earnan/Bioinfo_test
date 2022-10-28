#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   Gene_family.py
#         Author:   yujie
#    Description:   Gene_family.py
#        Version:   1.0
#           Time:   2022/10/11 10:49:39
#  Last Modified:   2022/10/11 10:49:39
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from humre import *  # 正则
from icecream import ic  # 打印
import argparse  # 命令行
import linecache  # 大文件行读取
import os  # 目录路径
import pretty_errors  # 错误提示
import re  # 正则
import sys
import time
import copy  # 深度拷贝
#import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   Gene_family.py\n\
#         Author:   yujie\n\
#    Description:   Gene_family.py\n\
#        Version:   1.0\n\
#           Time:   2022/10/11 10:50:32\n\
#  Last Modified:   2022/10/11 10:50:32\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   Gene_family.py\n\
Function:\n\
1.常规使用\n\
1.1 -i [ ] -o [ ] \n\
2.其他使用\n\
2.1 -i [ ] -o [ ] \n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\Bioinfo_test\202210\Gene_family.py\n\
Path: /share/nas1/yuj/script/chloroplast/assembly/Gene_family.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[help_information]')
args = parser.parse_args()


# 读取输入文件，获取最终分支id  调用一次即可
def read_nwk(innwk_path):
    with open(innwk_path, 'r') as innwk_handle:
        for line in innwk_handle:
            if line.startswith('# IDs of nodes:'):
                #print((re.search(r'[a-zA-Z]', line.split(':')[1])).group(0))
                list_branch = re.findall(r'[A-Za-z]+', line.split(':')[1])
                break
    return list_branch


# 先以一行为例 传进去一行  需要每次调用
def get_list_node(ori_line):
    ori_list_node = re.findall(r'[\w]+', ori_line)
    n = 0
    list_node = []
    list_node_edit = []
    for i in ori_list_node:
        if i.find('_') >= 0:
            n += 1
            list_node_edit.append(i.split('_')[0])
            list_node.append(str(n)+','+i)
    # ic(list_node)  # 形如 1,li_36
    return list_node_edit, list_node


if __name__ == '__main__':
    #################################################################
    # 格式化成2016-03-20 11: 45: 39形式
    begin_time = time.time()
    start_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('Start Time : {}'.format(start_time))
    #################################################################

    print('---------------------------00初始化------------------------------------')
    # 文件初始化
    innwk_path = "E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\202210\\resultfile.cafe"
    outtab_path = "E:\\OneDrive\\jshy信息部\\Script\\Bioinfo_test\\202210\\statistics.xls"
    # 函数初始化
    createVar = locals()

    print('---------------------------01获取最终分支,动态生成所需字典------------------------------------')
    list_branch = read_nwk(innwk_path)
    print(list_branch)
    for branch in list_branch:
        createVar['dict_'+branch+'_allinfo'] = {}  # 该物种字典,用于写入
        createVar['list_'+branch+'shrunk'] = []  # 该物种收缩列表
        createVar['list_'+branch+'expanded'] = []  # 该物种扩张列表
        createVar['list_'+branch+'unchanged'] = []  # 该物种不变列表
        createVar['dict_'+branch +
                  '_allinfo']['shrunk'] = createVar['list_'+branch+'shrunk']
        createVar['dict_'+branch +
                  '_allinfo']['expanded'] = createVar['list_'+branch+'expanded']
        createVar['dict_'+branch +
                  '_allinfo']['unchanged'] = createVar['list_'+branch+'unchanged']

    print('---------------------------02逐行读取文件进行提取------------------------------------')
    innwk_contents = linecache.getlines(innwk_path)
    innwk_contents = innwk_contents[10:]
    for line in innwk_contents:
        family_name = line.split('\t')[0]
        ori_line = line.split('\t')[1]
        list_node_edit, list_node = get_list_node(ori_line)

        # 获取index   需要每次调用,放主函数里
        list_node_index = []  # 列表
        num_node_index = 0  # 复用
        num_node_len = 0  # 复用
        for i in list_node:  # '6,_22'  '7,_22'
            num_node = i.split(',')[-1]  # _22
            createVar['dict_'+i+'_allinfo'] = {}
            createVar['dict_'+i+'_allinfo']['number'] = num_node.split('_')[-1]
            if num_node.startswith('_'):
                createVar['dict_'+i+'_allinfo']['flag'] = 'child node'
            else:
                createVar['dict_'+i+'_allinfo']['flag'] = 'leaf node'
            start_index = num_node_index+num_node_len  # 起始索引，每次更新到 剩余序列在整行中的绝对索引
            ori_index = ori_line[start_index:].index(
                num_node)  # 后半截序列里的节点id 后面要修改
            num_node_index = ori_index+num_node_index+num_node_len  # 修改成在整行中的绝对索引
            createVar['dict_'+i+'_allinfo']['index'] = num_node_index  # 绝对索引
            createVar['dict_'+i +
                      '_allinfo']['depth'] = ori_line[:num_node_index].count('(')-ori_line[:num_node_index].count(')')+1
            # 对应深度层级=前半截序列左括号个数 - 右括号个数 + 1
            list_node_index.append(num_node_index)
            num_node_len = len(num_node)  # 节点长度

        # 判断基因家族是否收缩扩张
        for i in list_branch:
            branch = i
            i = list_node[list_node_edit.index(i)]
            current_number = createVar['dict_'+i+'_allinfo']['number']
            current_depth = createVar['dict_'+i+'_allinfo']['depth']
            for j in list_node:
                if createVar['dict_'+j+'_allinfo']['depth'] == current_depth-1 and createVar['dict_'+j+'_allinfo']['flag'] == 'child node':
                    if current_number < createVar['dict_'+j+'_allinfo']['number']:
                        createVar['list_'+branch +
                                  'shrunk'].append(family_name)
                    elif current_number == createVar['dict_'+j+'_allinfo']['number']:
                        createVar['list_'+branch +
                                  'unchanged'].append(family_name)
                    elif current_number > createVar['dict_'+j+'_allinfo']['number']:
                        createVar['list_'+branch +
                                  'expanded'].append(family_name)
    # 写入文件
    with open(outtab_path, 'w') as out_handle:
        out_handle.write('node\tshrunk\texpanded\n')
        for branch in list_branch:
            out_handle.write(branch+'\t')
            s = ''
            for i in createVar['dict_'+branch+'_allinfo']['shrunk']:
                s += (i+',')
            out_handle.write(s.rstrip(',')+'\t')
            for i in createVar['dict_'+branch+'_allinfo']['expanded']:
                s += (i+',')
            out_handle.write(s.rstrip(',')+'\n')
    ###############################################################
    end_time = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
    print('End Time : {}'.format(end_time))
    print('Already Run {}s'.format(time.time()-begin_time))
    print('Done')
    ###############################################################
