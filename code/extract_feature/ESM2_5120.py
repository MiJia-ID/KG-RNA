import esm
import argparse
import torch
import numpy as np
import os

p = argparse.ArgumentParser()
p.add_argument('--device', type=str, default='3', help='')
p.add_argument('--output_path', type=str, default='/home/mijia/KG-RNA/code/data/esm2/', help='')
p.add_argument('--fasta_file', type=str, default='"/home/mijia/KG-RNA/code/data/Test_RNA.fasta"', help='Input fasta file')
p.add_argument('--esm2_model', type=str, default='/home/mijia/app/esm2/esm2_t48_15B_UR50D.pt', help='The path of esm2 model parameter.')
args = p.parse_args()

cuda_index = 3  # 这里可以设置为 0, 1, 2, ... 对应不同的 GPU

# 检查是否有足够的 GPU 设备
if torch.cuda.is_available() and cuda_index < torch.cuda.device_count():
    device = torch.device(f"cuda:{cuda_index}")
else:
    device = torch.device("cpu")
    print("CUDA is not available or the specified index is out of range. Falling back to CPU.")
# device = 'cpu'
esm_model, alphabet = esm.pretrained.load_model_and_alphabet_local(args.esm2_model)
batch_converter = alphabet.get_batch_converter()
esm_model = esm_model.to(device)
print(f"ESM2 运行设备: {next(esm_model.parameters()).device}")
esm_model.eval()

def get_esm(fasta_file,output_path):
    print(f"Input fasta file: {fasta_file}")
    print(f"Output path: {output_path}")
    ID_list = []
    seq_list = []
    with open(fasta_file, "r") as f:
        lines = f.readlines()
    for line in lines:
        line = line.strip()  # 移除行末尾的换行符和空白字符
        print(line)
        if line and line[0] == ">":  # 检查行是否非空且是否为ID行
            ID_list.append(line[1:])  # 去除ID行的大于号
            print("ID：",ID_list[-1])
        elif line:  # 检查行是否非空
            seq_list.append(" ".join(list(line)))  # 将非空行的序列加入序列列表中
            print("seq:", line)
            batch_labels, batch_strs, batch_tokens = batch_converter([('tmp', line)])
            batch_tokens = batch_tokens.to(device)
            try:
                with torch.no_grad():
                    results = esm_model(batch_tokens, repr_layers=[33], return_contacts=True)
                    token_representations = results["representations"][33][0].cpu().numpy().astype(np.float16)
                    esm_embed = token_representations[1:len(line) + 1]
                    # 打印特征的尺寸大小
                    print("特征尺寸大小:", esm_embed.shape)
                    print("len:", len(line))
                    np.save(os.path.join(output_path, ID_list[-1]), esm_embed)
            except Exception as e:
                # 捕获异常并打印异常信息，然后继续到下一个循环迭代
                print(e)

get_esm(args.fasta_file,args.output_path)

#CUDA_VISIBLE_DEVICES=2 nohup python ESM2_5120.py > ESM2.log 2>&1 &