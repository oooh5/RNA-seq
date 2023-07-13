# RNA-seq
1 SRA toolkit
从https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit下载SRA Toolkit
mac系统本地下载：MacOS 64 bit architecture
linux下载：CentOS Linux 64 bit architecture可本地下载后上传或者wget下载
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.5/sratoolkit.3.0.5-centos_linux64.tar.gz
tar zxvf sratoolkit.3.0.5-centos_linux64.tar.gz #解压

export PATH=$PATH:$./sratoolkit.3.0.5-centos_linux64/bin  #设置环境变量
nohup prefetch SRR12563916 &  
nohup sratoolkit-centos_linux64/bin/prefetch --option-file sra.txt &
nohup prefetch --option-file sra.txt &


#nohup命令
jobs#查看当前挂起的命令
ps#查看PID号码
top#查看服务器内所有挂起
kill PID号码#利用ps查看PID后用kill杀死命令
sra.txt文件是从GEO databse下载的



2 SRA转fastq
mkdir rawdata
mv SRR*/SRR*.sra ./rawdata
rm -r SRR*

#单个sra
fastq-dump --split-e --gzip xxx.sra -O ./
#批量sra， 'nohup' 需要一个单字命令及其参数 - 而不是 shell 循环结构。您必须使用:
nohup sh -c 'for i in *.sra; do echo $i; fastq-dump --split-e --gzip $i -O ./; done' > fastq_dump.log &

3 fastqc
mkdir fastqc
nohup fastqc -q -t 3 -o ./ SRR*.fastq.gz >qc.log & 
nohup multiqc ./*.zip -o ./ > ./multiqc.log &
#-q:静默运行，系统不会实时报告软件的运行情况
#-t:设置所使用的核数，20或者30
#-o:输出文件路径
#后面就是两个rawdata的fastq文件，可以同时检查正反向的数据

# or  使用管道符和shell命令来进行批量处理 
ls rawdata/*fq | while read id; do nohup fastqc -o ~/lailingling/fastqc/ -q -t 2 $id & done



4 cutadapter
pwd 
#/mnt/datadisk/renyuhan/lailingling/YTHDC1_seq_data
mkdir trim_galore
echo "SRR12563916" >> rawdata_name
echo "SRR12563917" >> rawdata_name
echo "SRR12563918" >> rawdata_name
echo "SRR12563919" >> rawdata_name
echo "SRR12563920" >> rawdata_name
echo "SRR12563921" >> rawdata_name
echo "SRR12563922" >> rawdata_name
echo "SRR12563923" >> rawdata_name
echo "SRR12563924" >> rawdata_name
echo "SRR12563925" >> rawdata_name
echo "SRR12563926" >> rawdata_name

vim  trim_galore.sh
###vim输入内容
rawdata=~/lailingling/YTHDC1_seq_data/rawdata
trim_galore=~/lailingling/YTHDC1_seq_data/trim_galore
cat rawdata_name | while read id
do
  trim_galore --phred33 --stringency 3 --length 36 --paired --fastqc --clip_R1 5 --three_prime_clip_R1 15 --clip_R2 5 --three_prime_clip_R2 15 -o ${trim_galore} ${rawdata}/${id}_1.fastq.gz ${rawdata}/${id}_2.fastq.gz
done
######运行
nohup sh trim_galore.sh >trim_galore.log &
#出现错误可以通过cat来查看错误记录
cat trim_galore.log
nohup multiqc *.zip -o ./trim_fastqc > ./multiqc_t.log &



#-q：设定Phred quality score阈值，默认为20。
#--phred33：选择-phred33或者-phred64，表示测序平台使用的Phred quality score。
#--adapter：输入adapter序列。也可以不输入，Trim Galore!会自动寻找可能性最高的平台对应的adapter。自动搜选的平台三个，也直接显式输入这三种平台，即--illumina、--nextera和--small_rna。
#--stringency：设定可以忍受的前后adapter重叠的碱基数，默认为1。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到。
#--length：设定输出reads长度阈值，小于设定值会被抛弃。
#--paired：对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个会被同样抛弃，而不管是否达到标准。
#--retain_unpaired：对于双端测序结果，一对reads中，如果一个read达到标准，但是对应的另一个要被抛弃，达到标准的read会被单独保存为一个文件。
#--gzip和--dont_gzip：清洗后的数据zip打包或者不打包。
#--output_dir：输入目录。需要提前建立目录，否则运行会报错。
#--trim-n : 移除read一端的reads。
#--clip_R1、 --three_prime_clip_R1、--clip_R2、--three_prime_clip_R2：R1 R2文件前后两段切除的碱基数目

5 STAR
###构建index
mkdir index_mus
cd index_mus
nohup wget -c https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz > genome.log&
nohup wget -c https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz >gtf.log&
nohup gunzip *.gz>unzip.log &

wget https://cloud.biohpc.swmed.edu/index.php/s/hisat2-210-Linux_x86_64/download
unzip download
echo 'export PATH=~/lailingling/YTHDC1_seq_data/index_mus/hisat2-2.2.1:$PATH' >> ~/lailingling/YTHDC1_seq_data/index_mus/.bashrc
source .bashrc


hisat2_extract_exons.py Mus_musculus.GRCm39.109.gtf > exons_Mus.txt
hisat2_extract_splice_sites.py Mus_musculus.GRCm39.109.gtf >splice_sites_Mus.txt
#利用hisat2提取外显子位置到exons_Mus.txt ，可变剪切位置到ss_Mus.txt
nohup hisat2-build -p 6 --ss splice_sites_Mus.txt --exon exons_Mus.txt Mus_musculus.GRCm39.dna.primary_assembly.fa GRCm39 > hisat2_bulid.log &
#利用hisat-build构建小麦转录组比对的索引文件， -p 8为 8线程，--ss 可变剪切文件，--exon 外显子文件，后接参考基因组，index文件前缀为GRCm39，，，，，需要二十多分钟

###比对
mkdir hisat2
echo "SRR12563916" >> rawdata_name.txt
echo "SRR12563917" >> rawdata_name.txt
echo "SRR12563918" >> rawdata_name.txt
echo "SRR12563919" >> rawdata_name.txt
echo "SRR12563920" >> rawdata_name.txt
echo "SRR12563921" >> rawdata_name.txt
echo "SRR12563922" >> rawdata_name.txt
echo "SRR12563923" >> rawdata_name.txt
echo "SRR12563924" >> rawdata_name.txt
echo "SRR12563925" >> rawdata_name.txt
echo "SRR12563926" >> rawdata_name.txt

cd ~/lailingling/YTHDC1_seq_data/hisat2
echo 'export PATH=~/lailingling/YTHDC1_seq_data/index_mus/hisat2-2.2.1:$PATH' >> ./.bashrc
source .bashrc
##用一下代码批量分析
#写一个bash脚本
vim hisat2.sh
##把下面的代码复制进去
cat rawdata_name.txt | while read id
do
    hisat2 -p 8 -x ~/lailingling/YTHDC1_seq_data/index/GRCm39 -1 ~/lailingling/YTHDC1_seq_data/trim_galore/trim_fqgz/${id}_1_val_1.fq.gz -2 ~/lailingling/YTHDC1_seq_data/trim_galore/trim_fqgz/${id}_2_val_2.fq.gz -S ~/lailingling/YTHDC1_seq_data/hisat2/${id%%_*}.sam > ${id%%_*}.log
done

##后台运行脚本
nohup bash hisat2.sh &

#-p 2: 指定线程数为2，表示使用两个处理器核心进行并行处理。
#-x ${hisat_index}: 指定hisat2索引的路径。
#-1 ${id}_1_val.fq.gz -2 ${id}_2_val.fq.gz: 指定输入的配对的FASTQ文件。
#2>${id%%_*}.log: 将标准错误输出重定向到名为${id%%_*}.log的日志文件中，${id%%_*}是去除文件名中第一个下划线及其后的部分，例如将文件名ERR1698194_1.fastq.gz转换为ERR1698194.log。
#|: 管道符号，用于将hisat2的输出传递给下一个命令。
#samtools sort -@ 2 -o ${id%%_*}.bam: 使用samtools sort对hisat2的输出进行排序并生成BAM文件，-@ 2指定线程数为2，-o ${id%%_*}.bam指定输出的BAM文件名，${id%%_*}同样是去除文件名中第一个下划线及其后的部分



samtools sort -@ 2 -o ${id%%_*}.bam >Hisat2Index.sh.log

vim sam2bam.sh
##把下面的代码复制进去
cat rawdata_name.txt | while read id
do
samtools view -@ 2 -bS ${id}.sam | samtools sort -@ 2 -o ${id%%_*}.bam > ${id%%_*}_samtools.log
done
#######

nohup bash sam2bam.sh &


6 定量
###each---可跑，为单个txt文件
mkdir featurecount__each
cd ~/lailingling/YTHDC1_seq_data/hisat2

for bam_file in *.bam; do
    output_file=~/lailingling/YTHDC1_seq_data/featurecounts_each/${bam_file%.bam}.txt
    nohup featureCounts -p -T 4 -a ~/lailingling/YTHDC1_seq_data/index_mus/Mus_musculus.GRCm39.109.gtf -o "$output_file" "$bam_file"
done > ~/lailingling/YTHDC1_seq_data/featurecounts_each/featurecounts.log 2>&1 &

multiqc *.summary -o ./ > ./multiqc_feature.log 2>&1

for feature_data in *.txt;do
less -S "$feature_data" | grep -v '#' | cut -f 1,7- > ${feature_data%.txt}_rawcounts.txt
done
head SRR16_counts.txt | column -t

###all----跑这个跑这个 别看上面那个 
mkdir featurecounts_all
cd featurecounts_all
gtf=~/lailingling/YTHDC1_seq_data/index_mus/Mus_musculus.GRCm39.109.gtf
inputdir=~/lailingling/YTHDC1_seq_data/hisat2
nohup featureCounts -T 6 -p -t exon -g gene_id -a $gtf -o all.id.txt $inputdir/*.bam > featurecounts3.log 2>&1 &

multiqc all.id.txt.summary -o ./
less -S all.id.txt |grep -v '#' |cut -f 1,7- |sed  's#/mnt/datadisk/renyuhan/lailingling/YTHDC1_seq_data/hisat2/##g' |sed 's#.bam##g' > raw_counts.txt
head raw_counts.txt  |column -t

#gtf=~/database/human/GRCh38.108/Homo_sapiens.GRCh38.108.gtf.gz：将GTF文件的路径保存在变量 gtf 中。请确保该路径是正确的。
#inputdir=$HOME/project/No4exam_batch/data/Mapping/Hisat2/：将输入文件目录的路径保存在变量 inputdir 中。请确保该路径是正确的。
#-T 3：指定线程数为3，你可以根据自己的需求进行调整。
#-p：启用并行模式，可以提高运行速度。
#-t exon：计数的特征类型为外显子。
#-g gene_id：使用基因ID进行计数。
#-a "$gtf"：指定GTF文件的路径。这里使用了双引号将变量 gtf 包裹起来，以防止路径中包含空格等特殊字符引发问题。
#-o all.id.txt：指定输出文件的名称为 all.id.txt，结果将保存在该文件中。
#"$inputdir"/*.sorted.bam：指定输入的BAM文件路径，使用通配符 *.sorted.bam 匹配该目录下的所有以 .sorted.bam 结尾的文件。
#> ~/lailingling/YTHDC1_seq_data/featurecounts/featurecounts.log 2>&1 的作用是将标准输出和标准错误都重定向到 featurecounts.log 文件中。这样做的目的是将命令的输出和错误信息都记录在同一个日志文件中，方便后续查看。如果发生错误或警告，它们将被写入 featurecounts.log 文件，以便进行排查和分析。


7 DEseq2
#10.10.140.182:8787      账号：renyuhan 密码：123

getwd()
setwd("~/lailingling/YTHDC1_seq_data/featurecounts_all/")


#install.packages("BiocManager")
#BiocManager::install("DEseq2")
library(DESeq2)
##导入counts
DC1_counts <- read.table("~/lailingling/YTHDC1_seq_data/featurecounts_all/raw_counts2.txt",row.names = 1 ,header=T,sep = '\t')  #若遇到报错“多字节字符串错”，加上参数 fileEncoding="UTF-16LE"即可
new_colnames <- c("f/f-1","f/f-2", "cKO-1","cKO-2","cKO-3", "wtRes-1","wtRes-2","wtRes-3", "W378A-1","W378A-2","W378A-3")  # 替换为您想要的新列名
colnames(DC1_counts) <- new_colnames

##构建分组信息
condition <- factor(c(rep("f/f",2), rep("cKO",3),rep("wtRes",3),rep("W378A",3))) 
coldata <- data.frame(row.names = colnames(DC1_counts), condition)
coldata    #显示coldata值,看看分组信息与真实数据是否一致，很关键！！！！

##构建dds矩阵：就是利用上面的counts.txt文件（即对象counts）和分组信息（即对象coldata）构建
dds <- DESeqDataSetFromMatrix(countData=DC1_counts, colData=coldata, design=~condition)


##Hierarchical clustering

#BiocManager::install('factoextra') 
library('factoextra')
vsd <- vst(dds,blind = TRUE)  #参数blind=TRUE是为了不把样本分组信息考虑在内——即以无偏的方式进行转换
sampleDists <- dist(t(assay(vsd))) #dist()计算样本间的距离
res1 <- hcut(sampleDists, k = 3, stand = FALSE, hc_method ="average" ) #hcut()函数进行层次聚类分析，其中sampleDists是距离矩阵，k = 2表示将样本分为2个簇，stand = FALSE表示不进行标准化，hc_method = "average"表示使用平均链接法进行层次聚类。
fviz_dend(res1,
          rect_fill = T,
          # 字体大小
          cex = 1,
          # 字体颜色
          color_labels_by_k=T,
          # 平行放置
          horiz=T)#rect_fill = T表示填充矩形，cex = 1设置字体大小为1，color_labels_by_k = T表示按簇进行标签颜色的区分，horiz = T表示水平放置树状图。








