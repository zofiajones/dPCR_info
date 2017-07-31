VCF0[0]="/home/zjones/exomes_1/RKO.snp.vcf.gz"
VCF0[1]="/home/zjones/exomes_1/HCT116.snp.vcf.gz"
VCF0[2]="/home/zjones/exomes_1/SW48.snp.vcf.gz"

VCF1[0]="/home/zjones/exomes_1/RKO.indel.vcf.gz"
VCF1[1]="/home/zjones/exomes_1/HCT116.indel.vcf.gz"
VCF1[2]="/home/zjones/exomes_1/SW48.indel.vcf.gz"

VCF2[0]="/home/zjones/exomes/hd701_tru_fb_1.vcf.gz"
VCF2[1]="/home/zjones/exomes/GTB_348.vcf.gz"
VCF2[2]="/home/zjones/Sickkids/HD701.gatk.snp.indel.vcf.gz"

export PATH=~/anaconda3/bin/:$PATH

rm flanking_snps.txt
while read p;do
array=($p)
#echo $p
#echo ${array[0]}":"${array[1]}"-"${array[1]}
j=0
rm check_1.txt
for i in "${VCF0[@]}";do
#echo $i
~/hb/tools/tabix $i ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150] | python parse_multi_WES_bed_2.py | awk '{print $6}' | grep chr > check.txt
~/hb/tools/tabix ${VCF1[$j]} ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150] | python parse_multi_WES_bed_2.py | awk '{print $6}' | grep chr >> check.txt
a=$( cat check.txt | wc -l )
#cat check.txt
if [[ $a == 0 ]];then
label[j]=0
else
label[j]=$( cat check.txt | grep chr | awk '{printf "%s,", $1}' )
cat check.txt >> check_1.txt
fi
let j++
done

echo -e ${array[0]}"\t"${array[1]} >> flanking_snps.txt
cat check_1.txt | sort | uniq | sed 's/@/\t/g' >> flanking_snps.txt

b=$(~/hb/tools/tabix ~/exomes/hd701_tru_fb_1.vcf.gz  ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v ${array[1]} | python parse_multi_bed.py | awk '{print $6}' | grep -v ${array[1]} | awk '{printf "%s," , $1 }' )
c=$(~/hb/tools/tabix ${VCF2[1]}  ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v ${array[1]} | python parse_multi_bed_1.py | awk '{print $6}' | grep chr | awk '{ printf "%s,", $1 }' )
d=$(~/hb/tools/tabix ${VCF2[2]}  ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v ${array[1]} | python parse_multi_WES_bed_3.py | awk '{print $6}' | awk '{printf "%s," , $1}' )
#| python parse_multi_WES_bed_3.py

if [[ -z $b  ]];then
b=","
fi
if [[ -z $c  ]];then
c=","
fi
if [[ -z $d  ]];then
d=","
fi


echo -e ${array[0]}"\t"${array[1]}"\t"${label[0]}"\t"${label[1]}"\t"${label[2]}"\t"$b"\t"$c"\t"$d

done < to_check.txt


while read p;do
array=($p)
#echo $p
#echo ${array[0]}":"${array[1]}"-"${array[1]}
j=0
rm check_1.txt
for i in "${VCF0[@]}";do
#echo $i
~/hb/tools/tabix $i ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150] | python parse_multi_WES_bed_2.py | awk '{print $6}' | grep chr > check.txt
~/hb/tools/tabix ${VCF1[$j]} ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150] | python parse_multi_WES_bed_2.py | awk '{print $6}' | grep chr  >> check.txt
a=$(cat check.txt | wc -l )
#cat check.txt
if [[ $a == 0 ]];then
label[j]=0
else
label[j]=$( cat check.txt | awk '{printf "%s,", $1}' )
cat check.txt >> check_1.txt
fi
let j++
done

echo -e ${array[0]}"\t"${array[1]} >> flanking_snps.txt
cat check_1.txt | sort | uniq | sed 's/@/\t/g' >> flanking_snps.txt

b=$(~/hb/tools/tabix ${VCF2[0]} ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v $[${array[1]}-1] | python parse_multi_bed.py | awk '{print $6}' | grep -v ${array[1]} | awk '{printf "%s,",$1}END{ print "" }' | grep chr  )
c=$(~/hb/tools/tabix ${VCF2[1]}  ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v $[${array[1]}-1] | python parse_multi_bed_1.py | awk '{print $6}' | awk '{printf "%s," , $1}END{ print "" }' | grep chr )
d=$(~/hb/tools/tabix ${VCF2[2]}  ${array[0]}":"$[${array[1]}-150]"-"$[${array[1]}+150]  | grep -v $[${array[1]}-1] | python parse_multi_WES_bed_3.py | awk '{print $6}' | awk '{ printf "%s," , $1}')

if [[ -z $b  ]];then
b=","
fi
if [[ -z $c  ]];then
c=","
fi
if [[ -z $d  ]];then
d=","
fi


echo -e ${array[0]}"\t"${array[1]}"\t"${label[0]}"\t"${label[1]}"\t"${label[2]}"\t"$b"\t"$c"\t"$d

done < to_check_1.txt

