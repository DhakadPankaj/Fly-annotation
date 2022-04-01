

for i in `cat sra.txt`
do 

cd $i
ls *.gz | xargs -n 1 gunzip -t 2>&1 | cut -f 2 -d: - | xargs -t -n 1 rm

if [ -f "$i"_1.fastq.gz ]
then 
echo "$i exists"
else
echo "download"
fi
cd ..
done
