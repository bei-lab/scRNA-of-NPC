module load java/1.8.0_60
java -Xmx5G -jar Aracne.jar -e path/sample/sample_matrix.txt  -o path/sample/output_sample --tfs path/sample/sample_TF.txt --pvalue 1E-8 --seed 1 \
--calculateThreshold

for i in {1..100}
do
java -Xmx5G -jar Aracne.jar -e path/sample/sample_matrix.txt  -o path/sample/output_sample --tfs path/sample/sample_TF.txt --pvalue 1E-8 --seed $i --threads 8
done


java -Xmx5G -jar Aracne.jar -o path/sample/output_sample --consolidate --consolidatepvalue 0.01 --threads 8