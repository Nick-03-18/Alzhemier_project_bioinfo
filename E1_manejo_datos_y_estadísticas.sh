#!/bin/bash
#Enfermedad de Alzheimer de Origen Tardío (LOAD): Análisis preliminar de estructura poblacional y frecuencia genotípica del gen APOE en el mundo
#Autores:J. Nicolas Avila M. & Sofia Lombana F.
#Datos: proyecto 1000 genomes fase 3 cromosoma 19 de 2504 individuos
VCF="ALL.chr19.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL="integrated_call_samples_v3.20130502.ALL.panel.txt"


#Creamos el directorio para colocar resultados de análisis
mkdir analysis


###Extraemos los SNPs de APOE a partir del VCF utilizando bcftools


bcftools view \
-r 19:45411941,19:45412079 \
-m2 -M2 -v snps \
$VCF \
-O z \
-o analysis/APOE.vcf.gz


##creamos un indice (.tbi) para el .vcf que construimos, -p indica que es un .vcf
tabix -p vcf analysis/APOE.vcf.gz
##Creamos el encabezado de lo que será la tabla resumen con frecuencia absoluta y relativa de los genotipos y las estadísticas poblacionales y se guarda como .tsv en la carpeta de resultados.
echo -e "SUPERPOP\tN\tE2E2\tE2E3\tE2E4\tE3E3\tE3E4\tE4E4\tE2E2_pct\tE2E3_pct\tE2E4_pct\tE3E3_pct\tE3E4_pct\tE4E4_pct\tPI\tHO\tHE\tFIS\tHWE_P" \
> analysis/APOE_summary.tsv


####Creamos un loop por superpoblaciones para calcular las frecuencias absoltuas y relativas de los genotipos y las estadisticas poblacionales de cada superpoblación


###Creamos el loop aplicado a las 5 superpoblaciones
for SUPER in AFR EUR EAS SAS AMR
do
##extraemos la muestra de cada superpoblación, con awk leemos el archivo de panel y hacemos que si la columna 3 de superpoblación es igual a la variable s que es la superpoblación que se esta trabajando entonces se imprima la columna 1 que contiene el sample (nombre de la meustra )y luego se guarde todo en un archivo con el nombre de la superpoblación .text
awk -v s="$SUPER" '$3==s {print $1}' $PANEL \
> analysis/${SUPER}.samples.txt

##Extraemos del .vcf general la información de los individuos de cada superpoblación y los guardamos en archivos separados
bcftools view \
-S analysis/${SUPER}.samples.txt \
analysis/APOE.vcf.gz \
-O z \
-o analysis/${SUPER}.vcf.gz

##creamos un indice (.tbi) para el .vcf que construimos, -p indica que es un .vcf
tabix -p vcf analysis/${SUPER}.vcf.gz

##construimos las frecuencias genotípicas

#extraemos los genotipos de cada variante
#extraemos información del vcf
bcftools query \
-f '[%GT\t]\n' \
analysis/${SUPER}.vcf.gz \
> analysis/${SUPER}.tmp

#en GENO queda guardado los conteos y porcentajes de los genotipos de las  superpoblaciones
#con awk se estan planteando unas indicaciones y acciones para armar las recuencias relativas y absolutas y estadisticas
#en el begin se ejecuta una vez antes de leer el archivo y inicializamos contadores, luego con NR toma cada linea y la divide usando taps en rs429358 y rs7412
#Ahora, con el END al final de leer el archivo se cuentan genotipos y porcentjes, vamos individuo por individuo separando los SNPs y poniendo reglas para detemrinar las combinaciones de estos que dan los alelos e2,e3 y e4 en cada copia
#luego rodenamos y definimos el genotipo
#luego contamos cada genotipo y sacamos el total por genotipo y global e imprimimos la tabla

GENO=$(awk '
BEGIN{
e2e2=0;e2e3=0;e2e4=0;
e3e3=0;e3e4=0;e4e4=0;
}
NR==1{split($0, rs429358, "\t")}
NR==2{split($0, rs7412, "\t")}

END{
for(i=1;i<=length(rs429358);i++){

split(rs429358[i], a, "[|/]")
split(rs7412[i], b, "[|/]")

if(a[1]==0 && b[1]==1){h1="E2"}
else if(a[1]==0 && b[1]==0){h1="E3"}
else if(a[1]==1 && b[1]==0){h1="E4"}

if(a[2]==0 && b[2]==1){h2="E2"}
else if(a[2]==0 && b[2]==0){h2="E3"}
else if(a[2]==1 && b[2]==0){h2="E4"}

if(h1 > h2){
tmp=h1
h1=h2
h2=tmp
}

g=h1 h2

if(g=="E2E2"){e2e2++}
else if(g=="E2E3"){e2e3++}
else if(g=="E2E4"){e2e4++}
else if(g=="E3E3"){e3e3++}
else if(g=="E3E4"){e3e4++}
else if(g=="E4E4"){e4e4++}

}

total=e2e2+e2e3+e2e4+e3e3+e3e4+e4e4

printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f",
total,
e2e2,e2e3,e2e4,e3e3,e3e4,e4e4,
(e2e2/total)*100,
(e2e3/total)*100,
(e2e4/total)*100,
(e3e3/total)*100,
(e3e4/total)*100,
(e4e4/total)*100)

}
' analysis/${SUPER}.tmp)

#luego en el loop sacamos las estadisticas por superpoblación
#calculamos el equilibrio de Hardey-Weinberg con vcftools
vcftools \
--gzvcf analysis/${SUPER}.vcf.gz \
--hardy \
--out analysis/${SUPER}
#luego calculamos el promedio de los p valor y se guarda en HWE, se lee con awk sin el header y hace que se sumen los p en la ultima columna y el numero de varinates
#luego en el end entonces hacemos que si hay variantes imprima la suma dividido en el n generando el promedio de los p valor de HWE
HWE=$(awk '
NR>1{
sum+=$NF
n++
}
END{
if(n>0) print sum/n
else print "NA"
}' analysis/${SUPER}.hwe)

#calculamos la diversidad nucleotidica con vcftools
vcftools \
--gzvcf analysis/${SUPER}.vcf.gz \
--site-pi \
--out analysis/${SUPER}

#luego calculamos el promedio de los pi entre los SNPs y se guarda en PI, se lee con awk sin el header y hace que se sumen los pi en la columna 3 y el numero de SNPs
#luego en el end entonces hacemos que si hay SNPs imprima la suma dividido en el n generando el promedio de los pi
PI=$(awk '
NR>1{
sum+=$3
n++
}
END{
if(n>0) print sum/n
else print "NA"
}' analysis/${SUPER}.sites.pi)

#Heterocigosidad
#calculamos heterocigosidad con vcftools
vcftools \
--gzvcf analysis/${SUPER}.vcf.gz \
--het \
--out analysis/${SUPER}_het

HO=$(awk '
NR>1{
ho = ($4 - $2)/$4
sum += ho
n++
}
END{
if(n>0) print sum/n
else print "NA"
}' analysis/${SUPER}_het.het)
#luego calculamos el promedio de las heterocigosidad esperadas entre los SNPs y se guarda en PI, se lee con awk sin el header y hace que se calcule la heterocigosidad esperadas ((total sitios menos homocigotos esperados) dividido en el total de sitios ) y el numero de SNPs
#luego en el end entonces hacemos que si hay SNPs imprima la suma dividido en el n generando el promedio de las heterocigosidades esperados
HE=$(awk '
NR>1{
he = ($4 - $3)/$4
sum += he
n++
}
END{
if(n>0) print sum/n
else print "NA"
}' analysis/${SUPER}_het.het)
#luego calculamos el promedio de las FIS entre los SNPs y se guarda en PI, se lee con awk sin el header y hace que se toma la FIS de la quinta columna  y el numero de SNPs
#luego en el end entonces hacemos que si hay SNPs imprima la suma dividido en el n generando el promedio de las FIS
FIS=$(awk '
NR>1{
sum += $5
n++
}
END{
if(n>0) print sum/n
else print "NA"
}' analysis/${SUPER}_het.het)

#escribimos los resutlados
echo -e "${SUPER}\t${GENO}\t${PI}\t${HO}\t${HE}\t${FIS}\t${HWE}" \
>> analysis/APOE_summary.tsv
#aquicerramos el loop
done

##Muestra la tabla final organizada

column -t analysis/APOE_summary.tsv
