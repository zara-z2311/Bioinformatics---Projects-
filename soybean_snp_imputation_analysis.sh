# Activate the project conda environment
Source ~/miniconda3/bin/activate
conda create -n impute -c conda-forge -c bioconda bcftools samtools minimac4 r-base r-tidyverse -y
conda activate impute
# To deactivate an active environment, use                                                                                                               
#     $ conda deactivate

#QC checks
bcftools view \ 
-m2 -M2 \ 
-v snps \ 
data/soybean_genotype.vcf \ 
-Oz \ 
-o qc/soybean_chr20.snps.biallelic.vcf.gz 
 
bcftools index qc/soybean_chr20.snps.biallelic.vcf.gz 
 
bcftools view \ 
-i 'F_MISSING < 0.1 && MAF > 0.05' \ 
qc/soybean_chr20.snps.biallelic.vcf.gz \ 
-Oz \ 
-o qc/soybean_chr20.qc.vcf.gz 
 
bcftools index -f qc/soybean_chr20.qc.vcf.gz

#Variants counts were checked before and after filtering
bcftools view -H data/soybean_genotype.vcf | wc -l 
 
bcftools view -H qc/soybean_chr20.qc.vcf.gz | wc -l
bcftools query -l qc/soybean_chr20.qc.vcf.gz > sample_lists/all_samples.txt

#The 80/20 split
Rscript -e ' 

set.seed(123) 
 
ids <- readLines("sample_lists/all_samples.txt") 
 
n <- length(ids) 
 
ref_idx <- sample( 
1:n, 
size = floor(0.8 * n) 
) 
 
writeLines( 
ids[ref_idx], 
"sample_lists/ref_samples.txt" 
) 
 
writeLines( 
ids[-ref_idx], 
"sample_lists/target_samples.txt" 
) 
,

#Reference and target VCF files were generated

bcftools view \ 
-S sample_lists/ref_samples.txt \ 
qc/soybean_chr20.qc.vcf.gz \ 
-Oz \ 
-o data/reference.vcf.gz 
 
bcftools view \ 
-S sample_lists/target_samples.txt \ 
qc/soybean_chr20.qc.vcf.gz \ 
-Oz \ 
-o data/target_dense.vcf.gz 


#Thinning each target SNP datasets
#To check if the file is chromosome 20 only (file not indexed)
bcftools view -h soybean_genotype.vcf.gz | grep "##contig"
###contig=<ID=20,length=47846026>

#Create 5% SNP density at random 
awk 'BEGIN{srand(123)} {if(rand() <= 0.05) print}' all_positions.txt > thin5_positions.txt
wc -l thin5_positions.txt
#14371 thin5_positions.txt

#Create 5% VCF files and index that  file
bcftools view -T thin5_positions.txt target_dense.vcf.gz -Oz -o target_thinned5.vcf.gz
bcftools index target_thinned5.vcf.gz

#Check final SNP counts
bcftools view -H target_thinned5.vcf.gz | wc -l
#14371 which is 5%

#Phasing with Beagle 
#First phase the reference 
java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
  gt=reference.vcf.gz \
  out=reference_phased \
  nthreads=4

#Index it
bcftools index reference_phased.vcf.gz

#Phase 5% target SNP
java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
  gt=target_thinned5.vcf.gz \
  out=target_thinned5_phased \
  nthreads=4

#Index it 
bcftools index target_thinned5_phased.vcf.gz

#Imputation with Beagle 
#Same command pattern for each of 5%, 10%, 30%, 50%. The gp=true flag emits genotype probabilities in the output for downstream evaluation. Each run was timed with /usr/bin/time -v.

#5% density
cd ~/work30

/usr/bin/time -v -o beagle_5_time.txt \
  java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
    gt=target_thinned5.vcf.gz \
    ref=reference_phased.vcf.gz \
    out=imputed_beagle_5 \
    nthreads=4 \
    gp=true

bcftools index -f imputed_beagle_5.vcf.gz

#10% density
/usr/bin/time -v -o beagle_10_time.txt \
  java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
    gt=target_thinned10.vcf.gz \
    ref=reference_phased.vcf.gz \
    out=imputed_beagle_10 \
    nthreads=4 \
    gp=true

bcftools index -f imputed_beagle_10.vcf.gz

#30% density

/usr/bin/time -v -o beagle_30_time.txt \
  java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
    gt=target_thinned30.vcf.gz \
    ref=reference_phased.vcf.gz \
    out=imputed_beagle_30 \
    nthreads=4 \
    gp=true

bcftools index -f imputed_beagle_30.vcf.gz

#50% density

/usr/bin/time -v -o beagle_50_time.txt \
  java -Xmx8g -jar ~/tools/beagle.27Feb25.75f.jar \
    gt=target_thinned50.vcf.gz \
    ref=reference_phased.vcf.gz \
    out=imputed_beagle_50 \
    nthreads=4 \
    gp=true

bcftools index -f imputed_beagle_50.vcf.gz

#Verify imputed outputs

for DENSITY in 5 10 30 50; do
  echo "=== Beagle ${DENSITY}% ==="
  echo "  Samples:  $(bcftools query -l imputed_beagle_${DENSITY}.vcf.gz | wc -l)"
  echo "  Variants: $(bcftools view -H imputed_beagle_${DENSITY}.vcf.gz | wc -l)"
done


#Imputation with Minimac4
#MINIMAC download for IMPUTATION 

mkdir -p ~/tools/minimac4
cd ~/tools/minimac4
wget https://github.com/statgen/Minimac4/releases/download/v4.1.6/minimac4-4.1.6-Linux-x86_64.sh

#Install

bash minimac4-4.1.6-Linux-x86_64.sh
ls
find ~/tools/minimac4 -name "minimac4*"
~/tools/minimac4/bin/minimac4 --help

#Imputation 
#First convert phased reference VCF to MSAV format
#Check header first 

bcftools view -h reference_phased.vcf.gz | grep "##contig"
###contig=<ID=20>

#Minimac4 reheader

bcftools view -h reference_phased.vcf.gz > reference_header.txt

sed 's/##contig=<ID=20>/##contig=<ID=20,length=47846026>/' reference_header.txt > reference_header.fixed.txt

bcftools reheader \
  -h reference_header.fixed.txt \
  -o reference_phased.fixed.vcf.gz \
  reference_phased.vcf.gz

bcftools index -f reference_phased.fixed.vcf.gz

bcftools view -h reference_phased.fixed.vcf.gz | grep "##contig"

#Same with the SNP target files
#Fix header for the 5% phased target file

bcftools view -h target_thinned5_phased.vcf.gz > target5_header.txt

sed 's/##contig=<ID=20>/##contig=<ID=20,length=47846026>/' \
target5_header.txt > target5_header.fixed.txt

bcftools reheader \
  -h target5_header.fixed.txt \
  -o target_thinned5_phased.fixed.vcf.gz \
  target_thinned5_phased.vcf.gz

bcftools index -f target_thinned5_phased.fixed.vcf.gz

#Check

bcftools view -h target_thinned5_phased.fixed.vcf.gz | grep "##contig"
###contig=<ID=20,length=47846026>

#Fix header for the 10% phased target file

bcftools view -h target_thinned10_phased.vcf.gz > target10_header.txt

sed 's/##contig=<ID=20>/##contig=<ID=20,length=47846026>/' \
target10_header.txt > target10_header.fixed.txt

bcftools reheader \
  -h target10_header.fixed.txt \
  -o target_thinned10_phased.fixed.vcf.gz \
  target_thinned10_phased.vcf.gz

bcftools index -f target_thinned10_phased.fixed.vcf.gz

#Check

bcftools view -h target_thinned10_phased.fixed.vcf.gz | grep "##contig"
###contig=<ID=20,length=47846026>

#Change the reference VCF to MSAV format 

~/tools/minimac4/bin/minimac4 \
  --compress-reference reference_phased.fixed.vcf.gz \
  -o reference_panel.msav

#Imputation with Minimac4
#Check the files first

ls -lh reference_panel.msav
ls -lh target_thinned5_phased.fixed.vcf.gz
ls -lh target_thinned10_phased.fixed.vcf.gz
ls -lh target_thinned30_phased.fixed.vcf.gz
ls -lh target_thinned50_phased.fixed.vcf.gz

#Minimac4 for 5% SNP dataset

~/tools/minimac4/bin/minimac4 \
  reference_panel.msav \
  target_thinned5_phased.fixed.vcf.gz \
  --all-typed-sites \
  --format GT,DS \
  --output minimac_5percent.vcf.gz \
  --output-format vcf.gz \
  --threads 4

#Index 5% output

bcftools index -f minimac_5percent.vcf.gz

#Minimac4 for 10% SNP dataset

~/tools/minimac4/bin/minimac4 \
  reference_panel.msav \
  target_thinned10_phased.fixed.vcf.gz \
  --all-typed-sites \
  --format GT,DS \
  --output minimac_10percent.vcf.gz \
  --output-format vcf.gz \
  --threads 4

#Index 10% output

bcftools index -f minimac_10percent.vcf.gz

#Check

bcftools view -H minimac_5percent.vcf.gz | wc -l
#544997

bcftools view -H minimac_10percent.vcf.gz | wc -l
#544997

bcftools view -H minimac_30percent.vcf.gz | wc -l
#544997

bcftools view -H minimac_50percent.vcf.gz | wc -l
#544997

#Data Analysis 
#Glimpse 2 concordance analysis

echo "20 data/target_dense.vcf.gz data/target_dense.vcf.gz results/target5.minimac.imputed.vcf.gz" \
> results/target5_minimac_concordance.txt
tools/GLIMPSE2_concordance_static \
  --gt-val \
  --gt-tar \
  --threads 2 \
  --bins 0 0.001 0.005 0.01 0.05 0.1 0.2 0.5 \
  --input results/target5_minimac_concordance.txt \
  --output results/target5_minimac_concordance

#R Analysis and visualisation
#Load package

library(tidyverse)

#Summary Tables:
#Genotype Error and Concordance Summary Tables 
#Find GLIMPSE2 genotype error output files

error_files <- list.files( 
"results", 
pattern = "concordance.error.grp.txt.gz$", 
full.names = TRUE 
)

#Read and combine all genotype error files

error_df <- map_dfr(error_files, function(file) { 
 
df <- read_tsv( 
file, 
comment = "#", 
col_names = TRUE, 
show_col_types = FALSE 
) 
 
df$file <- basename(file) 
 
df$thinning <- case_when( 
str_detect(df$file, "target5") ~ 5, 
str_detect(df$file, "target10") ~ 10, 
str_detect(df$file, "target30") ~ 30, 
str_detect(df$file, "target50") ~ 50 
) 
 
df$method <- case_when( 
str_detect(df$file, "beagle") ~ "Beagle", 
str_detect(df$file, "minimac") ~ "Minimac4" 
) 
 
return(df) 
})

#Create weighted mean genotype error summary

summary_df <- error_df %>% 
filter(n_genotypes > 0) %>% 
group_by(thinning, method) %>% 
summarise( 
mean_error = weighted.mean( 
error_rate, 
n_genotypes, 
na.rm = TRUE 
), 
.groups = "drop" 
)

#Save the genotype error summary table

write.csv( 
summary_df, 
"results/genotype_error_summary.csv", 
row.names = FALSE 
)

#Concordance was calculated as 100 minus genotype error

concordance_df <- summary_df 
 
concordance_df$concordance <- 100 - concordance_df$mean_error 
 
write.csv( 
concordance_df, 
"results/concordance_summary.csv", 
row.names = FALSE 
)

#R-squared summary table
#Locate GLIMPSE2 R-squared output files

rsq_files <- list.files( 
"results", 
pattern = "concordance.rsquare.grp.txt.gz$", 
full.names = TRUE 
)

#Read and combine all R-squared files

rsq_df <- map_dfr(rsq_files, function(file) { 
 
df <- read_tsv( 
file, 
comment = "#", 
col_names = TRUE, 
show_col_types = FALSE 
) 
 
df$file <- basename(file) 
 
df$thinning <- case_when( 
str_detect(df$file, "target5") ~ 5, 
str_detect(df$file, "target10") ~ 10, 
str_detect(df$file, "target30") ~ 30, 
str_detect(df$file, "target50") ~ 50 
) 
 
df$method <- case_when( 
str_detect(df$file, "beagle") ~ "Beagle", 
str_detect(df$file, "minimac") ~ "Minimac4" 
) 
 
return(df) 
})

#Create mean R-squared summary table

summary_rsq <- rsq_df %>% 
group_by(thinning, method) %>% 
summarise( 
mean_rsq = mean(rsquare, na.rm = TRUE), 
.groups = "drop" 
) 
 
write.csv( 
summary_rsq, 
"results/rsquare_summary.csv", 
row.names = FALSE 
)

#Allele frequency error summary table
#Save allele frequency error table for plotting

write.csv( 
error_df, 
"results/allele_frequency_error_table.csv", 
row.names = FALSE 
)

#Concordance across thinning levels

concordance <- read_csv("results/concordance_summary.csv") 
 
concordance$thinning <- factor( 
  concordance$thinning, 
  levels = c(5, 10, 30, 50) 
) 
 
p <- ggplot(concordance, 
            aes(x = thinning, 
                y = concordance, 
                color = method, 
                group = method)) + 
  geom_line(linewidth = 1.2) + 
  geom_point(size = 3) + 
  labs( 
    title = "Imputation concordance across thinning levels", 
    x = "Thinning level (%)", 
    y = "Genotype concordance", 
    color = "Method" 
  ) + 
  theme_minimal() + 
  theme( 
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), 
    axis.title.x = element_text(size = 17, face = "bold"), 
    axis.title.y = element_text(size = 17, face = "bold"), 
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14), 
    legend.title = element_text(size = 15, face = "bold"), 
    legend.text = element_text(size = 13) 
  ) 
 
print(p) 
 
ggsave( 
  "results/concordance_across_thinning_levels_updated.png", 
  p, 
  width = 8, 
  height = 6, 
  dpi = 300 
)

#R-squared across thinning levels

rsq <- read_csv("results/rsquare_summary.csv") 
 
rsq$thinning <- factor( 
rsq$thinning, 
levels = c(5, 10, 30, 50) 
) 
 
p <- ggplot(rsq, 
aes(x = thinning, 
y = mean_rsq, 
color = method, 
group = method)) + 
geom_line(linewidth = 1.2) + 
geom_point(size = 3) + 
labs( 
title = "Imputation r² across thinning levels", 
x = "Thinning level (%)", 
y = "Mean r²", 
color = "Method" 
) + 
theme_minimal() + 
theme( 
plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), 
axis.title.x = element_text(size = 17, face = "bold"), 
axis.title.y = element_text(size = 17, face = "bold"), 
axis.text.x = element_text(size = 14), 
axis.text.y = element_text(size = 14), 
legend.title = element_text(size = 15, face = "bold"), 
legend.text = element_text(size = 13) 
) 
 
print(p) 
 
ggsave( 
"results/R2_across_thinning_levels_updated.png", 
p, 
width = 8, 
height = 6, 
dpi = 300 
)

#Genotype error across thinning levels

error_df <- read_csv("results/genotype_error_summary.csv") 
 
error_df$thinning <- factor( 
error_df$thinning, 
levels = c(5, 10, 30, 50) 
) 
 
p <- ggplot(error_df, 
aes(x = thinning, 
y = mean_error, 
color = method, 
group = method)) + 
geom_line(linewidth = 1.2) + 
geom_point(size = 3) + 
labs( 
title = "Imputation genotype error across thinning levels", 
x = "Thinning level (%)", 
y = "Weighted genotype error", 
color = "Method" 
) + 
theme_minimal() + 
theme( 
plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), 
axis.title.x = element_text(size = 17, face = "bold"), 
axis.title.y = element_text(size = 17, face = "bold"), 
axis.text.x = element_text(size = 14), 
axis.text.y = element_text(size = 14), 
legend.title = element_text(size = 15, face = "bold"), 
legend.text = element_text(size = 13) 
) 
 
print(p) 
 
ggsave( 
"results/genotype_error_across_thinning_levels_updated.png", 
p, 
width = 8, 
height = 6, 
dpi = 300 
)

#Error across allele frequency bins

af_error <- read_csv("results/allele_frequency_error_table.csv") 
 
p <- ggplot(af_error, 
aes(x = mean_AF, 
y = error_rate, 
color = method, 
group = method)) + 
geom_line(linewidth = 1) + 
geom_point(size = 2) + 
facet_wrap(~ thinning) + 
labs( 
title = "Imputation error across allele frequency bins", 
x = "Mean allele frequency", 
y = "Error rate", 
color = "Method" 
) + 
theme_minimal() + 
theme( 
plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), 
axis.title.x = element_text(size = 17, face = "bold"), 
axis.title.y = element_text(size = 17, face = "bold"), 
axis.text.x = element_text(size = 13), 
axis.text.y = element_text(size = 13), 
strip.text = element_text(size = 14, face = "bold"), 
legend.title = element_text(size = 15, face = "bold"), 
legend.text = element_text(size = 13) 
) 
 
print(p) 
 
ggsave( 
"results/error_across_allele_frequency_bins_updated.png", 
p, 
width = 10, 
height = 8, 
dpi = 300 
)

#Concordance heatmap

df <- read_csv("results/concordance_summary.csv") 
 
df$thinning <- factor( 
df$thinning, 
levels = c(5, 10, 30, 50) 
) 
 
p <- ggplot(df, 
aes(x = thinning, 
y = method, 
fill = concordance)) + 
geom_tile(color = "white") + 
geom_text(aes(label = round(concordance, 3)), 
size = 5) + 
scale_fill_gradient( 
low = "lightblue", 
high = "darkblue" 
) + 
labs( 
title = "Genotype imputation concordance heatmap", 
x = "SNP retention (%)", 
y = "Imputation tool", 
fill = "Concordance" 
) + 
theme_bw() + 
theme( 
plot.title = element_text(hjust = 0.5, size = 22, face = "bold"), 
axis.title.x = element_text(size = 17, face = "bold"), 
axis.title.y = element_text(size = 17, face = "bold"), 
axis.text.x = element_text(size = 14), 
axis.text.y = element_text(size = 14), 
legend.title = element_text(size = 15, face = "bold"), 
legend.text = element_text(size = 13) 
) 
 
print(p) 
 
ggsave( 
"results/imputation_heatmap_new.png", 
p, 
width = 6, 
height = 4, 
dpi = 300 
)
