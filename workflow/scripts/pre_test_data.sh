mkdir -p ./test_fastqs

# Subset the first 400,000 lines (100k reads) of each FASTQ
# Note: Update the path '/path/to/your/liver_fastq/' to the real location of your files
for file in /projects/prjs1616/gatk_paper_demo/data/liver_fastq/*.fastq; do
    filename=$(basename "$file")
    echo "Subsetting $filename..."
    head -n 400000 "$file" > "./test_fastqs/$filename"
done
